import time
import sys
import traceback
import numpy as np
from collections import OrderedDict
from tqdm import tqdm

from . import packet as pkt
from . import conversion as conv
from . import constants as const
import gena.filter as flt
from . import pqt_io
import os


class _TqdmLogging(tqdm):
    """tqdm subclass that mirrors progress bar output to a log file."""

    def __init__(self, *args, log_file=None, **kwargs):
        self._log_file = log_file
        self._postfix_dict = {}
        super().__init__(*args, **kwargs)

    def set_postfix(self, ordered_dict=None, refresh=True, **kwargs):
        """Override to stash the raw dict before tqdm stringifies it."""
        pf = OrderedDict([] if ordered_dict is None else ordered_dict)
        for key in sorted(kwargs.keys()):
            pf[key] = kwargs[key]
        self._postfix_dict = pf
        super().set_postfix(ordered_dict=ordered_dict, refresh=refresh, **kwargs)

    def display(self, msg=None, pos=None):
        super().display(msg=msg, pos=pos)
        fd = self.format_dict
        if not fd:
            return
        postfix = self._postfix_dict
        step_parts = [f"{k}={v}" for k, v in postfix.items()]
        step_str = f" | {' | '.join(step_parts)}" if step_parts else ""
        elapsed = fd.get('elapsed', 0)
        remaining = fd.get('remaining', 0)
        eta_str = f"<{remaining:.0f}s" if remaining is not None and remaining >= 0 else ""
        line = (
            f"Chunk {fd.get('n', 0)}/{fd.get('total', 0)} "
            f"[{elapsed:.0f}s{eta_str}]"
            f"{step_str}\n"
        )
        if self._log_file is not None:
            self._log_file.write(line)
            self._log_file.flush()
        else:
            # verbose mode: print to terminal
            sys.__stdout__.write(line)
            sys.__stdout__.flush()

    def __del__(self, *args, **kwargs):
        try:
            super().__del__(*args, **kwargs)
        except AttributeError:
            pass


#from memory_profiler import profile
#@profile
#snid = 1: default is use orbit parameters
def decode(raw_file,output_folder='./',by_day=True,snid=1,verbose=False):
    params = const.init_param(snid)
    start_time_tot = time.time()
    directory = os.path.dirname(output_folder)
    if not os.path.exists(directory):
        os.makedirs(directory)

    # Redirect stdout/stderr to log files so they don't interrupt tqdm bar.
    # The progress bar itself writes to the real stderr (terminal).
    if not verbose:
        raw_basename = os.path.basename(raw_file)
        std_log = open(os.path.join(output_folder, f'{raw_basename}.stdout.log'), 'a')
        err_log = open(os.path.join(output_folder, f'{raw_basename}.stderr.log'), 'a')
        _orig_stdout = sys.stdout
        _orig_stderr = sys.stderr
        sys.stdout = std_log
        sys.stderr = err_log
    else:
        std_log = None
        err_log = None
        _orig_stdout = sys.stdout
        _orig_stderr = sys.stderr

    #if not use a single file instead of sorting by day
    #need to init all pqt files at first
    if not by_day:
        pqt_io.init_pqt(output_folder)

    #io_flag  = os.environ.get('GENA_IO_FLAG')
    #if io_flag == None:
    #     gio.gena_h5_init(output_folder)
    # if io_flag == 'pqt':
    #     #its a prefix _{pdf,l1a,l1b,l2a,l2b}.{pqt}
    #     directory = os.path.dirname(output_folder)
    #     if not os.path.exists(directory):
    #         os.makedirs(directory)
    #     gio.init_pqt(output_folder)
    #     #gio.init_pqt2(output_folder)
    # if io_flag == 'hdf':
    #     gio.init_hdf(output_folder)

    f = open(raw_file,'rb')
    fsize = f.seek(0,2)
    f.seek(0,0)
    cbytes = 2064*const.chunk_size

    chunk_size = np.min([cbytes,fsize])

    if np.mod(fsize,2064) == 0:
        print(f'file may contain {fsize/2064} packets')
    else:
        print(f'file may contain {fsize/2064} packets, but not aligned')
    if fsize < cbytes:
        print(f'Decode all packets in one chunk')
    else:
        print(f'Decode {const.chunk_size} packets in one chunk')

    total_chunks = int(np.ceil(fsize / chunk_size))

    read_raw_time = 0
    find_head_time = 0
    parse_pdf_time = 0
    parse_df0_time = 0
    save_pdf_time = 0
    save_l1_time = 0
    save_l1a_time = 0
    save_l1b_time = 0
    save_l2_time = 0

    chunk_start = 0
    total_packets = 0
    written_id = 0
    chunk_num = 0

    if not verbose:
        pbar = _TqdmLogging(total=total_chunks, desc="Decoding", unit="chunk",
                            file=sys.__stderr__, log_file=std_log)
    else:
        pbar = _TqdmLogging(total=total_chunks, desc="Decoding", unit="chunk",
                            file=open(os.devnull, 'w'), log_file=None)
    try:
        try:
            while(True):
                chunk_num += 1

                # Read chunk
                start_time = time.time()
                chunk_size = np.min((chunk_size,fsize-chunk_start))
                read_buffer = np.memmap(raw_file, offset = chunk_start, shape=(chunk_size), mode='r')
                end_time = time.time()

                read_raw_time += end_time - start_time
                pbar.set_postfix(step="Read", step_t=f"{end_time - start_time:.2f}s")

                # Find frame headers
                start_time = time.time()
                frame_index = pkt.find_frame_header_re(read_buffer)
                end_time = time.time()
                find_head_time += end_time - start_time

                pbar.set_postfix(step="FindHeaders", step_t=f"{end_time - start_time:.2f}s")

                # Read packets
                start_time = time.time()
                packets,last,residual = pkt.read_packet_opt(read_buffer,frame_index)
                # [MOD] validate packets: filter out structurally invalid ones (false headers with correct 2064 spacing but wrong structure)
                #valid_packets = []
                #for p in packets:
                #    mask = pkt.validate_packets(p)
                #    removed = len(p) - np.sum(mask)
                #    if removed > 0:
                #        pbar.write(f"  Removed {removed}/{len(p)} structurally invalid packets")
                #    if np.any(mask):
                #        valid_packets.append(p[mask])
                #packets = valid_packets

                #print("next start at:", last)
                len_packets  = np.sum([len(p) for p in packets])
                total_packets = total_packets + len_packets

                end_time = time.time()

                read_raw_time += end_time - start_time
                pbar.set_postfix(step="ReadPackets", step_t=f"{end_time - start_time:.2f}s")

                for p in packets:

                    # Parse PDF
                    start_time = time.time()
                    pdf = pkt.packets_2_pdf(p).drop_duplicates(subset=['Timestamp']).reset_index(drop=True)
                    #pdf = pkt.packets_2_pdf(p)
                    end_time = time.time()
                    parse_pdf_time += end_time - start_time
                    pbar.set_postfix(step="Parse PDF", step_t=f"{end_time - start_time:.2f}s")

                    # PDF data is always parsed because downstream event products
                    # need its attitude and engineering parameters.  This flag
                    # controls only whether the PDF product is written to disk.
                    if const.save_pdf == 1:
                        start_time = time.time()
                        if not by_day:
                            pqt_io.save_pqt(pdf,output_folder,'pdf')
                        else:
                            pqt_io.save_pqt_daily(pdf,output_folder,'pdf')
                        end_time = time.time()
                        save_pdf_time += end_time - start_time
                        pbar.set_postfix(step="Save PDF", step_t=f"{end_time - start_time:.2f}s")

                    if const.pdf_only == 0:
                        # Parse DF0
                        start_time = time.time()
                        df0 = pkt.packets_2_df0(p)
                        del p

                        end_time = time.time()
                        pbar.set_postfix(step="DF0", step_t=f"{end_time - start_time:.2f}s")

                        # L0 → L1
                        start_time = time.time()

                        # L0 → L1: Convert raw TDC to nanoseconds
                        df1_a = conv.l0_2_l1(df0,'a',snid).drop_duplicates().reset_index(drop=True)
                        df1_b = conv.l0_2_l1(df0,'b',snid).drop_duplicates().reset_index(drop=True)

                        end_time = time.time()
                        pbar.set_postfix(step="L0→L1", step_t=f"{end_time - start_time:.2f}s")

                        # Save L1 data
                        if (const.save_l1 == 1):
                            start_time = time.time()
                            if not by_day:
                                pqt_io.save_pqt(df1_a, output_folder, 'l1_a')
                                pqt_io.save_pqt(df1_b, output_folder, 'l1_b')
                            else:
                                pqt_io.save_pqt_daily(df1_a, output_folder, 'l1_a')
                                pqt_io.save_pqt_daily(df1_b, output_folder, 'l1_b')
                            end_time = time.time()
                            save_l1_time += end_time - start_time
                            pbar.set_postfix(step="Save L1", step_t=f"{end_time - start_time:.2f}s")
                        del df0

                        # L1 → L1A
                        start_time = time.time()
                        # L1 → L1A: Select 8-channel events and interpolate PDF parameters
                        df1a_a = conv.l1_2_l1a(df1_a, pdf).drop_duplicates().reset_index(drop=True)
                        df1a_b = conv.l1_2_l1a(df1_b, pdf).drop_duplicates().reset_index(drop=True)
                        end_time = time.time()
                        pbar.set_postfix(step="L1→L1A", step_t=f"{end_time - start_time:.2f}s")

                        # Save L1A data
                        if (const.save_l1a == 1):
                            start_time = time.time()
                            if not by_day:
                                pqt_io.save_pqt(df1a_a, output_folder, 'l1a_a')
                                pqt_io.save_pqt(df1a_b, output_folder, 'l1a_b')
                            else:
                                pqt_io.save_pqt_daily(df1a_a, output_folder, 'l1a_a')
                                pqt_io.save_pqt_daily(df1a_b, output_folder, 'l1a_b')
                            end_time = time.time()
                            save_l1a_time += end_time - start_time
                            pbar.set_postfix(step="Save L1A", step_t=f"{end_time - start_time:.2f}s")


                
                        if (const.repair == 1):
                            # L1A → L1B
                            start_time = time.time()
                            # L1A → L1B: Calculate TOF, Energy, positions
                            df1b_a = conv.l1a_2_l1b(df1a_a, params,prob='a',method=const.method).drop_duplicates().reset_index(drop=True)
                            df1b_b = conv.l1a_2_l1b(df1a_b, params,prob='b',method=const.method).drop_duplicates().reset_index(drop=True)
                            end_time = time.time()
                            pbar.set_postfix(step="L1A→L1B", step_t=f"{end_time - start_time:.2f}s")


                            # Save L1B data
                            if (const.save_l1b == 1):
                                start_time = time.time()
                                if not by_day:
                                    pqt_io.save_pqt(df1b_a, output_folder, 'l1b_a')
                                    pqt_io.save_pqt(df1b_b, output_folder, 'l1b_b')
                                else:
                                    pqt_io.save_pqt_daily(df1b_a, output_folder, 'l1b_a')
                                    pqt_io.save_pqt_daily(df1b_b, output_folder, 'l1b_b')
                                end_time = time.time()
                                save_l1b_time += end_time - start_time
                                pbar.set_postfix(step="Save L1B", step_t=f"{end_time - start_time:.2f}s")

                            # L1B → L2
                            start_time = time.time()
                            df2b_a = conv.l1_2_l2(df1b_a, prob='a').drop_duplicates().reset_index(drop=True)
                            df2b_b = conv.l1_2_l2(df1b_b, prob='b').drop_duplicates().reset_index(drop=True)
                            end_time = time.time()
                            pbar.set_postfix(step="L1B→L2", step_t=f"{end_time - start_time:.2f}s")
                
                            # Save L2 data
                            if (const.save_l2 == 1):
                                start_time = time.time()
                                if not by_day:
                                    pqt_io.save_pqt(df2b_a, output_folder, 'l2b_a')
                                    pqt_io.save_pqt(df2b_b, output_folder, 'l2b_b')
                                else:
                                    pqt_io.save_pqt_daily(df2b_a, output_folder, 'l2b_a')
                                    pqt_io.save_pqt_daily(df2b_b, output_folder, 'l2b_b')
                                end_time = time.time()
                                save_l2_time += end_time - start_time
                                pbar.set_postfix(step="Save L2", step_t=f"{end_time - start_time:.2f}s")

                            del  df1b_a, df1b_b, pdf, df2b_a, df2b_b
                        else:
                            #generate pure 1a l2 and soft l2 and strict l2
                
                            start_time = time.time()
                            df2_a = conv.l1_2_l2(df1a_a, params = params, prob='a').drop_duplicates().reset_index(drop=True)
                            df2_b = conv.l1_2_l2(df1a_b, params = params, prob='b').drop_duplicates().reset_index(drop=True)
                            end_time = time.time()
                            pbar.set_postfix(step="L1A -> L2", step_t=f"{end_time - start_time:.2f}s")
                
                            # L1A → L2S (soft quality selection)
                            start_time = time.time()
                            func_soft_l1_a = flt.soft_std(dir='all',params=params, prob='a')
                            func_soft_l1_b = flt.soft_std(dir='all',params=params, prob='b')
                            #need improve style
                            df2s_a = conv.l1_2_l2(df1a_a.loc[func_soft_l1_a(df1a_a)], params=params, prob='a').drop_duplicates().reset_index(drop=True)
                            df2s_b = conv.l1_2_l2(df1a_b.loc[func_soft_l1_b(df1a_b)], params=params, prob='b').drop_duplicates().reset_index(drop=True)
                
                            end_time = time.time()
                            pbar.set_postfix(step="L1A → L1S->L2S", step_t=f"{end_time - start_time:.2f}s")
 

                            # L1A → L2G (strict + spatial quality selection)
                            start_time = time.time()
                            func_gold_l1_a = flt.golden_std(dir='all',params=params, prob='a')
                            func_gold_l1_b = flt.golden_std(dir='all',params=params, prob='b')
                            #need improve style
                            df2g_a = conv.l1_2_l2(df1a_a.loc[func_gold_l1_a(df1a_a)], params=params, prob='a').drop_duplicates().reset_index(drop=True)
                            df2g_b = conv.l1_2_l2(df1a_b.loc[func_gold_l1_b(df1a_b)], params=params, prob='b').drop_duplicates().reset_index(drop=True)
                            # Golden products are the strict-filter result plus
                            # the spatial detector-window check.  Do not replace
                            # them with the softer L2S products here.
                            df2g_a = df2g_a.loc[flt.onboard_eff(df2g_a,c_size=[150,60])].reset_index(drop=True)
                            df2g_b = df2g_b.loc[flt.onboard_eff(df2g_b,c_size=[150,60])].reset_index(drop=True)
                
                            end_time = time.time()
                            pbar.set_postfix(step="L1A → L1G->L2G", step_t=f"{end_time - start_time:.2f}s")
 
                            # Save L2 data
                            if (const.save_l2 == 1):
                                start_time = time.time()
                                if not by_day:
                                    pqt_io.save_pqt(df2_a, output_folder, 'l2_a')
                                    pqt_io.save_pqt(df2_b, output_folder, 'l2_b')
                                    pqt_io.save_pqt(df2s_a, output_folder, 'l2s_a')
                                    pqt_io.save_pqt(df2s_b, output_folder, 'l2s_b')
                                    pqt_io.save_pqt(df2g_a, output_folder, 'l2g_a')
                                    pqt_io.save_pqt(df2g_b, output_folder, 'l2g_b')
                                else:
                                    pqt_io.save_pqt_daily(df2_a, output_folder, 'l2_a')
                                    pqt_io.save_pqt_daily(df2_b, output_folder, 'l2_b')
                                    pqt_io.save_pqt_daily(df2s_a, output_folder, 'l2s_a')
                                    pqt_io.save_pqt_daily(df2s_b, output_folder, 'l2s_b')
                                    pqt_io.save_pqt_daily(df2g_a, output_folder, 'l2g_a')
                                    pqt_io.save_pqt_daily(df2g_b, output_folder, 'l2g_b')
                                end_time = time.time()
                                save_l2_time += end_time - start_time
                                pbar.set_postfix(step="Save L2", step_t=f"{end_time - start_time:.2f}s")

                            del  df2_a, df2_b, df2s_a, df2s_b, df2g_a, df2g_b
                        del df1_a, df1_b, df1a_a, df1a_b, pdf
                        parse_df0_time += end_time - start_time

                        written_id +=1
                del packets
                pbar.set_postfix(step=f"Packets:{len_packets}")

                chunk_start = chunk_start + last
                pbar.update(1)

                if(chunk_start == fsize or chunk_start + residual == fsize):
                    break
        except:
            traceback.print_exc()
            if not verbose:
                traceback.print_exc(file=std_log)
                std_log.flush()
            sys.stderr.flush()
            raise
    finally:
        summary_parts = [
            '=='*10 + 'SUMMARY' + '=='*10,
            f'Num packets processed {total_packets}',
            f"Load RAW to memory costs: {read_raw_time:.2f} s",
            f"Searching frame index costs: {find_head_time:.2f} s",
            f"Parse and dump parameters costs: {parse_pdf_time:.2f} s",
            f"Parse and dump df0 costs: {parse_df0_time:.2f} s",
            f"Save PDF costs: {save_pdf_time:.2f} s",
            f"Save L1 costs: {save_l1_time:.2f} s",
            f"Save L1A costs: {save_l1a_time:.2f} s",
            f"Save L1B costs: {save_l1b_time:.2f} s",
            f"Save L2 costs: {save_l2_time:.2f} s",
        ]
        total_save = save_pdf_time + save_l1_time + save_l1a_time + save_l1b_time + save_l2_time
        summary_parts.append(f"Total save time: {total_save:.2f} s")
        end_time_tot = time.time()
        summary_parts.append(f"Total time costs: {end_time_tot - start_time_tot:.2f} s")

        for line in summary_parts:
            if not verbose:
                std_log.write(line + '\n')
            print(line, file=_orig_stdout)
        if not verbose:
            std_log.flush()

        pbar.close()
        if not verbose:
            sys.stdout = _orig_stdout
            sys.stderr = _orig_stderr
            std_log.close()
            err_log.close()

if __name__ == "__main__":
    rawfile = sys.argv[1]
    outfolder = sys.argv[2]
    decode(rawfile, outfolder)
