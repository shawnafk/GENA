import numpy as np
import pandas as pd
from . import structure as st
#bug!!!
'''
STEP 1: find header 0x(E2 25 0C 4E), decimal(226 37 12 28) position to locate each frame
input: 
    f: file object, open in binary mode
    or/
    buffer: bytes array f.read to memory to speed up searching (avoid many seek and read operations)
    buffer can be numpy memmap object, avoid acutally read to memory
return:
    an numpy array of index of each packet header

performance: (on WLS2,  AMD Ryzen 9 9950X 4.3GHz, 64GB DDR5 6400MT/s RAM) 
    for BIG file: 5153871936 bytes file, read cost 4.6 s 
    0. find_frame_header_seq_u: 2min 29s;
    1. find_frame_header_seq: 2min 38s (5.4 s in c++);
    2. find_frame_header_re: 1.37 s; (2.3s with memmap object)

    for SMALL file: 19727712 bytes file, read cost 9.12 ms
    0. find_frame_header_seq_u: 534 ms;
    1. find_frame_header_seq: 451 ms;
    2. find_frame_header_re: 3.78 ms;

    regular expression is the fastest way to find header, 
    fread have some implicit buffer mechanism,
    so there should be no difference between seq_u and seq,
    but seq_u is a little faster than seq when buffer is large.

    Just use regular expression to find header.
'''

'''
def size_of_file(f):
    f.seek(0, 2)
    size = f.tell()
    f.seek(0, 0)
    return size

def find_frame_header_seq_u(f):
    matches = []
    filesize = size_of_file(f)
    remaining_bytes = filesize
    while remaining_bytes > 3:
        headdata = f.read(4)
        if (headdata[0] == 226) and (headdata[1] == 37) and (headdata[2] == 12) and (headdata[3] == 78):
            idx = f.tell() - 4
            matches.append(idx)
            this_read = 4
        elif (headdata[1] == 226) and (headdata[2] == 37) and (headdata[3] == 12):
            f.seek(-3, 1)
            this_read = 1
        elif (headdata[2] == 226) and (headdata[3] == 37):
            f.seek(-2, 1)
            this_read = 2
        elif headdata[3] == 226:
            f.seek(-1, 1)
            this_read = 3
        else:
            this_read = 4
        remaining_bytes = remaining_bytes - this_read
    #add last as file size to keep the last frame when exam 2064 length by matches[1:] - headseek_num[:-1]
    matches.append(filesize)
    return np.array(matches)


#notes: change assignment to headdata by slicing buffer, extral cost, 620 ms -> 526 ms on 19727712 bytes;
def find_frame_header_seq(buffer):
    #set empty list to store header index
    matches = []
    #moving point to the start
    #set remaining bytes
    remaining_bytes = len(buffer)
    #frame header is a 4 bytes array: 0x- E2 25 0C 4E corresponding to 226, 37, 12 and 78
    #we read and calculate remaining bytes
    cur = 0
    while remaining_bytes > 3:
        if (buffer[cur] == 226) and (buffer[cur+1] == 37) and (buffer[cur+2] == 12) and (buffer[cur+3] == 78):
            matches.append(cur)
            this_read = 4
            cur += this_read
        elif (buffer[cur+1] == 226) and (buffer[cur+2] == 37) and (buffer[cur+3] == 12):
            this_read = 1
            cur += this_read
        elif (buffer[cur+2] == 226) and (buffer[cur+3] == 37):
            this_read = 2
            cur += this_read
        elif buffer[cur+3] == 226:
            this_read = 3
            cur += this_read
        else:
            this_read = 4
            cur += this_read
        remaining_bytes = remaining_bytes - this_read
    matches.append(len(buffer))
    return np.array(matches)
'''

#use regular expression to find header; fastest
import re
def find_frame_header_re(buffer, target_pattern = b'\xe2%\x0cN'):
    matches = [m.start() for m in re.finditer(target_pattern, buffer)]
    #add buffer end position to matches
    matches.append(len(buffer))
    return np.array(matches)

'''
Step 1: read each frame to an array in memory
input: 
    buffer/mmap: memory buffer or numpy memmap object
    matches: header index
return:
    an array of each frame, each frame is a predifined numpy structure

Performance:
    1. select_2064 using vectorized operation
    2. use np.empty instead of np.ones (initial memeory very costy, 7s for 2497030 packet_structs) to create a space to store each frame (4ms)
    3. real_packet still cost 4s v.s. 0.6s (C++) for 5153871936 bytes file
    4. using read_packet_opt to skip bad frames, cost 4us for 2497030 packets
'''


'''
#0. select_2064: select each frame with 2064 bytes length
#input: matches, header index
#return: idx, start index of each frame with length 2064

def select_2064(matches):
    #2064 bytes
    length = matches[1:] - matches[:-1]
    #qulified frame start location
    idx =  matches[:-1][length == 2064]
    return idx

def select_n2064(matches):
    #2064 bytes
    length = matches[1:] - matches[:-1]
    #qulified frame start location
    nidx =  matches[:-1][length != 2064]
    nlen = length[length != 2064]
    return nidx,nlen

#1. read_packet: read each frame to a numpy structure
def read_packet(buffer,idx):
    num_frame = len(idx)
    select_packets = np.empty(num_frame,dtype=st.packet_struct)
    for i in range(num_frame):
        select_packets[i] = np.frombuffer(buffer,dtype=st.packet_struct,offset= idx[i], count=1)
    return select_packets
'''

#given whole range, remove bad in between
#return start and length of each effective frame segments
#return last bad end
#calculate effective_frame_range from all matched header position

#buffer is numpy memmap file of the raw data; in chunk load mode, it is partial file (check cut file process, make sure use some redundency to avoid missing content)

#heads_idx are index of the b'\xe2%\x0cN' head, obtained from regex finding function;
def read_packet_opt(buffer,heads_idx):
    #!missing the last packet length, which should be buffer length - mathch[-1]
    #assuming there are N heads; len(heads_idx) = N;
    length = heads_idx[1:] - heads_idx[:-1]
    #for packet 1 to N - 1; get the indices of those length not equal to 2064
    bad_idx =  heads_idx[:-1][length != 2064]
    #for packet 1 to N - 1; get the actual length of those not equal to 2064 
    bad_len = length[length != 2064]
    
    # * if no bad packets, just load all data from buffer, start from the first heads and load N-1 packets
    if len(bad_idx) == 0:
        #total num readed is len match -1
        packets = [ 
            np.frombuffer(buffer, dtype=st.packet_struct, offset=heads_idx[0], count=len(heads_idx)-1)
        ]
        #? what are there two variables for, last and residual
        #*keep last header position
        last = heads_idx[-1]
        #*bytes not readed, for all good group, residual is 0
        residual = 0
    else:
    # *  if bad packets exits; exclude those bad and assemble those consecutively effive packet into eff_s, new len
    # * use np.frombuffer in group;
        eff_start = []
        eff_len = []
        #* make sure the first is included
        #! treat f as separated file; idx must start from 0, 
        #?so variable "last"  is meanningless?
        if bad_idx[0] != 0:
            #fixed by gemini and glm 5.1
            #eff_start.append(0)
            #eff_len.append(bad_idx[0])
            eff_start.append(heads_idx[0])
            eff_len.append(bad_idx[0] - heads_idx[0])
        for n in range(len(bad_idx)-1):
            #* look up to the next good start; if the first bad is not zero is not bad,
            #* if consecutive bad, will return a eff_len = 0
            eff_start.append(bad_idx[n] + bad_len[n])
            eff_len.append(bad_idx[n+1] - bad_idx[n] - bad_len[n]) 
        nonzero_idx = np.array(eff_len) != 0
        eff_start = np.array(eff_start)[nonzero_idx]
        eff_len = np.array(eff_len)[nonzero_idx]
        #if len(bad_idx == 1):
        #    #if only one found, can not be our bad, since we chose integer multiple of 2064
        #
        #    last = bad_idx[0] + bad_len[0]
        #else: 
        #    #last bad may due to out cut, so keep those bad
        #    last = bad_idx[-1] 
        #! heads_idx[-1] is the idx of the EOF
        #* if last bad start + bad length equal to the length of this file fragment; then the bad may be caused by cutting file, return last as next start point to read. 
        if bad_idx[-1] + bad_len[-1] == heads_idx[-1]:
            last = bad_idx[-1]
            residual = bad_len[-1] 
        else:
            last = bad_idx[-1] + bad_len[-1]
            residual = 0 
        #all >0
        #print(eff_len)
        #for first, if first frame is not bad, then start from 0
        #print('eff_len',eff_len)
        #for the last, if last frame start + len not equal EOF, then end with the last frame
        #print("last good is ",bad_idx[-1] + bad_len[-1])
        #print('wrong end',heads_idx[-1] + 2064)
        #if bad_idx[-1] + bad_len[-1] != heads_idx[-1] + 2064:
        #if bad_idx[-1] + bad_len[-1] != heads_idx[-1]:
        #    eff_start.append(bad_idx[-1] + bad_len[-1])
        #    eff_len.append(idx[-1] + 2064 - (bad_idx[-1] + bad_len[-1]))
        #last read from buffer
        #packets = np.concatenate(
        packets = [
                np.frombuffer(buffer,dtype=st.packet_struct,offset=eff_s, count=eff_l) 
                for eff_s, eff_l in 
                zip(
                    np.array(eff_start),
                    #how many packets in one group
                    (np.array(eff_len)/2064).astype(int)
                    )
            ]
    
    #flattened_packets = [item for sublist in packets for item in sublist]
    #return np.array(flattened_packets,dtype=st.packet_struct)
    #check packets idx
    return packets, last, residual 
    #return last, residual 

def eff_grp2(heads_idx):
    #!missing the last packet length, which should be buffer length - mathch[-1]
    #assuming there are N heads; len(heads_idx) = N;
    length = heads_idx[1:] - heads_idx[:-1]
    #for packet 1 to N - 1; get the indices of those length not equal to 2064
    bad_idx =  heads_idx[:-1][length != 2064]
    #for packet 1 to N - 1; get the actual length of those not equal to 2064 
    bad_len = length[length != 2064]
    
    # * if no bad packets, just load all data from buffer, start from the first heads and load N-1 packets
    if len(bad_idx) == 0:
        #total num readed is len match -1
        #packets = [ 
        #    np.frombuffer(buffer, dtype=st.packet_struct, offset=heads_idx[0], count=len(heads_idx)-1)
        #]
        #? what are there two variables for, last and residual
        #
        eff_start  =  0
        eff_len  =  len(heads_idx)
        last = heads_idx[-1]
        residual = 0
    else:
    # *  if bad packets exits; exclude those bad and assemble those consecutively effive packet into eff_s, new len
    # * use np.frombuffer in group;
        eff_start = []
        eff_len = []
        if bad_idx[0] != 0:
            eff_start.append(0)
            eff_len.append(bad_idx[0])
        for n in range(len(bad_idx)-1):
            eff_start.append(bad_idx[n] + bad_len[n])
            eff_len.append(bad_idx[n+1] - bad_idx[n] - bad_len[n]) 
        #if len(bad_idx == 1):
        #    #if only one found, can not be our bad, since we chose integer multiple of 2064
        #
        #    last = bad_idx[0] + bad_len[0]
        #else: 
        #    #last bad may due to out cut, so keep those bad
        #    last = bad_idx[-1] 
        if bad_idx[-1] + bad_len[-1] == heads_idx[-1]:
        #the bad is caused by cutting
            last = bad_idx[-1] 
            residual = bad_len[-1] 
        else: 
            last = bad_idx[-1] + bad_len[-1]
            residual = 0 
    return eff_start,eff_len

#!this function is too slow
def eff_idx(heads_idx):
    #good groups
    #head index of the first good one
    offsets = []
    #number of packets in each group 
    num = []
    #first one
    offsets.append(heads_idx[0])
    count = 0
    #from the second packet header
    i = 1
    while i < len(heads_idx):
        if (heads_idx[i] - heads_idx[i-1]) == 2064:
            count = count + 1
        else:
            #travel bad ones
            while (heads_idx[i+1]  - heads_idx[i]) !=2064 :
                i = i + 1
            offsets.append(heads_idx[i])
        i = i + 1
    return offsets,num


'''
Step 3: convert each frame to a human-readable measured physical data
contains a lot of individual functions to convert each field
'''
#bits_array is an array of uint8
#treat the 7th bit as sign bit not 8th
def signed_7bits(bits_array):
    bits_array = np.ascontiguousarray(bits_array)
    # 直接从字节数据创建有符号整数数组
    result = np.frombuffer(bits_array, dtype=np.int8)
    # 获取符号位（第7位）
    sign_bit = (result & 0x40) != 0
    # 清除符号位（第7位）
    result &= 0x3F
    # 应用符号
    result[sign_bit] *= -1
    return result

#in correct order
#[dim, num bytes should match the length of the dtype]
#4 bytes = u/int32 ...
def view_as(bytes_array,dtype):
    return np.ascontiguousarray(bytes_array).view(dtype).reshape(-1, 1)[:,0]

#------------------------------------------------------ new code ------------------------------------------------------
#note c++ cost only 0.6s something is wrong here
#using np u/int dtype directly, avoid int from bytes (not support vectorized naturally)
#packets alread np array, can do vec operation  
#p is packet array

def packets_2_pdf(p):
    #from keys almost no difference
    #pdf_dict = dict.fromkeys(keys)
    pdf = pd.DataFrame({col: pd.Series(dtype=dtype) for col, dtype in st.pdf_type.items()})

    #first_2bits
    #loc no need use start and end explicitly
    pdf['BM'] = (p[:]['seq1'] & 0b11).astype(np.uint8)
    # last_14bits 
    pdf['BX'] = (p[:]['seq1'] & 0b11111111111111).astype(np.uint16)
    #packet time
    #time_s is the seconds to 2021-01-01 00:00:00 BJT
    #2021-01-01 00:00:00 BJT is 2020-12-31 16:00:00 UTC
    #2020-12-31 16:00:00 UTC is 1609430400
    #so in UTC time zone, the base time is 1609430400
    #should check timestamp is an uuid or not
    #ms
    #pdf['Timestamp'] = (1609430400 + p[:]['time_s'] + p[:]['time_ms']*0.001).astype(np.float64)
    pdf['Timestamp'] = (1609430400 + p[:]['time_s'] + p[:]['time_ms']*0.001)

    #engineering parameter 1: voltage; temperature; Evt Number
    
    const1 = 50000/32768
    pdf['HV_N'] = p[:]['ep1']['HV_N'] * const1
    pdf['HV_P'] = p[:]['ep1']['HV_P'] * const1
    
    const2 = 10000/32768
    pdf['HV_C'] = p[:]['ep1']['HV_C'] * const2
    pdf['HV_G'] = p[:]['ep1']['HV_G'] * const2
    pdf['HV_A'] = p[:]['ep1']['HV_A'] * const2
    pdf['HV_B'] = p[:]['ep1']['HV_B'] * const2
    
    const3 = 50/128
    monitors =  p[:]['ep1']['monitors']
    pdf['A_N12'] = monitors[:,0] * const3
    
    const4 = 20/128
    pdf['A_N5'] = monitors[:,1] * const4

    const5 = 10/128
    pdf['A_3_3'] = monitors[:,2] * const5
    pdf['D_3_3'] = monitors[:,3] * const5

    pdf['HV_C_I'] = monitors[:,4] * const5
    pdf['HV_N_I'] = monitors[:,5] * const5
    
    # [MOD] guard np.log against invalid temps: temps==0 gives log(0)=-inf, temps>=64 makes log arg negative
    # Invalid readings produce NaN instead of RuntimeWarning
    temps_raw = p[:]['ep1']['temps'] * 0.078125
    log_arg = (temps_raw * 10000) / (5000 * (4.97 - temps_raw))
    log_arg = np.where(log_arg > 0, log_arg, np.nan)
    Temps = ((4100*298.15/(298.15*np.log(log_arg)+4100))-273.15)
    pdf['Temp_A'] = Temps[:,0].astype(np.float32)
    pdf['Temp_B'] = Temps[:,1].astype(np.float32)
    pdf['Temp_PCB'] = Temps[:,2].astype(np.float32)
    evt_num_tot = p[:]['ep1']['Evt_Num_Tot'].byteswap()
    evt_num_eff = p[:]['ep1']['Evt_Num_Eff'].byteswap()
    evt_num_aux = p[:]['ep1']['Evt_Num_Aux'].byteswap()
    pdf['Evt_Num_Tot'] = evt_num_tot.view(evt_num_tot.dtype.newbyteorder())
    pdf['Evt_Num_Eff'] = evt_num_eff.view(evt_num_eff.dtype.newbyteorder())
    pdf['Evt_Num_Aux'] = evt_num_aux.view(evt_num_aux.dtype.newbyteorder())
    
    #engineering parameter 2: attitude and position, reused frame
    #0: paras 1, q4, coords
    #1: paras 2, pointing 
    reuse_mask =  p[:]['ep2']['stat']
    
    # 复用帧第一组参数
    
    #q1, q2, q3, unit: q_dl = 2**-30
    q_dl = 2**-30
    pdf['MJ2000_q1'] = (view_as(p[:]['ep2']['flex'][:,0:4][:,::-1], np.int32) * q_dl  * (reuse_mask)).astype(np.float32)
    pdf['MJ2000_q2'] = (view_as(p[:]['ep2']['flex'][:,4:8][:,::-1], np.int32) * q_dl  * (reuse_mask)).astype(np.float32)
    pdf['MJ2000_q3'] = (view_as(p[:]['ep2']['flex'][:,8:12][:,::-1], np.int32) * q_dl  * (reuse_mask)).astype(np.float32)

    # 月心J2000惯性系卫星X\Y\Z方向位置: 补码，定点数，当量1m
    pdf['MJ2000_Star_X'] = view_as(p[:]['ep2']['flex'][:,12:16][:,::-1], np.int32) * (reuse_mask) 
    pdf['MJ2000_Star_Y'] = view_as(p[:]['ep2']['flex'][:,16:20][:,::-1], np.int32) * (reuse_mask)  
    pdf['MJ2000_Star_Z'] = view_as(p[:]['ep2']['flex'][:,20:24][:,::-1], np.int32) * (reuse_mask)   


    # 复用帧第二组参数
    #unit 0.01            
    #wrong: check
    #UTC_time = np.ascontiguousarray(p[:]['ep2']['flex'][:,0:6][:,::-1]).view(np.uint64).reshape(-1, 1)*0.01
    #arr = p[:]['ep2']['flex'][:,0:6]
    #result = (arr[:, 0].astype(np.uint64) * 281474976710656 +  # 256^5
    #          arr[:, 1].astype(np.uint64) * 1099511627776 +    # 256^4
    #          arr[:, 2].astype(np.uint64) * 4294967296 +       # 256^3
    #          arr[:, 3].astype(np.uint64) * 16777216 +         # 256^2
    #          arr[:, 4].astype(np.uint64) * 256 +              # 256^1
    #          arr[:, 5].astype(np.uint64))                     # 256^0
    #UTC_time = result * 0.01  * reuse_mask

    # 天线指向中心角:补码，定点数，当量0.00001度，范围0~+360度
    a_dl = 0.00001
    pdf['Ant_Ang'] = view_as(p[:]['ep2']['flex'][:,6:10][:,::-1], np.int32) * a_dl  *(1- reuse_mask)
    pdf['Ant_Azi'] = view_as(p[:]['ep2']['flex'][:,10:14][:,::-1], np.int32) * a_dl *(1- reuse_mask)
                    
    m_dl = 0.00001
    pdf['Star_L'] = view_as(p[:]['ep2']['flex'][:,14:18][:,::-1], np.int32) * m_dl  *(1- reuse_mask)
    pdf['Star_NZ'] = view_as(p[:]['ep2']['flex'][:,18:22][:,::-1], np.int32) * a_dl *(1- reuse_mask)
    pdf['Shadow'] = p[:]['ep2']['flex'][:,22] * (1-reuse_mask)
    
    #read first, then lookup table, avoid checking
    para_addr = p[:]['ep2']['args_addr']
    #four bytes
    paras = p[:]['ep2']['args_value']

    #select concerned parameters
    # [MOD] use st.pdf_type[name] instead of hardcoded np.uint16; uint16 causes FutureWarning when assigning int16 values from signed_7bits
    for name in ['Th_A_MCP','Th_A_X','Th_A_Y','Th_B_MCP','Th_B_X1','Th_B_X2','Th_B_Y1','Th_B_Y2','Th_B_G']:
        pdf[name] = np.zeros(len(p), dtype=st.pdf_type[name])
        
    pdf.loc[para_addr == 53,'Th_A_MCP'] = (view_as(paras[para_addr == 53,::-1],np.uint32) & 0b111111111111).astype(np.int16)
    pdf.loc[para_addr == 57,'Th_B_MCP'] = (view_as(paras[para_addr == 57,::-1],np.uint32) & 0b111111111111).astype(np.int16)
    
    pdf.loc[para_addr == 55,'Th_A_X']   = (view_as(paras[para_addr == 55,::-1],np.uint32) & 0b111111111111).astype(np.int16)
    pdf.loc[para_addr == 59,'Th_A_Y']   = (view_as(paras[para_addr == 59,::-1],np.uint32) & 0b111111111111).astype(np.int16)
    
    pdf.loc[para_addr == 142,'Th_B_X1'] = (signed_7bits(paras[para_addr == 142,1])).astype(np.int16)
    pdf.loc[para_addr == 142,'Th_B_X2'] = (signed_7bits(paras[para_addr == 142,3])).astype(np.int16)
    pdf.loc[para_addr == 146,'Th_B_Y1'] = (signed_7bits(paras[para_addr == 146,1])).astype(np.int16)
    pdf.loc[para_addr == 146,'Th_B_Y2'] = (signed_7bits(paras[para_addr == 146,3])).astype(np.int16)
    pdf.loc[para_addr == 150,'Th_B_G']  = view_as(paras[para_addr == 150][:,:1], np.int8)
    
    return pdf

# [MOD] validate_packets: detect structurally invalid packets that passed the 2064-byte spacing check
# but have wrong internal structure (e.g. false header from magic bytes appearing in event data)
def validate_packets(p, max_bx_jump=100):
    """
    Validate parsed packets and return boolean mask of good packets.

    Checks:
    1. Header bytes must match magic pattern (catches misalignment)
    2. BX sequence should be roughly monotonic (catches false headers)
    3. Timestamps should not decrease (catches garbage data)

    Args:
        p: numpy structured array of packets (from read_packet_opt)
        max_bx_jump: maximum allowed BX jump between consecutive packets
                     (default 100, for normal operation BX increments by 1)
    Returns:
        good_mask: boolean array, True for valid packets
    """
    n = len(p)
    if n == 0:
        return np.array([], dtype=bool)

    good = np.ones(n, dtype=bool)

    # Check 1: header bytes must be \xe2\x25\x0c\x4e
    header_valid = np.all(p[:]['header'] == [0xE2, 0x25, 0x0C, 0x4E], axis=1)
    good &= header_valid

    # Check 2: BX (14-bit counter) should be roughly monotonic
    bx = (p[:]['seq1'] & 0x3FFF).astype(np.int32)
    if n > 1:
        diff = np.diff(bx)
        # Allow wrap: if diff < -8000, it's likely a 14-bit wrap (16384)
        # Any other negative diff, or a huge positive jump, indicates a bad packet
        bad_jump = (diff < -100) | (diff > max_bx_jump)
        # Mark both neighbors of a bad jump as suspicious
        bad_seq = np.zeros(n, dtype=bool)
        bad_seq[:-1] |= bad_jump
        bad_seq[1:] |= bad_jump
        good &= ~bad_seq

    # Check 3: timestamp should not decrease
    if n > 1:
        time_diff = np.diff(p[:]['time_s'].astype(np.int64))
        bad_time = np.zeros(n, dtype=bool)
        bad_time[1:] |= (time_diff < 0)
        good &= ~bad_time

    return good

import uuid
#c++ parse packets cost: 4.09255 s
def packets_2_df0(p):
    df0 = pd.DataFrame({col: pd.Series(dtype=dtype) for col, dtype in st.df0_type.items()})
    num_packets = len(p)
    total_rows = num_packets * 62
    uuids = [str(uuid.uuid4()) for _ in range(num_packets)]
    df0['UUID'] = np.repeat(uuids, 62)
    
    df0['PID'] = np.tile(np.arange(1, 63, dtype=np.uint8), num_packets)

    df0['BM'] = np.repeat(p[:]['seq1'] & 0b11, 62)
    df0['BX'] = np.repeat(p[:]['seq1'] & 0b11111111111111, 62)
    df0['Timestamp'] = (np.repeat(1609430400 + p[:]['time_s'] + p[:]['time_ms']*0.001, 62)).astype(np.float64)
    
    p_dict_list = p[:]['events'][:] 
    # 1、探测器标志，取标志位，高位第1位
    df0['prob_flag'] = (p_dict_list['flag'] >> 7).flatten()
    # 1、事件有效性，取标志位，高位第2位
    df0['EvtLevel'] = ((p_dict_list['flag'] >> 6) & 0b1).flatten()
    # 1、探测器探测模式标记，取标志位，高位第3、4位
    #event_D_incid_m = (p_dict_list['flag'] >> 4) & 0b11
    # 1、探测器探测子模式标记，低4位
    #event_D_incid_z_m = (p_dict_list['flag'] >> 0) & 0b1111
    
    # 2、秒计数时间码，取时间码前2字节
    event_Times = p_dict_list['time_s']
    # 2、微秒计数时间码，取时间码后2字节？
    event_Timeus = p_dict_list['time_us']
    df0['TimeInternal'] = (event_Times + event_Timeus*0.00002).flatten().astype(np.float32)
    
    # 3、电子幅度，取低12位; similar to UTC which is 6 bytes
    df0['EP'] = ( ( (p_dict_list['se_amp'][..., 1].astype(np.uint16) << 8 ) | (p_dict_list['se_amp'][..., 2].astype(np.uint16) ) ) & 0b111111111111).flatten()
    
    #ch_fly_time_head =  p_dict_list['c0'] >> 13
    ch_fly_time_tail =  p_dict_list['c0'] & 0b1111111111111
    
    for i,name in enumerate(['ch_tof_0', 'ch_tof_1', 'ch_tof_2', 'ch_tof_3']):
        df0[name] = ch_fly_time_tail[...,i].flatten()
    
    
    #ch01_head =  p_dict_list['c1'] >> 13
    ch01_tail =  p_dict_list['c1'] & 0b1111111111111
    
    #ch03_head =  p_dict_list['c3'] >> 13
    ch03_tail =  p_dict_list['c3'] & 0b1111111111111
    
    #ch05_head =  p_dict_list['c5'] >> 13
    ch05_tail =  p_dict_list['c5'] & 0b1111111111111
    
    #ch07_head =  p_dict_list['c7'] >> 13
    ch07_tail =  p_dict_list['c7'] & 0b1111111111111
    
    df0['ch_1s'] = (ch01_tail[...,0].flatten()).astype(np.uint16)
    df0['ch_1e'] = (ch01_tail[...,1].flatten()).astype(np.uint16)
    df0['ch_3s'] = (ch03_tail[...,0].flatten()).astype(np.uint16)
    df0['ch_3e'] = (ch03_tail[...,1].flatten()).astype(np.uint16)
    df0['ch_5s'] = (ch05_tail[...,0].flatten()).astype(np.uint16)
    df0['ch_5e'] = (ch05_tail[...,1].flatten()).astype(np.uint16)
    df0['ch_7s'] = (ch07_tail[...,0].flatten()).astype(np.uint16)
    df0['ch_7e'] = (ch07_tail[...,1].flatten()).astype(np.uint16)
    
    return df0


if __name__ == "__main__":
    print('test')
