from common.config import get_configured_paths
import argparse
from multiprocessing import Pool, cpu_count
from functools import partial

def process_single_file(task, verbose=False):
    from gena.raw import decode

    raw_file, output_dir = task
    try:
        if not verbose:
            print(f"Processing: {raw_file}")
        decode.decode(raw_file, output_dir, verbose=verbose)
    except Exception as e:
        print(f"ERROR processing {raw_file}: {str(e)}")

def main():
    parser = argparse.ArgumentParser(
        description="Convert raw files to parquet (parallel processing)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Environment variables (set before starting Python; 0 disables and 1 enables a flag):
  GENA_CHUNK_SIZE      Packets decoded per chunk. Default: 100000.
  GENA_PDF_ONLY       Skip event conversion after parsing the PDF data. Default: 0.
  GENA_SAVE_PDF       Write the PDF product. Default: 1.
  GENA_SAVE_L1        Write L1 products. Default: 0.
  GENA_SAVE_L1A       Write L1A products. Default: 1.
  GENA_USE_REPAIR     Run L1A-to-L1B event repair. Default: 0.
  GENA_REPAIR_METHOD  Repair algorithm used when GENA_USE_REPAIR=1: 1 or 2. Default: 2.
  GENA_SAVE_L1B       Write repaired L1B products when repair is enabled. Default: 1.
  GENA_SAVE_L2        Write L2/L2S/L2G products. Default: 1.

Example:
  GENA_CHUNK_SIZE=200000 GENA_SAVE_L1=1 python gena/raw/raw_2_pqt.py -o OUTPUT RAW_FILE
""",
    )
    parser.add_argument("-o", "--outdir", type=str, help="Specify output directory")
    parser.add_argument("files", nargs="+", help="Raw input files")
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose mode: no pbar, all output to terminal")
    args = parser.parse_args()

    output_dir = args.outdir if args.outdir is not None else get_configured_paths()['pqt_base']

    # Prepare tasks
    tasks = [(f, output_dir) for f in args.files]
    
    # ✅ KEY FIX: auto use min(CPU cores, file count)
    num_workers = min(cpu_count(), len(tasks))
    print(f"Starting with {num_workers} workers...")

    with Pool(num_workers) as pool:
        pool.map(partial(process_single_file, verbose=args.verbose), tasks)

    print("\nAll files processed.")

if __name__ == "__main__":
    main()
