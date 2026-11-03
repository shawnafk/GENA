#!/usr/bin/env python3
"""Simple sequential CLI for GENA raw -> Parquet conversion."""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

if __package__ in (None, ""):
    # Allow direct execution as ./raw_2_pqt.py from gena/raw: make the repository
    # root importable before importing the package.
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    from common.config import get_configured_paths
    from gena.raw.decode import RawConfig, convert_raw
else:
    from common.config import get_configured_paths
    from .decode import RawConfig, convert_raw


def _products(value):
    return frozenset(item.strip() for item in value.split(',') if item.strip())


def main(argv=None):
    parser = argparse.ArgumentParser(
        description='Convert GENA raw files to Parquet (sequential, overwrite-based)')
    parser.add_argument('-o', '--outdir', type=str,
                        help='output directory; defaults to configured pqt_base')
    parser.add_argument('files', nargs='+', help='raw input files')
    parser.add_argument('--snid', type=int, choices=(1, 2), default=1)
    parser.add_argument('--chunk-packets', type=int, default=RawConfig.chunk_packets)
    parser.add_argument('--products', default=','.join(sorted(RawConfig().products)),
                        help='comma-separated: pdf,l1,l1a,l2,l2l,l2m,l2h,l2f')
    parser.add_argument('--no-day-split', action='store_true')
    parser.add_argument('--no-overwrite', action='store_true')
    parser.add_argument('--force', action='store_true')
    parser.add_argument('-v', '--verbose', action='store_true')
    parser.add_argument('-q', '--quiet', action='store_true',
                        help='disable the progress bar and per-file start message')
    args = parser.parse_args(argv)
    logging.basicConfig(
        level=logging.INFO if args.verbose else logging.WARNING,
        format='%(levelname)s %(name)s: %(message)s')
    output_dir = Path(args.outdir) if args.outdir else Path(get_configured_paths()['pqt_base'])
    config = RawConfig(
        snid=args.snid,
        chunk_packets=args.chunk_packets,
        by_day=not args.no_day_split,
        overwrite=not args.no_overwrite,
        force=args.force,
        products=_products(args.products),
        progress=not args.quiet,
    )
    from tqdm import tqdm
    failed = 0
    for raw_file in args.files:
        progress = None
        if not args.quiet:
            print(f'Processing {raw_file} ...', flush=True)
            progress = tqdm(total=None, unit='packet', desc=Path(raw_file).name,
                            leave=False)
        try:
            report = convert_raw(raw_file, output_dir, config=config,
                                 progress=progress)
            print(f'{raw_file}: {report.packet_count} packets, '
                  f'{sum(report.product_rows.values())} rows, '
                  f'{report.elapsed_seconds:.2f}s', flush=True)
        except Exception as exc:
            failed += 1
            print(f'ERROR {raw_file}: {exc}', flush=True)
        finally:
            if progress is not None:
                progress.close()
    return 1 if failed else 0


if __name__ == '__main__':
    raise SystemExit(main())
