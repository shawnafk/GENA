"""Raw packet conversion orchestration.

The implementation is intentionally kept in this module: framing and field
decoding live in :mod:`gena.raw.packet`, level physics lives in
:mod:`gena.raw.conversion`, and this module only coordinates passes, product
selection and I/O.
"""

from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import dataclass
import hashlib
import logging
from pathlib import Path
import time

import numpy as np
import pandas as pd
from tqdm import tqdm

from . import packet as pkt
from . import conversion as conv
from . import constants as const
from . import pqt_io
import gena.filter as flt


LOGGER = logging.getLogger(__name__)


def _normalize_probe(prob):
    """Return lowercase probe identity; decode owns public probe validation."""
    if not isinstance(prob, str) or prob.lower() not in ('a', 'b'):
        raise ValueError("prob must be 'a' or 'b'")
    return prob.lower()

ENGINEERING_COLUMNS = [
    'PacketID', 'Timestamp', 'HV_C', 'HV_G',
    'MJ2000_q1', 'MJ2000_q2', 'MJ2000_q3',
    'MJ2000_Star_X', 'MJ2000_Star_Y', 'MJ2000_Star_Z',
    'PosePresent',
]
EVENT_PRODUCTS = frozenset({"l1", "l1a", "l2", "l2s", "l2g"})
ALL_PRODUCTS = frozenset({"pdf"}) | EVENT_PRODUCTS

# L2G selection: apply the physical 90 x 30 mm carbon/start aperture. The end
# plane keeps the existing 150 x 60 mm engineering window.
L2G_START_APERTURE_MM = (90, 30)
L2G_END_WINDOW_MM = (150, 60)


@dataclass(frozen=True)
class RawConfig:
    """Explicit configuration for one raw conversion."""

    snid: int = 1
    chunk_packets: int = 4096
    by_day: bool = True
    overwrite: bool = True
    force: bool = False
    products: frozenset[str] = frozenset({"pdf", "l1a", "l2", "l2s", "l2g"})
    progress: bool = True

    def validate(self):
        if self.snid not in (1, 2):
            raise ValueError("snid must be 1 (orbit) or 2 (ground)")
        if not isinstance(self.chunk_packets, int) or self.chunk_packets < 1:
            raise ValueError("chunk_packets must be a positive integer")
        unknown = set(self.products) - ALL_PRODUCTS
        if unknown:
            raise ValueError(f"unknown products: {sorted(unknown)}")
        if not self.products:
            raise ValueError("at least one product must be requested")


@dataclass
class RunReport:
    source: Path
    source_bytes: int
    source_sha256: str
    snid: int
    packet_count: int
    duplicate_packets: int
    truncated_bytes: int
    product_rows: dict[str, int]
    product_paths: dict[str, Path]
    elapsed_seconds: float


def _sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def _iter_unique_packets(raw_file, config, seen, progress=None, stats=None):
    """Yield (packet array, PacketID array, duplicate count) in file order."""
    next_id = 0
    for packets in pkt.iter_packets(raw_file, config.chunk_packets, stats=stats):
        for packet in packets:
            uuids = pkt.packet_uuids(packet)
            keep_list = []
            for uuid_value in uuids:
                accepted = uuid_value not in seen
                keep_list.append(accepted)
                if accepted:
                    seen.add(uuid_value)
            keep = np.asarray(keep_list, dtype=bool)
            count = int(keep.sum())
            if progress is not None:
                progress.update(len(uuids))
            if count == 0:
                yield packet[keep], np.empty(0, dtype=np.uint32), len(uuids)
                continue
            ids = np.arange(next_id, next_id + count, dtype=np.uint32)
            next_id += count
            yield packet[keep], ids, len(uuids) - count


def _engineering_frame(pdf):
    columns = [name for name in ENGINEERING_COLUMNS if name in pdf.columns]
    return pdf[columns].copy()


def _scan_engineering(raw_file, config, writer, progress=None):
    """First pass: validate/deduplicate packets, write PDF, collect index."""
    seen = set()
    parts = []
    packet_count = 0
    duplicate_packets = 0
    stats = {}
    for packet, ids, duplicates in _iter_unique_packets(raw_file, config, seen, progress, stats=stats):
        packet_count += len(ids)
        duplicate_packets += duplicates
        if len(ids) == 0:
            continue
        pdf = pkt.packets_2_pdf(packet)
        pdf['PacketID'] = ids
        if 'pdf' in config.products:
            writer.write(pdf, 'pdf')
        parts.append(_engineering_frame(pdf))
    engineering = (pd.concat(parts, ignore_index=True)
                   if parts else pd.DataFrame(columns=ENGINEERING_COLUMNS))
    return engineering, packet_count, duplicate_packets, int(stats.get('truncated_bytes', 0))


def _quality_masks(l1a, l2, params, probe):
    soft = flt.level2_soft(l1a, params=params, prob=probe)
    gold = flt.level2_tight(
        l1a, l2, params=params, prob=probe,
        c_size=L2G_START_APERTURE_MM,
        m_size=L2G_END_WINDOW_MM,
    )
    return np.asarray(soft, dtype=bool), np.asarray(gold, dtype=bool)


def _convert_events(raw_file, config, engineering, writer, progress=None):
    """Second pass: decode events and write L1/L1A/L2 products."""
    params = const.init_param(config.snid)
    seen = set()
    counts = Counter()
    packet_count = 0
    duplicate_packets = 0
    want_l1 = 'l1' in config.products
    want_l1a = 'l1a' in config.products or bool(config.products & {'l2', 'l2s', 'l2g'})
    want_l2 = bool(config.products & {'l2', 'l2s', 'l2g'})

    for packet, ids, duplicates in _iter_unique_packets(raw_file, config, seen, progress):
        packet_count += len(ids)
        duplicate_packets += duplicates
        if len(ids) == 0:
            continue
        df0 = pkt.packets_2_df0(packet)
        repeated = np.repeat(ids, 62)
        pid = np.tile(np.arange(1, 63, dtype=np.uint64), len(ids))
        df0['PacketID'] = repeated
        df0['EventID'] = (repeated.astype(np.uint64) << 6) | (pid - 1)
        # L0 expands all 62 fixed slots. ``is_empty_slot`` is true only when
        # the whole 32-byte event slot is zero; such padding has flag=0, so
        # without this drop l0_2_l1() would classify it as a probe-A event and
        # write an all-zero L1 row. Drop structural padding before probe split.
        df0 = df0.loc[~df0['is_empty_slot']].reset_index(drop=True)
        if df0.empty:
            continue
        for probe in ('a', 'b'):
            probe = _normalize_probe(probe)
            l1 = conv.l0_2_l1(df0, probe, config.snid).reset_index(drop=True)
            if want_l1:
                writer.write(l1, f'l1_{probe}')
                counts[f'l1_{probe}'] += len(l1)
            if not want_l1a:
                continue
            l1a = conv.l1_2_l1a(l1, engineering)
            if len(l1a) == 0:
                continue
            if 'l1a' in config.products:
                writer.write(l1a, f'l1a_{probe}')
                counts[f'l1a_{probe}'] += len(l1a)
            if not want_l2:
                continue
            l2 = conv.l1a_2_l2(l1a, params=params, prob=probe)
            if 'l2' in config.products:
                writer.write(l2, f'l2_{probe}')
                counts[f'l2_{probe}'] += len(l2)
            soft, gold = _quality_masks(l1a, l2, params, probe)
            if 'l2s' in config.products:
                # l2.loc[soft] is already the selected subset. Mutate that small
                # frame in place; never copy the full L2 table here.
                l2s = l2.loc[soft]
                l2s = conv.l2s_with_xy_tof(l2s)
                writer.write(l2s, f'l2s_{probe}')
                counts[f'l2s_{probe}'] += len(l2s)
            if 'l2g' in config.products:
                writer.write(l2.loc[gold], f'l2g_{probe}')
                counts[f'l2g_{probe}'] += int(np.count_nonzero(gold))
    return packet_count, duplicate_packets, counts


def convert_raw(raw_file, output_folder, config=None, progress=None):
    """Convert one raw file into the configured product set.

    The function runs a packet-level engineering pass, builds a global
    :class:`conversion.EngineeringIndex`, then runs the event pass. Existing
    products are replaced atomically; no append or partial Parquet file is
    published.
    """
    config = config or RawConfig()
    config.validate()
    source = Path(raw_file).expanduser().resolve()
    target = Path(output_folder).expanduser().resolve()
    if not source.is_file():
        raise FileNotFoundError(source)
    target.mkdir(parents=True, exist_ok=True)
    LOGGER.info('converting %s -> %s (snid=%d, products=%s)',
                source, target, config.snid, ','.join(sorted(config.products)))
    before = source.stat()
    source_sha256 = _sha256_file(source)
    start = time.perf_counter()
    writer = pqt_io.ProductWriter(
        target,
        by_day=config.by_day,
        overwrite=config.overwrite,
        force=config.force,
        metadata={
            'gena_source_sha256': source_sha256,
            'gena_snid': config.snid,
            'gena_schema_version': pqt_io.SCHEMA_VERSION,
        },
    )
    try:
        engineering, packet_count, duplicates, truncated_bytes = _scan_engineering(
            source, config, writer, progress)
        LOGGER.info('engineering pass: %d unique packets, %d duplicate packets, %d truncated bytes',
                    packet_count, duplicates, truncated_bytes)
        if duplicates:
            LOGGER.warning('%d duplicate packet(s) skipped', duplicates)
        if truncated_bytes:
            LOGGER.warning('%d truncated byte(s) at end of source', truncated_bytes)
        if config.products & EVENT_PRODUCTS:
            event_engineering = conv.EngineeringIndex(engineering)
            event_packets, event_duplicates, _ = _convert_events(
                source, config, event_engineering, writer, progress)
            if event_packets != packet_count or event_duplicates != duplicates:
                raise RuntimeError(
                    'source changed between engineering and event passes')
        after = source.stat()
        if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
            raise RuntimeError('source changed during conversion')
        rows, paths = writer.finish()
        LOGGER.info('wrote %d product file(s), %d total rows',
                    len(rows), sum(rows.values()))
    except BaseException:
        writer.abort()
        raise
    elapsed = time.perf_counter() - start
    return RunReport(
        source=source,
        source_bytes=before.st_size,
        source_sha256=source_sha256,
        snid=config.snid,
        packet_count=packet_count,
        duplicate_packets=duplicates,
        truncated_bytes=truncated_bytes,
        product_rows=rows,
        product_paths=paths,
        elapsed_seconds=elapsed,
    )


def _parse_products(value):
    if not value:
        return RawConfig().products
    return frozenset(item.strip() for item in value.split(',') if item.strip())


def main(argv=None):
    parser = argparse.ArgumentParser(description='Convert GENA raw telemetry to Parquet.')
    parser.add_argument('raw_file')
    parser.add_argument('output_folder')
    parser.add_argument('--snid', type=int, choices=(1, 2), default=1)
    parser.add_argument('--chunk-packets', type=int, default=4096)
    parser.add_argument('--products', default=','.join(sorted(RawConfig().products)),
                        help='comma-separated: pdf,l1,l1a,l2,l2s,l2g')
    parser.add_argument('--no-day-split', action='store_true')
    parser.add_argument('--no-overwrite', action='store_true')
    parser.add_argument('--force', action='store_true')
    parser.add_argument('--quiet', action='store_true')
    args = parser.parse_args(argv)
    logging.basicConfig(
        level=logging.ERROR if args.quiet else logging.INFO,
        format='%(levelname)s %(name)s: %(message)s')
    config = RawConfig(
        snid=args.snid,
        chunk_packets=args.chunk_packets,
        by_day=not args.no_day_split,
        overwrite=not args.no_overwrite,
        force=args.force,
        products=_parse_products(args.products),
        progress=not args.quiet,
    )
    progress = None
    if config.progress:
        progress = tqdm(total=None, unit='packet', desc='raw', leave=False)
    try:
        report = convert_raw(args.raw_file, args.output_folder, config=config,
                             progress=progress)
    except Exception:
        LOGGER.exception('conversion failed: %s', args.raw_file)
        return 1
    finally:
        if progress is not None:
            progress.close()
    print(f'source: {report.source}')
    print(f'packets: {report.packet_count} (duplicates {report.duplicate_packets})')
    for name in sorted(report.product_rows):
        print(f'{name}: {report.product_rows[name]}')
    print(f'elapsed: {report.elapsed_seconds:.2f} s')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
