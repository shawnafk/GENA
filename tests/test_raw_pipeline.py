import json
from collections import Counter
from pathlib import Path

import numpy as np
import pytest

from gena.raw import packet, structure
from gena.raw.pipeline import iter_packets


def _packets(count=1):
    packets = np.zeros(count, dtype=structure.packet_struct)
    packets['header'] = [0xE2, 0x25]
    packets['seq1'] = 0x0C4E
    packets['seq2'] = np.arange(count)
    packets['dataLen'] = packet.PACKET_DATA_LENGTH
    packets['time_s'] = 121513700 + np.arange(count)
    packets['time_ms'] = 197
    raw = packets.view(np.uint8).reshape(count, packet.PACKET_SIZE)
    packets['verification_sum'] = raw[:, 8:-2].sum(axis=1, dtype=np.uint64) & 0xFFFF
    return packets


def test_iter_packets_keeps_last_packet_and_split_headers(tmp_path):
    packets = _packets(3)
    payload = b'prefix' + packets[:1].tobytes() + packets[1:].tobytes() + b'tail'
    source = tmp_path / 'sample.raw'
    source.write_bytes(payload)
    stats = Counter()
    rejected = (tmp_path / 'rejected.jsonl').open('w')
    try:
        arrays = list(iter_packets(source, packet.PACKET_SIZE, stats, rejected))
    finally:
        rejected.close()
    result = np.concatenate([part for part, _ in arrays], dtype=structure.packet_struct)
    assert result.tobytes() == packets.tobytes()
    assert [offsets.tolist() for _, offsets in arrays] == [[6], [6 + packet.PACKET_SIZE],
                                                            [6 + 2 * packet.PACKET_SIZE]]
    assert stats['verified_packets'] == 3
    assert stats['source_bytes'] == len(payload)


def test_iter_packets_reports_rejected_and_truncated_candidates(tmp_path):
    packets = _packets(2)
    corrupt = bytearray(packets[0].tobytes())
    corrupt[100] ^= 1
    truncated = packets[1].tobytes()[:100]
    source = tmp_path / 'bad.raw'
    source.write_bytes(bytes(corrupt) + truncated)
    stats = Counter()
    rejected_path = tmp_path / 'rejected.jsonl'
    with rejected_path.open('w') as rejected:
        arrays = list(iter_packets(source, packet.PACKET_SIZE, stats, rejected))
    assert arrays == []
    rows = [json.loads(line) for line in rejected_path.read_text().splitlines()]
    assert any('checksum' in row['reasons'] for row in rows)
    assert any('truncated_tail' in row['reasons'] for row in rows)
    assert stats['verified_packets'] == 0


def test_iter_packets_does_not_split_on_magic_inside_valid_payload(tmp_path):
    packets = _packets(2)
    raw = packets.view(np.uint8).reshape(2, packet.PACKET_SIZE)
    raw[0, 100:104] = np.frombuffer(packet.PACKET_MAGIC, dtype=np.uint8)
    raw[0, 8:-2]  # checksum is recomputed after the payload edit
    packets['verification_sum'] = raw[:, 8:-2].sum(axis=1, dtype=np.uint64) & 0xFFFF
    source = tmp_path / 'payload_magic.raw'
    source.write_bytes(packets.tobytes())
    stats = Counter()
    with (tmp_path / 'rejected.jsonl').open('w') as rejected:
        arrays = list(iter_packets(source, packet.PACKET_SIZE * 2, stats, rejected))
    result = np.concatenate([part for part, _ in arrays], dtype=structure.packet_struct)
    assert result.tobytes() == packets.tobytes()
    assert stats['verified_packets'] == 2


def test_frozen_validation_sample_matches_expected_manifest(tmp_path):
    from gena.raw.pipeline import process_raw

    root = Path(__file__).parents[1] / 'docs/raw_audit_20261027/validation_samples_v1'
    expected = json.loads((root / 'expected.json').read_text())
    manifest = process_raw(root / 'sample.raw', tmp_path / 'package', chunk_packets=1,
                           by_day=False)
    audit = manifest['packet_audit']
    assert audit['source_sha256'] == expected['input_sha256']
    assert audit['source_bytes'] == expected['input_bytes']
    for key, value in expected['packet_audit'].items():
        assert audit[key] == value
    assert {probe: dict(values) for probe, values in manifest['funnel'].items()} == expected['funnel']
    assert {name: info['rows'] for name, info in manifest['files'].items()} == expected['products']
    assert manifest['profile']['profile_id'] == expected['profile_id']


def test_frozen_processing_refuses_to_overwrite_a_package(tmp_path):
    from gena.raw.pipeline import process_raw

    root = Path(__file__).parents[1] / 'docs/raw_audit_20261027/validation_samples_v1'
    output = tmp_path / 'package'
    process_raw(root / 'sample.raw', output, chunk_packets=2, by_day=False)
    with pytest.raises(FileExistsError):
        process_raw(root / 'sample.raw', output, chunk_packets=2, by_day=False)
