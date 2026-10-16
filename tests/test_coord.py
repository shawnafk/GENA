from datetime import datetime, timezone

import numpy as np
import pytest

from common.coord import AGSM_ANGLE_DEG, coord, transform_vectors



def test_gsm_to_agsm_aligns_nominal_tail_with_negative_x():
    angle = np.radians(AGSM_ANGLE_DEG)
    tail_gsm = np.array([-np.cos(angle), np.sin(angle), 0.0])

    tail_agsm = transform_vectors(tail_gsm, source='GSM', target='AGSM')

    assert tail_agsm == pytest.approx([-1.0, 0.0, 0.0], abs=1e-15)


@pytest.mark.parametrize(
    "vectors",
    [
        np.array([1.0, 2.0, 3.0]),
        np.array([[1.0, 2.0, 3.0], [-4.0, 5.0, -6.0]]),
    ],
)
def test_agsm_helpers_round_trip(vectors):
    assert transform_vectors(transform_vectors(vectors, source='GSM', target='AGSM'), source='AGSM', target='GSM') == pytest.approx(vectors)


def test_coord_converts_between_spacepy_frames_and_agsm():
    dates = [
        datetime(2024, 1, 1, tzinfo=timezone.utc),
        datetime(2024, 1, 2, tzinfo=timezone.utc),
    ]
    gsm = np.array([[10.0, 1.0, -2.0], [-20.0, 3.0, 4.0]])

    position = coord(dates, gsm, "GSM")
    position.convert("AGSM")
    expected_agsm = transform_vectors(gsm, source='GSM', target='AGSM')
    assert position["AGSM"] == pytest.approx(expected_agsm)

    restored = coord(dates, position["AGSM"], "AGSM")
    restored.convert("GSE")

    direct = coord(dates, gsm, "GSM")
    direct.convert("GSE")
    assert restored["GSE"] == pytest.approx(direct["GSE"])


def test_agsm_helpers_reject_invalid_shapes():
    with pytest.raises(ValueError, match=r"expected \(3,\) or \(N, 3\)"):
        transform_vectors(np.ones((3, 2)), source='GSM', target='AGSM')


def test_probe_mounting_matches_gena_spacecraft_y_signs():
    source_boresight = np.array([0.0, 0.0, -1.0])

    probe_a = transform_vectors(source_boresight, source='A', target='QQ')
    probe_b = transform_vectors(source_boresight, source='B', target='QQ')

    assert probe_a[1] < 0.0
    assert probe_b[1] > 0.0
    assert probe_a[2] == pytest.approx(probe_b[2])
