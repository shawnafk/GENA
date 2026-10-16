from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import numpy as np
import pandas as pd
import pytest


SCRIPT = Path(__file__).resolve().parents[1] / "1.Paper_BBF" / "ena_cg_fov_peak" / "run_conditional_cg_duration_sweep.py"
SPEC = spec_from_file_location("run_conditional_cg_duration_sweep", SCRIPT)
MODULE = module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def test_centered_window_preserves_duration_and_center_away_from_edges():
    start = pd.Timestamp("2025-01-01T00:00:00")
    end = start + pd.Timedelta(hours=2)
    center = (start + pd.Timedelta(hours=1)).timestamp()
    window = MODULE.centered_window_bounds(
        interval_start=start,
        interval_end=end,
        center_unix=center,
        duration_seconds=1800.0,
    )

    assert np.isclose(window["end_unix"] - window["start_unix"], 1800.0)
    assert np.isclose(0.5 * (window["start_unix"] + window["end_unix"]), center)


def test_centered_window_clamps_to_interval_without_shortening():
    start = pd.Timestamp("2025-01-01T00:00:00")
    end = start + pd.Timedelta(hours=1)
    window = MODULE.centered_window_bounds(
        interval_start=start,
        interval_end=end,
        center_unix=(start + pd.Timedelta(minutes=2)).timestamp(),
        duration_seconds=1800.0,
    )

    assert window["start_unix"] == start.timestamp()
    assert np.isclose(window["end_unix"] - window["start_unix"], 1800.0)


def test_three_segments_are_gap_free_and_cover_parent_window():
    segments = MODULE.split_window_bounds(1000.0, 2800.0, n_segments=3)

    assert len(segments) == 3
    assert segments[0]["start_unix"] == 1000.0
    assert segments[-1]["end_unix"] == 2800.0
    assert all(np.isclose(segment["window_seconds"], 600.0) for segment in segments)
    assert all(segments[index]["end_unix"] == segments[index + 1]["start_unix"] for index in range(2))


def test_invalid_window_requests_are_rejected():
    with pytest.raises(ValueError):
        MODULE.split_window_bounds(2.0, 1.0)
    with pytest.raises(ValueError):
        MODULE.centered_window_bounds(
            interval_start="2025-01-01T00:00:00",
            interval_end="2025-01-01T00:10:00",
            center_unix=pd.Timestamp("2025-01-01T00:05:00").timestamp(),
            duration_seconds=900.0,
        )
