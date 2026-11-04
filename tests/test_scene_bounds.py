"""Pose-coverage boundaries for the vectorized 3-D scene."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from scipy.spatial.transform import Rotation

from plots.v3d import SceneContext

# 2025-11-13T02:23:06.046 has no exact float64 Unix-second representation, so a
# datetime64[ns] round-trip lands one float64 unit (about 0.24 us) past it.
BOUNDARY = 1763000586.046


def _context(timestamps) -> SceneContext:
    timestamps = np.asarray(timestamps, dtype=float)
    return SceneContext(
        pose_timestamps=timestamps,
        pose_positions=np.column_stack(
            [
                np.linspace(1.0, 2.0, len(timestamps)),
                np.zeros(len(timestamps)),
                np.zeros(len(timestamps)),
            ]
        ),
        pose_rotations=Rotation.identity(len(timestamps)),
        cache_timestamps=timestamps,
        orbit_timestamps=timestamps,
        satellite_orbits={},
        space_weather={},
        space_weather_source="none",
        include_field_lines=False,
        include_neutral_sheet=False,
        include_reference_planes=False,
    )


def test_scene_time_on_last_pose_sample_stays_inside_coverage():
    context = _context([BOUNDARY - 60.0, BOUNDARY])

    _, timestamps = context._normalize_times(
        pd.to_datetime([BOUNDARY], unit="s", utc=True)
    )

    assert timestamps[0] <= context.pose_timestamps[-1]


@pytest.mark.parametrize("offset", [-61.0, 1.0])
def test_scene_time_outside_pose_coverage_is_rejected(offset):
    context = _context([BOUNDARY - 60.0, BOUNDARY])

    with pytest.raises(ValueError, match="outside loaded pose coverage"):
        context._normalize_times(
            pd.to_datetime([BOUNDARY + offset], unit="s", utc=True)
        )


def test_timerange_limits_on_pose_samples_keep_both_endpoints():
    context = _context([BOUNDARY - 60.0, BOUNDARY])

    figure = context.make_timerange([BOUNDARY - 60.0, BOUNDARY], max_frames=2)

    assert len(figure.frames) == 2
