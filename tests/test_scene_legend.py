"""Legend entries for the probe FOV cones and satellite body axes."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.spatial.transform import Rotation

from plots.v3d import SceneContext

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


def test_scene_legend_names_probe_fovs_and_body_axes():
    context = _context([BOUNDARY - 60.0, BOUNDARY])

    figure = context.make_timerange([BOUNDARY - 60.0, BOUNDARY], max_frames=2)

    names = {trace.name for trace in figure.data if trace.showlegend}
    assert {'A probe', 'B probe', 'sat +X', 'sat +Y', 'sat +Z'} <= names
