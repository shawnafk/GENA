"""Smoke tests for the precomputed-frame animation path."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from plots.ani import AnimationSpec, gena_fov_frames, render_fov


def _event_frame(count: int = 8) -> pd.DataFrame:
    times = pd.date_range("2025-12-04T00:00:00", periods=count, freq="1min")
    return pd.DataFrame(
        {
            "Datetime": times,
            "Timestamp": times.view("int64") / 1e9,
            "Energy": np.linspace(10.0, 90.0, count),
            "PX": np.linspace(-10.0, 10.0, count),
            "PY": np.linspace(-5.0, 5.0, count),
        }
    )


def test_animation_spec_validation():
    with pytest.raises(ValueError):
        AnimationSpec(fps=0)
    with pytest.raises(ValueError):
        AnimationSpec(dpi=0)
    with pytest.raises(ValueError):
        AnimationSpec(figsize=(10, 0))


def test_frames_cap_and_global_scale():
    frames = gena_fov_frames(_event_frame(10), window="1min", max_frames=4)
    assert frames is not None
    assert len(frames) == 4
    assert frames.energy_counts.shape[0] == 4
    assert frames.fov_counts.shape[0] == 4
    assert frames.cmax >= 1.0
    assert frames.nummax >= 1.0


def test_render_writes_mp4_and_poster(tmp_path: Path):
    frames = gena_fov_frames(_event_frame(5), window="1min")
    assert frames is not None
    target = tmp_path / "fov.mp4"
    progress: list[tuple[int, int]] = []
    result = render_fov(
        frames,
        target,
        spec=AnimationSpec(fps=2, dpi=40, figsize=(4, 3), poster=True),
        progress=lambda i, n: progress.append((i, n)),
    )
    assert result == target
    assert target.exists() and target.stat().st_size > 0
    assert target.with_suffix(".jpg").exists()
    assert progress and progress[-1] == (len(frames), len(frames))


def test_render_cancel_leaves_no_output(tmp_path: Path):
    frames = gena_fov_frames(_event_frame(5), window="1min")
    assert frames is not None
    target = tmp_path / "cancelled.mp4"
    result = render_fov(
        frames,
        target,
        spec=AnimationSpec(fps=2, dpi=40, figsize=(4, 3), poster=False),
        cancel=lambda: True,
    )
    assert result is None
    assert not target.exists()
    assert not target.with_suffix(".jpg").exists()
