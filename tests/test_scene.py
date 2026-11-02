"""Contract tests for the shared Plotly scene layer."""

from __future__ import annotations

from pathlib import Path

import plotly.graph_objects as go
import pytest

from plots.scene import SceneBundle, SceneSpec, add_target, target_point, to_html, write_html


def test_target_points():
    assert target_point("earth") == (0.0, 0.0, 0.0)
    assert target_point("mid_tail") == (-30.0, 0.0, 0.0)
    assert target_point("far_tail") == (-60.0, 0.0, 0.0)
    assert target_point("other") is None
    assert target_point("unknown") is None


def test_add_target_marks_and_planes():
    figure = go.Figure()
    add_target(figure, "mid_tail")
    marker = [trace for trace in figure.data if trace.name == "mid_tail target"]
    assert len(marker) == 1
    assert len(figure.data) >= 2  # marker + target plane


def test_add_target_ignores_other():
    figure = go.Figure()
    add_target(figure, "other")
    assert len(figure.data) == 0


def test_to_html_uses_shared_plotly():
    figure = go.Figure(go.Scatter3d(x=[0], y=[0], z=[0], mode="markers"))
    document = to_html(figure, include_plotlyjs=False, plotlyjs_url="/static/plotly.min.js")
    assert '/static/plotly.min.js' in document
    assert 'id="gena-scene"' in document
    assert "cdn.plot.ly" not in document


def test_write_html_is_atomic(tmp_path: Path):
    figure = go.Figure(go.Scatter3d(x=[0], y=[0], z=[0], mode="markers"))
    target = tmp_path / "scene.html"
    result = write_html(figure, target)
    assert result == target
    assert target.exists()
    assert "gena-scene" in target.read_text(encoding="utf-8")


def test_scene_spec_validation():
    with pytest.raises(ValueError):
        SceneSpec(trange=("2025-12-04", "2025-12-05"), target_mode="lunar")
    with pytest.raises(ValueError):
        SceneSpec(trange=("2025-12-04", "2025-12-05"), quality="ultra")


def test_scene_bundle_reports_weather_and_omissions():
    class DummyContext:
        cache_timestamps = [0.0, 1.0, 2.0]
        space_weather = {"BZ_GSM": [0.0, 0.0, 0.0]}
        space_weather_source = "provided"
        omitted_satellites = ("THA",)

        def make_scene(self, time):
            return go.Figure()

    bundle = SceneBundle(
        context=DummyContext(),
        spec=SceneSpec(trange=("2025-12-04", "2025-12-05"), target_mode="far_tail"),
    )
    assert bundle.omitted_satellites == ("THA",)
    assert bundle.space_weather_status == "provided"
    figure = bundle.make_scene("2025-12-04T12:00:00")
    assert any(trace.name == "far_tail target" for trace in figure.data)


def test_scene_bundle_marks_loader_weather_source():
    class DummyContext:
        space_weather = {"BZ_GSM": [0.0]}
        space_weather_source = "load_sw"
        omitted_satellites = ()

    bundle = SceneBundle(
        context=DummyContext(),
        spec=SceneSpec(trange=("2025-12-04", "2025-12-05")),
    )
    assert bundle.space_weather_status == "load_sw"
