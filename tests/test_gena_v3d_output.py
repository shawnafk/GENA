import importlib

import numpy as np
import pandas as pd
import pytest

from plots import v3d
from plots.basic import gena as vis

gena_load = importlib.import_module("load.load")


def _configure_pqt(monkeypatch, path):
    configured = lambda: {"pqt_base": str(path)}
    monkeypatch.setattr(gena_load, "get_configured_paths", configured)
    monkeypatch.setattr(v3d, "get_configured_paths", configured)


def _write_pdf(path):
    times = pd.date_range("2026-04-07T01:00:00Z", periods=4, freq="10min")
    pd.DataFrame(
        {
            "Timestamp": [value.timestamp() for value in times],
            "MJ2000_Star_X": np.linspace(60_000_000.0, 61_000_000.0, 4),
            "MJ2000_Star_Y": np.linspace(10_000_000.0, 11_000_000.0, 4),
            "MJ2000_Star_Z": np.linspace(5_000_000.0, 5_500_000.0, 4),
            "MJ2000_q1": np.full(4, 0.10),
            "MJ2000_q2": np.full(4, 0.12),
            "MJ2000_q3": np.full(4, 0.08),
        }
    ).to_parquet(path, index=False)


def test_dump_geo_has_animation_controls(tmp_path, monkeypatch):
    pqt_dir = tmp_path / "pqt"
    output_dir = tmp_path / "results"
    pqt_dir.mkdir()
    _configure_pqt(monkeypatch, pqt_dir)
    _write_pdf(pqt_dir / "20260407_pdf.pqt")

    result = v3d.dumpGeo(
        dates="2026-04-07",
        output_dir=output_dir,
        max_frames=3,
        use_average_sw_params=False,
        include_field_lines=False,
        include_neutral_sheet=False,
        include_plotlyjs=False,
    )

    path = result["20260407"]
    html = path.read_text(encoding="utf-8")
    assert path == output_dir / "20260407" / "3d.html"
    assert "Plotly.newPlot" in html
    assert "Plotly.addFrames" in html
    assert "Play" in html
    assert "Pause" in html
    assert "X GSM (R_E)" in html


def test_no_date_arguments_discovers_every_local_pdf_day(tmp_path, monkeypatch):
    _configure_pqt(monkeypatch, tmp_path)
    for day in ("20260407", "20260427"):
        (tmp_path / f"{day}_pdf.pqt").touch()
    calls = []

    def fake_output_one_day(day, **kwargs):
        calls.append(day)
        return tmp_path / day / "3d.html"

    monkeypatch.setattr(v3d, "_output_one_day", fake_output_one_day)
    result = v3d.dumpGeo(
        output_dir=tmp_path / "results",
    )

    assert calls == ["20260407", "20260427"]
    assert list(result) == ["20260407", "20260427"]


def test_vis_does_not_duplicate_animation_or_3d_implementations():
    assert not hasattr(vis, "AnimationDemo3")
    assert not hasattr(vis, "SceneContext")
    assert not hasattr(vis, "dumpGeo")


def test_scene_rejects_nonphysical_satellite_coordinates():
    dates = pd.date_range("2024-06-26", periods=2, freq="10min").tolist()
    with pytest.raises(ValueError, match="expects GSM positions in R_E"):
        v3d.scene(
            dates,
            (np.zeros((2, 3)), None),
            [np.full((3, 2), 1e11)],
            None,
            None,
        )


def test_scene_keeps_spacecraft_metres_and_satellite_earth_radii_distinct(
    monkeypatch,
):
    class IdentityCoord:
        def __init__(self, dates, values, system, rotation):
            self.values = np.asarray(values, dtype=float)
            self.system = system
            self.converted = {}

        def convert(self, system):
            self.converted[system] = self.values

        def __getitem__(self, system):
            return self.converted.get(system, self.values)

    class Boundary:
        @staticmethod
        def surface_xyz():
            zeros = np.zeros((2, 2))
            return zeros, zeros, zeros

    def moon_in_km(times):
        count = 1 if np.isscalar(times) else len(times)
        moon_x_km = 10 * v3d.RE / 1000
        return np.tile([moon_x_km, 0.0, 0.0], (count, 1))

    monkeypatch.setattr(v3d, "coord", IdentityCoord)
    monkeypatch.setattr(v3d.epms, "moon", moon_in_km)
    satellite_re = np.array([[2.0], [3.0], [4.0]])
    trajectory_re = np.array(
        [
            [1.5, 2.0, 2.5],
            [2.5, 3.0, 3.5],
            [3.5, 4.0, 4.5],
        ]
    )

    figure = v3d.scene(
        [pd.Timestamp("2026-04-07T01:50:00Z")],
        (np.array([[v3d.RE, 0.0, 0.0]]), None),
        [satellite_re],
        Boundary(),
        Boundary(),
        sat_trajectories=[trajectory_re],
        satellite_names=["THA"],
    )

    trajectory = next(trace for trace in figure.data if trace.name == "THA trajectory")
    current = next(
        trace for trace in figure.data if trace.name == "THA current position"
    )
    spacecraft = next(
        trace
        for trace in figure.data
        if trace.type == "scatter3d"
        and trace.mode == "markers"
        and trace.name is None
        and trace.marker.symbol == "circle"
    )
    np.testing.assert_allclose(trajectory.x, [1.5, 2.0, 2.5])
    assert current.x[0] == 2.0
    assert current.y[0] == 3.0
    assert current.z[0] == 4.0
    assert current.marker.size == 6
    assert spacecraft.x[0] == 11.0
    assert any(trace.name == "THA direction" for trace in figure.data)
    assert any(trace.name == "GSM Z=0 reference" for trace in figure.data)
    assert any(trace.name == "GSM Y=0 reference" for trace in figure.data)


def test_scene_context_loads_once_for_repeated_instantaneous_scenes(
    tmp_path, monkeypatch
):
    pqt_dir = tmp_path / "pqt"
    pqt_dir.mkdir()
    _configure_pqt(monkeypatch, pqt_dir)
    _write_pdf(pqt_dir / "20260407_pdf.pqt")
    calls = []
    real_load = v3d.load_gena

    def counted_load(*args, **kwargs):
        calls.append((args, kwargs))
        return real_load(*args, **kwargs)

    monkeypatch.setattr(v3d, "load_gena", counted_load)
    context = v3d.prepare_scene(
        ["2026-04-07T01:00:00Z", "2026-04-07T01:30:00Z"],
        satellites=(),
        use_space_weather=False,
        include_field_lines=False,
        include_neutral_sheet=False,
    )

    first = context.make_scene("2026-04-07T01:05:00Z")
    second = v3d.make_scene(
        "2026-04-07T01:15:00Z",
        context=context,
    )
    animation = context.make_timerange(max_frames=3)

    assert len(calls) == 1
    assert len(first.frames) == 0
    assert len(second.frames) == 0
    assert first.layout.title.text == "time: 2026-04-07 01:05:00"
    assert len(animation.frames) == 3
    assert animation.layout.updatemenus[0].buttons[0].label == "Play"


def test_scene_context_preloads_satellite_trajectory_with_one_hour_padding(
    tmp_path,
    monkeypatch,
):
    pqt_dir = tmp_path / "pqt"
    pqt_dir.mkdir()
    _configure_pqt(monkeypatch, pqt_dir)
    _write_pdf(pqt_dir / "20260407_pdf.pqt")
    calls = []

    def fake_satellites(timestamps, satellite, *, update):
        timestamps = np.asarray(timestamps, dtype=float)
        calls.append((timestamps.copy(), satellite, update))
        return np.column_stack(
            [
                np.linspace(2.0, 2.5, len(timestamps)),
                np.linspace(3.0, 3.5, len(timestamps)),
                np.linspace(4.0, 4.5, len(timestamps)),
            ]
        )

    monkeypatch.setattr(v3d.epms, "satellites", fake_satellites)
    context = v3d.prepare_scene(
        ["2026-04-07T01:00:00Z", "2026-04-07T01:30:00Z"],
        satellites=("THA",),
        use_space_weather=False,
        include_field_lines=False,
        include_neutral_sheet=False,
    )

    assert len(calls) == 1
    timestamps, satellite, update = calls[0]
    assert satellite == "THA"
    assert update is False
    assert timestamps[0] == pd.Timestamp("2026-04-07T00:00:00Z").timestamp()
    assert timestamps[-1] == pd.Timestamp("2026-04-07T02:30:00Z").timestamp()
    np.testing.assert_array_equal(context.orbit_timestamps, timestamps)


def test_scene_context_caches_time_dependent_t96_field_lines(
    tmp_path,
    monkeypatch,
):
    pqt_dir = tmp_path / "pqt"
    pqt_dir.mkdir()
    _configure_pqt(monkeypatch, pqt_dir)
    _write_pdf(pqt_dir / "20260407_pdf.pqt")
    calls = []

    class FakeT96:
        def __init__(self, ut, pressure, dst, by, bz, Nlat):
            calls.append((ut, pressure, dst, by, bz, Nlat))
            self.ut = ut

        def trace_lines_y0(self, *, footprint_latitude_range, rlim):
            assert footprint_latitude_range == (70.0, 82.0)
            return [
                np.array(
                    [
                        [-20.0, -10.0, 0.0],
                        [0.0, 0.0, 0.0],
                        [self.ut % 10, self.ut % 10 + 1, self.ut % 10 + 2],
                    ]
                )
            ]

    monkeypatch.setattr(v3d, "T96_fieldlines", FakeT96)
    context = v3d.prepare_scene(
        ["2026-04-07T01:00:00Z", "2026-04-07T01:30:00Z"],
        satellites=(),
        use_space_weather=False,
        field_line_nlat=3,
        field_line_cadence="1h",
        include_neutral_sheet=False,
    )

    first = context.make_scene("2026-04-07T01:05:00Z")
    context.make_scene("2026-04-07T01:05:00Z")
    animation = context.make_timerange(max_frames=3)

    assert len(calls) == 2
    assert calls[0][1:] == (5.0, -40.0, 0.0, 0.0, 3)
    assert any(trace.name == "T96 field lines" for trace in first.data)
    assert all(
        any(trace.name == "T96 field lines" for trace in frame.data)
        for frame in animation.frames
    )


def test_scene_context_caches_and_annotates_t96_neutral_sheet(
    tmp_path,
    monkeypatch,
):
    pqt_dir = tmp_path / "pqt"
    pqt_dir.mkdir()
    _configure_pqt(monkeypatch, pqt_dir)
    _write_pdf(pqt_dir / "20260407_pdf.pqt")
    calls = []

    class FakeNeutralSheet:
        def __init__(self, ut, pressure, dst, by, bz, *, half_thickness):
            calls.append((ut, pressure, dst, by, bz, half_thickness))
            self.ut = ut
            self.parameters = np.array([pressure, dst, by, bz])
            self.dipole_tilt_degrees = 4.25
            self.half_thickness = half_thickness

        @staticmethod
        def surface_xyz():
            x, y = np.meshgrid([-20.0, -10.0], [-5.0, 5.0])
            return x, y, np.zeros_like(x)

        def thickness_xyz(self):
            x, y, z = self.surface_xyz()
            return (x, y, z - self.half_thickness), (x, y, z + self.half_thickness)

        @staticmethod
        def z_at(x, y):
            return 0.5

    monkeypatch.setattr(v3d, "T96NeutralSheet", FakeNeutralSheet)
    context = v3d.prepare_scene(
        ["2026-04-07T01:00:00Z", "2026-04-07T01:30:00Z"],
        satellites=(),
        use_space_weather=False,
        include_field_lines=False,
        neutral_sheet_half_thickness=2.5,
    )

    first = context.make_scene("2026-04-07T01:05:00Z")
    context.make_scene("2026-04-07T01:05:00Z")
    animation = context.make_timerange(max_frames=3)

    assert len(calls) == 2
    assert all(call[1:] == (5.0, -40.0, 0.0, 0.0, 2.5) for call in calls)
    assert any(trace.name == "T96 Bx=0 neutral sheet" for trace in first.data)
    assert "Dipole tilt ψ: +4.25°" in first.layout.annotations[0].text
    assert all(
        any(trace.name == "T96 Bx=0 neutral sheet" for trace in frame.data)
        and "T96+IGRF Bx=0" in frame.layout.annotations[0].text
        for frame in animation.frames
    )
