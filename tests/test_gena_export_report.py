from pathlib import Path
import importlib

import matplotlib
import numpy as np
import pandas as pd
import pytest
from PIL import Image
from openpyxl import load_workbook

matplotlib.use("Agg")

from gena.frame import ParameterDataFrame, ParticleDataFrame
from gena.post import export

gena_load = importlib.import_module("load.load")


def _write(path: Path, data):
    pd.DataFrame(data).to_parquet(path, index=False)


def _configure_pqt(monkeypatch, path):
    monkeypatch.setattr(
        gena_load,
        "get_configured_paths",
        lambda: {"pqt_base": str(path)},
    )


def _make_daily_products(folder: Path, day="20260407", prob="a"):
    times = pd.date_range("2026-04-07T01:00:00Z", periods=30, freq="min")
    timestamps = np.array([value.timestamp() for value in times])
    engineering = {
        "Timestamp": timestamps,
        "HV_A": np.full(30, 2850.0),
        "Th_A_MCP": np.full(30, 500.0),
        "Th_A_X": np.full(30, 200.0),
        "Th_A_Y": np.full(30, 210.0),
        "Temp_A": np.linspace(5.0, 8.0, 30),
        "HV_B": np.full(30, 2850.0),
        "Th_B_MCP": np.full(30, 500.0),
        "Th_B_X1": np.full(30, 100.0),
        "Th_B_X2": np.full(30, 100.0),
        "Th_B_Y1": np.full(30, 105.0),
        "Th_B_Y2": np.full(30, 105.0),
        "Th_B_G": np.full(30, 10.0),
        "Temp_B": np.linspace(5.0, 8.0, 30),
        "MJ2000_Star_X": np.linspace(60_000_000.0, 61_000_000.0, 30),
        "MJ2000_Star_Y": np.linspace(10_000_000.0, 11_000_000.0, 30),
        "MJ2000_Star_Z": np.linspace(5_000_000.0, 5_500_000.0, 30),
        "MJ2000_q1": np.full(30, 0.10),
        "MJ2000_q2": np.full(30, 0.12),
        "MJ2000_q3": np.full(30, 0.08),
    }
    _write(folder / f"{day}_pdf.pqt", engineering)

    l1a = {"UUID": [f"event-{i}" for i in range(30)], "Timestamp": timestamps}
    for offset, name in enumerate(("TX_1s", "TX_2s", "TY_1s", "TY_2s")):
        l1a[name] = 20.0 + offset + np.arange(30) % 10
    for offset, name in enumerate(("TX_1e", "TX_2e", "TY_1e", "TY_2e")):
        l1a[name] = 100.0 + offset + np.arange(30) % 20
    _write(folder / f"{day}_l1a_{prob}.pqt", l1a)

    l2 = {
        "UUID": [f"event-{i}" for i in range(30)],
        "Timestamp": timestamps,
        "Energy": np.geomspace(0.6, 80.0, 30),
        "X1": np.linspace(-40.0, 40.0, 30),
        "Y1": np.linspace(-15.0, 15.0, 30),
        "X2": np.linspace(-35.0, 35.0, 30),
        "Y2": np.linspace(-12.0, 12.0, 30),
        "PX": np.linspace(-20.0, 20.0, 30),
        "PY": np.linspace(-10.0, 10.0, 30),
    }
    for product in ("l2", "l2s", "l2g"):
        _write(folder / f"{day}_{product}_{prob}.pqt", l2)


def test_normalize_report_times_infers_daily_range_and_exact_interval():
    single = export.normalize_report_times("2026-04-07")
    assert len(single) == 1
    assert single[0][0] == pd.Timestamp("2026-04-07T00:00:00")
    assert single[0][1] == pd.Timestamp("2026-04-07T23:59:59.999999999")

    daily = export.normalize_report_times(["2026-04-07", "2026-04-09"])
    assert [period[0].strftime("%Y%m%d") for period in daily] == [
        "20260407",
        "20260408",
        "20260409",
    ]

    exact = export.normalize_report_times(
        ["2026-04-07T23:50:00Z", "2026-04-08T00:10:00Z"]
    )
    assert exact == (
        (
            pd.Timestamp("2026-04-07T23:50:00"),
            pd.Timestamp("2026-04-08T00:10:00"),
        ),
    )

    typed_daily = export.normalize_report_times(
        pd.date_range("2026-04-07", "2026-04-09", freq="D")
    )
    assert [period[0].strftime("%Y%m%d") for period in typed_daily] == [
        "20260407",
        "20260408",
        "20260409",
    ]


def test_load_products_clips_one_exact_interval(tmp_path, monkeypatch):
    _configure_pqt(monkeypatch, tmp_path)
    _make_daily_products(tmp_path)

    products = export.load_products(
        ["2026-04-07T01:05:00Z", "2026-04-07T01:10:00Z"],
        "a",
    )

    assert {name: len(frame) for name, frame in products.items()} == {
        "pdf": 6,
        "l1a": 6,
        "l2": 6,
        "l2s": 6,
        "l2g": 6,
    }


def test_report_time_input_validation(monkeypatch, tmp_path):
    captured = []

    def fake_load_products(period, prob, **kwargs):
        captured.append((period[0].strftime("%Y%m%d"), prob))
        return {key: pd.DataFrame() for key in ("pdf", "l1a", "l2", "l2s", "l2g")}

    monkeypatch.setattr(export, "load_products", fake_load_products)
    monkeypatch.setattr(export, "batch_plot", lambda *args, **kwargs: {})
    monkeypatch.setattr(export, "save_figures", lambda *args, **kwargs: {})

    export.report(
        "A",
        ["2026-04-07", "2026-04-09"],
        output_dir=tmp_path,
    )
    assert captured == [
        ("20260407", "a"),
        ("20260408", "a"),
        ("20260409", "a"),
    ]

    captured.clear()
    export.report(
        "B",
        ["2026-04-27", "2026-04-07", "2026-04-27"],
        output_dir=tmp_path,
    )
    assert captured == [("20260407", "b"), ("20260427", "b")]

    with pytest.raises(ValueError, match="prob"):
        export.report("C", "2026-04-07")


@pytest.mark.parametrize("prob", ["a", "b"])
def test_build_observation_summary_splits_and_classifies_pointing(monkeypatch, prob):
    dates = pd.to_datetime(
        [
            "2026-04-07T00:00:00Z",
            "2026-04-07T00:01:00Z",
            "2026-04-07T00:10:00Z",
            "2026-04-07T02:00:00Z",
        ]
    )
    pdf = ParameterDataFrame({"Datetime": pd.to_datetime([value.timestamp() for value in dates], unit="s"), "Timestamp": [value.timestamp() for value in dates]})
    ena_dates = pd.to_datetime(
        [
            "2026-04-07T00:00:30Z",
            "2026-04-07T00:00:45Z",
            "2026-04-07T00:10:00Z",
            "2026-04-07T02:00:00Z",
        ]
    )
    ena = ParticleDataFrame({"Datetime": pd.to_datetime([value.timestamp() for value in ena_dates], unit="s"), "Timestamp": [value.timestamp() for value in ena_dates]})
    probes = []

    def fake_pointing_geometry(frame, *, probe, coord, forward):
        probes.append(probe)
        start = frame["Datetime"].iloc[0]
        angle = 0.0 if start.minute == 0 and start.hour == 0 else (30.0 if start.hour == 0 else 150.0)
        angle_rad = np.deg2rad(angle)
        return pd.DataFrame(
            {
                "Timestamp": frame["Timestamp"].to_numpy(),
                "Datetime": frame["Datetime"].to_numpy(),
                "spacecraft_x": np.full(len(frame), -60.0),
                "spacecraft_y": np.zeros(len(frame)),
                "spacecraft_z": np.ones(len(frame)),
                "look_x": np.full(len(frame), np.cos(angle_rad)),
                "look_y": np.full(len(frame), np.sin(angle_rad)),
                "look_z": np.zeros(len(frame)),
                "xz_x": np.zeros(len(frame)),
                "xz_z": np.zeros(len(frame)),
            }
        )

    monkeypatch.setattr(export.gut, "pointing_geometry", fake_pointing_geometry)
    result = export.build_observation_summary(pdf, ena, prob)

    assert result["target_mode"].tolist() == ["earth", "mid_tail", "far_tail"]
    np.testing.assert_allclose(result["sun_earth_angle_deg"], [0.0, 30.0, 150.0])
    assert result["parent_l_level_index"].tolist() == [1, 1, 2]
    assert result["sub_observation_index"].tolist() == [1, 2, 1]
    assert result["ena_total_counts"].tolist() == [2, 1, 1]
    assert result["ena_product"].unique().tolist() == [f"l2g_{prob}"]
    assert probes == [prob.upper()] * 3


def test_observation_summary_uses_report_time_rules_and_writes_workbook(
    monkeypatch, tmp_path
):
    data_dir = tmp_path / "data"
    output_dir = tmp_path / "reports"
    data_dir.mkdir()
    _make_daily_products(data_dir)

    def fake_pointing_geometry(frame, **kwargs):
        angle_rad = np.deg2rad(150.0)
        return pd.DataFrame(
            {
                "Timestamp": frame["Timestamp"].to_numpy(),
                "Datetime": frame["Datetime"].to_numpy(),
                "spacecraft_x": np.full(len(frame), -60.0),
                "spacecraft_y": np.zeros(len(frame)),
                "spacecraft_z": np.ones(len(frame)),
                "look_x": np.full(len(frame), np.cos(angle_rad)),
                "look_y": np.full(len(frame), np.sin(angle_rad)),
                "look_z": np.zeros(len(frame)),
                "xz_x": np.full(len(frame), -60.0),
                "xz_z": np.zeros(len(frame)),
            }
        )

    monkeypatch.setattr(export.gut, "pointing_geometry", fake_pointing_geometry)
    _configure_pqt(monkeypatch, data_dir)
    result = export.observation_summary(
        "A",
        ["2026-04-07T01:05:00Z", "2026-04-07T01:10:00Z"],
        output_dir=output_dir,
    )

    output_path = output_dir / "observation_summary.xlsx"
    assert result.attrs["output_path"] == output_path.resolve()
    assert result["start_time"].iloc[0] == pd.Timestamp("2026-04-07T01:05:00")
    assert result["end_time"].iloc[0] == pd.Timestamp("2026-04-07T01:10:00")
    assert result["target_mode"].tolist() == ["far_tail"]
    assert output_path.is_file()
    workbook = load_workbook(output_path, read_only=False)
    sheet = workbook["observations"]
    assert sheet.freeze_panes == "A2"
    assert sheet.auto_filter.ref == sheet.dimensions
    assert sheet.cell(1, 1).value == "observation_index"


def test_report_automatic_time_keeps_exact_interval_as_one_artifact_group(
    monkeypatch, tmp_path
):
    captured = {}
    empty_products = {
        key: pd.DataFrame()
        for key in ("pdf", "l1a", "l2", "l2s", "l2g")
    }

    def fake_load_products(period, prob, **kwargs):
        captured["period"] = period
        return empty_products

    monkeypatch.setattr(export, "load_products", fake_load_products)
    monkeypatch.setattr(export, "batch_plot", lambda *args, **kwargs: {})
    monkeypatch.setattr(export, "save_figures", lambda *args, **kwargs: {})

    result = export.report(
        "A",
        ["2026-04-07T23:50:00Z", "2026-04-08T00:10:00Z"],
        output_dir=tmp_path,
    )

    assert captured["period"] == (
        pd.Timestamp("2026-04-07T23:50:00"),
        pd.Timestamp("2026-04-08T00:10:00"),
    )
    assert list(result) == ["20260407T235000Z_20260408T001000Z"]


def test_report_delegates_geometry_for_one_day(monkeypatch, tmp_path):
    captured = {}

    def fake_save_geometry(period, **kwargs):
        captured.update(period=period, **kwargs)
        return kwargs["output_path"]

    empty_products = {
        key: pd.DataFrame()
        for key in ("pdf", "l1a", "l2", "l2s", "l2g")
    }
    monkeypatch.setattr(export, "load_products", lambda *args, **kwargs: empty_products)
    monkeypatch.setattr(export, "batch_plot", lambda *args, **kwargs: {})
    monkeypatch.setattr(export, "save_figures", lambda *args, **kwargs: {})
    monkeypatch.setattr(export, "_save_geometry", fake_save_geometry)

    result = export.report(
        "A",
        "2026-04-07",
        output_dir=tmp_path,
        include_geometry=True,
        geometry_kwargs={
            "max_frames": 12,
            "field_line_nlat": 7,
            "field_line_cadence": "2h",
            "update_ephemeris": True,
        },
    )

    assert result["20260407"]["geometry"] == tmp_path / "20260407" / "3d.html"
    assert captured["period"] == (
        pd.Timestamp("2026-04-07T00:00:00"),
        pd.Timestamp("2026-04-07T23:59:59.999999999"),
    )
    assert captured["geometry_options"]["max_frames"] == 12
    assert captured["geometry_options"]["field_line_nlat"] == 7
    assert captured["geometry_options"]["field_line_cadence"] == "2h"
    assert captured["geometry_options"]["update_ephemeris"] is True


def test_report_exact_interval_geometry_uses_one_preloaded_scene(monkeypatch, tmp_path):
    captured = {}
    empty_products = {
        key: pd.DataFrame()
        for key in ("pdf", "l1a", "l2", "l2s", "l2g")
    }

    class FakeFigure:
        def write_html(self, path, **kwargs):
            captured["write_html"] = kwargs
            Path(path).write_text("<html></html>")

    class FakeContext:
        def make_timerange(self, **kwargs):
            captured["make_timerange"] = kwargs
            return FakeFigure()

    def fake_prepare_scene(trange, **kwargs):
        captured["trange"] = trange
        captured["prepare_scene"] = kwargs
        return FakeContext()

    monkeypatch.setattr(export, "load_products", lambda *args, **kwargs: empty_products)
    monkeypatch.setattr(export, "batch_plot", lambda *args, **kwargs: {})
    monkeypatch.setattr(export, "save_figures", lambda *args, **kwargs: {})
    monkeypatch.setattr(export.gena_v3d, "prepare_scene", fake_prepare_scene)

    result = export.report(
        "A",
        ["2026-04-07T23:50:00Z", "2026-04-08T00:10:00Z"],
        output_dir=tmp_path,
        include_geometry=True,
        geometry_kwargs={"max_frames": 12, "field_line_nlat": 7},
    )

    key = "20260407T235000Z_20260408T001000Z"
    assert captured["trange"] == [
        pd.Timestamp("2026-04-07T23:50:00"),
        pd.Timestamp("2026-04-08T00:10:00"),
    ]
    assert captured["prepare_scene"]["field_line_nlat"] == 7
    assert captured["make_timerange"] == {"max_frames": 12}
    assert captured["write_html"]["auto_play"] is False
    assert result[key]["geometry"] == tmp_path / key / "3d.html"
    assert result[key]["geometry"].is_file()


def test_report_delegates_animation_to_gena_ani(monkeypatch, tmp_path):
    captured = {}
    products = {
        "pdf": pd.DataFrame(),
        "l1a": pd.DataFrame(),
        "l2": pd.DataFrame({"Energy": [10.0]}),
        "l2s": pd.DataFrame(),
        "l2g": pd.DataFrame(),
    }

    def fake_dump_fov(l2, output_path, **kwargs):
        captured.update(l2=l2, output_path=output_path, **kwargs)
        return output_path

    monkeypatch.setattr(export, "load_products", lambda *args, **kwargs: products)
    monkeypatch.setattr(export, "batch_plot", lambda *args, **kwargs: {})
    monkeypatch.setattr(export, "save_figures", lambda *args, **kwargs: {})
    monkeypatch.setattr(export.gena_ani, "dumpFov", fake_dump_fov)

    result = export.report(
        "B",
        "2026-04-07",
        output_dir=tmp_path,
        include_animation=True,
        window="2min",
    )

    assert captured["l2"] is products["l2"]
    assert captured["output_path"] == tmp_path / "20260407" / "b_fov.mp4"
    assert captured["window"] == "2min"
    assert result["20260407"]["fov_animation"] == captured["output_path"]


@pytest.mark.parametrize("prob", ["a", "b"])
def test_report_generates_independent_300_dpi_jpegs(tmp_path, monkeypatch, prob):
    data_dir = tmp_path / "data"
    output_dir = tmp_path / "reports"
    data_dir.mkdir()
    _configure_pqt(monkeypatch, data_dir)
    _make_daily_products(data_dir, prob=prob)

    result = export.report(
        prob,
        "2026-04-07",
        output_dir=output_dir,
    )

    artifacts = result["20260407"]
    assert set(artifacts) == {
        "summary",
        "levels_count",
        "parameters",
        "orbit_pointing",
        "start_signal_distribution",
        "end_signal_distribution",
        "joint_signal_distribution",
        "hist_E",
        "hist_Layer1",
        "hist_Layer2",
        "hist_FOV",
    }
    assert all(path.is_file() and path.stat().st_size > 0 for path in artifacts.values())
    assert all(path.suffix == ".jpg" for path in artifacts.values())
    assert not list(output_dir.rglob("*.pdf"))
    for path in artifacts.values():
        with Image.open(path) as image:
            assert image.format == "JPEG"
            assert image.info["dpi"] == pytest.approx((300, 300), rel=0.01)
