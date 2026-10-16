"""Behavior checks for the design section 11 changes, using local data only."""

import inspect

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from common import config
from common.frame import TimeSeriesDataFrame, split_time_intervals
from common.time import Time
from gena.frame import ParticleDataFrame
from gena import filter as gf
from plots import savefig


@pytest.fixture
def runtime_paths(monkeypatch, tmp_path):
    monkeypatch.setattr(config, "_RUNTIME_FILE_BASE", None)
    monkeypatch.setattr(config, "_RUNTIME_OUTPUT_DIR", None)
    monkeypatch.setattr(config, "_CONFIG_CACHE", {})
    monkeypatch.setenv("FILE_BASE", str(tmp_path / "absent-data"))
    monkeypatch.delenv("GENA_OUTPUT_DIR", raising=False)
    return tmp_path


def test_runtime_configuration_is_used_by_later_loaders(runtime_paths):
    base = runtime_paths / "data"
    base.mkdir()
    output = runtime_paths / "run"
    configured = config.configure(base, output_dir=output)
    assert config.get_configured_paths() == configured
    assert config.resolve_output_path("figures/a.png") == output / "figures/a.png"
    assert config.resolve_output_path(runtime_paths / "override.png") == runtime_paths / "override.png"
    returned = config.get_configured_paths()
    returned["pqt_base"] = "modified by caller"
    assert config.get_configured_paths()["pqt_base"] == str(base / "GENA_PQT")


def test_environment_changes_do_not_reuse_stale_cache(runtime_paths, monkeypatch):
    first = config.get_configured_paths(mkdir=False, validate=False)
    monkeypatch.setenv("FILE_BASE", str(runtime_paths / "second"))
    second = config.get_configured_paths(mkdir=False, validate=False)
    assert first["file_base"] != second["file_base"]
    assert not (runtime_paths / "second").exists()


def test_savefig_uses_runtime_output_without_requiring_input_data(runtime_paths):
    output = runtime_paths / "run"
    config.configure(output_dir=output, mkdir=False, validate=False)
    assert not output.exists()
    figure, axis = plt.subplots()
    axis.plot([0, 1], [0, 1])
    destination = savefig(figure, "figures/a.png")
    assert destination == output / "figures/a.png"
    assert destination.is_file()
    assert not (runtime_paths / "absent-data").exists()


def test_report_default_uses_runtime_output(runtime_paths, monkeypatch):
    from gena.post import export

    config.configure(output_dir=runtime_paths / "run", mkdir=False, validate=False)
    monkeypatch.setattr(export, "load_products", lambda *a, **k: {})
    monkeypatch.setattr(export, "batch_plot", lambda *a, **k: {})
    destinations = []
    def save(figures, prob, *, output_dir):
        destinations.append(output_dir)
        return {}
    monkeypatch.setattr(export, "save_figures", save)
    export.report("A", "2026-04-07")
    assert destinations == [runtime_paths / "run/reports/20260407"]


@pytest.mark.parametrize("prob", ["a", "b"])
def test_probe_rules_keep_mask_semantics_and_frame_metadata(prob):
    frame = ParticleDataFrame({"Datetime": pd.to_datetime([1., 2., 3.], unit="s"), "Timestamp": [1., 2., 3.], "Energy": [2., 10., 100.],
                          "TX_1s": [5., 15., 30.]})
    frame.attrs["probe"] = prob
    params = {"A_X_SUM": 10., "B_X_SUM": 20.}
    rule = gf.check("length", dir="X", edge="1", channel="s", params=params,
                    prob=prob, re_err=0) & gf.energy_eff(Emin=1, Emax=50)
    expected = (frame["TX_1s"] < params[f"{prob.upper()}_X_SUM"]) & (frame["Energy"] < 50)
    np.testing.assert_array_equal(rule(frame), expected)
    selected = frame.select(rule)
    assert isinstance(selected, ParticleDataFrame)
    assert selected.attrs == frame.attrs
    assert selected["Timestamp"].tolist() == frame.loc[expected, "Timestamp"].tolist()


def test_filter_direct_keyword_and_deferred_calls_agree():
    frame = ParticleDataFrame({"Datetime": pd.to_datetime([1., 2., 3.], unit="s"), "Timestamp": [1., 2., 3.], "Energy": [2., 10., 100.]})
    direct = gf.energy_eff(frame, Emin=4, Emax=50)
    pd.testing.assert_series_equal(direct, gf.energy_eff(df=frame, Emin=4, Emax=50))
    pd.testing.assert_series_equal(direct, gf.energy_eff(Emin=4, Emax=50)(frame))
    np.testing.assert_array_equal((~gf.energy_eff(Emin=4, Emax=50))(frame), ~direct)
    assert list(inspect.signature(gf.energy_eff).parameters) == ["df", "Emin", "Emax"]


def test_time_helpers_preserve_explicit_numeric_meaning_and_utc():
    assert Time([2451545]).datetime[0] == pd.Timestamp("1970-01-29T08:59:05")
    assert Time(2451545, format="jd").datetime == pd.Timestamp("2000-01-01T12:00:00")
    assert Time(["2026-09-17T08:00:00+08:00"]).datetime[0] == pd.Timestamp("2026-09-17")
    assert len(split_time_intervals(pd.to_datetime([0., 10., 21.], unit="s"), max_gap="10s")) == 2
