"""Offline, mocked boundary tests for the ``load`` facade.

These tests exercise the loader plumbing with fake backend payloads, so they
do not depend on local data or network access. They assert the *correct*
behavior of the fixed loader boundaries; the still-open performance issue
(``columns=`` is not pushed down for event products) is labelled explicitly.

Run with::

    python -m pytest tests/test_load_audit.py -v
"""

from __future__ import annotations

import importlib
import sys
import types

import numpy as np
import pandas as pd
import pytest

import load
import load.api.spedas as spedas_api
import load.loaders.gena as gena_loader
import load.loaders.indices as indices_loader
import load.loaders.spedas as spedas_loader

# ``load.load`` the module and ``load.load`` the function share a name, so the
# package attribute points at the function; fetch the module from sys.modules.
load_module = importlib.import_module("load.load")

TRANGE = ["2024-01-01", "2024-01-02"]


# ---------------------------------------------------------------------------
# Index dict shape
# ---------------------------------------------------------------------------
def _fake_indices_record(day: str):
    times = pd.DatetimeIndex([pd.Timestamp(day)])
    record = {
        "x": times,
        "y": np.array([-10.0]),
        "units": "nT",
        "attrs": {"source": "kyoto"},
    }
    layout = {"dst": {"kind": "scalar", "time_axis": 0, "columns": ["dst"]}}
    return {"dst": record}, layout


def test_single_interval_index_dict_is_flat(monkeypatch):
    calls = []

    def fake_read_indices(trange, product=None, *, prefix=None, download="local_first", timeout=30):
        calls.append((trange, product))
        return _fake_indices_record("2024-01-01")

    monkeypatch.setattr(indices_loader, "read_indices", fake_read_indices)
    record = load.load_indices(trange=TRANGE, product="dst", return_type="dict")
    assert set(record) == {"x", "y", "units", "attrs"}
    assert len(calls) == 1


# ---------------------------------------------------------------------------
# MMS level and datatype
# ---------------------------------------------------------------------------
def test_read_spedas_forwards_mms_level(monkeypatch):
    captured = {}

    def fake_read_mms(trange, *, product, variables, probe, prefix, download, level="l2", **kwargs):
        captured["level"] = level
        return {}, {}

    monkeypatch.setattr(spedas_api, "_read_mms", fake_read_mms)
    spedas_api.read_spedas("mms", TRANGE, product="fpi", probe=1, level="ql")
    assert captured["level"] == "ql"


def test_mms_provenance_records_requested_level(monkeypatch):
    monkeypatch.setattr(
        spedas_loader, "read_spedas",
        lambda *args, **kwargs: ({}, {}),
    )
    frame = spedas_loader.load_spedas(
        "mms", TRANGE, "fpi", probe=1, level="ql", download=False,
    )
    assert frame.attrs["level"] == "ql"


def test_mms_fpi_energy_spectrum_requests_dis_dist(monkeypatch, tmp_path):
    captured = {}

    fake_mms = types.ModuleType("pyspedas.projects.mms")

    def mms_load_fpi(**kwargs):
        captured.update(kwargs)
        return {}

    fake_mms.mms_load_fpi = mms_load_fpi
    fake_mms.mms_load_hpca = lambda **kwargs: {}
    fake_config = types.ModuleType("pyspedas.projects.mms.mms_config")
    fake_config.CONFIG = {}
    fake_projects = types.ModuleType("pyspedas.projects")
    fake_projects.mms = fake_mms

    monkeypatch.setitem(sys.modules, "pyspedas.projects", fake_projects)
    monkeypatch.setitem(sys.modules, "pyspedas.projects.mms", fake_mms)
    monkeypatch.setitem(sys.modules, "pyspedas.projects.mms.mms_config", fake_config)
    monkeypatch.setenv("MMS_DATA_DIR", str(tmp_path))

    spedas_api.read_spedas(
        "mms", TRANGE, product="fpi",
        variables=["energyspectr_omni_fast"], probe=1,
        prefix=str(tmp_path), download=False,
    )
    assert captured["datatype"] == ["dis-dist"]
    assert captured["level"] == "l2"


# ---------------------------------------------------------------------------
# load_sw per-variable fallback
# ---------------------------------------------------------------------------
def _sw_variable(value: float) -> dict:
    times = pd.DatetimeIndex([
        "2023-12-31T22:00:00",
        "2024-01-01T00:00:00",
        "2024-01-01T02:00:00",
        "2024-01-01T04:00:00",
    ])
    return {"x": times, "y": np.full(len(times), value, dtype=float)}


def _install_sw_fakes(monkeypatch, *, include_sym_h: bool):
    complete = {
        "Pressure": _sw_variable(5.0),
        "BY_GSM": _sw_variable(2.0),
        "BZ_GSM": _sw_variable(-3.0),
        "Beta": _sw_variable(0.4),
        "Mgs_mach_num": _sw_variable(6.5),
        "SYM_H": _sw_variable(-20.0),
    }
    if not include_sym_h:
        complete.pop("SYM_H")
    monkeypatch.setattr(load_module, "load_spedas", lambda *args, **kwargs: complete)
    monkeypatch.setattr(
        load_module, "load_indices",
        lambda *args, **kwargs: _sw_variable(-15.0),
    )


def test_load_sw_keeps_real_values_when_complete(monkeypatch):
    _install_sw_fakes(monkeypatch, include_sym_h=True)
    values = load_module.load_sw(["2024-01-01T00:30:00", "2024-01-01T01:30:00"])
    assert values["Pressure"][0] == pytest.approx(5.0)
    assert values["SYM_H"][0] == pytest.approx(-20.0)


def test_load_sw_falls_back_per_variable(monkeypatch):
    _install_sw_fakes(monkeypatch, include_sym_h=False)
    values = load_module.load_sw(["2024-01-01T00:30:00", "2024-01-01T01:30:00"])
    # Only the unusable variable keeps its documented constant.
    assert values["Pressure"][0] == pytest.approx(5.0)
    assert values["BY_GSM"][0] == pytest.approx(2.0)
    assert values["Mgs_mach_num"][0] == pytest.approx(6.5)
    assert values["SYM_H"][0] == pytest.approx(0.0)
    assert values["dst"][0] == pytest.approx(-15.0)


# ---------------------------------------------------------------------------
# Empty backend results are labelled
# ---------------------------------------------------------------------------
def test_empty_backend_result_is_labelled(monkeypatch):
    monkeypatch.setattr(spedas_loader, "read_spedas", lambda *args, **kwargs: ({}, {}))
    frame = spedas_loader.load_spedas(
        "themis", TRANGE, "esa", probe="b", download=False,
    )
    assert len(frame) == 0
    assert list(frame.columns) == ["Datetime", "Timestamp"]
    assert frame.attrs["variables"] == {}
    assert "empty_reason" in frame.attrs
    assert "themis esa" in frame.attrs["empty_reason"]


# ---------------------------------------------------------------------------
# Public exports and GENA frame selection
# ---------------------------------------------------------------------------
def test_assemble_is_exported():
    from load import assemble

    assert callable(assemble)


def _write_engineering_parquet(path) -> None:
    pd.DataFrame(
        {
            "Timestamp": [0.0, 1.0],
            "MJ2000_Star_X": [1.0, 2.0],
            "MJ2000_Star_Y": [3.0, 4.0],
            "MJ2000_Star_Z": [5.0, 6.0],
            "MJ2000_q1": [0.0, 0.1],
            "MJ2000_q2": [0.0, 0.1],
            "MJ2000_q3": [0.0, 0.1],
        }
    ).to_parquet(path, engine="pyarrow")


def test_conventional_pdf_name_returns_parameter_frame(tmp_path):
    path = tmp_path / "20250101_pdf.pqt"
    _write_engineering_parquet(path)
    frame = gena_loader.load_gena(file=path)
    assert frame.__class__.__name__ == "ParameterDataFrame"
    assert hasattr(frame, "pos_at")


def test_engineering_columns_determine_frame_type(tmp_path):
    path = tmp_path / "qq_eng.pqt"
    _write_engineering_parquet(path)
    frame = gena_loader.load_gena(file=path)
    assert frame.__class__.__name__ == "ParameterDataFrame"
    assert hasattr(frame, "pos_at")


def test_gena_provenance_is_attached(tmp_path):
    path = tmp_path / "20250101_l2s_a.pqt"
    pd.DataFrame({"Timestamp": [0.0], "X1": [1.0], "Tof": [10.0]}).to_parquet(
        path, engine="pyarrow"
    )
    frame = gena_loader.load_gena(file=path)
    assert frame.__class__.__name__ == "ParticleDataFrame"
    assert frame.attrs["source"] == "gena"
    assert frame.attrs["product"] == "l2s_a"
    assert frame.attrs["level"] == "l2s"
    assert frame.attrs["probe"] == "a"


def test_columns_are_not_pushed_down_for_event_products(tmp_path, monkeypatch):
    """Known performance issue: event products always read every column."""
    path = tmp_path / "20250101_l2_a.pqt"
    pd.DataFrame({"Timestamp": [0.0, 1.0], "Energy": [5.0, 6.0]}).to_parquet(
        path, engine="pyarrow"
    )
    seen = {}
    original = gena_loader.pq.read_table

    def spy(file_path, *, filters=None, columns=None):
        seen["columns"] = columns
        return original(file_path, filters=filters, columns=columns)

    monkeypatch.setattr(gena_loader.pq, "read_table", spy)
    frame = gena_loader.load_gena(file=path, columns=["Energy"])
    assert list(frame.columns) == ["Timestamp", "Energy", "Datetime"]
    assert seen["columns"] is None
