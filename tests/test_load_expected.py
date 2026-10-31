"""Expected-behavior contract tests for the loaders.

These started as ``xfail`` entries for the reported defects; after the fixes
they are ordinary regression tests. They are intentionally redundant with
``tests/test_load_real_cases.py`` (which also exercises the happy path): this
file is the short list that must never regress.

All cases read the local data root and never download; cases skip themselves
when the expected local file is absent.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import load.api.spedas as spedas_api
import load.loaders.gena as gena_loader
from common.config import resolve_paths
from load import load, load_pos, load_sw

TRANGE = ["20251204", "20251205"]
TIMES = ["2025-12-04T00:00:00", "2025-12-04T00:05:00"]

_PATHS = resolve_paths()
PQT_BASE = Path(_PATHS["pqt_base"])
SPEDAS_BASE = Path(_PATHS["spedas_base"])


def _require(path: Path) -> Path:
    if not Path(path).exists():
        pytest.skip(f"local data missing: {path}")
    return Path(path)


def _require_state_file() -> Path:
    matches = sorted(
        SPEDAS_BASE.glob("themis_data/thb/l1/state/2025/thb_l1_state_20251204.cdf")
    )
    if not matches:
        pytest.skip("THEMIS state CDF for 2025-12-04 is missing")
    return matches[0]


# ---------------------------------------------------------------------------
# THEMIS position
# ---------------------------------------------------------------------------
def test_themis_state_product_loads():
    _require_state_file()
    frame = load("themis", TRANGE, "state", probe="b", download=False)
    assert len(frame) > 0
    assert "thb_pos_gsm_x" in frame.columns
    assert frame.attrs["level"] == "l1"


def test_load_pos_themis_returns_positions():
    _require_state_file()
    positions = load_pos("themis", TIMES, "b")
    assert np.asarray(positions).shape == (2, 3)
    assert np.isfinite(positions).all()


# ---------------------------------------------------------------------------
# MMS: an empty local result is labelled
# ---------------------------------------------------------------------------
def test_empty_mms_result_is_distinguishable():
    local = list(SPEDAS_BASE.glob("mms_data/mms1/fpi/**/mms1_fpi_*_2025120*.cdf"))
    if local:
        pytest.skip("local MMS FPI files exist for this window")
    frame = load("mms", TRANGE, "fpi", download=False)
    assert len(frame) > 0 or "empty_reason" in frame.attrs


def test_mms_level_is_forwarded(monkeypatch):
    captured = {}

    def fake_read_mms(trange, *, product, variables, probe, base, download, level="l2", **kwargs):
        captured["level"] = level
        return {}, {}

    monkeypatch.setattr(spedas_api, "_read_mms", fake_read_mms)
    spedas_api.read_spedas("mms", TRANGE, product="fpi", probe=1, level="ql")
    assert captured["level"] == "ql"


# ---------------------------------------------------------------------------
# load_sw: real OMNI values survive a single unusable variable
# ---------------------------------------------------------------------------
def test_load_sw_keeps_real_omni_values():
    _require(SPEDAS_BASE / "omni_data/hro_1min/2025/omni_hro_1min_20251201_v01.cdf")
    _require(SPEDAS_BASE / "indices/kyoto/hour_dst/202512.json")
    values = load_sw(TIMES[:1] + ["2025-12-04T06:00:00"], download=False)
    fallback = np.array([1.0, 0.0, 0.0, 1.0, 5.0, 0.0])
    actual = np.array([
        values["Pressure"][0], values["BY_GSM"][0], values["BZ_GSM"][0],
        values["Beta"][0], values["Mgs_mach_num"][0], values["SYM_H"][0],
    ])
    assert not np.allclose(actual, fallback)


# ---------------------------------------------------------------------------
# GENA: provenance and frame type
# ---------------------------------------------------------------------------
def test_gena_frame_carries_provenance():
    _require(PQT_BASE / "20251204_l2s_a.pqt")
    frame = load("gena", TRANGE, "l2s_a")
    assert frame.attrs["source"] == "gena"
    assert frame.attrs["product"] == "l2s_a"


def test_gena_frame_type_follows_data(tmp_path):
    path = tmp_path / "qq_eng.pqt"
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
    frame = gena_loader.load_gena(file=path)
    assert hasattr(frame, "pos_at")


# ---------------------------------------------------------------------------
# Documented facade exports
