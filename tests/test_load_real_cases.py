"""Real, offline integration cases for the public ``load`` facade.

Every case here uses the configured local data root (``FILE_BASE``) and the
2025-12-04/05 window present on this machine. Nothing downloads: calls pass
``download=False`` or use loaders that only read local Parquet/CDF files.
Cases skip themselves when the expected local file is absent.

Run with::

    python -m pytest tests/test_load_real_cases.py -v
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import load as load_pkg
from common.config import resolve_paths
from common.interp import f_interp_vars
from load import load, load_pos, load_spedas, load_sw

TRANGE = ["20251204", "20251205"]
TRANGE_NEXT_DAY = ["20251204", "20251206"]

_PATHS = resolve_paths()
PQT_BASE = Path(_PATHS["pqt_base"])
SPEDAS_BASE = Path(_PATHS["spedas_base"])


def _require(path: Path) -> Path:
    if not Path(path).exists():
        pytest.skip(f"local data missing: {path}")
    return Path(path)


# ---------------------------------------------------------------------------
# GENA: product loading and inclusive-midnight semantics
# ---------------------------------------------------------------------------
def test_gena_two_date_endpoints_cover_only_the_first_day():
    _require(PQT_BASE / "20251204_l2s_a.pqt")
    frame = load("gena", TRANGE, "l2s_a")
    assert len(frame) > 10_000
    assert frame["Datetime"].min().date() == pd.Timestamp("2025-12-04").date()
    # Time("20251205") is midnight, so the inclusive end holds one instant.
    assert frame["Datetime"].max() < pd.Timestamp("2025-12-05")


def test_gena_three_endpoints_cover_both_days():
    _require(PQT_BASE / "20251205_l2s_a.pqt")
    frame = load("gena", TRANGE_NEXT_DAY, "l2s_a")
    assert frame["Datetime"].max() >= pd.Timestamp("2025-12-05")
    assert len(frame) > len(load("gena", TRANGE, "l2s_a"))


def test_gena_frames_carry_provenance():
    _require(PQT_BASE / "20251204_l2s_a.pqt")
    frame = load("gena", TRANGE, "l2s_a")
    assert frame.attrs["source"] == "gena"
    assert frame.attrs["product"] == "l2s_a"
    assert frame.attrs["level"] == "l2s"
    assert frame.attrs["probe"] == "a"
    assert frame.__class__.__name__ == "ParticleDataFrame"
    assert {"Energy", "Tof", "Spec", "X1", "Y1", "MJ2000_Star_X"} <= set(frame.columns)


# ---------------------------------------------------------------------------
# THEMIS position: state product now loads from the local CDF
# ---------------------------------------------------------------------------
def _require_state_file() -> Path:
    matches = sorted(
        SPEDAS_BASE.glob("themis_data/thb/l1/state/2025/thb_l1_state_20251204.cdf")
    )
    if not matches:
        pytest.skip("THEMIS state CDF for 2025-12-04 is missing")
    return matches[0]


def test_themis_state_product_loads_and_positions_work():
    _require_state_file()
    frame = load("themis", TRANGE, "state", probe="b", download=False)
    assert len(frame) > 0
    assert "thb_pos_gsm_x" in frame.columns
    assert frame.attrs["level"] == "l1"

    times = ["2025-12-04T00:00:00", "2025-12-04T00:05:00"]
    positions = load_pos("themis", times, "b")
    assert positions.shape == (2, 3)
    assert np.isfinite(positions).all()
    # thb_l1_state positions are kilometres in GSM; the Moon is ~3.8e5 km away.
    assert np.all(np.abs(positions) > 1_000)


def test_themis_esa_loads_offline():
    _require(SPEDAS_BASE / "themis_data/thb/l2/esa/2025/thb_l2_esa_20251204_v01.cdf")
    frame = load("themis", TRANGE, "esa", probe="b", download=False)
    assert len(frame) > 1_000
    assert frame.attrs["source"] == "themis"
    assert frame.attrs["probe"] == "b"
    assert frame.attrs["level"] == "l2"
    assert frame["Datetime"].max() <= pd.Timestamp("2025-12-05")


# ---------------------------------------------------------------------------
# MMS: a missing local result is labelled instead of looking successful
# ---------------------------------------------------------------------------
def _mms_december_files():
    return list(SPEDAS_BASE.glob("mms_data/mms1/fpi/**/mms1_fpi_*_2025120*.cdf"))


def test_mms_fpi_without_local_files_is_labelled():
    if _mms_december_files():
        pytest.skip("local MMS FPI files exist for this window; behavior differs")
    frame = load("mms", TRANGE, "fpi", download=False)
    assert len(frame) == 0
    assert list(frame.columns) == ["Datetime", "Timestamp"]
    assert frame.attrs["product"] == "fpi"
    assert frame.attrs["level"] == "l2"
    assert "empty_reason" in frame.attrs
    assert "mms fpi" in frame.attrs["empty_reason"]


def test_mms_energy_spectrum_request_is_labelled_empty():
    if _mms_december_files():
        pytest.skip("local MMS FPI files exist for this window; behavior differs")
    frame = load(
        "mms", TRANGE, "fpi",
        variables=["energyspectr_omni_fast"], download=False,
    )
    assert len(frame) == 0
    assert "energyspectr_omni_fast" not in frame.columns
    assert "empty_reason" in frame.attrs


# ---------------------------------------------------------------------------
# OMNI / indices: local products load
# ---------------------------------------------------------------------------
def test_omni_loads_from_local_cdf():
    _require(SPEDAS_BASE / "omni_data/hro_1min/2025/omni_hro_1min_20251201_v01.cdf")
    frame = load("omni", TRANGE, download=False)
    assert len(frame) > 1_400
    assert list(frame.columns) == [
        "Datetime", "Timestamp", "BZ_GSM", "BY_GSM",
        "flow_speed", "proton_density", "Pressure",
    ]
    assert frame.attrs["product"] == "hro"


def test_dst_loads_from_local_cache():
    _require(SPEDAS_BASE / "indices/kyoto/hour_dst/202512.json")
    frame = load("indices", TRANGE, "dst", download=False)
    assert 20 <= len(frame) <= 24
    stamps = pd.DatetimeIndex(frame["Datetime"])
    assert (stamps.minute == 29).all()
    assert (stamps.second == 30).all()


# ---------------------------------------------------------------------------
# load_sw: a single unusable OMNI variable no longer resets the others
# ---------------------------------------------------------------------------
def test_load_sw_keeps_real_omni_values():
    _require(SPEDAS_BASE / "omni_data/hro_1min/2025/omni_hro_1min_20251201_v01.cdf")
    _require(SPEDAS_BASE / "indices/kyoto/hour_dst/202512.json")
    values = load_sw(
        ["2025-12-04T00:00:00", "2025-12-04T06:00:00"],
        download=False,
    )
    fallback = np.array([1.0, 0.0, 0.0, 1.0, 5.0, 0.0])
    actual = np.array([
        values["Pressure"][0], values["BY_GSM"][0], values["BZ_GSM"][0],
        values["Beta"][0], values["Mgs_mach_num"][0], values["SYM_H"][0],
    ])
    assert not np.allclose(actual, fallback), (
        "every OMNI variable was replaced by its fallback constant"
    )
    assert np.isfinite(values["dst"]).all()
    assert abs(values["dst"][0]) > 1.0


def test_all_nan_sym_h_still_breaks_plain_interpolation():
    """Diagnostic: load_sw now isolates this failure per variable."""
    _require(SPEDAS_BASE / "omni_data/hro_1min/2025/omni_hro_1min_20251201_v01.cdf")
    data = load_spedas(
        "omni", ["2025-12-03T23:00:00", "2025-12-04T07:00:00"],
        product="hro", variables=["SYM_H", "Pressure"],
        return_type="dict", download=False,
    )
    assert np.isfinite(data["Pressure"]["y"]).any()
    assert np.isnan(data["SYM_H"]["y"]).all()
    with pytest.raises(ValueError):
        f_interp_vars(data["SYM_H"])


# ---------------------------------------------------------------------------
# Public surface
# ---------------------------------------------------------------------------
def test_load_loc_is_not_a_public_name():
    assert not hasattr(load_pkg, "load_loc")
    assert hasattr(load_pkg, "load_pos")
