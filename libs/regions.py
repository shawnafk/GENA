#!/usr/bin/env python3
"""Cross-mission plasma-region and BBF/PSBL high-speed-flow classifiers.

The implementation is intentionally conservative for the Queqiao-2 / lunar
distance BBF study: broad plasma-region context is separated from CPS BBF
identification, and PSBL beams are reported as candidates unless directional
ion information is available.
"""

from __future__ import annotations

from collections.abc import Callable

import numpy as np
import pandas as pd

from libs import utils as cu
from common.frame import TimeSeriesDataFrame


# Region codes are shared by classifier output and plotting. Keeping this map
# fixed makes the same region use the same color in every event and probe.
REGION_LABELS = {
    0: "Unclassified",
    1: "Solar wind",
    2: "Magnetosheath",
    3: "Mantle",
    4: "Lobe",
    5: "PSBL",
    6: "Plasma sheet",
}



def _require_timeseries(frame, name):
    if not isinstance(frame, TimeSeriesDataFrame):
        raise TypeError(f"{name} must be a TimeSeriesDataFrame")


def _find_column(frame: pd.DataFrame, *suffixes: str) -> str:
    matches = [
        str(column)
        for column in frame.columns
        if any(str(column).endswith(suffix) for suffix in suffixes)
    ]
    if len(matches) != 1:
        raise KeyError(f"Expected one column ending in {suffixes}, found {matches}")
    return matches[0]


def _vector_columns(frame: pd.DataFrame, *root_suffixes: str) -> list[str]:
    x_column = _find_column(frame, *(f"{suffix}_x" for suffix in root_suffixes))
    root = x_column[:-2]
    columns = [f"{root}_{axis}" for axis in ("x", "y", "z")]
    missing = [column for column in columns if column not in frame]
    if missing:
        raise KeyError(f"Missing vector columns: {missing}")
    return columns


def _column_frame(frame: TimeSeriesDataFrame, columns) -> TimeSeriesDataFrame:
    """Copy selected science columns while preserving time and provenance."""
    data = {"Datetime": frame["Datetime"], "Timestamp": frame["Timestamp"]}
    data.update({column: frame[column].to_numpy(copy=False) for column in columns})
    result = TimeSeriesDataFrame(data)
    result.attrs = dict(frame.attrs)
    return result


def _resample_frames(*frames, tres="4s", max_gap="15s"):
    """Synchronize mission frames to one regular grid over common coverage."""
    if not frames:
        return ()
    for index, frame in enumerate(frames):
        _require_timeseries(frame, f"frames[{index}]")
        if frame.empty:
            raise ValueError("Cannot classify an empty time series")

    starts = [frame["Datetime"].min() for frame in frames]
    ends = [frame["Datetime"].max() for frame in frames]
    start, end = max(starts), min(ends)
    if start > end:
        raise ValueError("Time series do not have overlapping coverage")

    reference = frames[0].resample(tres=tres, max_gap=max_gap)
    dates = reference["Datetime"]
    reference = reference.select((dates >= start) & (dates <= end))
    if reference.empty:
        raise ValueError("No synchronized samples fall inside common coverage")
    return (reference, *(frame.resample_to(reference, max_gap=max_gap) for frame in frames[1:]))


def _evaluate_position(pos_gsm, timestamps):
    position = np.asarray(pos_gsm(timestamps), dtype=float)
    if position.shape == (3, len(timestamps)):
        position = position.T
    if position.shape != (len(timestamps), 3):
        raise ValueError(
            f"pos_gsm must return shape ({len(timestamps)}, 3), got {position.shape}"
        )
    return position





def classify_liuzzo(
    esa: TimeSeriesDataFrame,
    fgm: TimeSeriesDataFrame,
    pos_gsm: Callable[[np.ndarray], np.ndarray],
    *,
    tres="4s",
    max_gap="15s",
    columns=None,
) -> TimeSeriesDataFrame:
    """Classify Liuzzo lobe, plasma-sheet, and magnetosheath-like regions.

    ``columns`` maps ``density`` and ``magnetic`` to an exact scalar column and
    three XYZ column names, respectively. Without it, THEMIS names are resolved.
    Supply density in cm^-3, magnetic field in GSM nT, and GSM position in R_E.
    Column mapping enables other missions; it does not validate scientific
    applicability of the original thresholds to their observation conditions.
    """

    _require_timeseries(esa, "esa")
    _require_timeseries(fgm, "fgm")
    if columns is None:
        n_key = _find_column(esa, "_peir_density", "_peif_density")
        b_columns = _vector_columns(fgm, "_fgs_gsm", "_fgl_gsm")
    else:
        n_key = columns["density"]
        b_columns = list(columns["magnetic"])
        if len(b_columns) != 3:
            raise ValueError("magnetic must contain three GSM XYZ column names")
    esa, fgm = _resample_frames(
        _column_frame(esa, [n_key]),
        _column_frame(fgm, b_columns),
        tres=tres,
        max_gap=max_gap,
    )
    time = esa["Datetime"]
    t = cu.Time(time).to_timestamp()
    pos = _evaluate_position(pos_gsm, t)
    b = fgm[b_columns].to_numpy(dtype=float)
    bmag = np.linalg.norm(b, axis=1)
    phase = np.degrees(np.arctan2(pos[:, 1], pos[:, 0])) % 360.0
    near_midnight = np.abs(phase - 180.0) <= 20.0

    bavg = {}
    for minutes in (30, 60, 90, 120):
        avg = (
            pd.DataFrame(b, index=time)
            .rolling(f"{minutes}min", center=True)
            .mean()
            .to_numpy()
        )
        bavg[minutes] = (avg, np.linalg.norm(avg, axis=1))

    lobe = near_midnight.copy()
    for minutes, threshold in ((30, 0.93), (60, 0.92), (90, 0.91), (120, 0.90)):
        avg, avgmag = bavg[minutes]
        lobe &= (np.abs(avg[:, 0]) > threshold * avgmag) & (np.abs(bmag - avgmag) <= 0.05 * avgmag)
    plasma_sheet = near_midnight & (bavg[30][1] < 0.75 * bavg[120][1])
    magnetosheath = ~near_midnight & (esa[n_key].to_numpy(dtype=float) >= 1.0)
    region = np.zeros(len(time), dtype=int)
    region[lobe], region[plasma_sheet], region[magnetosheath] = 4, 6, 2
    result = TimeSeriesDataFrame(
        {
            "Datetime": time, "Timestamp": pd.DatetimeIndex(time).as_unit("ns").asi8.astype(float) / 1e9,
            "region_code": region,
            "lobe": lobe,
            "plasma_sheet": plasma_sheet,
            "magnetosheath": magnetosheath,
        }
    )
    result.attrs = {
        "source": esa.attrs.get("source"),
        "classifier": "Liuzzo 2022",
        "probe": esa.attrs.get("probe"),
    }
    return result


def classify_kiehas(
    esa: TimeSeriesDataFrame, *, tres="4s", max_gap="15s", columns=None
) -> TimeSeriesDataFrame:
    """Kiehas 2018 plasma-sheet gate: Ti > 300 eV and n < 0.5 cm^-3.

    Pass ``columns={"density": name, "temperature": name}`` for any mission's
    ion moments. Temperature must be scalar ion temperature in eV, not just
    a parallel or perpendicular component. Default resolution uses THEMIS names.
    """

    _require_timeseries(esa, "esa")
    if columns is None:
        n_key = _find_column(esa, "_peir_density", "_peif_density")
        t_key = n_key.replace("_density", "_avgtemp")
    else:
        n_key, t_key = columns["density"], columns["temperature"]
    if t_key not in esa:
        raise KeyError(f"Missing ion-temperature column {t_key!r}")
    esa, = _resample_frames(
        _column_frame(esa, [n_key, t_key]), tres=tres, max_gap=max_gap
    )
    time = esa["Datetime"]
    n = esa[n_key].to_numpy(dtype=float)
    ti = esa[t_key].to_numpy(dtype=float)
    plasma_sheet = np.isfinite(n) & np.isfinite(ti) & (ti > 300.0) & (n < 0.5)
    result = TimeSeriesDataFrame({"Datetime": time, "Timestamp": pd.DatetimeIndex(time).as_unit("ns").asi8.astype(float) / 1e9, "plasma_sheet": plasma_sheet})
    result.attrs = {
        "source": esa.attrs.get("source"),
        "classifier": "Kiehas 2018",
        "probe": esa.attrs.get("probe"),
    }
    return result


def classify_runov(
    esa: TimeSeriesDataFrame,
    sst: TimeSeriesDataFrame,
    pos_gsm: Callable[[np.ndarray], np.ndarray],
    *,
    tres="4s",
    max_gap="15s",
    columns=None,
) -> TimeSeriesDataFrame:
    """Runov sample gate; ``pos_gsm`` evaluates Unix-second timestamps.

    Explicit ``columns`` maps ion_density, ion_temperature, electron_density,
    electron_temperature in the first frame, plus hot_density, hot_temperature
    in the second frame. Densities are cm^-3, temperatures eV, positions GSM R_E.
    All six physically comparable moments are required; mapping variable names
    alone does not make another instrument's energy coverage equivalent to SST.
    Default resolution uses THEMIS ESA/SST names.
    """

    _require_timeseries(esa, "esa")
    _require_timeseries(sst, "sst")
    if columns is None:
        ni_key = _find_column(esa, "_peir_density", "_peif_density")
        prefix = ni_key.split("_", 1)[0]
        ti_key = ni_key.replace("_density", "_avgtemp")
        ne_key, te_key = f"{prefix}_peer_density", f"{prefix}_peer_avgtemp"
        ns_key, ts_key = f"{prefix}_psif_density", f"{prefix}_psif_avgtemp"
    else:
        ni_key, ti_key = columns["ion_density"], columns["ion_temperature"]
        ne_key, te_key = columns["electron_density"], columns["electron_temperature"]
        ns_key, ts_key = columns["hot_density"], columns["hot_temperature"]
    missing = [key for key in (ti_key, ne_key, te_key) if key not in esa]
    missing += [key for key in (ns_key, ts_key) if key not in sst]
    if missing:
        raise KeyError(f"Missing Runov input columns: {missing}")
    esa, sst = _resample_frames(
        _column_frame(esa, [ni_key, ti_key, ne_key, te_key]),
        _column_frame(sst, [ns_key, ts_key]),
        tres=tres,
        max_gap=max_gap,
    )
    time = esa["Datetime"]
    t = cu.Time(time).to_timestamp()
    pos = _evaluate_position(pos_gsm, t)
    ni = esa[ni_key].to_numpy(dtype=float)
    ti = esa[ti_key].to_numpy(dtype=float)
    ne = esa[ne_key].to_numpy(dtype=float)
    te = esa[te_key].to_numpy(dtype=float)
    ns = sst[ns_key].to_numpy(dtype=float)
    ts = sst[ts_key].to_numpy(dtype=float)
    ti_esa_sst = np.divide(
        ni * ti + ns * ts,
        ni + ns,
        out=np.full_like(ni, np.nan),
        where=np.isfinite(ni + ns) & ((ni + ns) != 0),
    )
    ti_ratio = np.divide(
        ti_esa_sst,
        ti,
        out=np.full_like(ti_esa_sst, np.nan),
        where=np.isfinite(ti) & (ti != 0),
    )
    pressure_gate = (ni > 0.05) & (ni * ti > 10.0) & (ne * te > 1.0)
    y_gate = np.abs(pos[:, 1]) < 20.0
    plasma_sheet = y_gate & pressure_gate & (ti_ratio < 1.25)
    result = TimeSeriesDataFrame(
        {
            "Datetime": time, "Timestamp": pd.DatetimeIndex(time).as_unit("ns").asi8.astype(float) / 1e9,
            "plasma_sheet": plasma_sheet,
            "y_agsm_gate": y_gate,
            "pressure_gate": pressure_gate,
            "ti_esa_sst_over_esa": ti_ratio,
        }
    )
    result.attrs = {
        "source": esa.attrs.get("source"),
        "classifier": "Runov 2023",
        "probe": esa.attrs.get("probe"),
    }
    return result


def classify_plasma_sheet_flows(
    df: pd.DataFrame,
    *,
    esa_l1: TimeSeriesDataFrame | None = None,
    esa_l2: TimeSeriesDataFrame | None = None,
) -> pd.DataFrame:
    """Classify CPS BBFs and PSBL high-speed-flow candidates inside plasma sheet.

    The function does not promote 1-D omni spectra to VDF evidence. If
    loaded THEMIS L1 and L2 ESA frames are supplied, Grigorenko-style
    support columns are added.
    """

    out = df.copy()
    plasma = out["plasma_sheet_like"].to_numpy(dtype=bool)
    beta_i = out["beta_i"].to_numpy(dtype=float)
    beta_xy = out["beta_xy"].to_numpy(dtype=float)
    b_nt = out["b_nt"].to_numpy(dtype=float)
    bx_nt = out["bx_nt"].to_numpy(dtype=float)
    vx = out["vx_kms"].to_numpy(dtype=float)
    v = out["v_kms"].to_numpy(dtype=float)
    vpar_abs = np.abs(out["v_parallel_kms"].to_numpy(dtype=float))
    vperp = out["v_perp_kms"].to_numpy(dtype=float)
    vperp_x_abs = np.abs(out["v_perp_x_kms"].to_numpy(dtype=float))

    bx_over_b = np.divide(np.abs(bx_nt), b_nt, out=np.full_like(b_nt, np.nan), where=b_nt > 0)
    central_or_outer_ps = (np.abs(bx_nt) <= 10.0) | (bx_over_b < 0.8)
    high_speed = (np.abs(vx) >= 400.0) | (v >= 400.0)
    convective_fast = plasma & (vperp_x_abs >= 250.0) & (vperp >= 250.0)
    convective_bbf = convective_fast & (np.abs(vx) >= 400.0) & (vperp >= vpar_abs)

    raj_like_bbf = plasma & (np.abs(vx) >= 400.0) & (vperp >= 250.0) & (beta_xy >= 0.5)
    cao_cps_bbf = (
        plasma
        & convective_bbf
        & (beta_i >= 0.5)
        & central_or_outer_ps
        & out["field_direction_reliable"].to_numpy(dtype=bool)
    )
    cps_bbf = cao_cps_bbf & (vperp >= vpar_abs)

    psbl_candidate = plasma & high_speed & ~cps_bbf
    field_aligned_candidate = psbl_candidate & (
        (vpar_abs >= 400.0) | (vpar_abs > vperp) | (beta_i < 0.5) | (bx_over_b >= 0.8)
    )
    ambiguous_inner_outer = plasma & ~cps_bbf & ~field_aligned_candidate & (esa_l1 is None)

    out["fast_flow_sample"] = high_speed
    out["convective_fast_flow_sample"] = convective_fast
    out["convective_bbf_sample"] = convective_bbf
    out["raj_like_bbf_sample"] = raj_like_bbf
    out["cao_cps_bbf_sample"] = cao_cps_bbf
    out["bbf_sample"] = cps_bbf
    out["cps_bbf_sample"] = cps_bbf
    out["parallel_fast_psbl_sample"] = field_aligned_candidate
    out["psbl_high_speed_flow_candidate"] = psbl_candidate
    out["ambiguous_inner_psbl_outer_ps"] = ambiguous_inner_outer

    flow_class = np.full(len(out), "non_plasma_sheet_or_quiet", dtype=object)
    flow_class[plasma & ~high_speed] = "quiet_plasma_sheet_like"
    flow_class[psbl_candidate] = "PSBL_high_speed_flow_candidate"
    flow_class[cps_bbf] = "CPS_BBF"
    out["flow_class"] = flow_class

    out["vdf_check_available"] = False
    out["grigorenko_field_aligned_ratio_gt2kev"] = np.nan
    out["grigorenko_vhigh_parallel_x_proxy_kms"] = np.nan
    out["grigorenko_nhigh_proxy"] = np.nan
    out["grigorenko_psbl_beam_support"] = False

    if esa_l1 is not None and esa_l2 is not None:
        out = _attach_grigorenko_reduced_esa_check(out, esa_l1, esa_l2)

    out["bbf_event"] = out["cps_bbf_sample"]
    return out


def _attach_grigorenko_reduced_esa_check(
    df: pd.DataFrame,
    esa_l1: TimeSeriesDataFrame,
    esa_l2: TimeSeriesDataFrame,
) -> pd.DataFrame:
    """Attach Grigorenko metrics from already-loaded THEMIS ESA data."""

    if df.empty:
        return df

    _require_timeseries(esa_l1, "esa_l1")
    _require_timeseries(esa_l2, "esa_l2")

    metrics = _grigorenko_metrics(esa_l1, esa_l2)
    if len(metrics["time"]) == 0:
        return df
    target = cu.Time(df["time"]).to_timestamp()
    nearest = np.searchsorted(metrics["time"], target)
    nearest = np.clip(nearest, 1, len(metrics["time"]) - 1)
    left = nearest - 1
    nearest = np.where(
        np.abs(metrics["time"][left] - target) <= np.abs(metrics["time"][nearest] - target),
        left,
        nearest,
    )
    available = np.abs(metrics["time"][nearest] - target) <= 2.0
    ratio = np.where(available, metrics["ratio"][nearest], np.nan)
    vhigh_x = np.where(available, metrics["vhigh_parallel_x"][nearest], np.nan)
    nhigh = np.where(available, metrics["nhigh_proxy"][nearest], np.nan)

    df = df.copy()
    df["vdf_check_available"] = available
    df["grigorenko_field_aligned_ratio_gt2kev"] = ratio
    df["grigorenko_vhigh_parallel_x_proxy_kms"] = vhigh_x
    df["grigorenko_nhigh_proxy"] = nhigh
    df["grigorenko_psbl_beam_support"] = (
        available
        & df["plasma_sheet_like"].to_numpy(dtype=bool)
        & (nhigh >= 0.01)
        & ((ratio >= 3.0) | (np.abs(vhigh_x) >= 600.0))
    )
    return df


def _grigorenko_metrics(
    esa_l1: TimeSeriesDataFrame, esa_l2: TimeSeriesDataFrame
):
    """Calculate THEMIS-only reduced-ESA metrics from canonical mission frames."""
    from load import esa_mode2_angles, themis_decompress
    _require_timeseries(esa_l1, "esa_l1")
    _require_timeseries(esa_l2, "esa_l2")
    l1_key = _find_column(esa_l1, "_eir_050x24")
    prefix = l1_key.split("_", 1)[0]
    eflux_key = f"{prefix}_peir_en_eflux"
    energy_key = f"{eflux_key}_energy"
    magnetic_columns = _vector_columns(esa_l2, "_peir_magf")
    missing = [key for key in (eflux_key, energy_key) if key not in esa_l2]
    if missing:
        raise KeyError(f"Missing Grigorenko input columns: {missing}")

    l1_valid = esa_l1[l1_key].map(
        lambda value: isinstance(value, (np.ndarray, list, tuple))
    )
    l1 = esa_l1.select(l1_valid)
    l2_valid = esa_l2[eflux_key].map(
        lambda value: isinstance(value, (np.ndarray, list, tuple))
    ) & esa_l2[energy_key].map(
        lambda value: isinstance(value, (np.ndarray, list, tuple))
    )
    l2_reference = esa_l2.select(l2_valid)
    if l1.empty or l2_reference.empty:
        empty = np.array([], dtype=float)
        return {
            "time": empty,
            "ratio": empty,
            "vhigh_parallel_x": empty,
            "nhigh_proxy": empty,
        }
    l2 = esa_l2.resample_to(l2_reference, max_gap="10s")

    l1_time = cu.Time(l1["Datetime"]).to_timestamp()
    l2_time = cu.Time(l2["Datetime"]).to_timestamp()
    compressed = np.stack(l1[l1_key].to_numpy())
    energy = np.stack(l2[energy_key].to_numpy()).astype(float)
    omni = np.stack(l2[eflux_key].to_numpy()).astype(float)
    magnetic = l2[magnetic_columns].to_numpy(dtype=float)

    look, solid_angle = esa_mode2_angles()
    counts_all = themis_decompress(compressed).transpose(0, 2, 1)  # time, energy, angle
    nearest_l2 = np.searchsorted(l2_time, l1_time)
    nearest_l2 = np.clip(nearest_l2, 1, len(l2_time) - 1)
    left = nearest_l2 - 1
    choose_left = np.abs(l2_time[left] - l1_time) <= np.abs(l2_time[nearest_l2] - l1_time)
    nearest_l2 = np.where(choose_left, left, nearest_l2)
    matched = np.abs(l2_time[nearest_l2] - l1_time) <= 2.0

    ratio = np.full(len(l1_time), np.nan)
    vhigh_x = np.full(len(l1_time), np.nan)
    nhigh = np.full(len(l1_time), np.nan)
    for i in np.where(matched)[0]:
        j = nearest_l2[i]
        energies = energy[j, :24]
        if not np.isfinite(energies).any():
            continue
        high = energies >= 2000.0
        if not np.any(high):
            continue
        counts = counts_all[i].astype(float)
        mean_counts = np.sum(counts * solid_angle, axis=1) / np.sum(solid_angle, axis=1)
        directional = np.divide(
            omni[j, :24, None] * counts,
            mean_counts[:, None],
            out=np.full_like(counts, np.nan, dtype=float),
            where=mean_counts[:, None] > 0,
        )
        b_norm = np.linalg.norm(magnetic[j])
        if not np.isfinite(b_norm) or b_norm <= 0:
            continue
        b_hat = magnetic[j] / b_norm
        cos_pa = np.einsum("eai,i->ea", look, b_hat)
        fa = np.abs(cos_pa) >= np.cos(np.deg2rad(30.0))
        perp = np.abs(cos_pa) <= np.sin(np.deg2rad(30.0))
        field_flux = np.nanmean(np.where(high[:, None] & fa, directional, np.nan))
        perp_flux = np.nanmean(np.where(high[:, None] & perp, directional, np.nan))
        ratio[i] = field_flux / perp_flux if np.isfinite(perp_flux) and perp_flux > 0 else np.nan
        speed = np.sqrt(2 * energies * 1.602176634e-19 / 1.67262192369e-27) / 1e3
        v_parallel = speed[:, None] * cos_pa
        weighted_vpar = np.nanmean(np.where(high[:, None] & fa, directional * v_parallel, np.nan))
        weighted_flux = np.nanmean(np.where(high[:, None] & fa, directional, np.nan))
        vhigh_x[i] = weighted_vpar * b_hat[0] / weighted_flux if np.isfinite(weighted_flux) and weighted_flux > 0 else np.nan
        nhigh[i] = np.nanmean(np.where(high[:, None], directional, np.nan)) / 1e5
    return {"time": l1_time, "ratio": ratio, "vhigh_parallel_x": vhigh_x, "nhigh_proxy": nhigh}
