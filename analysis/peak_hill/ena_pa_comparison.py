#!/usr/bin/env python3
"""Compare GENA ENA counts with THEMIS-B ion pitch-angle structure.

The script answers one question: between the 2025-11-03/04 ENA peak and the
2025-11-04/05 ENA "hill", does the ion pitch-angle distribution at THEMIS-B
become more field-aligned, and can that concentration explain the ENA count
contrast when the omnidirectional ion spectra are similar?

Products
--------
* GENA-A ``l2f_a`` event counts and proxy ENA energy spectra.
* THEMIS-B ESA ``peir`` ion spectra, SST ``psif`` high-energy ion spectra.
* THEMIS-B ESD ``peif`` three-dimensional ion distributions and FGM.

The ESD pitch-angle reduction follows the project's validated recipe: ESD
``look`` vectors are already particle-flight directions (cross-checked against
``peif_velocity_dsl``), and the matching ESA ``peif_magf`` gives B in DSL.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import urllib.request
import warnings

import matplotlib

matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
import pandas as pd
from cdflib import CDF

from common.config import get_configured_paths
from load import load
from gena.post import statistic as stat
from gena.post.project import event_vectors
from libs.plasma import sigma_HH

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parent
DATA_DIR = ROOT / "data"
OUT_DIR = ROOT / "out"
PROBE = "b"
PREFIX = f"th{PROBE}"
GENA_PRODUCT = "l2f_a"

TRANGE = ("2025-11-03T15:00:00", "2025-11-06T00:00:00")
#: The two intervals supplied by the user.  They overlap on 2025-11-04
#: 09:00--23:00, so two tighter sub-windows isolate the ENA peak and hill.
WINDOWS = {
    "Peak (T1)": ("2025-11-03T15:00:00", "2025-11-04T23:00:00"),
    "Hill (T2)": ("2025-11-04T09:00:00", "2025-11-05T23:00:00"),
}
FOCUS = {
    "Peak burst": ("2025-11-03T19:00:00", "2025-11-04T00:00:00"),
    "Hill core": ("2025-11-04T09:00:00", "2025-11-04T23:00:00"),
}

ENA_EDGES = np.logspace(np.log10(0.5), np.log10(100.0), 16)
ENERGY_EV = (4_000.0, 24_600.0)          # ESA ion overlap band for PA work
PA_EDGES = np.arange(0.0, 181.0, 15.0)
PA_CENTERS = (PA_EDGES[:-1] + PA_EDGES[1:]) / 2.0
COUNT_FOV_DEG = (20.0, 10.0)             # GENA effective FOV used by the project
ESD_URL = "https://themis.ssl.berkeley.edu/data/themis/{p}/l2/esd/{y}/{p}_l2_esa_peif_{d}_v01.cdf"


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------
def _nearest(source, target):
    source = np.asarray(source, dtype=float)
    target = np.asarray(target, dtype=float)
    right = np.clip(np.searchsorted(target, source), 1, len(target) - 1)
    left = right - 1
    index = np.where(np.abs(source - target[left]) <= np.abs(source - target[right]), left, right)
    return index, np.abs(source - target[index])


def _stack_array_column(frame, column, width):
    """Stack an object array column into a fixed-width matrix with NaN padding."""
    stacked = np.full((len(frame), width), np.nan)
    mask = np.zeros(len(frame), dtype=bool)
    for row, item in enumerate(frame[column].to_numpy()):
        value = np.atleast_1d(np.asarray(item, dtype=float))
        if value.size == width:
            stacked[row] = value
            mask[row] = True
    return stacked, mask


def _window(frame, start, end):
    return frame.select(
        (frame["Datetime"] >= pd.Timestamp(start)) & (frame["Datetime"] <= pd.Timestamp(end))
    )


def ensure_esd(day: pd.Timestamp) -> Path:
    """Return the local ESD CDF path, downloading it when absent."""
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    path = DATA_DIR / f"{PREFIX}_l2_esa_peif_{day:%Y%m%d}_v01.cdf"
    if path.exists():
        return path
    url = ESD_URL.format(p=PREFIX, y=f"{day:%Y}", d=f"{day:%Y%m%d}")
    print(f"downloading {url}")
    urllib.request.urlretrieve(url, path)
    return path


# ---------------------------------------------------------------------------
# GENA
# ---------------------------------------------------------------------------
def load_ena():
    columns = [
        "Timestamp", "X1", "X2", "Y1", "Y2", "Tof", "Energy", "Spec",
        "MJ2000_q1", "MJ2000_q2", "MJ2000_q3",
        "MJ2000_Star_X", "MJ2000_Star_Y", "MJ2000_Star_Z",
    ]
    return load("gena", TRANGE, GENA_PRODUCT, columns=columns)


def ena_direction_angle(events, fgm):
    """Angle between the GENA particle-flight direction and THEMIS-B B (GSM)."""
    quaternion = events[["MJ2000_q1", "MJ2000_q2", "MJ2000_q3"]].to_numpy(dtype=float)
    valid = np.all(np.isfinite(quaternion), axis=1) & (np.sum(quaternion**2, axis=1) <= 1.0 + 1e-8)
    valid &= np.all(np.isfinite(events[["X1", "X2", "Y1", "Y2", "Tof"]].to_numpy(dtype=float)), axis=1)
    vectors = event_vectors(events.loc[valid].reset_index(drop=True), coord="GSM", probe="A")
    direction = vectors[["DX", "DY", "DZ"]].to_numpy(dtype=float)
    field = fgm[[f"{PREFIX}_fgs_gsm_x", f"{PREFIX}_fgs_gsm_y", f"{PREFIX}_fgs_gsm_z"]].to_numpy(dtype=float)
    field_time = fgm["Timestamp"].to_numpy(dtype=float)
    index, _ = _nearest(vectors["Timestamp"].to_numpy(dtype=float), field_time)
    local = field[index]
    norm = np.linalg.norm(local, axis=1)
    cosine = np.einsum("ij,ij->i", direction, local) / norm
    vectors["LOS_PA"] = np.rad2deg(np.arccos(np.clip(cosine, -1.0, 1.0)))
    vectors["COS"] = cosine
    return vectors


def ena_metrics(events):
    rows = {}
    per_minute = events.set_index("Datetime").resample("1min").size()
    for name, (start, end) in {**WINDOWS, **FOCUS}.items():
        selection = _window(events, start, end)
        seconds = (pd.Timestamp(end) - pd.Timestamp(start)).total_seconds()
        energy, flux, _ = stat.flux_e(selection, ENA_EDGES, exposure=seconds)
        rows[name] = {
            "start": start, "end": end, "seconds": seconds,
            "events": int(len(selection)),
            "counts_per_min": float(len(selection) / seconds * 60.0),
            "energy_kev": np.asarray(energy, dtype=float),
            "flux": np.asarray(flux, dtype=float),
            "band_rate": {
                f"{low:g}-{high:g} keV": float(len(selection[selection["Energy"].between(low, high)]) / seconds)
                for low, high in ((2.0, 10.0), (10.0, 30.0))
            },
        }
    return rows, per_minute


# ---------------------------------------------------------------------------
# THEMIS ESA / SST / ESD
# ---------------------------------------------------------------------------
def load_esa():
    esa = load("themis", TRANGE, "esa", probe=PROBE, download=False)
    fgm = load("themis", TRANGE, "fgm", probe=PROBE, variables=["fgs_gsm"], download=False)
    energy, keep = _stack_array_column(esa, f"{PREFIX}_peir_en_eflux_energy", 32)
    eflux, keep_flux = _stack_array_column(esa, f"{PREFIX}_peir_en_eflux", 32)
    keep = keep & keep_flux
    eflux = np.where(eflux <= -1e30, np.nan, eflux)
    energy = np.where(energy <= -1e30, np.nan, energy)
    times = pd.DatetimeIndex(esa["Datetime"])[keep]
    return times, energy[keep], eflux[keep], esa, fgm


def load_sst():
    sst = load("themis", TRANGE, "sst", probe=PROBE, download=False)
    energy, keep = _stack_array_column(sst, f"{PREFIX}_psif_en_eflux_energy", 16)
    eflux, keep_flux = _stack_array_column(sst, f"{PREFIX}_psif_en_eflux", 16)
    keep = keep & keep_flux
    eflux = np.where(eflux <= -1e30, np.nan, eflux)
    energy = np.where(energy <= -1e30, np.nan, energy)
    return pd.DatetimeIndex(sst["Datetime"])[keep], energy[keep], eflux[keep], keep


def load_esd_records(days):
    """Read quality-zero ESD records with their matching DSL magnetic field."""
    records = []
    for day in days:
        esd_path = ensure_esd(pd.Timestamp(day))
        base = Path(get_configured_paths()["spedas_base"]) / "themis_data" / PREFIX / "l2" / "esa"
        esa_path = sorted((base / f"{day:%Y}").glob(f"{PREFIX}_l2_esa_{day:%Y%m%d}_v*.cdf"))
        if not esa_path:
            raise FileNotFoundError(f"THEMIS ESA L2 file is missing for {day:%Y-%m-%d}")
        with CDF(str(esd_path)) as cdf:
            arrays = {name: np.asarray(cdf.varget(name)) for name in
                      ("time", "valid", "data_quality", "bins", "eflux", "en_ind", "an_ind", "phi_offset")}
            tables = {name: np.asarray(cdf.varget(name), dtype=float) for name in
                      ("energy", "denergy", "theta", "phi", "domega")}
        with CDF(str(esa_path[-1])) as cdf:
            b_time = np.asarray(cdf.varget(f"{PREFIX}_peif_time"), dtype=float)
            b_dsl = np.asarray(cdf.varget(f"{PREFIX}_peif_magf"), dtype=float)
            velocity_dsl = np.asarray(cdf.varget(f"{PREFIX}_peif_velocity_dsl"), dtype=float)
        index, offset = _nearest(arrays["time"].astype(float), b_time)
        for record, stamp in enumerate(arrays["time"].astype(float)):
            if arrays["valid"][record] != 1 or arrays["data_quality"][record] != 0 or offset[record] > 10.0:
                continue
            field = b_dsl[index[record]]
            if not np.all(np.isfinite(field)) or np.linalg.norm(field) == 0:
                continue
            en_mode = int(np.clip(arrays["en_ind"][record], 0, tables["energy"].shape[1] - 1))
            an_mode = int(np.clip(arrays["an_ind"][record], 0, tables["theta"].shape[2] - 1))
            theta = np.deg2rad(tables["theta"][:, :, an_mode])
            phi = np.deg2rad(tables["phi"][:, :, an_mode] + arrays["phi_offset"][record])
            look = np.stack((np.cos(theta) * np.cos(phi), np.cos(theta) * np.sin(phi), np.sin(theta)), axis=-1)
            records.append({
                "time": pd.Timestamp(stamp, unit="s"),
                "eflux": arrays["eflux"][record].astype(float),
                "bins": arrays["bins"][record].astype(bool),
                "energy_ev": tables["energy"][:, en_mode],
                "denergy_ev": tables["denergy"][:, en_mode],
                "look": look,
                "domega": tables["domega"][:, :, an_mode],
                "b_dsl": field,
                "velocity_dsl": velocity_dsl[index[record]],
            })
    return records


def reduce_pa(record):
    """Reduce one 3-D ESD distribution to J(E, PA) and the omnidirectional J(E)."""
    energy = record["energy_ev"]
    particle_flux = record["eflux"] / energy[:, None]
    b_hat = record["b_dsl"] / np.linalg.norm(record["b_dsl"])
    pitch = np.rad2deg(np.arccos(np.clip(np.einsum("eai,i->ea", record["look"], b_hat), -1.0, 1.0)))
    pa_index = np.digitize(pitch, PA_EDGES) - 1
    usable = (
        record["bins"] & np.isfinite(particle_flux) & (particle_flux > 0)
        & np.isfinite(record["domega"]) & (record["domega"] > 0)
        & (pa_index >= 0) & (pa_index < len(PA_CENTERS))
    )
    energy_index, _ = np.nonzero(usable)
    numerator = np.zeros((len(energy), len(PA_CENTERS)))
    denominator = np.zeros_like(numerator)
    np.add.at(numerator, (energy_index, pa_index[usable]), (particle_flux * record["domega"])[usable])
    np.add.at(denominator, (energy_index, pa_index[usable]), record["domega"][usable])
    j_pa = np.divide(numerator, denominator, out=np.full_like(numerator, np.nan), where=denominator > 0)
    iso_num = np.sum(np.where(usable, particle_flux * record["domega"], 0.0), axis=1)
    iso_den = np.sum(np.where(usable, record["domega"], 0.0), axis=1)
    j_iso = np.divide(iso_num, iso_den, out=np.full(len(energy), np.nan), where=iso_den > 0)
    weight = np.where(usable, particle_flux * record["domega"], 0.0)
    flow = np.sum(record["look"] * weight[..., None], axis=(0, 1))
    velocity = record["velocity_dsl"]
    denominator = np.linalg.norm(flow) * np.linalg.norm(velocity)
    cosine = float(np.dot(flow, velocity) / denominator) if denominator > 0 else np.nan
    return {"j_pa": j_pa, "j_iso": j_iso, "pitch": pitch, "usable": usable,
            "particle_flux": particle_flux, "look_velocity_cosine": cosine}


def fov_grid(look_pitch_deg, full_width_deg=COUNT_FOV_DEG, step_deg=2.0):
    """Gnomonic grid for a rectangular FOV centred on a pitch angle."""
    nx = max(1, int(np.ceil(full_width_deg[0] / step_deg)))
    ny = max(1, int(np.ceil(full_width_deg[1] / step_deg)))
    x_edges = np.linspace(-full_width_deg[0] / 2, full_width_deg[0] / 2, nx + 1)
    y_edges = np.linspace(-full_width_deg[1] / 2, full_width_deg[1] / 2, ny + 1)
    x = np.deg2rad((x_edges[:-1] + x_edges[1:]) / 2)
    y = np.deg2rad((y_edges[:-1] + y_edges[1:]) / 2)
    xx, yy = np.meshgrid(x, y)
    alpha = np.deg2rad(look_pitch_deg)
    bore = np.array([np.sin(alpha), 0.0, np.cos(alpha)])
    ex = np.array([np.cos(alpha), 0.0, -np.sin(alpha)])
    ey = np.array([0.0, 1.0, 0.0])
    direction = bore + np.tan(xx)[..., None] * ex + np.tan(yy)[..., None] * ey
    direction /= np.linalg.norm(direction, axis=-1, keepdims=True)
    pitch = np.rad2deg(np.arccos(np.clip(direction[..., 2], -1.0, 1.0)))
    dx, dy = np.diff(np.deg2rad(x_edges))[0], np.diff(np.deg2rad(y_edges))[0]
    jacobian = (1 / np.cos(xx) ** 2) * (1 / np.cos(yy) ** 2) / (1 + np.tan(xx) ** 2 + np.tan(yy) ** 2) ** 1.5
    return pitch, jacobian * dx * dy


def directional_flux(record, reduced, look_pitch_deg):
    """Angle-averaged J(E) over the GENA FOV centred on ``look_pitch_deg``."""
    pitch, solid_angle = fov_grid(look_pitch_deg)
    band = (record["energy_ev"] >= ENERGY_EV[0]) & (record["energy_ev"] <= ENERGY_EV[1])
    values = np.empty((np.count_nonzero(band),) + pitch.shape)
    for out_index, in_index in enumerate(np.flatnonzero(band)):
        valid = np.isfinite(reduced["j_pa"][in_index]) & (reduced["j_pa"][in_index] > 0)
        values[out_index] = (
            np.interp(pitch, PA_CENTERS[valid], reduced["j_pa"][in_index, valid])
            if valid.sum() >= 2 else reduced["j_iso"][in_index]
        )
    return np.nansum(values * solid_angle[None, :, :], axis=(1, 2)) / solid_angle.sum(), band


def predict_ena(energy_ev, j_fov, band):
    """Relative ENA production for a directional ion flux (arbitrary units)."""
    energy = energy_ev[band] / 1_000.0
    return float(np.nansum(j_fov * sigma_HH(energy)))


def _energy_edges(energy):
    """Logarithmic cell edges from channel centres."""
    energy = np.asarray(energy, dtype=float)
    edges = np.empty(len(energy) + 1)
    edges[1:-1] = np.sqrt(energy[:-1] * energy[1:])
    edges[0] = energy[0] ** 2 / edges[1]
    edges[-1] = energy[-1] ** 2 / edges[-2]
    return edges


def gain_spectrum(j_pa, j_iso, energy_ev, look_pitch_deg):
    """FOV-averaged directional flux divided by the omnidirectional flux, per energy."""
    band = (energy_ev >= ENERGY_EV[0]) & (energy_ev <= ENERGY_EV[1])
    pitch, solid_angle = fov_grid(look_pitch_deg)
    gain = np.full(len(energy_ev), np.nan)
    for index in np.flatnonzero(band):
        valid = np.isfinite(j_pa[index]) & (j_pa[index] > 0)
        if valid.sum() < 2:
            gain[index] = 1.0
            continue
        values = np.interp(pitch, PA_CENTERS[valid], j_pa[index, valid])
        gain[index] = np.nansum(values * solid_angle) / solid_angle.sum() / j_iso[index]
    return gain


def window_analysis(rows, peir_time, peir_energy, peir_eflux):
    """Per-window pitch-angle metrics, LOS gains and a relative ENA prediction."""
    out = {}
    alpha_grid = np.arange(0.0, 180.1, 5.0)
    selected = {**WINDOWS, **FOCUS}
    for name, (start, end) in selected.items():
        start, end = pd.Timestamp(start), pd.Timestamp(end)
        subset = [row for row in rows if start <= row["time"] <= end]
        if not subset:
            continue
        j_pa = np.nanmean([row["reduced"]["j_pa"] for row in subset], axis=0)
        j_iso = np.nanmean([row["reduced"]["j_iso"] for row in subset], axis=0)
        energy = subset[0]["record"]["energy_ev"]
        energy = energy[np.argsort(energy)]
        order = np.argsort(subset[0]["record"]["energy_ev"])
        j_pa = j_pa[order]
        j_iso = j_iso[order]

        mask = (peir_time >= start) & (peir_time <= end)
        peir_j = np.nanmedian(peir_eflux[mask] / peir_energy[mask], axis=0) if mask.any() else np.full(32, np.nan)
        peir_e = np.nanmedian(peir_energy[mask], axis=0) if mask.any() else np.full(32, np.nan)
        good = np.isfinite(peir_e) & np.isfinite(peir_j) & (peir_j > 0) & (peir_e > 0)
        order_e = np.argsort(peir_e[good])
        peir_j = np.exp(np.interp(np.log(energy), np.log(peir_e[good][order_e]),
                                  np.log(peir_j[good][order_e]), left=np.nan, right=np.nan))

        los_pa = np.array([row["los_pa"] for row in subset])
        los_median = float(np.nanmedian(los_pa))
        gain = np.array([gain_spectrum(j_pa, j_iso, energy, angle) for angle in alpha_grid])

        band = (energy >= ENERGY_EV[0]) & (energy <= ENERGY_EV[1])
        sigma = sigma_HH(energy[band] / 1_000.0)
        weight = np.nan_to_num(peir_j[band]) * sigma

        # Evaluate the window-mean gain curve at every record's own LOS angle so
        # the result follows the full LOS distribution rather than its median.
        def _per_record_gain(angles):
            gains = np.empty((len(angles), np.count_nonzero(band)))
            for index in range(len(angles)):
                gains[index] = [
                    np.interp(angles[index], alpha_grid, gain[:, channel])
                    for channel in np.flatnonzero(band)
                ]
            return np.nan_to_num(gains)

        direct = _per_record_gain(los_pa)
        mirror = _per_record_gain(180.0 - los_pa)
        predicted = float(np.nanmean(direct @ weight)) if len(direct) else np.nan
        predicted_mirror = float(np.nanmean(mirror @ weight)) if len(mirror) else np.nan
        predicted_iso = float(np.nansum(weight))
        gain_los = np.nanmean(direct, axis=0) if len(direct) else np.ones(np.count_nonzero(band))
        gain_mirror = np.nanmean(mirror, axis=0) if len(mirror) else np.ones(np.count_nonzero(band))

        par = j_pa[:, (PA_CENTERS < 30) | (PA_CENTERS > 150)]
        perp = j_pa[:, (PA_CENTERS > 60) & (PA_CENTERS < 120)]
        anisotropy = (np.nanmean(par, axis=1) - np.nanmean(perp, axis=1)) \
            / (np.nanmean(par, axis=1) + np.nanmean(perp, axis=1))
        solid = np.sin(np.deg2rad(PA_CENTERS))[None, :]
        aligned = np.nansum(j_pa * solid * ((PA_CENTERS < 30) | (PA_CENTERS > 150))[None, :], axis=1) \
            / np.nansum(j_pa * solid, axis=1)

        out[name] = {
            "records": len(subset),
            "energy_ev": energy,
            "j_pa": j_pa,
            "j_iso": j_iso,
            "peir_j": peir_j,
            "anisotropy": anisotropy,
            "aligned_fraction": aligned,
            "los_pa_median": los_median,
            "los_pa_p16": float(np.nanpercentile(los_pa, 16)),
            "los_pa_p84": float(np.nanpercentile(los_pa, 84)),
            "alpha_grid": alpha_grid,
            "gain_alpha": gain,
            "gain_los": gain_los,
            "gain_mirror": gain_mirror,
            "gain_band": band,
            "predicted": predicted,
            "predicted_mirror": predicted_mirror,
            "predicted_iso": predicted_iso,
            "los_velocity_cosine": float(np.nanmedian([row["reduced"]["look_velocity_cosine"] for row in subset])),
        }
    return out


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------
def plot_overview(ena, per_minute, esa, sst):
    times, energy, eflux = esa[0], esa[1], esa[2]
    sst_times, sst_energy, sst_eflux = sst[0], sst[1], sst[2]
    figure, axes = plt.subplots(3, 1, figsize=(11.5, 10.4), constrained_layout=True)
    axis = axes[0]
    axis.plot(per_minute.index, per_minute.to_numpy(), lw=0.7, color="C0")
    for (name, (start, end)), color in zip(WINDOWS.items(), ("C0", "C3")):
        axis.axvspan(pd.Timestamp(start), pd.Timestamp(end), color=color, alpha=0.10, label=name)
    for name, (start, end) in FOCUS.items():
        axis.axvspan(pd.Timestamp(start), pd.Timestamp(end), color="k", alpha=0.05, hatch="//",
                     label=name)
    axis.set(ylabel="GENA-A counts / min", title="GENA-A l2f event rate (1-min bins)")
    axis.legend(loc="upper right", fontsize=7)
    axis.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d %H:%M"))

    axis = axes[1]
    for name, color in zip(WINDOWS, ("C0", "C3")):
        item = ena[name]
        axis.plot(item["energy_kev"], item["flux"], "o-", ms=3, color=color, label=name)
    axis.set(xscale="log", yscale="log", xlabel="ENA energy (keV)",
             ylabel=r"proxy ENA flux (counts s$^{-1}$ keV$^{-1}$)",
             title="GENA-A ENA energy spectrum (geometry factor = 1)")
    axis.legend(fontsize=8)

    axis = axes[2]
    for name, color in zip(WINDOWS, ("C0", "C3")):
        start, end = WINDOWS[name]
        mask = (times >= pd.Timestamp(start)) & (times <= pd.Timestamp(end))
        e = np.nanmedian(energy[mask], axis=0)
        j = np.nanmedian(eflux[mask] / energy[mask], axis=0)
        axis.plot(e, j, "o-", ms=3, color=color, label=f"{name} ESA peir")
        smask = (sst_times >= pd.Timestamp(start)) & (sst_times <= pd.Timestamp(end))
        axis.plot(np.nanmedian(sst_energy[smask], axis=0), np.nanmedian(sst_eflux[smask] / sst_energy[smask], axis=0),
                  "s--", ms=3, color=color, alpha=0.7, label=f"{name} SST psif")
    axis.set(xscale="log", yscale="log", xlabel="ion energy (eV)",
             ylabel=r"$j_{ion}$ (cm$^{-2}$ s$^{-1}$ sr$^{-1}$ eV$^{-1}$)",
             title=f"THEMIS-{PROBE.upper()} ion spectra (ESA reduced + SST)")
    axis.legend(fontsize=7, ncol=2)
    figure.savefig(OUT_DIR / "fig1_overview.png", dpi=170)
    plt.close(figure)


def plot_pa(analysis):
    names = list(FOCUS)
    figure, axes = plt.subplots(2, 2, figsize=(11.5, 8.2), constrained_layout=True)
    values = np.concatenate([
        analysis[name]["j_pa"][(analysis[name]["energy_ev"] >= ENERGY_EV[0])
                               & (analysis[name]["energy_ev"] <= ENERGY_EV[1])].ravel()
        for name in names
    ])
    values = values[np.isfinite(values) & (values > 0)]
    norm = LogNorm(vmin=np.percentile(values, 5), vmax=np.percentile(values, 99))
    mesh = None
    for column, name in enumerate(names):
        item = analysis[name]
        band = (item["energy_ev"] >= ENERGY_EV[0]) & (item["energy_ev"] <= ENERGY_EV[1])
        axis = axes[0, column]
        mesh = axis.pcolormesh(PA_EDGES, _energy_edges(item["energy_ev"][band]) / 1e3,
                               item["j_pa"][band], shading="flat", cmap="turbo", norm=norm)
        axis.set(xlabel="pitch angle (deg)", ylabel="energy (keV)", xlim=(0, 180),
                 title=f"{name}: {item['records']} ESD records, "
                       f"LOS PA={item['los_pa_median']:.0f}° "
                       f"[{item['los_pa_p16']:.0f}-{item['los_pa_p84']:.0f}°]")
        axis.set_xticks(np.arange(0, 181, 30))
        axis = axes[1, column]
        energy = item["energy_ev"] / 1e3
        for (low, high), color in zip(((4, 6), (8, 12), (15, 24.6)), ("C0", "C1", "C2")):
            subset = (energy >= low) & (energy < high)
            profile = np.nanmean(item["j_pa"][subset] / item["j_iso"][subset][:, None], axis=0)
            axis.plot(PA_CENTERS, profile, "o-", color=color, label=f"{low:g}-{high:g} keV")
        axis.axhline(1.0, color="k", lw=0.5)
        axis.set(xlabel="pitch angle (deg)", ylabel=r"$J(\alpha)/\langle J\rangle$", xlim=(0, 180))
        axis.legend(fontsize=7)
    figure.colorbar(mesh, ax=axes[0, :], label=r"$j_{ion}$ (cm$^{-2}$ s$^{-1}$ sr$^{-1}$ eV$^{-1}$)",
                    pad=0.01, shrink=0.9)
    figure.suptitle(f"THEMIS-{PROBE.upper()} ESD ion pitch-angle distribution (15° bins)", fontsize=11)
    figure.savefig(OUT_DIR / "fig2_esd_pa.png", dpi=170)
    plt.close(figure)


def plot_geometry(rows, per_minute):
    """Time series of GENA rate, ENA/B geometry and the per-record PA anisotropy."""
    frame = pd.DataFrame([{"time": row["time"], "los_pa": row["los_pa"]} for row in rows])
    anisotropy = []
    for row in rows:
        energy = row["record"]["energy_ev"]
        band = (energy >= ENERGY_EV[0]) & (energy <= ENERGY_EV[1])
        j_pa = row["reduced"]["j_pa"][band]
        par = np.nanmean(j_pa[:, (PA_CENTERS < 30) | (PA_CENTERS > 150)], axis=1)
        perp = np.nanmean(j_pa[:, (PA_CENTERS > 60) & (PA_CENTERS < 120)], axis=1)
        anisotropy.append(np.nanmean((par - perp) / (par + perp)))
    frame["anisotropy"] = anisotropy

    figure, axes = plt.subplots(3, 1, figsize=(11.5, 8.6), constrained_layout=True, sharex=True)
    axis = axes[0]
    axis.plot(per_minute.index, per_minute.to_numpy(), lw=0.6, color="C0")
    axis.set(ylabel="GENA counts/min", title="ENA rate, GENA line-of-sight geometry and ion anisotropy")
    axis = axes[1]
    axis.plot(frame["time"], frame["los_pa"], "o", ms=3, color="C3", alpha=0.6)
    axis.plot(frame.set_index("time")["los_pa"].resample("30min").median(), color="k", lw=1.2)
    axis.set(ylabel=r"$\angle(v_{ENA}, B_{GSM})$ [deg]", ylim=(0, 180))
    axis = axes[2]
    axis.plot(frame["time"], frame["anisotropy"], "o", ms=3, color="C2", alpha=0.45)
    smoothed = frame.set_index("time")["anisotropy"].rolling("3h", min_periods=5).median()
    axis.plot(smoothed.index, smoothed.to_numpy(), color="C2", lw=1.8, label="3 h rolling median")
    axis.axhline(0.0, color="k", lw=0.5)
    axis.legend(fontsize=7, loc="lower right")
    axis.set(ylabel="PA anisotropy A (4-25 keV)", xlabel="time (UTC)")
    for axis in axes:
        for (name, (start, end)), color in zip(WINDOWS.items(), ("C0", "C3")):
            axis.axvspan(pd.Timestamp(start), pd.Timestamp(end), color=color, alpha=0.07, label=name)
        for start, end in FOCUS.values():
            axis.axvspan(pd.Timestamp(start), pd.Timestamp(end), color="k", alpha=0.04, hatch="//")
        axis.set_xlim(pd.Timestamp(TRANGE[0]), pd.Timestamp(TRANGE[1]))
    axis.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d %H:%M"))
    axes[0].legend(fontsize=7, loc="upper right")
    figure.savefig(OUT_DIR / "fig3_geometry.png", dpi=170)
    plt.close(figure)
    return frame


def plot_field_direction(rows, fgm):
    """Where do PA 0 deg and 180 deg point in GSM during the two focus windows?"""
    field = fgm[["%s_fgs_gsm_x" % PREFIX, "%s_fgs_gsm_y" % PREFIX, "%s_fgs_gsm_z" % PREFIX]].to_numpy(dtype=float)
    stamp = fgm["Timestamp"].to_numpy(dtype=float)
    figure, axes = plt.subplots(1, 2, figsize=(11.0, 4.6), constrained_layout=True)
    for name, color in zip(FOCUS, ("C0", "C3")):
        start, end = pd.Timestamp(FOCUS[name][0]), pd.Timestamp(FOCUS[name][1])
        times = [row["time"].timestamp() for row in rows if start <= row["time"] <= end]
        index, _ = _nearest(times, stamp)
        vector = field[index]
        vector = vector / np.linalg.norm(vector, axis=1)[:, None]
        axes[0].scatter(vector[:, 0], vector[:, 1], s=16, color=color, alpha=0.5, label=name)
        axes[1].scatter(vector[:, 0], vector[:, 2], s=16, color=color, alpha=0.5, label=name)
        mean = vector.mean(axis=0)
        axes[0].arrow(0, 0, mean[0], mean[1], color=color, width=0.01, length_includes_head=True)
        axes[1].arrow(0, 0, mean[0], mean[2], color=color, width=0.01, length_includes_head=True)
    for axis in axes:
        axis.axhline(0, color="k", lw=0.4)
        axis.axvline(0, color="k", lw=0.4)
        axis.set_xlim(-1.1, 1.1)
        axis.set_ylim(-1.1, 1.1)
        axis.set_aspect("equal")
        axis.set_xlabel(r"$B_X$ / |B|  (GSM, +X toward Earth)")
    axes[0].set_ylabel(r"$B_Y$ / |B|")
    axes[1].set_ylabel(r"$B_Z$ / |B|")
    axes[0].annotate("tailward (-X)", (-1.0, 0.05), fontsize=8, ha="right")
    axes[0].annotate("earthward (+X)", (1.0, 0.05), fontsize=8, ha="left")
    axes[0].set_title("GSM X-Y: 0 deg points along +B")
    axes[1].set_title("GSM X-Z")
    axes[0].legend(fontsize=7, loc="lower right")
    figure.savefig(OUT_DIR / "fig6_field_direction.png", dpi=170)
    plt.close(figure)


def plot_ceiling(analysis):
    """Show that even the most favourable line of sight cannot close the ENA gap."""
    names = list(FOCUS)
    figure, axes = plt.subplots(1, 2, figsize=(11.5, 4.4), constrained_layout=True,
                                width_ratios=(1.25, 1.0))
    axis = axes[0]
    for name, color in zip(names, ("C0", "C3")):
        item = analysis[name]
        band = item["gain_band"]
        profile = np.nanmean(item["j_pa"][band] / item["j_iso"][band][:, None], axis=0)
        axis.plot(PA_CENTERS, profile, "o-", color=color,
                  label=f"{name} (LOS {item['los_pa_median']:.0f}°)")
        axis.axvline(item["los_pa_median"], color=color, ls="--", lw=1)
        axis.axvspan(item["los_pa_median"] - 10, item["los_pa_median"] + 10, color=color, alpha=0.08)
    axis.axhline(1.0, color="k", lw=0.5)
    axis.set(xlabel="pitch angle (deg)", ylabel=r"$J(\alpha)/\langle J\rangle$ (4-25 keV)",
             xlim=(0, 180), title="Narrow beam vs. flux along the GENA line of sight")
    axis.legend(fontsize=7)

    axis = axes[1]
    peak, hill = names
    # How much could PA at most *favour the hill* relative to the peak?
    max_ratio = np.nanmax(analysis[hill]["gain_alpha"][:, analysis[hill]["gain_band"]]) \
        / np.nanmax(analysis[peak]["gain_alpha"][:, analysis[peak]["gain_band"]])
    los_ratio = np.nanmean(analysis[peak]["gain_los"]) / np.nanmean(analysis[hill]["gain_los"])
    labels = ["ENA counts\n(observed)", "PA ceiling\n(best LOS for hill)",
              "PA at actual\nGENA LOS"]
    values = [23.2, max_ratio, los_ratio]
    axis.bar(np.arange(3), values, color=["C0", "C2", "C3"])
    axis.set_xticks(np.arange(3), labels, fontsize=8)
    axis.set_yscale("log")
    axis.axhline(1.0, color="k", lw=0.8)
    axis.set_ylabel("peak / hill ratio")
    axis.set_title("Even the best geometry gives only ~1.7x")
    for index, value in enumerate(values):
        axis.text(index, value * 1.08, f"{value:.2f}", ha="center", fontsize=9)
    figure.savefig(OUT_DIR / "fig5_pa_ceiling.png", dpi=170)
    plt.close(figure)
    return {"max_gain_ratio": max_ratio, "los_gain_ratio": los_ratio}


def plot_summary(ena, analysis):
    peak, hill = list(WINDOWS)
    focus_peak, focus_hill = list(FOCUS)
    labels, values = [], []
    labels.append("GENA counts/min\n(burst/core)")
    values.append(ena[focus_peak]["counts_per_min"] / ena[focus_hill]["counts_per_min"])
    for band, low, high in (("2-10 keV", 2.0, 10.0), ("10-30 keV", 10.0, 30.0)):
        labels.append(f"ENA {band}\n(burst/core)")
        values.append(ena[focus_peak]["band_rate"][f"{low:g}-{high:g} keV"]
                      / ena[focus_hill]["band_rate"][f"{low:g}-{high:g} keV"])
    labels.append("ion production\n$\\int j\\sigma\\,dE$")
    values.append(analysis[focus_peak]["predicted_iso"] / analysis[focus_hill]["predicted_iso"])
    labels.append("ion + PA gain\n(predicted ENA)")
    values.append(analysis[focus_peak]["predicted"] / analysis[focus_hill]["predicted"])
    labels.append("PA gain only\n(at GENA LOS)")
    values.append(np.nanmean(analysis[focus_peak]["gain_los"])
                  / np.nanmean(analysis[focus_hill]["gain_los"]))

    figure, axis = plt.subplots(figsize=(8.6, 4.6), constrained_layout=True)
    colors = ["C0", "C0", "C0", "C1", "C2", "C3"]
    axis.bar(np.arange(len(values)), values, color=colors)
    axis.set_xticks(np.arange(len(values)), labels, fontsize=8)
    axis.set_yscale("log")
    axis.axhline(1.0, color="k", lw=0.8)
    axis.set_ylabel("ratio peak / hill")
    axis.set_title("Decomposition of the ENA contrast (burst 11-03 19:00 / core 11-04 09:00)")
    for index, value in enumerate(values):
        axis.text(index, value * 1.1, f"{value:.2f}", ha="center", fontsize=9)
    figure.savefig(OUT_DIR / "fig4_summary.png", dpi=170)
    plt.close(figure)
    return labels, values


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--no-figures", action="store_true")
    arguments = parser.parse_args()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    events = load_ena()
    fgm = load("themis", TRANGE, "fgm", probe=PROBE, variables=["fgs_gsm"], download=False)
    los = ena_direction_angle(events, fgm)
    ena, per_minute = ena_metrics(events)
    esa = load_esa()
    sst = load_sst()

    start, end = pd.Timestamp(TRANGE[0]), pd.Timestamp(TRANGE[1])
    days = pd.date_range(start.floor("D"), (end - pd.Timedelta(seconds=1)).floor("D"), freq="D")
    records = load_esd_records(days)
    los_time = los["Datetime"].astype("int64").to_numpy(dtype=float) / 1e9
    order = np.argsort(los_time)
    los_time, los_pa = los_time[order], los["LOS_PA"].to_numpy(dtype=float)[order]
    rows = []
    for record in records:
        if not (start <= record["time"] <= end):
            continue
        stamp = record["time"].timestamp()
        nearest = int(np.clip(np.searchsorted(los_time, stamp), 0, len(los_time) - 1))
        rows.append({
            "time": record["time"],
            "record": record,
            "reduced": reduce_pa(record),
            "los_pa": float(los_pa[nearest]) if abs(los_time[nearest] - stamp) <= 1200.0 else np.nan,
        })
    print(f"ESD records: {len(rows)}")
    analysis = window_analysis(rows, esa[0], esa[1], esa[2])

    summary = {
        "ena": {name: {"events": item["events"], "seconds": item["seconds"],
                       "counts_per_min": item["counts_per_min"]} for name, item in ena.items()},
        "pa": {name: {key: value for key, value in item.items()
                      if key not in ("j_pa", "energy_ev", "j_iso", "anisotropy", "aligned_fraction",
                                     "alpha_grid", "gain_alpha", "gain_los", "gain_mirror",
                                     "peir_j", "gain_band")}
               for name, item in analysis.items()},
    }
    for name, item in analysis.items():
        band = (item["energy_ev"] >= ENERGY_EV[0]) & (item["energy_ev"] <= ENERGY_EV[1])
        summary["pa"][name]["anisotropy_mean"] = float(np.nanmean(item["anisotropy"][band]))
        summary["pa"][name]["aligned_fraction_mean"] = float(np.nanmean(item["aligned_fraction"][band]))
        summary["pa"][name]["anisotropy_by_band"] = {
            f"{low:g}-{high:g} keV": float(np.nanmean(
                item["anisotropy"][(item["energy_ev"] >= low * 1e3) & (item["energy_ev"] < high * 1e3)]))
            for low, high in ((4, 6), (6, 10), (10, 15), (15, 25))
        }
        summary["pa"][name]["gain_los_mean"] = float(np.nanmean(item["gain_los"]))
        summary["pa"][name]["gain_mirror_mean"] = float(np.nanmean(item["gain_mirror"]))
    (OUT_DIR / "summary.json").write_text(json.dumps(summary, indent=2, default=float))
    print(json.dumps(summary, indent=2, default=float))

    if not arguments.no_figures:
        plot_overview(ena, per_minute, esa, sst)
        plot_pa(analysis)
        plot_geometry(rows, per_minute)
        print(list(zip(*plot_summary(ena, analysis))))
        print(plot_ceiling(analysis))
        plot_field_direction(rows, fgm)


if __name__ == "__main__":
    main()
