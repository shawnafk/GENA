#!/usr/bin/env python3
"""Reproduce the pitch-angle decomposition for THEMIS-A (near-Earth control).

THEMIS-B sits near the Moon (GSE X ~ -55 R_E), which is very far from the
likely ENA source in the near-Earth plasma sheet.  THEMIS-A provides a second,
near-Earth ion sample for the same intervals.  It has ESD coverage only through
2025-11-04, so the two focus windows are used and the given T1/T2 pair is
truncated at 2025-11-04 23:00.

The script imports the main comparison module, overrides its module-level
probe configuration, and writes ``out/tha_crosscheck.json``.
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("ena_pa", ROOT / "ena_pa_comparison.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)

module.PROBE = "a"
module.PREFIX = "tha"
module.DATA_DIR = ROOT / "data_tha"
module.TRANGE = ("2025-11-03T15:00:00", "2025-11-05T00:00:00")
module.WINDOWS = {
    "Peak (T1)": ("2025-11-03T15:00:00", "2025-11-04T09:00:00"),
    "Hill (T2)": ("2025-11-04T09:00:00", "2025-11-04T23:00:00"),
}
module.FOCUS = {
    "Peak burst": ("2025-11-03T19:00:00", "2025-11-04T00:00:00"),
    "Hill core": ("2025-11-04T09:00:00", "2025-11-04T23:00:00"),
}


def main():
    events = module.load_ena()
    fgm = module.load_data("themis", module.TRANGE, "fgm", probe="a",
                      variables=["fgs_gsm"], download=False)
    los = module.ena_direction_angle(events, fgm)
    ena, _ = module.ena_metrics(events)

    days = [pd.Timestamp("2025-11-03"), pd.Timestamp("2025-11-04")]
    records = module.load_esd_records(days)
    los_time = los["Datetime"].astype("int64").to_numpy(dtype=float) / 1e9
    los_pa = los["LOS_PA"].to_numpy(dtype=float)
    order = np.argsort(los_time)
    los_time, los_pa = los_time[order], los_pa[order]
    rows = []
    for record in records:
        stamp = record["time"].timestamp()
        nearest = int(np.clip(np.searchsorted(los_time, stamp), 0, len(los_time) - 1))
        rows.append({
            "time": record["time"],
            "record": record,
            "reduced": module.reduce_pa(record),
            "los_pa": float(los_pa[nearest]) if abs(los_time[nearest] - stamp) <= 1200.0 else np.nan,
        })
    times, energy, eflux, _, _ = module.load_esa()
    analysis = module.window_analysis(rows, times, energy, eflux)

    summary = {"probe": "A", "records": len(rows), "windows": {}}
    for name, item in analysis.items():
        band = item["gain_band"]
        summary["windows"][name] = {
            "records": item["records"],
            "counts_per_min": ena[name]["counts_per_min"] if name in ena else None,
            "anisotropy_mean": float(np.nanmean(item["anisotropy"][band])),
            "aligned_fraction_mean": float(np.nanmean(item["aligned_fraction"][band])),
            "gain_los_mean": float(np.nanmean(item["gain_los"])),
            "gain_mirror_mean": float(np.nanmean(item["gain_mirror"])),
            "predicted": item["predicted"],
            "predicted_iso": item["predicted_iso"],
        }
    burst, core = "Peak burst", "Hill core"
    summary["ratio_burst_core"] = {
        "counts": ena[burst]["counts_per_min"] / ena[core]["counts_per_min"],
        "ion_production": analysis[burst]["predicted_iso"] / analysis[core]["predicted_iso"],
        "ion_and_pa": analysis[burst]["predicted"] / analysis[core]["predicted"],
        "pa_gain": float(np.nanmean(analysis[burst]["gain_los"]) / np.nanmean(analysis[core]["gain_los"])),
    }
    (ROOT / "out" / "tha_crosscheck.json").write_text(json.dumps(summary, indent=2, default=float))
    print(json.dumps(summary, indent=2, default=float))


if __name__ == "__main__":
    main()
