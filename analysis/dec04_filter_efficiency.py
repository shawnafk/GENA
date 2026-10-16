"""Audit GENA sorting losses around 2025-12-04 from raw telemetry and products."""

from __future__ import annotations

from pathlib import Path

import cdflib
import matplotlib

matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from gena import filter as gf
from gena.constants import init_param
from gena.raw.structure import packet_struct


BASE = Path("/Users/shawn/ENA_BASE")
PQT = BASE / "GENA_PQT"
RAW = next((BASE / "GENA_RAW").glob("021.251204-251209*.raw"))
OUT = BASE / "GENA_RES" / "ion_rejection_20251204"
START = pd.Timestamp("2025-12-04 08:25:00", tz="UTC")
STOP = pd.Timestamp("2025-12-06 00:00:00", tz="UTC")
EDGES = pd.date_range(START, STOP, freq="5min")
EDGE_SECONDS = EDGES.to_numpy(dtype="datetime64[ns]").astype(np.int64) / 1e9
EVENTS = {
    "E1": ("2025-12-04 14:10", "2025-12-04 15:40"),
    "E2": ("2025-12-04 17:20", "2025-12-04 18:30"),
    "E3": ("2025-12-04 21:20", "2025-12-05 01:10"),
}
CONTROLS = {
    "E1": [("2025-12-04 12:00", "2025-12-04 14:10"), ("2025-12-04 15:40", "2025-12-04 17:20")],
    "E2": [("2025-12-04 16:30", "2025-12-04 17:20"), ("2025-12-04 18:30", "2025-12-04 19:10")],
    "E3": [("2025-12-04 20:00", "2025-12-04 21:20"), ("2025-12-05 01:10", "2025-12-05 02:00")],
}
COUNT_KEYS = ("downlinked", "triggered", "triggered_12", "four_start", "eight_channel", "soft", "tof_consistent", "golden")
END_KEYS = tuple(f"end_{count}" for count in range(5))


def _hist_packet_counts(packet_times: np.ndarray, per_packet: np.ndarray) -> np.ndarray:
    bins = np.searchsorted(EDGE_SECONDS, packet_times, side="right") - 1
    inside = (bins >= 0) & (bins < len(EDGES) - 1)
    return np.bincount(bins[inside], weights=per_packet[inside], minlength=len(EDGES) - 1).astype(np.int64)


def raw_counts() -> dict[str, dict[str, np.ndarray]]:
    packets = np.memmap(RAW, dtype=packet_struct, mode="r")
    result = {probe: {key: np.zeros(len(EDGES) - 1, dtype=np.int64) for key in (*COUNT_KEYS[:5], *END_KEYS)} for probe in "ab"}
    result["padding"] = np.zeros(len(EDGES) - 1, dtype=np.int64)
    start_s = START.timestamp()
    stop_s = STOP.timestamp()
    chunk_size = 10000
    for left in range(0, len(packets), chunk_size):
        packet = packets[left:left + chunk_size]
        times = 1609430400 + packet["time_s"].astype(np.int64)
        in_range = (times >= start_s) & (times < stop_s)
        if not in_range.any():
            continue
        packet = packet[in_range]
        times = times[in_range]
        if not np.all((packet["header"][:, 0] == 226) & (packet["header"][:, 1] == 37)):
            raise ValueError(f"Unaligned raw packets near index {left}")
        events = packet["events"]
        flags = events["flag"]
        present = flags != 0
        any_channel = np.zeros(flags.shape, dtype=bool)
        start = np.ones(flags.shape, dtype=bool)
        all_channels = np.ones(flags.shape, dtype=bool)
        end_hits = np.zeros(flags.shape, dtype=np.uint8)
        for name in ("c1", "c3", "c5", "c7"):
            channels = (events[name] & 0x1FFF) > 94
            any_channel |= channels.any(axis=2)
            start &= channels[:, :, 0]
            all_channels &= channels.all(axis=2)
            end_hits += channels[:, :, 1]
        if (any_channel & ~present).any():
            raise ValueError(f"An empty event slot has a position TDC signal near packet {left}")
        any_channel_12 = any_channel | ((events["c0"] & 0x1FFF) > 94).any(axis=2)
        result["padding"] += _hist_packet_counts(times, (~present).sum(axis=1))
        for probe, bit in (("a", 0), ("b", 1)):
            selected = present & ((flags >> 7) == bit)
            for key, mask in (("downlinked", selected), ("triggered", selected & any_channel),
                              ("triggered_12", selected & any_channel_12),
                              ("four_start", selected & start),
                              ("eight_channel", selected & all_channels)):
                result[probe][key] += _hist_packet_counts(times, mask.sum(axis=1))
            for count, key in enumerate(END_KEYS):
                result[probe][key] += _hist_packet_counts(
                    times, (selected & start & (end_hits == count)).sum(axis=1)
                )
    return result


def _event_hist(frame: pd.DataFrame, mask=None) -> np.ndarray:
    timestamp = frame["Timestamp"].to_numpy(dtype=float)
    if mask is not None:
        timestamp = timestamp[np.asarray(mask, dtype=bool)]
    return np.histogram(timestamp, bins=EDGE_SECONDS)[0]


def product_counts(result: dict[str, dict[str, np.ndarray]]) -> None:
    params = init_param(1)
    for probe in "ab":
        pieces = {product: [] for product in ("l1a", "l2", "l2s", "l2g")}
        for day in ("20251204", "20251205"):
            for product in pieces:
                pieces[product].append(pd.read_parquet(PQT / f"{day}_{product}_{probe}.pqt"))
        frames = {key: pd.concat(value, ignore_index=True) for key, value in pieces.items()}
        l1a, l2 = frames["l1a"], frames["l2"]
        if len(l1a) != len(l2) or not l1a[["UUID", "PID"]].equals(l2[["UUID", "PID"]]):
            raise ValueError(f"L1A and L2 rows do not align for probe {probe}")
        soft = np.asarray(gf.soft_std(params=params, prob=probe)(l1a), dtype=bool)
        tof = np.asarray(gf.tofxy_eff(l1a), dtype=bool)
        strict = np.asarray(gf.golden_std(params=params, prob=probe)(l1a), dtype=bool)
        if not np.array_equal(strict, soft & tof):
            raise ValueError(f"Current golden rule disagrees with soft + TOF for probe {probe}")
        spatial = np.asarray(gf.onboard_eff(l2, c_size=[150, 60]), dtype=bool)
        if not np.array_equal(_event_hist(l1a), result[probe]["eight_channel"]):
            raise ValueError(f"Raw 8-channel counts do not match L1A for probe {probe}")
        l2s, l2g = frames["l2s"], frames["l2g"]
        soft_ids = pd.MultiIndex.from_frame(l2s[["UUID", "PID"]])
        golden_ids = pd.MultiIndex.from_frame(l2g[["UUID", "PID"]])
        archived_golden = soft_ids.isin(golden_ids)
        archived_spatial = np.asarray(gf.onboard_eff(l2s, c_size=[150, 60]), dtype=bool)
        if not np.array_equal(archived_golden, archived_spatial):
            raise ValueError(f"Archived L2G is not the spatial subset of L2S for probe {probe}")
        tof_on_soft = ((l2s["Tof_X"] / l2s["Tof_Y"] > 0.937)
                       & (l2s["Tof_X"] / l2s["Tof_Y"] < 1.117)
                       & (l2s["Tof_X"] > 0) & (l2s["Tof_Y"] > 0))
        result[probe]["soft"] = _event_hist(l2s)
        result[probe]["tof_consistent"] = _event_hist(l2s, tof_on_soft)
        result[probe]["golden"] = _event_hist(l2g)
        result[probe]["current_soft"] = _event_hist(l1a, soft)
        result[probe]["current_golden"] = _event_hist(l1a, strict & spatial)


def themis_density() -> pd.DataFrame:
    frames = []
    for day in ("20251204", "20251205"):
        path = BASE / "SPEDAS" / "themis_data" / "thb" / "l2" / "esa" / "2025" / f"thb_l2_esa_{day}_v01.cdf"
        cdf = cdflib.CDF(str(path))
        frames.append(pd.DataFrame({"Date": pd.to_datetime(cdf.varget("thb_peir_time"), unit="s", utc=True),
                                    "thb_density_cm3": cdf.varget("thb_peir_density")}))
    density = pd.concat(frames, ignore_index=True).set_index("Date").sort_index()
    density = density.where(density["thb_density_cm3"].between(0, 100))
    return density.resample("5min", origin=START).median()


def build_table(counts: dict[str, dict[str, np.ndarray]]) -> pd.DataFrame:
    table = pd.DataFrame({"start_utc": EDGES[:-1], "end_utc": EDGES[1:]})
    for probe in "ab":
        for key in COUNT_KEYS:
            table[f"{probe}_{key}"] = counts[probe][key]
        for key in END_KEYS:
            table[f"{probe}_{key}"] = counts[probe][key]
        table[f"{probe}_current_soft"] = counts[probe]["current_soft"]
        table[f"{probe}_current_golden"] = counts[probe]["current_golden"]
        n = table[f"{probe}_downlinked"].to_numpy()
        triggered = table[f"{probe}_triggered"].to_numpy()
        triggered_12 = table[f"{probe}_triggered_12"].to_numpy()
        eight = table[f"{probe}_eight_channel"].to_numpy()
        soft = table[f"{probe}_soft"].to_numpy()
        tof = table[f"{probe}_tof_consistent"].to_numpy()
        golden = table[f"{probe}_golden"].to_numpy()
        for key, numerator, denominator in (
            ("reject_8", n - eight, n),
            ("reject_soft", eight - soft, eight),
            ("reject_tof", soft - tof, soft),
            ("reject_spatial", soft - golden, soft),
            ("reject_total", n - golden, n),
            ("accept_8_of_four", eight, table[f"{probe}_four_start"].to_numpy()),
            ("accept_final_end", table[f"{probe}_end_4"].to_numpy(),
             (table[f"{probe}_end_3"] + table[f"{probe}_end_4"]).to_numpy()),
            ("accept_l2g_of_triggered", golden, triggered),
            ("accept_current_golden_of_triggered", table[f"{probe}_current_golden"].to_numpy(), triggered),
            ("accept_l2g_of_triggered_12", golden, triggered_12),
            ("accept_current_golden_of_triggered_12", table[f"{probe}_current_golden"].to_numpy(), triggered_12),
        ):
            table[f"{probe}_{key}"] = np.divide(numerator, denominator,
                                                  out=np.full(len(n), np.nan), where=denominator > 0)
    table["padding_slots"] = counts["padding"]
    density = themis_density().reindex(EDGES[:-1])
    table["thb_density_cm3"] = density["thb_density_cm3"].to_numpy()
    return table


def _window_mask(table: pd.DataFrame, windows: list[tuple[str, str]]) -> np.ndarray:
    center = table["start_utc"] + pd.Timedelta(minutes=2.5)
    selected = np.zeros(len(table), dtype=bool)
    for start, stop in windows:
        selected |= (center >= pd.Timestamp(start, tz="UTC")) & (center < pd.Timestamp(stop, tz="UTC"))
    return selected


def summarize(table: pd.DataFrame) -> pd.DataFrame:
    rows = []
    windows = [("whole", "all", np.ones(len(table), dtype=bool))]
    for name, interval in EVENTS.items():
        windows.extend(((name, "event", _window_mask(table, [interval])),
                        (name, "control", _window_mask(table, CONTROLS[name]))))
    for name, kind, mask in windows:
        subset = table.loc[mask]
        for probe in "ab":
            totals = {key: int(subset[f"{probe}_{key}"].sum()) for key in COUNT_KEYS}
            totals.update({key: int(subset[f"{probe}_{key}"].sum()) for key in END_KEYS})
            n, triggered, triggered_12, eight, soft, tof, gold = (totals[key] for key in
                                                    ("downlinked", "triggered", "triggered_12", "eight_channel",
                                                     "soft", "tof_consistent", "golden"))
            current_gold = int(subset[f"{probe}_current_golden"].sum())
            rows.append({"window": name, "kind": kind, "probe": probe, "bins": len(subset), **totals,
                         "current_soft": int(subset[f"{probe}_current_soft"].sum()),
                         "current_golden": current_gold,
                         "trigger_rate_min": triggered / (5 * len(subset)) if len(subset) else np.nan,
                         "trigger_12_rate_min": triggered_12 / (5 * len(subset)) if len(subset) else np.nan,
                         "accept_l2g_of_triggered": gold / triggered if triggered else np.nan,
                         "accept_current_golden_of_triggered": current_gold / triggered if triggered else np.nan,
                         "accept_l2g_of_triggered_12": gold / triggered_12 if triggered_12 else np.nan,
                         "accept_current_golden_of_triggered_12": current_gold / triggered_12 if triggered_12 else np.nan,
                         "reject_8": (n - eight) / n if n else np.nan,
                         "reject_soft": (eight - soft) / eight if eight else np.nan,
                         "reject_tof": (soft - tof) / soft if soft else np.nan,
                         "reject_spatial": (soft - gold) / soft if soft else np.nan,
                         "reject_total": (n - gold) / n if n else np.nan,
                         "accept_8_of_four": eight / totals["four_start"] if totals["four_start"] else np.nan,
                         "accept_final_end": totals["end_4"] / (totals["end_3"] + totals["end_4"])
                         if totals["end_3"] + totals["end_4"] else np.nan,
                         "thb_density_median_cm3": subset["thb_density_cm3"].median()})
    return pd.DataFrame(rows)


def plot(table: pd.DataFrame) -> None:
    time = table["start_utc"].dt.tz_localize(None) + pd.Timedelta(minutes=2.5)
    fig, axes = plt.subplots(4, 1, figsize=(14, 10), sharex=True,
                             gridspec_kw={"height_ratios": [1.15, 1.05, 1, 0.85]})
    colors = {"downlinked": "#75808c", "eight_channel": "#df8c32", "golden": "#0b847c"}
    for key, label in (("downlinked", "Nonpadding telemetry records"), ("eight_channel", "8-channel"),
                       ("golden", "Archived L2G")):
        axes[0].plot(time, table[f"a_{key}"] / 5, color=colors[key], lw=1.3, label=label)
    axes[0].plot(time, table["b_eight_channel"] / 5, color="#a64d79", lw=1, label="B 8-channel")
    axes[0].set_ylabel("GENA-A count / min")
    axes[0].set_yscale("log")
    axes[0].legend(loc="upper right", ncol=3, fontsize=9)

    axes[1].plot(time, 100 * (1 - table["a_reject_8"]), color="#b94b3a", lw=1.2,
                 label="A: 8-channel / nonpadding records")
    axes[1].plot(time, 100 * table["a_accept_8_of_four"], color="#df8c32", lw=1.2,
                 label="A: 8-channel / 4-start")
    axes[1].plot(time, 100 * (1 - table["b_reject_8"]), color="#a64d79", lw=1,
                 label="B: 8-channel / nonpadding records")
    axes[1].set_yscale("log")
    axes[1].set_ylabel("8-channel acceptance (%)")
    axes[1].legend(loc="upper right", ncol=2, fontsize=9)

    axes[2].plot(time, 100 * (1 - table["a_reject_soft"]), color="#0b847c", lw=1.2,
                 label="A: archived L2S / 8-channel")
    axes[2].plot(time, 100 * (1 - table["a_reject_tof"]), color="#4173ab", lw=1.1,
                 label="A: current X/Y TOF criterion / L2S")
    axes[2].set_ylabel("After 8 channels: accepted (%)")
    axes[2].set_ylim(0, 105)
    axes[2].legend(loc="upper right", fontsize=9)

    axes[3].plot(time, table["thb_density_cm3"], color="#416c9c", lw=1.2)
    axes[3].set_yscale("log")
    axes[3].set_ylabel(r"THEMIS-B ion density (cm$^{-3}$)")
    axes[3].set_xlabel("2025-12-04 to 05 UTC; 5-min bins")
    for ax in axes:
        ax.grid(alpha=0.2)
        for name, (start, stop) in EVENTS.items():
            ax.axvspan(pd.Timestamp(start), pd.Timestamp(stop), color="#e6a63b", alpha=0.13)
        ax.set_xlim(time.iloc[0], time.iloc[-1])
    for name, (start, stop) in EVENTS.items():
        midpoint = pd.Timestamp(start) + (pd.Timestamp(stop) - pd.Timestamp(start)) / 2
        axes[0].text(midpoint, 0.97, name, transform=axes[0].get_xaxis_transform(),
                     ha="center", va="top", fontsize=9)
    axes[-1].xaxis.set_major_locator(mdates.HourLocator(interval=4))
    axes[-1].xaxis.set_major_formatter(mdates.DateFormatter("%d %H:%M"))
    fig.suptitle("GENA sorting losses around the 2025-12-04 ENA intervals", fontsize=14)
    fig.text(0.01, 0.01, "Nonpadding records have flag != 0. Archived L2G equals L2S + spatial cut. "
             "THEMIS-B density is remote, not QQ-2 in-situ.", fontsize=8)
    fig.tight_layout(rect=(0, 0.025, 1, 0.97))
    fig.savefig(OUT / "filter_efficiency_timeseries.png", dpi=180)
    fig.savefig(OUT / "filter_efficiency_timeseries.pdf")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
    selected = table["a_four_start"] >= 100
    event_bins = _window_mask(table, list(EVENTS.values()))
    for event, color, label in ((False, "#81919c", "Other times"),
                                (True, "#0b847c", "ENA intervals")):
        subset = selected & (event_bins == event)
        axes[0].scatter(table.loc[subset, "a_four_start"] / 5,
                        100 * table.loc[subset, "a_accept_8_of_four"],
                        s=14, alpha=0.55, color=color, label=label)
        with_density = subset & table["thb_density_cm3"].notna()
        axes[1].scatter(table.loc[with_density, "thb_density_cm3"],
                        100 * table.loc[with_density, "a_accept_8_of_four"],
                        s=14, alpha=0.55, color=color, label=label)
    axes[0].set_xscale("log")
    axes[0].set_xlabel("GENA-A 4-start records / min")
    axes[0].set_ylabel("8-channel / 4-start (%)")
    axes[1].set_xscale("log")
    axes[1].set_xlabel(r"THEMIS-B ion density (cm$^{-3}$; remote)")
    axes[1].set_ylabel("8-channel / 4-start (%)")
    for ax in axes:
        ax.grid(alpha=0.2)
    axes[0].legend(frameon=False)
    fig.tight_layout()
    fig.savefig(OUT / "filter_efficiency_scatter.png", dpi=180)
    plt.close(fig)


def plot_event_control(summary: pd.DataFrame) -> None:
    data = summary[(summary["probe"] == "a") & (summary["window"] != "whole")].copy()
    labels = list(EVENTS)
    x = np.arange(len(labels))
    fig, axes = plt.subplots(1, 4, figsize=(17, 4.3))
    for axis, numerator, denominator, ylabel in (
        (axes[0], "four_start", "downlinked", "4-start / nonpadding records (%)"),
        (axes[1], "eight_channel", "four_start", "8-channel / 4-start (%)"),
        (axes[2], "end_4", "near_complete", "4-end / (3-end + 4-end) (%)"),
        (axes[3], "soft", "eight_channel", "Archived L2S / 8-channel (%)"),
    ):
        for kind, offset, color in (("event", -0.19, "#0b847c"), ("control", 0.19, "#acb3b9")):
            part = data[data["kind"] == kind].set_index("window").loc[labels]
            denom = (part["end_3"] + part["end_4"]) if denominator == "near_complete" else part[denominator]
            values = 100 * part[numerator].to_numpy() / denom.to_numpy()
            axis.bar(x + offset, values, width=0.34, color=color, label=kind.capitalize())
            for location, value in zip(x + offset, values):
                axis.text(location, value, f"{value:.2f}", ha="center", va="bottom", fontsize=8)
        axis.set_xticks(x, labels)
        axis.set_ylabel(ylabel)
        axis.set_ylim(0, axis.get_ylim()[1] * 1.12)
        axis.grid(axis="y", alpha=0.2)
    axes[0].legend(frameon=False)
    fig.suptitle("GENA-A sorting acceptance: ENA intervals vs adjacent controls")
    fig.tight_layout()
    fig.savefig(OUT / "event_control_efficiency.png", dpi=180)
    plt.close(fig)


def plot_end_multiplicity(summary: pd.DataFrame) -> None:
    data = summary[(summary["probe"] == "a") & (summary["window"] != "whole")]
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.8), sharey=True)
    x = np.arange(5)
    for axis, name in zip(axes, EVENTS):
        for kind, offset, color in (("event", -0.19, "#0b847c"), ("control", 0.19, "#acb3b9")):
            row = data[(data["window"] == name) & (data["kind"] == kind)].iloc[0]
            values = 100 * row[list(END_KEYS)].to_numpy(dtype=float) / row["four_start"]
            axis.bar(x + offset, values, width=0.34, color=color, label=kind.capitalize())
        axis.set_title(name)
        axis.set_xticks(x, [str(value) for value in x])
        axis.set_xlabel("Number of end channels fired")
        axis.grid(axis="y", alpha=0.2)
    axes[0].set_ylabel("Among 4-start records (%)")
    axes[0].legend(frameon=False)
    fig.suptitle("GENA-A end-channel pairing conditional on all four start channels")
    fig.tight_layout()
    fig.savefig(OUT / "end_channel_pairing.png", dpi=180)
    plt.close(fig)


def neighbor_segments(table: pd.DataFrame) -> pd.DataFrame:
    slots = table["padding_slots"] + table["a_downlinked"] + table["b_downlinked"]
    if (slots % 62 != 0).any():
        raise ValueError("Event slots do not form complete 62-slot packets")
    rows = []
    for name, interval in EVENTS.items():
        for phase, windows in (("Before", [CONTROLS[name][0]]),
                               ("Enhanced", [interval]),
                               ("After", [CONTROLS[name][1]])):
            selected = _window_mask(table, windows)
            part = table.loc[selected]
            packets = slots.loc[selected].sum() / 62
            triggered = part["a_triggered"].sum()
            eight = part["a_eight_channel"].sum()
            final = part["a_current_golden"].sum()
            rows.append({"window": name, "phase": phase, "start_utc": windows[0][0],
                         "end_utc": windows[0][1], "bins": len(part), "packets": int(packets),
                         "triggered": int(triggered), "eight_channel": int(eight),
                         "current_golden": int(final), "packet_rate_s": packets / (300 * len(part)),
                         "trigger_per_packet": triggered / packets,
                         "trigger_rate_s": triggered / (300 * len(part)),
                         "eight_pct": 100 * eight / triggered,
                         "final_of_eight_pct": 100 * final / eight,
                         "final_pct": 100 * final / triggered})
    return pd.DataFrame(rows)


def plot_trigger_efficiency(table: pd.DataFrame, segments: pd.DataFrame) -> None:
    time = table["start_utc"].dt.tz_localize(None) + pd.Timedelta(minutes=2.5)
    fig, axes = plt.subplots(3, 1, figsize=(13, 8), sharex=True)
    axes[0].plot(time, table["a_triggered"] / 300, color="#b45536", lw=1.25)
    axes[0].set_ylabel("A downlinked triggers / s")
    eight_yield = np.divide(table["a_eight_channel"], table["a_triggered"],
                            out=np.full(len(table), np.nan), where=table["a_triggered"] > 0)
    axes[1].plot(time, 100 * eight_yield, color="#cf842d", lw=1.25)
    axes[1].set_yscale("log")
    axes[1].set_ylabel("8-channel / triggered (%)")
    axes[2].plot(time, 100 * table["a_accept_current_golden_of_triggered"],
                 color="#087e79", lw=1.3)
    axes[2].set_yscale("log")
    axes[2].set_ylabel("Current final / triggered (%)")
    axes[2].set_xlabel("2025-12-04 to 05 UTC; 5-min bins")
    for ax in axes:
        ax.grid(alpha=0.2)
        for name, (start, stop) in EVENTS.items():
            ax.axvspan(pd.Timestamp(start), pd.Timestamp(stop), color="#e6a63b", alpha=0.14)
        ax.set_xlim(time.iloc[0], time.iloc[-1])
    for name, (start, stop) in EVENTS.items():
        midpoint = pd.Timestamp(start) + (pd.Timestamp(stop) - pd.Timestamp(start)) / 2
        axes[0].text(midpoint, 0.97, name, transform=axes[0].get_xaxis_transform(),
                     ha="center", va="top", fontsize=9)
    axes[-1].xaxis.set_major_locator(mdates.HourLocator(interval=4))
    axes[-1].xaxis.set_major_formatter(mdates.DateFormatter("%d %H:%M"))
    fig.suptitle("GENA-A: downlinked trigger records and current filter yield")
    fig.text(0.01, 0.01, "Triggered = >=1 of 8 position TDC channels; final = golden_std + onboard_eff. "
             "Rates describe downlinked records, not all front-end triggers.", fontsize=8)
    fig.tight_layout(rect=(0, 0.03, 1, 0.97))
    fig.savefig(OUT / "trigger_l2g_timeseries.png", dpi=180)
    plt.close(fig)

    x = np.arange(len(EVENTS))
    fig, axes = plt.subplots(2, 3, figsize=(15, 7.3))
    metrics = (
        ("packet_rate_s", "Downlinked packets / s", "{:.1f}"),
        ("trigger_per_packet", "A triggered / packet", "{:.1f}"),
        ("trigger_rate_s", "A downlinked triggers / s", "{:.0f}"),
        ("eight_pct", "8-channel / triggered (%)", "{:.3f}"),
        ("final_of_eight_pct", "Current final / 8-channel (%)", "{:.1f}"),
        ("final_pct", "Current final / triggered (%)", "{:.3f}"),
    )
    for ax, (column, ylabel, number_format) in zip(axes.flat, metrics):
        for phase, offset, color in (("Before", -0.25, "#9aa6ad"),
                                     ("Enhanced", 0, "#087e79"),
                                     ("After", 0.25, "#d38a3e")):
            values = segments.loc[segments["phase"] == phase].set_index("window").loc[list(EVENTS), column]
            bars = ax.bar(x + offset, values, width=0.23, color=color, label=phase)
            ax.bar_label(bars, labels=[number_format.format(value) for value in values],
                         padding=2, fontsize=7)
        ax.set_xticks(x, list(EVENTS))
        ax.set_ylabel(ylabel)
        ax.set_ylim(0, segments[column].max() * 1.18)
        ax.grid(axis="y", alpha=0.2)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.945),
               ncol=3, frameon=False)
    fig.suptitle("GENA-A: ENA-enhanced intervals and neighboring time windows")
    fig.text(0.01, 0.01, "Before/After are temporal neighbors, not verified plasma-quiet intervals. "
             "Packets/s x A triggers/packet = A downlinked triggers/s.", fontsize=8)
    fig.tight_layout(rect=(0, 0.035, 1, 0.895))
    fig.savefig(OUT / "trigger_l2g_event_control.png", dpi=180)
    plt.close(fig)


def plot_trigger_rate_relation(table: pd.DataFrame) -> pd.DataFrame:
    data = table.loc[table["a_triggered"] > 0].copy()
    data["trigger_rate_s"] = data["a_triggered"] / 300
    data["final_pct"] = 100 * data["a_current_golden"] / data["a_triggered"]
    data["quartile"] = pd.qcut(data["trigger_rate_s"], 4, labels=["Q1", "Q2", "Q3", "Q4"])
    event = _window_mask(data, list(EVENTS.values()))
    control = _window_mask(data, [interval for windows in CONTROLS.values() for interval in windows])

    rows = []
    for name, part in data.groupby("quartile", observed=True):
        rows.append({"quartile": name, "bins": len(part),
                     "trigger_rate_min_s": part["trigger_rate_s"].min(),
                     "trigger_rate_max_s": part["trigger_rate_s"].max(),
                     "trigger_rate_median_s": part["trigger_rate_s"].median(),
                     "final_pct_median": part["final_pct"].median(),
                     "final_pct_pooled": 100 * part["a_current_golden"].sum() / part["a_triggered"].sum()})
    quartiles = pd.DataFrame(rows)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.7), gridspec_kw={"width_ratios": [1.65, 1]})
    for selected, label, color in ((~(event | control), "Other times", "#919ca5"),
                                   (control, "Adjacent controls", "#d38a3e"),
                                   (event, "Selected intervals", "#087e79")):
        part = data.loc[selected & (data["final_pct"] > 0)]
        axes[0].scatter(part["trigger_rate_s"], part["final_pct"], s=19, alpha=0.7,
                        color=color, edgecolors="none", label=label)
    axes[0].set_yscale("log")
    axes[0].set_xlabel("A downlinked triggers / s")
    axes[0].set_ylabel("Current final / triggered (%)")
    axes[0].legend(frameon=False, fontsize=8)
    axes[0].grid(alpha=0.2)

    x = np.arange(len(quartiles))
    axes[1].bar(x, quartiles["final_pct_pooled"], width=0.65, color="#087e79", alpha=0.85,
                label="Pooled fraction")
    axes[1].scatter(x, quartiles["final_pct_median"], color="#263943", s=35, zorder=3,
                    label="Median 5-min fraction")
    axes[1].set_xticks(x, quartiles["quartile"])
    axes[1].set_xlabel("Quartile of A downlinked trigger rate")
    axes[1].set_ylabel("Current final / triggered (%)")
    axes[1].legend(frameon=False, fontsize=8)
    axes[1].grid(axis="y", alpha=0.2)
    fig.suptitle("GENA-A: recorded trigger rate vs current final yield (5-min bins)")
    fig.text(0.01, 0.01, "Three bins with zero final events are omitted from the log scatter; all bins enter the quartiles. "
             "Rate is limited by telemetry coverage and capacity.", fontsize=8)
    fig.tight_layout(rect=(0, 0.03, 1, 0.96))
    fig.savefig(OUT / "trigger_l2g_rate_relation.png", dpi=180)
    plt.close(fig)
    return quartiles


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    counts = raw_counts()
    product_counts(counts)
    table = build_table(counts)
    summary = summarize(table)
    table.to_csv(OUT / "five_minute_counts.csv", index=False)
    summary.to_csv(OUT / "event_control_summary.csv", index=False)
    plot_end_multiplicity(summary)
    segments = neighbor_segments(table)
    plot_trigger_efficiency(table, segments)
    segments.to_csv(OUT / "a_filter_neighbor_segments.csv", index=False)
    quartiles = plot_trigger_rate_relation(table)
    table[["start_utc", "end_utc", "a_triggered", "a_eight_channel", "a_current_soft",
           "a_current_golden", "a_accept_current_golden_of_triggered"]].to_csv(
               OUT / "a_filter_five_minute.csv", index=False)
    summary.loc[summary["probe"] == "a", ["window", "kind", "bins", "triggered",
             "eight_channel", "current_soft", "current_golden", "trigger_rate_min",
             "accept_current_golden_of_triggered"]].to_csv(OUT / "a_filter_event_control.csv", index=False)
    quartiles.to_csv(OUT / "a_trigger_rate_quartiles.csv", index=False)
    focused = summary.loc[summary["probe"] == "a", ["window", "kind", "triggered",
                       "eight_channel", "current_golden", "trigger_rate_min",
                       "accept_current_golden_of_triggered"]]
    print(focused.to_string(index=False, float_format=lambda x: f"{x:.6f}"))
    print(f"Saved to {OUT}")


if __name__ == "__main__":
    main()
