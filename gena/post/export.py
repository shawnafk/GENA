"""Export reusable GENA post-processing products and observation reports."""

from __future__ import annotations

import datetime as dt
from pathlib import Path
import re

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common.config import resolve_output_path
from common.time import Time
from libs.utils import normalize_probe
from plots.vis import savefig
from plots import ani as gena_ani
from gena import utils as gut
from plots import v3d as gena_v3d
from plots.basic import gena as gena_vis
from plots.advance import env
from gena.frame import ParameterDataFrame
from load import load_gena

from .project import get_frame, get_ray, get_spacecraft_pos

__all__ = [
    "batch_plot",
    "build_observation_summary",
    "dump",
    "load_products",
    "normalize_report_times",
    "observation_summary",
    "report",
    "save_figures",
]

_REPORT_DPI = 300
_OBSERVATION_SUMMARY_COLUMNS = [
    "observation_index",
    "parent_l_level_index",
    "sub_observation_index",
    "probe",
    "start_time",
    "end_time",
    "duration_seconds",
    "spacecraft_x_mean",
    "spacecraft_y_mean",
    "spacecraft_z_mean",
    "pointing_x_min",
    "pointing_x_max",
    "pointing_x_mean",
    "pointing_z_min",
    "pointing_z_max",
    "pointing_z_mean",
    "sun_earth_angle_deg",
    "target_mode",
    "aiming",
    "ena_product",
    "ena_total_counts",
    "ena_counts_per_second",
]
_FIGURE_FILENAMES = {
    "summary": "{prob}_summary.jpg",
    "levels_count": "{prob}_levels_count.jpg",
    "parameters": "{prob}_parameters.jpg",
    "orbit_pointing": "{prob}_orbit_pointing.jpg",
    "environment": "{prob}_environment.jpg",
    "start_signal_distribution": "{prob}_start_signal_distribution.jpg",
    "end_signal_distribution": "{prob}_end_signal_distribution.jpg",
    "joint_signal_distribution": "{prob}_joint_signal_distribution.jpg",
    "hist_E": "{prob}_hist_E.jpg",
    "hist_Layer1": "{prob}_hist_Layer1.jpg",
    "hist_Layer2": "{prob}_hist_Layer2.jpg",
    "hist_FOV": "{prob}_hist_FOV.jpg",
}


def _day_period(day):
    start = Time(day).datetime.floor("D")
    return start, start + pd.Timedelta(days=1) - pd.Timedelta(nanoseconds=1)


def _has_explicit_clock(value):
    """Return whether one input explicitly represents sub-day time precision."""
    if isinstance(value, str):
        return bool(re.search(r"(?:[Tt]|\s)\d{1,2}:\d{2}", value.strip()))
    if isinstance(value, bytes):
        return _has_explicit_clock(value.decode())
    if isinstance(value, dt.datetime):
        return any(
            (value.hour, value.minute, value.second, value.microsecond)
        )
    if isinstance(value, dt.date):
        return False
    if isinstance(value, np.datetime64):
        return value != value.astype("datetime64[D]")
    if isinstance(value, (int, float, np.integer, np.floating)):
        return True
    return False


def normalize_report_times(time):
    """Normalize report intervals, assuming all input times represent UTC.

    A scalar or one-element input selects one complete UTC day. A two-element
    input containing an explicit clock time selects one exact interval and is
    not split at midnight. Two date-only endpoints expand to an inclusive
    sequence of daily intervals. Longer date-only sequences select discrete
    UTC days.
    """
    scalar_types = (
        str,
        bytes,
        dt.date,
        pd.Timestamp,
        np.datetime64,
        int,
        float,
        np.integer,
        np.floating,
    )
    if isinstance(time, scalar_types):
        requested = [time]
    else:
        try:
            requested = list(time)
        except TypeError:
            requested = [time]
    if not requested:
        raise ValueError("time must contain at least one value")
    if len(requested) == 1:
        return (_day_period(requested[0]),)
    if len(requested) == 2 and any(_has_explicit_clock(value) for value in requested):
        start, stop = Time(requested).datetime
        if start > stop:
            raise ValueError("report time-range start must not be later than end")
        return ((start, stop),)
    if any(_has_explicit_clock(value) for value in requested):
        raise ValueError("clock-precision report time must contain exactly [start, end]")
    if len(requested) == 2:
        start, stop = Time(requested).datetime
        if start > stop:
            raise ValueError("day range start must not be later than end")
        days = pd.date_range(start.floor("D"), stop.floor("D"), freq="D", tz="UTC")
    else:
        days = Time(requested).floor("D").tz_localize("UTC").drop_duplicates().sort_values()
    return tuple(_day_period(day) for day in days)


def build_observation_summary(
    pdf,
    ena,
    prob,
    *,
    large_gap="1h",
    small_gap="5min",
    coord="GSM",
    ena_product="l2g",
):
    """Build one row per continuous GENA observation segment.

    Engineering samples are first split at ``large_gap`` and then at
    ``small_gap``. Target mode is classified from the angle between the mean
    probe ``-Z`` look direction and the GSM ``+X`` Sun-Earth axis: below 10
    degrees is Earth, 15--40 degrees is mid-tail, and above 120 degrees is
    far-tail. Angles outside those intervals are reported as unknown.
    """
    prob = normalize_probe(prob)
    large_gap = pd.Timedelta(large_gap)
    small_gap = pd.Timedelta(small_gap)
    if len(pdf) == 0:
        result = pd.DataFrame(columns=_OBSERVATION_SUMMARY_COLUMNS)
        result.attrs.update(
            coordinate=str(coord).upper(),
            units="R_E",
        )
        return result

    _, _, parent_frames = pdf.split(large_gap.total_seconds())
    ena_times = ena["Datetime"] if len(ena) else pd.DatetimeIndex([])
    rows = []
    observation_index = 1
    for parent_index, parent in enumerate(parent_frames, start=1):
        _, _, subframes = parent.split(small_gap.total_seconds())
        for sub_index, subframe in enumerate(subframes, start=1):
            if len(subframe) == 0:
                continue
            sub_times = subframe["Datetime"]
            start_time = sub_times.iloc[0]
            end_time = sub_times.iloc[-1]
            duration = float((end_time - start_time).total_seconds())
            geometry = gut.pointing_geometry(
                subframe, probe=prob.upper(), coord=coord, forward=True
            )

            px = geometry["xz_x"] if "xz_x" in geometry else np.array([])
            pz = geometry["xz_z"] if "xz_z" in geometry else np.array([])
            px_mean = np.nanmean(px)
            pz_mean = np.nanmean(pz)
            look_columns = ["look_x", "look_y", "look_z"]
            if all(column in geometry for column in look_columns):
                look = geometry.loc[:, look_columns].to_numpy(dtype=float)
                look = look[np.isfinite(look).all(axis=1)]
            else:
                look = np.empty((0, 3))
            if len(look):
                mean_look = np.mean(look, axis=0)
                mean_look_norm = np.linalg.norm(mean_look)
                sun_earth_angle = (
                    float(
                        np.degrees(
                            np.arccos(
                                np.clip(mean_look[0] / mean_look_norm, -1.0, 1.0)
                            )
                        )
                    )
                    if mean_look_norm > 0
                    else np.nan
                )
            else:
                sun_earth_angle = np.nan

            if sun_earth_angle < 10:
                target_mode = "earth"
            elif 15 <= sun_earth_angle <= 40:
                target_mode = "mid_tail"
            elif sun_earth_angle > 120:
                target_mode = "far_tail"
            else:
                target_mode = "unknown"

            in_interval = (
                (ena_times >= start_time) & (ena_times <= end_time)
                if len(ena_times)
                else np.zeros(0, dtype=bool)
            )
            ena_count = int(np.count_nonzero(in_interval))
            count_rate = ena_count / duration if duration > 0 else np.nan
            rows.append(
                {
                    "observation_index": observation_index,
                    "parent_l_level_index": parent_index,
                    "sub_observation_index": sub_index,
                    "probe": prob.upper(),
                    "start_time": start_time,
                    "end_time": end_time,
                    "duration_seconds": duration,
                    "spacecraft_x_mean": np.nanmean(geometry["spacecraft_x"]),
                    "spacecraft_y_mean": np.nanmean(geometry["spacecraft_y"]),
                    "spacecraft_z_mean": np.nanmean(geometry["spacecraft_z"]),
                    "pointing_x_min": np.nanmin(px),
                    "pointing_x_max": np.nanmax(px),
                    "pointing_x_mean": px_mean,
                    "pointing_z_min": np.nanmin(pz),
                    "pointing_z_max": np.nanmax(pz),
                    "pointing_z_mean": pz_mean,
                    "sun_earth_angle_deg": sun_earth_angle,
                    "target_mode": target_mode,
                    "aiming": target_mode,
                    "ena_product": f"{ena_product}_{prob}",
                    "ena_total_counts": ena_count,
                    "ena_counts_per_second": count_rate,
                }
            )
            observation_index += 1
    result = pd.DataFrame(rows, columns=_OBSERVATION_SUMMARY_COLUMNS)
    result.attrs.update(
        coordinate=str(coord).upper(),
        units="R_E",
        large_gap=str(large_gap),
        small_gap=str(small_gap),
    )
    return result


def _save_observation_summary(summary, output_path):
    from openpyxl import load_workbook
    from openpyxl.styles import Font, PatternFill
    from openpyxl.utils import get_column_letter

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
        summary.to_excel(writer, sheet_name="observations", index=False)

    workbook = load_workbook(output_path)
    sheet = workbook["observations"]
    sheet.freeze_panes = "A2"
    sheet.auto_filter.ref = sheet.dimensions
    header_fill = PatternFill("solid", fgColor="1F4E78")
    for cell in sheet[1]:
        cell.fill = header_fill
        cell.font = Font(color="FFFFFF", bold=True)

    palette = ("E2F0D9", "DDEBF7", "FCE4D6", "E4DFEC", "FFF2CC")
    parent_column = next(
        (cell.column for cell in sheet[1] if cell.value == "parent_l_level_index"),
        None,
    )
    if parent_column is not None:
        for row in range(2, sheet.max_row + 1):
            parent = sheet.cell(row=row, column=parent_column).value
            if parent is None:
                continue
            fill = PatternFill("solid", fgColor=palette[(int(parent) - 1) % len(palette)])
            for column in range(1, sheet.max_column + 1):
                sheet.cell(row=row, column=column).fill = fill

    headers = {cell.value: cell.column for cell in sheet[1]}
    for name in ("start_time", "end_time"):
        column = headers.get(name)
        if column is not None:
            for row in range(2, sheet.max_row + 1):
                sheet.cell(row=row, column=column).number_format = "yyyy-mm-dd hh:mm:ss"
    for column_cells in sheet.columns:
        values = [str(cell.value) for cell in column_cells if cell.value is not None]
        width = min(max([len(value) for value in values] + [10]) + 2, 32)
        sheet.column_dimensions[get_column_letter(column_cells[0].column)].width = width
    workbook.save(output_path)
    return output_path


def observation_summary(
    prob,
    time,
    *,
    output_dir="reports",
    large_gap="1h",
    small_gap="5min",
    coord="GSM",
    ena_product="l2g",
):
    """Generate ``observation_summary.xlsx`` using :func:`report` time rules.

    ``time`` accepts one UTC day, an inclusive date-only range, one exact
    clock-precision interval, or a sequence of discrete UTC days. The returned
    DataFrame stores the workbook path in ``attrs['output_path']``.
    Relative output_dir values are resolved under the configured result root.
    """
    prob = normalize_probe(prob)
    periods = normalize_report_times(time)
    output_dir = resolve_output_path(output_dir)
    intervals = [list(period) for period in periods]
    pdf = load_gena(trange=intervals, product="pdf")
    ena = load_gena(
        trange=intervals,
        product=f"{ena_product}_{prob}",
    )
    summary = build_observation_summary(
        pdf,
        ena,
        prob,
        large_gap=large_gap,
        small_gap=small_gap,
        coord=coord,
        ena_product=ena_product,
    )
    output_path = _save_observation_summary(
        summary, output_dir / "observation_summary.xlsx"
    )
    summary.attrs["output_path"] = output_path
    return summary



def dump( l2df, prob="A", coord="GSM", *, extra_columns=("Energy",), normalize_ray=True, output_path=None,):
    """Build a tabular L2 geometry product from a GENA L2 DataFrame."""
    prob = prob.upper()
    pos = get_spacecraft_pos(l2df, coord=coord)
    ray_obj = get_ray(l2df, prob=prob, coord=coord)
    ray = np.asarray(ray_obj[coord], dtype=float)
    if normalize_ray:
        norm = np.linalg.norm(ray, axis=1)
        ray = ray / np.where(norm == 0, np.nan, norm)[:, np.newaxis]

    frame_xyz = get_frame(l2df, prob=prob, coord=coord)
    out = pd.DataFrame(
        {
            "Timestamp": l2df.loc[:, "Timestamp"].to_numpy(),
            f"QQ_X_{coord}": pos[coord][:, 0],
            f"QQ_Y_{coord}": pos[coord][:, 1],
            f"QQ_Z_{coord}": pos[coord][:, 2],
            f"XAXIS_X_{coord}": frame_xyz[0][coord][:, 0],
            f"XAXIS_Y_{coord}": frame_xyz[0][coord][:, 1],
            f"XAXIS_Z_{coord}": frame_xyz[0][coord][:, 2],
            f"YAXIS_X_{coord}": frame_xyz[1][coord][:, 0],
            f"YAXIS_Y_{coord}": frame_xyz[1][coord][:, 1],
            f"YAXIS_Z_{coord}": frame_xyz[1][coord][:, 2],
            f"ZAXIS_X_{coord}": frame_xyz[2][coord][:, 0],
            f"ZAXIS_Y_{coord}": frame_xyz[2][coord][:, 1],
            f"ZAXIS_Z_{coord}": frame_xyz[2][coord][:, 2],
            f"ENA_X_{coord}": ray[:, 0],
            f"ENA_Y_{coord}": ray[:, 1],
            f"ENA_Z_{coord}": ray[:, 2],
        }
    )
    for column in extra_columns:
        if column in l2df.columns:
            out[column] = l2df.loc[:, column].to_numpy()
    if output_path is not None:
        destination = resolve_output_path(output_path)
        destination.parent.mkdir(parents=True, exist_ok=True)
        out.to_csv(destination, index=False)
    return out


def load_products(period, prob):
    """Load all products needed by one daily or exact-interval report."""
    interval = list(period)
    product_names = {
        "pdf": "pdf",
        "l1a": f"l1a_{prob}",
        "l2": f"l2_{prob}",
        "l2s": f"l2s_{prob}",
        "l2g": f"l2g_{prob}",
    }
    return {
        name: load_gena(trange=interval, product=product)
        for name, product in product_names.items()
    }


def batch_plot(
    period,
    prob,
    products,
    *,
    period_label,
    daily=False,
    window="5min",
    max_engineering_samples=50_000,
    environment_sources=None,
    environment_cols=("dst", "al", "flow_speed", "BZ_GSM", "BY_GSM"),
):
    """Build every available static figure for one report interval."""
    figures = {}

    fig, ax = plt.subplots(figsize=(8.27, 11.69))
    ax.axis("off")
    ax.text(
        0.5,
        0.92,
        "GENA Daily Observation Report" if daily else "GENA Observation Report",
        ha="center",
        va="center",
        fontsize=20,
        weight="bold",
    )
    ax.text(
        0.5,
        0.865,
        f"{period_label} · Probe {prob.upper()}",
        ha="center",
        va="center",
        fontsize=14,
    )
    labels = [
        "Engineering (PDF)",
        "L1A events",
        "L2 events",
        "L2 soft events",
        "L2 golden events",
    ]
    keys = ["pdf", "l1a", "l2", "l2s", "l2g"]
    rows = []
    for label, key in zip(labels, keys):
        frame = products[key]
        coverage = "—"
        if not frame.empty:
            coverage = (
                f"{frame['Datetime'].iloc[0]:%H:%M:%S}–"
                f"{frame['Datetime'].iloc[-1]:%H:%M:%S}"
            )
        rows.append([label, f"{len(frame):,}", coverage])
    table = ax.table(
        cellText=rows,
        colLabels=["Product", "Rows", "UTC coverage"],
        colWidths=[0.42, 0.18, 0.32],
        cellLoc="left",
        loc="center",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(11)
    table.scale(1, 1.6)
    ax.text(
        0.06,
        0.14,
        "Products are loaded from GENA daily Parquet files and clipped to the requested report interval. Counts are event rows, not exposure-normalized flux.",
        ha="left",
        va="top",
        fontsize=9,
        wrap=True,
    )
    figures["summary"] = fig

    level_keys = ["l1a", "l2", "l2s", "l2g"]
    fig, ax = gena_vis.counts(
        [products[key] for key in level_keys],
        labels=["L1A", "L2", "L2 soft", "L2 golden"],
        window=window,
        vmax=None,
    )
    ax.set_title("Data-level event counts")
    ax.set_xlabel("Time (UTC)")
    ax.grid(True, alpha=0.25)
    figures["levels_count"] = fig

    engineering = products["pdf"]
    if not engineering.empty:
        if len(engineering) > max_engineering_samples:
            stride = int(np.ceil(len(engineering) / max_engineering_samples))
            engineering = ParameterDataFrame(
                engineering.iloc[::stride].reset_index(drop=True)
            )
        fig, _ = gena_vis.parameters(engineering, prob)
        fig.suptitle(f"Probe {prob.upper()} engineering parameters")
        figures["parameters"] = fig

        pose_columns = {
            "Timestamp", "MJ2000_Star_X", "MJ2000_Star_Y", "MJ2000_Star_Z",
            "MJ2000_q1", "MJ2000_q2", "MJ2000_q3",
        }
        if pose_columns.issubset(engineering.columns):
            fig, _ = gut.orbit(engineering, probe=prob, coord="GSM", xlim=period)
            fig.suptitle(f"Probe {prob.upper()} QQ orbit and FOV-plane intersections")
            figures["orbit_pointing"] = fig

    if environment_sources:
        fig, _ = env(
            *environment_sources,
            cols=environment_cols,
            xlim=period,
        )
        fig.suptitle("Solar-wind and geomagnetic environment")
        figures["environment"] = fig

    l1a = products["l1a"]
    for artifact, plotter, title in (
        (
            "start_signal_distribution",
            gena_vis.start_signal,
            "start-signal distribution",
        ),
        (
            "end_signal_distribution",
            gena_vis.end_signal,
            "end-signal distribution",
        ),
    ):
        if l1a.empty:
            break
        fig, axes = plotter(l1a, win=window)
        fig.suptitle(f"Probe {prob.upper()} {title}")
        locator = mdates.AutoDateLocator(minticks=4, maxticks=8)
        axes[-1].xaxis.set_major_locator(locator)
        axes[-1].xaxis.set_major_formatter(
            mdates.ConciseDateFormatter(locator)
        )
        figures[artifact] = fig

    joint_columns = {
        "TX_1s", "TX_2s", "TY_1s", "TY_2s",
        "TX_1e", "TX_2e", "TY_1e", "TY_2e",
    }
    if not l1a.empty and joint_columns.issubset(l1a.columns):
        fig, _ = gena_vis.joint(l1a)
        fig.suptitle(f"Probe {prob.upper()} joint signal distribution")
        figures["joint_signal_distribution"] = fig

    l2 = products["l2"]
    if not l2.empty:
        fig, _ = gena_vis.energy(l2)
        fig.suptitle(f"Probe {prob.upper()} energy distribution")
        figures["hist_E"] = fig

        for layer in (1, 2):
            fig, _ = gena_vis.layer(l2, layer=layer)
            fig.suptitle(f"Probe {prob.upper()} detector layer {layer}")
            figures[f"hist_Layer{layer}"] = fig

        fig, _ = gena_vis.fov(l2)
        fig.suptitle(f"Probe {prob.upper()} field of view")
        figures["hist_FOV"] = fig
    return figures


def save_figures(figures, prob, *, output_dir, dpi=_REPORT_DPI):
    """Save a static-figure batch directly into the supplied report directory."""
    report_dir = Path(output_dir)
    artifacts = {}
    for artifact, figure in figures.items():
        path = report_dir / _FIGURE_FILENAMES[artifact].format(prob=prob)
        artifacts[artifact] = savefig(
            figure,
            path,
            dpi=dpi,
            tight=True,
        )
    return artifacts


def _save_geometry(period, *, output_path, geometry_options):
    options = dict(geometry_options)
    max_frames = options.pop("max_frames", 48)
    include_plotlyjs = options.pop("include_plotlyjs", True)
    overwrite = options.pop("overwrite", True)
    output_path = Path(output_path)
    if output_path.exists() and not overwrite:
        return output_path
    context = gena_v3d.prepare_scene(
        list(period),
        **options,
    )
    figure = context.make_timerange(max_frames=max_frames)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = output_path.with_name(f".{output_path.name}.tmp")
    try:
        figure.write_html(
            temporary,
            include_plotlyjs=include_plotlyjs,
            full_html=True,
            auto_play=False,
        )
        temporary.replace(output_path)
    finally:
        temporary.unlink(missing_ok=True)
    return output_path


def report(
    prob,
    time,
    *,
    output_dir="reports",
    include_animation=False,
    include_geometry=False,
    window="5min",
    max_engineering_samples=50_000,
    geometry_kwargs=None,
    environment_sources=None,
    environment_cols=("dst", "al", "flow_speed", "BZ_GSM", "BY_GSM"),
):
    """Generate GENA observation artifacts for automatic time input.

    ``time`` automatically distinguishes a single UTC day, a two-endpoint date
    range, one exact clock-precision interval, and a sequence of discrete days.
    Static plots are independent 300 dpi JPEG files. Set
    ``include_animation=True`` to add an MP4 generated by :mod:`plots.ani`, and
    ``include_geometry=True`` to add a Plotly HTML scene generated by
    :mod:`plots.v3d`. Pass already-loaded OMNI/index frames through
    ``environment_sources`` to add an environment figure without implicit I/O.
    Relative output_dir values are resolved under the configured result root.
    """
    prob = normalize_probe(prob)
    periods = normalize_report_times(time)

    if not isinstance(max_engineering_samples, int) or max_engineering_samples < 1:
        raise ValueError("max_engineering_samples must be a positive integer")
    try:
        pd.Timedelta(window)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"invalid report window: {window!r}") from exc

    output_dir = resolve_output_path(output_dir)

    geometry_options = {}
    if include_geometry:
        geometry_options = dict(geometry_kwargs or {})
        reserved = {"trange", "output_dir"}
        duplicated = reserved.intersection(geometry_options)
        if duplicated:
            names = ", ".join(sorted(duplicated))
            raise ValueError(
                f"geometry_kwargs must not override report-owned options: {names}"
            )

    results = {}
    for period in periods:
        start, stop = period
        daily = (
            start == start.floor("D")
            and stop == start + pd.Timedelta(days=1) - pd.Timedelta(nanoseconds=1)
        )
        if daily:
            key = start.strftime("%Y%m%d")
            label = start.strftime("%Y-%m-%d UTC")
        else:
            key = f"{start:%Y%m%dT%H%M%SZ}_{stop:%Y%m%dT%H%M%SZ}"
            if start.date() == stop.date():
                label = f"{start:%Y-%m-%d %H:%M:%S}–{stop:%H:%M:%S} UTC"
            else:
                label = f"{start:%Y-%m-%d %H:%M:%S}–{stop:%Y-%m-%d %H:%M:%S} UTC"
        report_dir = output_dir / key
        products = load_products(period, prob)
        figures = batch_plot(
            period,
            prob,
            products,
            period_label=label,
            daily=daily,
            window=window,
            max_engineering_samples=max_engineering_samples,
            environment_sources=environment_sources,
            environment_cols=environment_cols,
        )
        artifacts = save_figures(
            figures,
            prob,
            output_dir=report_dir,
        )
        if include_animation and not products["l2"].empty:
            animation_path = report_dir / f"{prob}_fov.mp4"
            animation_path = gena_ani.dumpFov(
                products["l2"],
                animation_path,
                window=window,
            )
            if animation_path is not None:
                artifacts["fov_animation"] = animation_path
        if include_geometry:
            artifacts["geometry"] = _save_geometry(
                period,
                output_path=report_dir / "3d.html",
                geometry_options=geometry_options,
            )
        results[key] = artifacts
    return results
