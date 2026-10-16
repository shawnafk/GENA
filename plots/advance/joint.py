"""Cross-mission plotting of preloaded data; no loading or synchronization."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from inspect import signature
from numbers import Integral
from typing import Callable

import matplotlib.dates as mdates
import numpy as np
import pandas as pd

from common.frame import TimeSeriesDataFrame
from libs.regions import REGION_LABELS
from common.time import Time
from plots.vis import axes_for_columns as _axes_for_columns, column_groups as _column_groups

__all__ = ["env", "Panel", "overview", "plot_region_code"]


REGION_COLORS = {
    0: "#d4d4d4",  # neutral gray: no reliable classification
    1: "#f2c94c",  # yellow: solar wind
    2: "#56a3d9",  # blue: magnetosheath
    3: "#8b5fbf",  # purple: mantle
    4: "#6bbf59",  # green: magnetotail lobe
    5: "#f2994a",  # orange: plasma sheet boundary layer
    6: "#d94b4b",  # red: plasma sheet
}


def plot_region_code(
    ax, result: TimeSeriesDataFrame, label: str | None = None, *, legend: bool = False
):
    """Plot a categorical magnetospheric-region time strip.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axis receiving the one-row region strip.
    result : TimeSeriesDataFrame
        Classification result containing canonical ``Datetime``/``Timestamp``
        columns and an integer ``region_code`` column from 0 through 6. The
        frame returned by :func:`classify_liuzzo` can be passed directly.
    label : str, optional
        Text shown at the center of the y axis, for example ``"ARTEMIS 1"``.
        When omitted, the y tick is left blank.
    legend : bool, default False
        Add a legend for the region codes present in *result*.

    Returns
    -------
    matplotlib.collections.QuadMesh
        Artist created by ``Axes.pcolormesh``.

    Notes
    -----
    Region codes are categorical. ``ListedColormap`` and ``BoundaryNorm`` are
    used instead of a continuous colormap. The mesh is rasterized so a long,
    high-cadence interval does not produce an unnecessarily large PDF.
    """
    import matplotlib.dates as mdates
    from matplotlib.colors import BoundaryNorm, ListedColormap
    from matplotlib.patches import Patch

    # Sort time and code together. THEMIS daily products can contain non-monotonic records, while
    # pcolormesh requires ordered x coordinates.
    if not isinstance(result, TimeSeriesDataFrame):
        raise TypeError("result must be a TimeSeriesDataFrame")
    if result.empty:
        raise ValueError("Cannot plot an empty region time series")
    time = pd.DatetimeIndex(result["Datetime"])
    code = np.asarray(result["region_code"], dtype=int)
    order = np.argsort(time.asi8)
    time = time[order]
    code = code[order]

    # Classification timestamps are cell centers. Convert them to Matplotlib
    # date numbers and use neighboring midpoints as pcolormesh cell edges.
    x = mdates.date2num(time.to_pydatetime())
    if len(x) == 1:
        x_edges = np.array([x[0] - 0.5 / 86400.0, x[0] + 0.5 / 86400.0])
    else:
        middle = (x[:-1] + x[1:]) / 2.0
        x_edges = np.r_[x[0] - (middle[0] - x[0]), middle,
                        x[-1] + (x[-1] - middle[-1])]

    # Half-integer boundaries assign every integer region code to one fixed
    # bin: [-0.5, 0.5) is code 0, [0.5, 1.5) is code 1, and so forth.
    cmap = ListedColormap([REGION_COLORS[i] for i in range(7)])
    norm = BoundaryNorm(np.arange(-0.5, 7.5, 1.0), cmap.N)
    mesh = ax.pcolormesh(
        x_edges,
        [0.0, 1.0],
        code[None, :],
        cmap=cmap,
        norm=norm,
        shading="flat",
        rasterized=True,
    )

    ax.set_xlim(x_edges[0], x_edges[-1])
    ax.set_ylim(0.0, 1.0)
    ax.set_yticks([0.5], [label] if label else [""])

    if legend:
        # Only show categories that occur in the plotted interval.
        handles = [
            Patch(facecolor=REGION_COLORS[i], edgecolor="none", label=REGION_LABELS[i])
            for i in np.unique(code)
        ]
        ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.35),
                  ncol=len(handles), frameon=False)

    return mesh


@dataclass
class Panel:
    """One plotting call and its arguments in a joint time-series figure.

    ``plotter`` must accept ``ax`` as a keyword; ``f`` is passed when explicitly
    present in its signature. For a multi-panel product (e.g. MMS FPI), set
    ``nrows`` to its required axes count. Arguments contain already-loaded data.
    Use only time-axis plots in overview; spatial/energy-axis plots stay standalone.
    """

    plotter: Callable
    args: tuple = ()
    kwargs: dict = field(default_factory=dict)
    nrows: int = 1

    def __post_init__(self):
        if not callable(self.plotter):
            raise TypeError("plotter must be callable")
        if isinstance(self.nrows, bool) or not isinstance(self.nrows, Integral) or self.nrows < 1:
            raise ValueError("nrows must be a positive integer")
        self.args = tuple(self.args)
        self.kwargs = dict(self.kwargs)
        if {"f", "ax"} & self.kwargs.keys():
            raise ValueError("overview owns f and ax; do not put them in Panel.kwargs")


def overview(panels, *, trange=None, f=None, ax=None, figsize=None):
    """Compose mission plotting calls on a shared UTC time axis.

    Returns ``(figure, axes)`` with a flat list of data axes (excluding colorbars).
    Data retain their native timestamps. ``trange`` controls display limits only;
    no loading, resampling, or derived science calculation is performed here.
    Caller-supplied axes must belong to one figure. Each Panel's nrows must match
    the product function's layout; plotters must reuse the supplied axes.
    """
    panels = list(panels)
    if not panels or any(not isinstance(panel, Panel) for panel in panels):
        raise ValueError("panels must be a nonempty sequence of Panel objects")
    limits = None
    if trange is not None:
        limits = Time(trange).datetime
        if isinstance(limits, pd.Timestamp) or len(limits) != 2:
            raise ValueError("trange must contain exactly two time values")
        if any(pd.isna(value) for value in limits) or limits[0] >= limits[1]:
            raise ValueError("trange must have finite times with start before end")

    # Validate call signatures before creating or modifying a figure.
    calls = []
    for panel in panels:
        sig = signature(panel.plotter)
        keywords = dict(panel.kwargs, ax=None)
        if "f" in sig.parameters:
            keywords["f"] = None
        sig.bind(*panel.args, **keywords)
        calls.append(keywords)
    if ax is not None:
        supplied = list(np.asarray(ax, dtype=object).reshape(-1))
        if supplied:
            owner = f if f is not None else supplied[0].get_figure()
            if any(axis.get_figure() is not owner for axis in supplied):
                raise ValueError("all axes must belong to the supplied figure")
    created_figure = f is None and ax is None
    f, axes = _axes_for_columns(
        sum(panel.nrows for panel in panels), f=f, ax=ax, figsize=figsize
    )
    try:
        for axis in axes[1:]:
            if not axis.get_shared_x_axes().joined(axes[0], axis):
                axis.sharex(axes[0])
        offset = 0
        for panel, keywords in zip(panels, calls):
            selected = axes[offset:offset + panel.nrows]
            keywords["ax"] = selected[0] if panel.nrows == 1 else selected
            if "f" in keywords:
                keywords["f"] = f
            panel.plotter(*panel.args, **keywords)
            offset += panel.nrows
        if limits is not None:
            axes[-1].set_xlim(limits)
        for axis in axes[:-1]:
            axis.set_xlabel("")
            axis.tick_params(axis="x", labelbottom=False)
        locator = mdates.AutoDateLocator(minticks=4, maxticks=10)
        axes[-1].xaxis.set_major_locator(locator)
        axes[-1].xaxis.set_major_formatter(mdates.ConciseDateFormatter(locator))
        axes[-1].set_xlabel("UTC")
    except Exception:
        if created_figure:
            from matplotlib import pyplot as plt
            plt.close(f)
        raise
    return f, axes


_ENV_LABELS = {
    "flow_speed": r"$V_{SW}$ (km/s)",
    "proton_density": r"$n_{SW}$ (cm$^{-3}$)",
    "Pressure": r"$P_{dyn}$ (nPa)",
    "BZ_GSM": r"$B_z^{GSM}$ (nT)",
    "BY_GSM": r"$B_y^{GSM}$ (nT)",
    "dst": "Dst (nT)",
    "ae": "AE (nT)",
    "al": "AL (nT)",
    "au": "AU (nT)",
    "kp": "Kp",
}
_ENV_COLORS = {"ae": "tab:red", "al": "tab:green", "au": "tab:blue"}

def _source_series(source, name):
    if source is None:
        return None
    if isinstance(source, pd.DataFrame) and name in source.columns:
        return source["Datetime"], source[name].to_numpy()
    if isinstance(source, Mapping) and name in source:
        value = source[name]
        if isinstance(value, Mapping) and "x" in value and "y" in value:
            return pd.to_datetime(value["x"], utc=True).tz_convert(None), np.asarray(value["y"])
    return None


def env(
    *sources,
    cols=("dst", "al", "flow_speed", "BZ_GSM", "BY_GSM"),
    f=None,
    ax=None,
    labels=None,
    ylabels=None,
    colors=None,
    figsize=None,
    xlim=None,
    grid=True,
):
    """Plot selected solar-wind and geomagnetic variables on dynamic axes.

    Data loading remains the caller's responsibility. Each entry in ``cols``
    creates one panel. A tuple/list groups several variables on the same panel,
    for example ``cols=["dst", ("ae", "al", "au"), "flow_speed"]``.
    Every variable is resolved from the first supplied DataFrame (or legacy
    ``{"x", "y"}`` variable dictionary) that contains it.
    """
    groups = _column_groups(cols)
    f, axes = _axes_for_columns(
        len(groups), f=f, ax=ax, figsize=figsize, sharex=True
    )
    labels = dict(labels or {})
    ylabels = dict(ylabels or {})
    colors = dict(colors or {})

    for axis, group in zip(axes, groups):
        plotted = []
        for name in group:
            series = next(
                (candidate for source in sources if (candidate := _source_series(source, name)) is not None),
                None,
            )
            if series is None:
                continue
            times, values = series
            axis.plot(
                times,
                values,
                color=colors.get(name, _ENV_COLORS.get(name, "black")),
                label=labels.get(name, name),
            )
            plotted.append(name)
        if not plotted:
            axis.text(
                0.5,
                0.5,
                f"No data: {', '.join(group)}",
                ha="center",
                va="center",
                transform=axis.transAxes,
            )
        elif len(group) > 1:
            axis.legend(loc="upper right", ncol=min(3, len(plotted)))
        default_ylabel = (
            _ENV_LABELS.get(group[0], group[0])
            if len(group) == 1
            else "/".join(labels.get(name, name.upper()) for name in group)
        )
        axis.set_ylabel(ylabels.get(group[0], default_ylabel))
        if any(name in {"BZ_GSM", "BY_GSM", "dst", "ae", "al", "au"} for name in group):
            axis.axhline(0, color="0.5", linestyle=":", linewidth=0.8)
        if xlim is not None:
            axis.set_xlim(xlim)
        if grid:
            axis.grid(alpha=0.2, linewidth=0.5)

    locator = mdates.AutoDateLocator(minticks=4, maxticks=10)
    axes[-1].xaxis.set_major_locator(locator)
    axes[-1].xaxis.set_major_formatter(mdates.ConciseDateFormatter(locator))
    axes[-1].set_xlabel("UTC")
    f.autofmt_xdate()
    return f, axes[0] if len(axes) == 1 else axes
