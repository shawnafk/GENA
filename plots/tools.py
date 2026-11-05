import numpy as np

from matplotlib import pyplot as plt
import matplotlib.dates as mdates
import matplotlib.ticker as mticker
import matplotlib.transforms as mtransforms
from matplotlib.patches import Rectangle
from matplotlib.axes import Axes
from lunardate import LunarDate
import datetime

import functools
from collections.abc import Sequence
import inspect
from pathlib import Path
from common.config import resolve_output_path

from libs import ephemeris as epms
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

__all__ = [
    "add_bar_right",
    "add_bar_top",
    "auto_loc",
    "axes_for_columns",
    "column_groups",
    "draw_rectangle",
    "format_lunar_date",
    "format_pos_date",
    "pass_axes",
    "savefig",
]


def savefig( figure, path, *, dpi=300, tight=True, close=True, facecolor="None", **kwargs,):
    """Save a Matplotlib figure with shared repository defaults.

    Parent directories are created automatically, the format is inferred from
    the filename unless supplied in ``kwargs``, and the output path is returned.
    Relative filenames use common.config's runtime result directory.
    """
    destination = resolve_output_path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    if tight:
        figure.tight_layout()
    save_kwargs = {
        "dpi": dpi,
        "bbox_inches": "tight" if tight else None,
        "facecolor": facecolor,
    }
    save_kwargs.update(kwargs)
    figure.savefig(destination, **save_kwargs)
    if close:
        plt.close(figure)
    return destination


def draw_rectangle(ax, width, height, *, center=(0, 0), **kwargs):
    """Draw a detector-style rectangle centered on *center*.

    The helper replaces the former repository-external ``myplot`` dependency.
    Extra keyword arguments are forwarded to ``matplotlib.patches.Rectangle``.
    """
    style = {"fill": False, "edgecolor": "black"}
    style.update(kwargs)
    x0 = center[0] - width / 2
    y0 = center[1] - height / 2
    patch = Rectangle((x0, y0), width, height, **style)
    ax.add_patch(patch)
    return patch

def pass_axes(subplot_spec=None):
    """Supply figure/axes for a fixed-layout plotting function.

    ``subplot_spec`` is (nrows, ncols, sharex, sharey), defaulting to one
    subplot. Reuse supplied axes, create axes in a supplied figure, or create
    both when neither is supplied. Explicit None is treated as omission.
    The wrapped function receives one Axes for a single subplot, otherwise a
    flat list in row-major order. Existing axes must match the requested count
    and belong to the supplied figure (or to one inferred common figure).

    Plotters must declare f and ax parameters. Prefer keyword arguments at
    call sites; existing positional calls remain supported.
    """
    nrows, ncols, sharex, sharey = (1, 1, False, False) if subplot_spec is None else subplot_spec
    if any(isinstance(n, (bool, np.bool_)) or not isinstance(n, (int, np.integer))
           or n < 1 for n in (nrows, ncols)):
        raise ValueError("subplot rows and columns must be positive integers")
    count = nrows * ncols

    def decorator(func):
        signature = inspect.signature(func)
        if not {"f", "ax"}.issubset(signature.parameters):
            raise TypeError("pass_axes requires explicit f and ax parameters")

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            bound = signature.bind(*args, **kwargs)
            bound.apply_defaults()
            f, ax = bound.arguments["f"], bound.arguments["ax"]
            if ax is None:
                if f is None:
                    f = plt.figure()
                ax = f.subplots(nrows, ncols, sharex=sharex, sharey=sharey,
                                squeeze=False)
            axes = list(np.asarray(ax, dtype=object).reshape(-1))
            if len(axes) != count:
                raise ValueError(f"Expected {count} axes, got {len(axes)}")
            if not all(isinstance(axis, Axes) for axis in axes):
                raise TypeError("ax must contain Matplotlib Axes objects")
            if f is None:
                f = axes[0].get_figure()
            if any(axis.get_figure() is not f for axis in axes):
                raise ValueError("all axes must belong to the same supplied or inferred figure")
            bound.arguments.update(f=f, ax=axes[0] if count == 1 else axes)
            return func(*bound.args, **bound.kwargs)

        return wrapper
    return decorator


def axes_for_columns(n_panels, *, f=None, ax=None, figsize=None, sharex=True):
    """Create or validate a dynamic vertical axes layout."""
    if n_panels < 1:
        raise ValueError("at least one panel is required")
    if ax is None:
        size = figsize or (9, max(2.2 * n_panels, 3.0))
        if f is None:
            f, created = plt.subplots(
                n_panels, 1, figsize=size, sharex=sharex, squeeze=False
            )
        else:
            created = f.subplots(n_panels, 1, sharex=sharex, squeeze=False)
        axes = list(np.asarray(created, dtype=object).reshape(-1))
    else:
        axes = list(np.asarray(ax, dtype=object).reshape(-1))
        if len(axes) != n_panels:
            raise ValueError(
                f"cols requires {n_panels} axes, but {len(axes)} were supplied"
            )
        if f is None:
            f = axes[0].get_figure()
    return f, axes


def column_groups(cols):
    if isinstance(cols, str):
        cols = [cols]
    groups = []
    for item in cols:
        if isinstance(item, str):
            groups.append((item,))
        elif isinstance(item, Sequence) and not isinstance(item, (str, bytes)):
            group = tuple(str(name) for name in item)
            if not group:
                raise ValueError("cols cannot contain an empty group")
            groups.append(group)
        else:
            raise TypeError("each cols item must be a column name or a sequence")
    if not groups:
        raise ValueError("cols must contain at least one column")
    return groups


#! this is for plotting function defined like f(args, f=None,ax=None, **kwargs)
#! it will creates require f and ax, by assign subplots_kw or defaults
#! it also take f and ax from outside, if you want to merge multiple plots

#matplotlib figure utils
def add_bar_top( ax, mesh, *, label='counts', fraction=0.1, pad=0.05, aspect=None, shrink=1.0, labelsize=None, tick_labelsize=10, labelpad=4, label_location="left",):
    """Attach a horizontal colorbar above *ax*.

    The colorbar is positioned relative to *ax* without resizing it.  With the
    default ``shrink=1``, the colorbar has exactly the same width as *ax*.
    ``aspect`` remains accepted for backward compatibility; the inset geometry
    is controlled directly by ``fraction`` instead.
    """
    if fraction <= 0:
        raise ValueError("fraction must be positive")
    if not 0 < shrink <= 1:
        raise ValueError("shrink must be in the interval (0, 1]")

    # ``figure.colorbar(..., aspect=...)`` can shorten a horizontal colorbar
    # because its aspect is resolved in physical figure units.  Anchor an
    # inset axis directly to the parent instead so the two widths stay locked.
    _ = aspect  # Retained in the public signature for existing callers.
    x_offset = (1.0 - shrink) / 2.0
    cax = inset_axes(
        ax,
        width=f"{100.0 * shrink:g}%",
        height=f"{100.0 * fraction:g}%",
        loc="lower left",
        bbox_to_anchor=(x_offset, 1.0 + pad, 1.0, 1.0),
        bbox_transform=ax.transAxes,
        borderpad=0,
    )

    cbar = ax.figure.colorbar(
        mesh,
        cax=cax,
        orientation="horizontal",
        location="top",
    )
    cbar.set_label(
        label,
        fontsize=labelsize,
        labelpad=labelpad,
        loc=label_location,
    )
    cbar.ax.tick_params(labelsize=tick_labelsize)
    return cbar

def add_bar_right( ax, mesh, *, label="counts", width="1.8%", gap=0.015, labelpad=6, tick_labelsize=10,):
    """在坐标轴右侧添加 colorbar，不改变原坐标轴宽度。"""
    cax = inset_axes(
        ax,
        width=width,
        height="100%",
        loc="lower left",
        bbox_to_anchor=(1.0 + gap, 0.0, 1.0, 1.0),
        bbox_transform=ax.transAxes,
        borderpad=0,
    )

    cbar = ax.figure.colorbar(
        mesh,
        cax=cax,
        orientation="vertical",
    )
    cbar.set_label(label, labelpad=labelpad)
    cbar.ax.tick_params(labelsize=tick_labelsize)
    return cbar

def auto_loc(ax):
    locator = mdates.AutoDateLocator(minticks=4, maxticks=12)
    minor_locator = mticker.AutoMinorLocator(n=4)
    ax.xaxis.set_minor_locator(minor_locator)
    ax.xaxis.set_major_locator(locator)
    
#format x axis (t axis)
def format_lunar_date(ax,xdates):
    lunar_dates = []
    for date in xdates:
        lunar = LunarDate.fromSolarDate(date.year, date.month, date.day)
        lunar_str = f"{lunar.year}-{lunar.month}-{lunar.day}"
        lunar_dates.append(lunar_str)
    def format_func(x, pos = None):
        # 将Matplotlib数值转为公历日期
        gregorian = mdates.num2date(x).strftime('%Y-%m-%d  %H:%M')
        # 查找对应的农历日期
        idx = int(round(x - mdates.date2num(xdates.iloc[0])))  # 获取索引
        if 0 <= idx < len(lunar_dates):
            return f"{gregorian}\n{lunar_dates[idx]}"
        else:
            return gregorian
    ax.xaxis.set_major_formatter(format_func)


def _sat_coord_formatter(x, pos_xyz):
    """Format one time tick using a precomputed satellite position.

    The caller supplies ``pos_xyz`` so this function only handles presentation.
    """
    dt = mdates.num2date(x, tz=datetime.timezone.utc)

    try:
        if pos_xyz is None:
            raise ValueError("Satellite position is unavailable")
        pos_xyz = np.asarray(pos_xyz, dtype=float)
        if pos_xyz.shape != (3,) or not np.all(np.isfinite(pos_xyz)):
            raise ValueError("Satellite position must contain three finite values")
        x_str = f"{pos_xyz[0]:.1f}"
        y_str = f"{pos_xyz[1]:.1f}"
        z_str = f"{pos_xyz[2]:.1f}"
    except (TypeError, ValueError):
        x_str, y_str, z_str = "N/A", "N/A", "N/A"

    time_str = dt.strftime('%d %H:%M')

    return f"{time_str}\n{x_str}\n{y_str}\n{z_str}"


def format_pos_date(ax, sat_name='MMS1', *, fontsize=12, max_ticks=8):
    """Show UTC time and Earth-centred GSM XYZ in Earth radii on an axis.

    Positions come from ``load.load_pos``. Names include THA--THE,
    ATMS1/ATMS2 (THEMIS B/C), MMS1--MMS4, GOES-N, and QUEQIAO/QQ/GENA.
    Queqiao's Moon-relative MJ2000 metres are translated to Earth-centred
    ECI2000 and then converted to GSM; other missions return GSM kilometres.

    Each tick position is cached in this formatter closure.  Zooming introduces
    new tick positions, which are interpolated once and then reused on redraw.
    Time and position values use four compact, aligned rows.  The first and last
    visible tick labels extend into the axes so they do not collide with the row
    headings or run outside the figure.  ``fontsize`` controls both the values
    and headings; ``max_ticks`` limits crowding on narrow multi-panel figures.
    """
    from load import load_pos
    from common.constants import RE
    from common.coord import coord
    from common.frame import read_variable

    name = str(sat_name).strip().upper()
    if name in {'QUEQIAO', 'QQ', 'GENA'}:
        mission, satellite_id = 'queqiao', None
    elif name in {'ATMS1', 'ATMS2'}:
        mission, satellite_id = 'themis', {'ATMS1': 'b', 'ATMS2': 'c'}[name]
    elif name in {'THA', 'THB', 'THC', 'THD', 'THE'}:
        mission, satellite_id = 'themis', name[-1].lower()
    elif name.startswith('MMS') and name[3:] in {'1', '2', '3', '4'}:
        mission, satellite_id = 'mms', name[3:]
    elif name.startswith('GOES') and name[4:].lstrip('-').isdigit():
        mission, satellite_id = 'goes', name[4:].lstrip('-')
    else:
        raise ValueError(f'Unsupported satellite name: {sat_name!r}')

    auto_loc(ax)
    ax.xaxis.set_major_locator(
        mdates.AutoDateLocator(minticks=4, maxticks=max_ticks)
    )
    position_cache = {}

    def format_func(x, _tick_index=None):
        if x not in position_cache:
            try:
                dt = mdates.num2date(x, tz=datetime.timezone.utc)
                pad = datetime.timedelta(hours=1)
                frame = load_pos(
                    mission, [dt - pad, dt + pad], satellite_id, coord='gsm',
                )
                if mission == 'queqiao':
                    columns = ['MJ2000_Star_X', 'MJ2000_Star_Y', 'MJ2000_Star_Z']
                    positions = frame.resample_at([dt])[columns].to_numpy(dtype=float)
                    earth_centred = positions + np.asarray(epms.moon([dt])) * 1000.0
                    orbit = coord([dt], earth_centred, 'ECI2000', 0)
                    orbit.convert('GSM')
                    positions = np.asarray(orbit['GSM'], dtype=float) / RE
                else:
                    variable = next(iter(frame.attrs['variables']), None)
                    if variable is None:
                        raise ValueError("No satellite position samples available")
                    positions = read_variable(frame.resample_at([dt]), variable)['values']
                    positions = positions * 1000.0 / RE
                if positions.shape != (1, 3):
                    raise ValueError(
                        "Satellite ephemeris must return one XYZ row per tick"
                    )
                position_cache[x] = positions[0]
            except (KeyError, OSError, TypeError, ValueError):
                position_cache[x] = None
        return _sat_coord_formatter(x, position_cache[x])

    ax.xaxis.set_major_formatter(mticker.FuncFormatter(format_func))
    tick_fontsize = fontsize
    tick_labels = ax.get_xticklabels()
    plt.setp(
        tick_labels,
        rotation=0,
        fontsize=tick_fontsize,
        ha='center',
        va='top',
        linespacing=1.2,
    )

    # Repeated calls should update, rather than duplicate, the row headings.
    header_gid = 'satellite-position-row-heading'
    for artist in list(ax.texts):
        if artist.get_gid() == header_gid:
            artist.remove()

    if tick_labels:
        # Use one multiline Text with the exact transform and typography of a
        # tick label.  Separate annotations drift because glyph metrics and
        # multiline line spacing are not equivalent to hand-tuned offsets.
        reference_label = tick_labels[0]
        heading_transform = mtransforms.blended_transform_factory(
            ax.transAxes,
            reference_label.get_transform(),
        )
        heading_transform = mtransforms.offset_copy(
            heading_transform,
            fig=ax.figure,
            x=-8,
            units='points',
        )
        heading = ax.text(
            0,
            reference_label.get_position()[1],
            'T\nX\nY\nZ',
            transform=heading_transform,
            ha='right',
            va=reference_label.get_va(),
            fontsize=tick_fontsize,
            linespacing=1.2,
            clip_on=False,
        )
        heading.set_gid(header_gid)

    def align_edge_ticklabels(_axis=None):
        labels = [label for label in ax.get_xticklabels() if label.get_visible()]
        plt.setp(labels, ha='center')
        if labels:
            labels[0].set_ha('left')
            labels[-1].set_ha('right')

    previous_callback = getattr(ax, '_format_pos_date_xlim_callback', None)
    if previous_callback is not None:
        ax.callbacks.disconnect(previous_callback)
    ax._format_pos_date_xlim_callback = ax.callbacks.connect(
        'xlim_changed', align_edge_ticklabels
    )
    align_edge_ticklabels()
