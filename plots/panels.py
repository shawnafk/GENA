"""Declarative single-axis panels built from the generic ``tplot`` primitive.

Many mission "plot functions" are only a variable lookup plus a ylabel, a log
scale, a zero line or a legend. :class:`PanelSpec` captures that metadata so
those wrappers do not need to exist as separate functions.
"""

from __future__ import annotations

from dataclasses import dataclass

from plots import tools
from plots.timeseries import tplot

__all__ = ["PanelSpec", "plot_panel", "plot_panels"]


@dataclass(frozen=True)
class PanelSpec:
    """Declarative description of one time-series axis.

    ``keys`` are passed to :func:`plots.tplot`; multiple keys share the axis.
    ``key_kwargs`` may be empty or contain one dict per key, allowing a scalar
    to be black while a following vector keeps its component colours.
    """

    keys: tuple[str, ...]
    ylabel: str | None = None
    yscale: str | None = None
    zero: bool = False
    legend: bool = False
    key_kwargs: tuple[dict, ...] = ()
    label: str | None = None

    def __post_init__(self):
        if not self.keys:
            raise ValueError("PanelSpec.keys must contain at least one variable")
        if self.yscale not in (None, "linear", "log"):
            raise ValueError("yscale must be None, 'linear' or 'log'")
        if self.key_kwargs and len(self.key_kwargs) != len(self.keys):
            raise ValueError("key_kwargs must be empty or match the number of keys")


@tools.pass_axes()
def plot_panel(data, spec: PanelSpec, *, trange=None, prefix=None, f=None, ax=None):
    """Draw one :class:`PanelSpec` on a single axis."""
    if not isinstance(spec, PanelSpec):
        raise TypeError("spec must be a PanelSpec")
    for index, key in enumerate(spec.keys):
        kwargs = dict(spec.key_kwargs[index]) if spec.key_kwargs else {}
        tplot(
            data,
            key,
            trange=trange,
            prefix=prefix,
            f=f,
            ax=ax,
            label=spec.label if len(spec.keys) == 1 else None,
            **kwargs,
        )
    if spec.yscale is not None:
        ax.set_yscale(spec.yscale)
    if spec.zero:
        ax.axhline(0, ls="--", color="k")
    if spec.ylabel:
        ax.set_ylabel(spec.ylabel)
    if spec.legend and ax.get_legend_handles_labels()[0]:
        ax.legend(frameon=False)
    return f, ax


def plot_panels(data, specs, *, trange=None, prefix=None, f=None, ax=None):
    """Draw several panel specs on a shared UTC axis."""
    specs = list(specs)
    if not specs:
        raise ValueError("specs must contain at least one PanelSpec")
    f, axes = tools.axes_for_columns(len(specs), f=f, ax=ax, sharex=True)
    for axis, spec in zip(axes, specs):
        plot_panel(data, spec, trange=trange, prefix=prefix, f=f, ax=axis)
    axes[-1].set_xlabel("Time (UTC)")
    f.autofmt_xdate()
    return f, axes
