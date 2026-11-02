"""Repository-level scientific plots and cross-mission composition.

Mission products live in the flat modules ``plots.gena``, ``plots.themis``,
``plots.mms`` and ``plots.omni``. Shared plotting primitives live in
``plots.tools``, generic time-series plotting in ``plots.timeseries`` and
cross-mission composition in ``plots.joint``. Submodules are imported lazily
so importing :mod:`plots` does not load every mission.

The old ``plots.basic``, ``plots.advance`` and ``plots.vis`` paths are gone.
Use the flat module names and :mod:`plots.tools` instead.
"""

from importlib import import_module

__all__ = [
    "Panel",
    "ani",
    "env",
    "gena",
    "joint",
    "labels",
    "mms",
    "omni",
    "overview",
    "panels",
    "plot_region_code",
    "savefig",
    "scene",
    "themis",
    "theme",
    "timeseries",
    "tools",
    "tplot",
    "v3d",
]

_MODULES = {
    "ani": "ani",
    "gena": "gena",
    "joint": "joint",
    "labels": "labels",
    "mms": "mms",
    "omni": "omni",
    "panels": "panels",
    "scene": "scene",
    "themis": "themis",
    "theme": "theme",
    "timeseries": "timeseries",
    "tools": "tools",
    "v3d": "v3d",
}

_REMOVED = {
    "advance": "plots.advance was removed; use plots.joint and the flat mission modules",
    "basic": "plots.basic was removed; use plots.gena, plots.themis, plots.mms, plots.omni",
    "vis": "plots.vis was removed; use plots.tools",
}


def __getattr__(name):
    if name in _REMOVED:
        raise AttributeError(_REMOVED[name])
    if name in _MODULES:
        value = import_module(f"{__name__}.{_MODULES[name]}")
    elif name == "savefig":
        value = import_module(f"{__name__}.tools").savefig
    elif name == "tplot":
        value = import_module(f"{__name__}.timeseries").tplot
    elif name in {"Panel", "overview", "env", "plot_region_code"}:
        value = getattr(import_module(f"{__name__}.joint"), name)
    else:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    globals()[name] = value
    return value
