"""Repository-level scientific plots and cross-mission composition.

Use ``plots.basic`` for mission products and ``plots.advance`` for joint plots.
Submodules are loaded on demand so importing plots does not load every mission.
"""

from importlib import import_module

__all__ = ["basic", "advance", "vis", "savefig", "tplot", "Panel", "overview", "env", "plot_region_code"]


def __getattr__(name):
    if name in ("basic", "advance", "vis"):
        value = import_module(f"{__name__}.{name}")
    elif name == "savefig":
        value = import_module(f"{__name__}.vis").savefig
    elif name == "tplot":
        value = import_module(f"{__name__}.basic").tplot
    elif name in ("Panel", "overview", "env", "plot_region_code"):
        value = getattr(import_module(f"{__name__}.advance"), name)
    else:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    globals()[name] = value
    return value
