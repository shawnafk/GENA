"""Mission-specific product plots sharing the plots.vis primitives."""

from importlib import import_module

__all__ = ["gena", "themis", "mms", "omni", "tplot"]


def __getattr__(name):
    if name not in __all__:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    if name == "tplot":
        value = import_module(f"{__name__}.timeseries").tplot
        globals()[name] = value
        return value
    module = import_module(f"{__name__}.{name}")
    globals()[name] = module
    return module
