"""Lazy public entries for the shared GENA foundation modules.

Table types and operations live in :mod:`common.frame`; time normalization,
coordinates and geometry, configuration, and constants remain in their owning
submodules. Only the explicit mapping below is re-exported here.
"""

import importlib
import sys

# ---------------------------------------------------------------------------
# Module-level __getattr__ for lazy re-exports.
# Sub-modules that do not yet exist are imported on first access without
# raising noise when ``from common import ...`` is evaluated.
# ---------------------------------------------------------------------------

_SYMBOL_TO_MODULE: dict[str, str] = {
    # Constants are lazy too: importing array/time types must not require Astropy.
    "SPEED_C": "common.constants",
    "MASS_E": "common.constants",
    "MASS_H": "common.constants",
    "MASS_O": "common.constants",
    "RE": "common.constants",
    "RM": "common.constants",
    "X_MCP": "common.constants",
    "Y_MCP": "common.constants",

    # Only the main coordinate data type is re-exported at package level.
    # Array transforms are imported explicitly from common.coord.
    "Coordinates": "common.coord",
    # Time class, interval tools, and miscellaneous/interpolation helpers
    "Time": "common.time",
    "timer": "common.decos",
    "logger": "common.decos",
    "repeat": "common.decos",
    "mid_ary": "common.interp",
    "i1d": "common.interp",
    "i_rot": "common.interp",
    "i_pos": "common.interp",
    "f_interp_vars": "common.interp",
    "clean_nan": "common.interp",
    # common.frame — shared mission DataFrame types and operations
    "TimeIndexedDataFrame": "common.frame",
    "TimeSeriesDataFrame": "common.frame",
    "concat": "common.frame",
    "read_variable": "common.frame",
    "split_time_intervals": "common.frame",
}


def __getattr__(name):
    """Lazy-import the sub-module that owns *name* and return the symbol."""
    if name in _SYMBOL_TO_MODULE:
        mod_name = _SYMBOL_TO_MODULE[name]
        try:
            mod = importlib.import_module(mod_name)
        except ImportError:
            raise AttributeError(
                f"'{__name__}' has no attribute '{name}'. "
                f"The sub-module '{mod_name}' has not been created yet."
            ) from None
        attr = getattr(mod, name)
        # Cache on this module so subsequent lookups skip __getattr__
        setattr(sys.modules[__name__], name, attr)
        return attr

    raise AttributeError(f"module '{__name__}' has no attribute '{name}'")


def __dir__():
    """List all known re-exported symbols for tab-completion / introspection."""
    return sorted(set(_SYMBOL_TO_MODULE.keys()) | {"SPEED_C", "MASS_E", "MASS_H",
                   "MASS_O", "RE", "RM", "X_MCP", "Y_MCP"})


__all__ = sorted(set(_SYMBOL_TO_MODULE.keys()) | {"SPEED_C", "MASS_E", "MASS_H",
                "MASS_O", "RE", "RM", "X_MCP", "Y_MCP"})
