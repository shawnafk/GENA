#config runtime for GENA CODE

"""Runtime data and output paths; importing this module performs no I/O."""
import os
from pathlib import Path

__all__ = ["resolve_paths", "configure", "get_configured_paths", "resolve_output_path"]

_RUNTIME_FILE_BASE = None
_RUNTIME_OUTPUT_DIR = None
_CONFIG_CACHE = {}


def resolve_paths(file_base=None, *, output_dir=None):
    """Resolve paths without validation or directory creation.

    Explicit arguments take precedence over runtime configuration, then the
    environment. GENA_RAW/PQT, SPEDAS and GENA_RES are data-root-relative;
    GENA_OUTPUT_DIR (or output_dir) selects an independent result directory.
    """
    base = file_base if file_base is not None else _RUNTIME_FILE_BASE
    base = Path(base or os.environ.get("FILE_BASE", "/Volumes/Rsync/ENA/")).expanduser().resolve()
    paths = {"file_base": str(base)}
    for name, key in (("GENA_RAW", "raw"), ("GENA_PQT", "pqt"),
                      ("SPEDAS", "spedas"), ("GENA_RES", "res")):
        component = os.environ.get(name, f"/{name}/")
        paths[name.lower()] = component
        paths[f"{key}_base"] = str((base / component.lstrip("/")).resolve())
    output = output_dir if output_dir is not None else _RUNTIME_OUTPUT_DIR
    if output is None:
        output = os.environ.get("GENA_OUTPUT_DIR")
    if output is not None:
        paths["res_base"] = str(Path(output).expanduser().resolve())
    return paths


def _prepare_paths(paths, *, mkdir, validate):
    base = Path(paths["file_base"])
    if validate and not base.is_dir():
        raise ValueError(f"FILE_BASE path '{base}' is not a directory. Set FILE_BASE or pass file_base=.")
    if mkdir:
        for key in ("raw_base", "pqt_base", "spedas_base", "res_base"):
            Path(paths[key]).mkdir(parents=True, exist_ok=True)
    return dict(paths)


def configure(file_base=None, *, output_dir=None, mkdir=True, validate=True):
    """Set process-wide paths, optionally validating and creating directories.

    For output-only work use configure(output_dir=..., mkdir=False,
    validate=False); writers create their own parent directories on demand.
    """
    global _RUNTIME_FILE_BASE, _RUNTIME_OUTPUT_DIR
    paths = _prepare_paths(resolve_paths(file_base, output_dir=output_dir),
                           mkdir=mkdir, validate=validate)
    if file_base is not None:
        _RUNTIME_FILE_BASE = paths["file_base"]
    if output_dir is not None:
        _RUNTIME_OUTPUT_DIR = paths["res_base"]
    _CONFIG_CACHE.clear()
    return paths


def get_configured_paths(file_base=None, *, mkdir=True, validate=True, refresh=False):
    """Validate resolved paths once; environment/runtime changes refresh the cache."""
    paths = resolve_paths(file_base)
    key = (tuple(paths.items()), mkdir, validate)
    if refresh or key not in _CONFIG_CACHE:
        _CONFIG_CACHE[key] = _prepare_paths(paths, mkdir=mkdir, validate=validate)
    return dict(_CONFIG_CACHE[key])


def resolve_output_path(path=""):
    """Resolve relative output paths under res_base; absolute paths override it.

    No directories are created and no input-data path is validated here.
    """
    destination = Path(path).expanduser()
    if not destination.is_absolute():
        destination = Path(resolve_paths()["res_base"]) / destination
    return destination.resolve()


#constants should be stored in some higher format, a class
#contains methods like unit conversion, derivation, etc. 
"""
Shared physical and geometric constants for the GENA codebase.

All physical constants are derived from scipy and astropy reference values.
Geometric constants (X_MCP, Y_MCP) define the GENA detector geometry.

Constants defined here:
    SPEED_C   : speed of light in vacuum (m/s)
    MASS_E    : electron rest mass energy (eV)
    MASS_H    : proton rest mass energy (eV)
    MASS_O    : oxygen rest mass energy (eV)
    mass_H_kg : proton rest mass (kg)
    e_charge  : elementary charge (C)
    RE        : Earth equatorial radius (m)
    RM        : Moon radius (m)
    X_MCP     : MCP detector width (mm)
    Y_MCP     : MCP detector height (mm)
"""

from scipy import constants as sci_const
from astropy import constants as astro_const

# ── Physical Constants ──────────────────────────────────────────────

# Speed of light in vacuum (m/s)
# Derived from scipy.constants.speed_of_light
SPEED_C: float = sci_const.speed_of_light  # 299792458 m/s

# Electron rest mass energy (eV)
# Standard value from NIST CODATA
MASS_E: float = 510998.918  # eV

# Proton rest mass energy (eV)
# Proton/electron mass ratio = 1837.15267343
MASS_H: float = MASS_E * 1837.15267343  # eV

# Proton rest mass (kg)
# Standard value from CODATA 2018
mass_H_kg: float = 1.6726e-27  # kg

# Elementary charge (C, J/eV)
# Standard value from CODATA 2018
e_charge: float = 1.6022e-19  # C (coulombs)

# Oxygen rest mass energy (eV)
# Oxygen/proton mass ratio = 15.8695
MASS_O: float = MASS_H * 15.8695  # eV

# ── Astronomical Constants ──────────────────────────────────────────

# Earth equatorial radius (m)
# Derived from astropy.constants.R_earth.value
# NOTE: astropy value is 6378100.0 m (not the WGS-84 equatorial radius 6378137.0 m)
RE: float = astro_const.R_earth.value  # 6378100.0 m

# Mean lunar radius (m)
# Approximated as 27% of Earth's radius
RM: float = 0.27 * RE  # ~1722087 m


