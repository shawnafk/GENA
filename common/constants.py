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


