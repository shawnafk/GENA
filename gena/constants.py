import numpy as np
from scipy.spatial.transform import Rotation

from common.constants import (
    MASS_E,
    MASS_H,
    MASS_O,
    RE,
    RM,
    SPEED_C,
    e_charge,
    mass_H_kg,
)

# ── GENA Detector Geometry ───────────────────────────────────────────────
# MCP detector width (mm)
# GENA geometric size — microchannel plate detector X dimension
#in mm
X_MCP: int = 145  # mm
X_CF: int = 90  # mm

# MCP detector height (mm)
# GENA geometric size — microchannel plate detector Y dimension
Y_MCP: int = 55  
Y_CF: int = 30  

# Effective full-field-of-view geometric factor reported for GENA.  Science
# product functions require callers to pass this value explicitly; defining it
# here does not make it an implicit calibration default.
GENA_GEOMETRIC_FACTOR_CM2_SR: float = 19.2

D_int: int = 157 
D_ext: int = 512 
D_all: int = 669 

WX_int: int = 241 
WY_int: int = 59

# in degree
THETA_X = np.rad2deg(np.abs(np.arctan((WX_int - X_CF)/2/D_int)))
THETA_Y = np.rad2deg(np.abs(np.arctan((WY_int - Y_CF)/2/D_int)))


#trustable is > 4
def ENA_PRE(E_AFT):
    E_PRE = np.ones_like(E_AFT)
    E_spx = 4
    mask_low  = E_AFT < E_spx
    mask_high = E_AFT >= E_spx
    E_PRE[mask_low] =  (E_AFT[mask_low] + 0.02)/0.873
    E_PRE[mask_high] =  ((E_AFT[mask_high] + 2.95) - 0.976)/0.961
    return E_PRE


# Relative energy efficiency from the fitted curves in the calibration figure.
# Both functions are normalized to 10 keV.
def eff_H(E):
    """Return H relative efficiency for the valid 1–60 keV range."""
    E = np.asarray(E, dtype=float)
    E = np.clip(E, 1, 60)

    low_energy = (
        63.31 * (1 - np.exp(-E / 0.56))
        + 12.76 * (1 - np.exp(-E / 8.59))
        - 52.96
    )
    high_energy = -3.41 * np.exp(E / 43.52) + 28.73
    efficiency = np.where(E <= 25, low_energy, high_energy)

    efficiency_at_10keV = (
        63.31 * (1 - np.exp(-10 / 0.56))
        + 12.76 * (1 - np.exp(-10 / 8.59))
        - 52.96
    )
    return efficiency / efficiency_at_10keV


def eff_O(E):
    """Return O relative efficiency for the valid 5–60 keV range."""
    E = np.asarray(E, dtype=float)
    E = np.clip(E, 5, 60)

    efficiency = -37.66 * np.exp(-E / 62.53) + 39.81
    efficiency_at_10keV = -37.66 * np.exp(-10 / 62.53) + 39.81
    return efficiency / efficiency_at_10keV


def eff_E(E, spec="H"):
    """Return relative energy efficiency normalized at 10 keV."""
    if spec == "O" or spec == 2:
        return eff_O(E)
    return eff_H(E)

def _rotation_matrix(axis: str, degrees: float) -> np.ndarray:
    """Return a 3x3 single-axis rotation matrix without loading coordinate backends."""
    return Rotation.from_euler(axis, degrees, degrees=True).as_matrix()


#: Probe A boresight: -5° rotation of [0,0,-1] about X.
va0: np.ndarray = _rotation_matrix('x', -5) @ np.array([0.0, 0.0, -1.0])

#: Probe B boresight: +5° rotation of [0,0,-1] about X.
vb0: np.ndarray = _rotation_matrix('x', 5) @ np.array([0.0, 0.0, -1.0])

# ---------------------------------------------------------------------------
# Corner FOV vectors (instrument frame)
# ---------------------------------------------------------------------------

# Probe A — four quadrant corners
#: Probe A, top-left quadrant.
va1: np.ndarray = _rotation_matrix('y', 22.5) @ _rotation_matrix('x', -5) @ va0
#: Probe A, top-right quadrant.
va2: np.ndarray = _rotation_matrix('y', 22.5) @ _rotation_matrix('x', 5) @ va0
#: Probe A, bottom-right quadrant.
va3: np.ndarray = _rotation_matrix('y', -22.5) @ _rotation_matrix('x', 5) @ va0
#: Probe A, bottom-left quadrant.
va4: np.ndarray = _rotation_matrix('y', -22.5) @ _rotation_matrix('x', -5) @ va0

# Probe B — four quadrant corners
#: Probe B, top-left quadrant.
vb1: np.ndarray = _rotation_matrix('y', 22.5) @ _rotation_matrix('x', -5) @ vb0
#: Probe B, top-right quadrant.
vb2: np.ndarray = _rotation_matrix('y', 22.5) @ _rotation_matrix('x', 5) @ vb0
#: Probe B, bottom-right quadrant.
vb3: np.ndarray = _rotation_matrix('y', -22.5) @ _rotation_matrix('x', 5) @ vb0
#: Probe B, bottom-left quadrant.
vb4: np.ndarray = _rotation_matrix('y', -22.5) @ _rotation_matrix('x', -5) @ vb0

# ---------------------------------------------------------------------------
# Full FOV grids  (instrument frame)
# ---------------------------------------------------------------------------

#: Probe A FOV set: 11 × 19 = 209 vectors spanning ±5° in X and ±22.5° in Y.
fov_set_a: list[np.ndarray] = []
for x in np.linspace(-5.0, 5.0, 11):
    for y in np.linspace(-22.5, 22.5, 19):
        fov_set_a.append(_rotation_matrix('y', y) @ _rotation_matrix('x', x) @ va0)

#: Probe B FOV set: 11 × 19 = 209 vectors spanning ±5° in X and ±22.5° in Y.
fov_set_b: list[np.ndarray] = []
for x in np.linspace(-5.0, 5.0, 11):
    for y in np.linspace(-22.5, 22.5, 19):
        fov_set_b.append(_rotation_matrix('y', y) @ _rotation_matrix('x', x) @ vb0)




#params func1, calculate K and B for linear transformation from ns to mm
def K_B(params,probe,dir):
    L = params[probe+f'_{dir}L']
    R = params[probe+f'_{dir}R']
    D = params['D'+dir]
    K = D/(R - L)
    B = - D*(R+L)/(R - L)/2
    return K,B

#params func2, cast all values in dict to float32, to save memory and speed up computation
def cast_dict_to_float32(params_dict):
    """
    convert default float to float32
    """
    return {
        k: np.float32(v) if isinstance(v, (int, float, np.number)) else v 
        for k, v in params_dict.items()
    }

def init_param(SN_ID):
    """
    Initialize ENA detector parameters.

    Args:
        SN_ID: 1 for orbit, 2 for ground

    Returns:
        Dictionary containing all detector parameters.
    """
    params = {}

    params['DX'] = 145
    params['DY'] = 55
    # Probe A constants
    params['A_XL']  = -54.91726
    params['A_XR'] = 50.44742
    params['A_YL']  = -16.68598
    params['A_YR'] = 41.62719

    # Probe B constants
    params['B_XL']  = -51.18360
    params['B_XR'] = 53.40479

    if SN_ID == 1:
        params['B_YL']  = -16.68598
        params['B_YR'] = 41.62719
    elif SN_ID == 2:
        params['B_YL']  = -34.71570
        params['B_YR'] = 23.59747

    params['A_KX'],params['A_BX'] = K_B(params,'A','X')
    params['A_KY'],params['A_BY'] = K_B(params,'A','Y')
    
    params['B_KX'],params['B_BX'] = K_B(params,'B','X')
    params['B_KY'],params['B_BY'] = K_B(params,'B','Y')
    
    
    #read from real data distribution
    params['A_X_SUM'] = 50.0
    params['A_Y_SUM'] = 38.5


    params['B_X_SUM'] = 55.
    params['B_Y_SUM'] = 43.

    # Cast before the closure is created so the callable and the public mapping
    # share exactly the same dictionary.  Mutating a calibration coefficient in
    # the returned mapping now also updates the conversion performed below.
    params = cast_dict_to_float32(params)

    def ns_2_xy(df, flag, prob, prefix=''):
        """Convert TDC nanoseconds to X/Y coordinates."""
        p = prob.upper()
        KX, BX = params[f'{p}_KX'], params[f'{p}_BX']
        KY, BY = params[f'{p}_KY'], params[f'{p}_BY']
        if flag == 0:
            x = KX * (df[f'{prefix}TX_1s'] - df[f'{prefix}TX_2s']) + BX
            y = KY * (df[f'{prefix}TY_1s'] - df[f'{prefix}TY_2s']) + BY
        elif flag == 1:
            x = KX * (df[f'{prefix}TX_1e'] - df[f'{prefix}TX_2e']) + BX
            y = KY * (df[f'{prefix}TY_1e'] - df[f'{prefix}TY_2e']) + BY
        else:
            raise ValueError("flag must be 0 (carbon) or 1 (mcp)")
        return x, y
    params['ns_2_xy'] = ns_2_xy
    return params
