import numpy as np

MU_0 = 4 * np.pi * 1e-7   # T·m/A
M_P = 1.6726219e-27       # kg
K_B = 1.380649e-23        # J/K


def alfven_velocity(B, n, m=1):
    """Alfvén velocity in m/s.

    Parameters
    ----------
    B : array_like
        Magnetic field magnitude in Gauss (G).
    n : array_like
        Number density in cm⁻³.
    m : float
        Ion mass in proton masses (default 1).

    Returns
    -------
    v_A : ndarray
        Alfvén velocity in m/s.
    """
    B_T = np.asarray(B) * 1e-4        # G → T
    n_m3 = np.asarray(n) * 1e6        # cm⁻³ → m⁻³
    rho = n_m3 * m * M_P
    with np.errstate(divide="ignore", invalid="ignore"):
        v_A = B_T / np.sqrt(MU_0 * rho)
    return np.where(np.isfinite(v_A), v_A, 0.0)


def beta(n, T, B):
    """Plasma beta (thermal / magnetic pressure).

    Parameters
    ----------
    n : array_like
        Number density in cm⁻³.
    T : array_like
        Temperature in eV.
    B : array_like
        Magnetic field magnitude in Gauss.

    Returns
    -------
    beta : ndarray
        Plasma beta (dimensionless).
    """
    n_m3 = np.asarray(n) * 1e6
    T_K = np.asarray(T) * 11604.5     # eV → K (1 eV ≈ 11604.5 K)
    P_th = n_m3 * K_B * T_K
    B_T = np.asarray(B) * 1e-4
    P_B = B_T**2 / (2 * MU_0)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(P_B > 0, P_th / P_B, np.inf)


"""
Shared physics models for ENA simulation across the GENA codebase.

Consolidates three heavily duplicated physics functions:
    sigma_HH      — H⁺ + H → H + H⁺ charge-exchange cross-section
    density_atom  — Hydrogen geocorona density model
    f_maxwellian  — Drift Maxwellian velocity distribution

All formulas preserve exact numerical behavior from their canonical sources.
References:
    - Equation (2) from the provided PDF (sigma_HH)
    - Chamberlain (1963) or Bishop (1999) style geocorona model (density_atom)
    - Standard kinetic theory Maxwellian (f_maxwellian)
"""

import numpy as np
from common.constants import mass_H_kg, e_charge


def sigma_HH(E):
    """
    Charge exchange cross-section for H⁺ + H → H + H⁺.

    Formula:
        sigma(E) = (a1 - a2 * ln(E))² * (1 - exp(-a3 / E))^4.5 × 1e-16 cm²

    where:
        a1 = 4.15, a2 = 0.531, a3 = 67.3

    Source: Equation (2) in the provided PDF.
    Canonical source: merge/map_T.py:36.

    Parameters
    ----------
    E : float or ndarray
        Energy in keV. Values ≤ 1e-5 are clamped to avoid log(0).

    Returns
    -------
    float or ndarray
        Cross-section in cm².
    """
    # Prevent log(0) errors with a small epsilon (canonical from map_T.py)
    E = np.maximum(E, 1e-5)
    a1 = 4.15
    a2 = 0.531
    a3 = 67.3
    return (a1 - a2 * np.log(E)) ** 2 * (1 - np.exp(-a3 / E)) ** 4.5 * 1e-16


def density_atom(x, y, z):
    """
    Neutral hydrogen (geocorona) density model.

    Formula:
        n_H(r) = 25 × (10 / max(r, 1))³ × 10⁶  [m⁻³]

    where r = sqrt(x² + y² + z²) in Earth radii (RE).

    This is a simplified spherically-symmetric geocorona model
    approximating a Chamberlain exosphere with scale height ~10 RE.

    Source: bbf_ena_simple.py:44 (canonical).
    References: Chamberlain (1963), Bishop (1999).

    Parameters
    ----------
    x, y, z : float or ndarray
        Position in Earth radii (RE), geocentric GSM coordinates.

    Returns
    -------
    float or ndarray
        Neutral hydrogen density in m⁻³.
    """
    r = np.sqrt(x**2 + y**2 + z**2)
    r = np.maximum(r, 1.0)
    return 25 * (10 / r) ** 3 * 1e6


def f_maxwellian(vx, vy, vz, Ux, Uy, Uz, Ti_eV):
    """
    Local drift Maxwellian velocity distribution.

    Formula:
        f(v) = (1 / (v_th √π)³) × exp(-|v - U|² / v_th²)

    where:
        v_th = √(2 × Ti_eV × e_charge / mass_H_kg)

    This is a 3D Maxwellian drifting at bulk velocity (Ux, Uy, Uz)
    with thermal energy Ti_eV.

    Source: bbf_ena_simple.py:50-56 (canonical).

    Parameters
    ----------
    vx, vy, vz : float or ndarray
        Particle velocity components (m/s).
    Ux, Uy, Uz : float or ndarray
        Bulk drift velocity components (m/s).
    Ti_eV : float
        Ion thermal energy (eV).

    Returns
    -------
    float or ndarray
        Distribution function value (s³/m⁶).
    """
    v_th = np.sqrt(2 * Ti_eV * e_charge / mass_H_kg)
    coeff = 1.0 / (v_th * np.sqrt(np.pi)) ** 3
    exponent = -((vx - Ux) ** 2 + (vy - Uy) ** 2 + (vz - Uz) ** 2) / v_th ** 2
    return coeff * np.exp(exponent)


def alfven_mach(B_nT, n_cm3, v_km_s):
    """Alfvén Mach number for B [nT], n [cm^-3] and v [km/s]."""
    return np.asarray(v_km_s) * 1000.0 / alfven_velocity(
        np.asarray(B_nT) * 1e-5,
        n_cm3,
        1,
    )


def plasma_beta(B_nT, n_cm3, T_eV):
    """Plasma beta for B [nT], n [cm^-3] and T [eV]."""
    return beta(n_cm3, T_eV, np.asarray(B_nT) * 1e-5)
