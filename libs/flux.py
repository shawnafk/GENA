"""Maxwellian-core and Kappa-tail ion flux models for ENA drift studies."""

import numpy as np
from scipy.optimize import least_squares
from scipy.special import gammaln

from common.constants import e_charge, mass_H_kg, RE

from collections.abc import Sequence

from libs.plasma import sigma_HH
'''
introduce the unit system used here:
eV for energy,
cm^-3 for density,
km/s for speed, 
particle fluxes in cm^-2 s^-1 sr^-1 eV^-1.
'''

_PARTICLE_FLUX_PREFACTOR = (
    1.0e2 * np.sqrt(e_charge / (2.0 * mass_H_kg)) / np.pi**1.5
)
_DRIFT_MU, _DRIFT_WEIGHT = np.polynomial.legendre.leggauss(32)


def _validate_energy_eV(E_eV):
    energy = np.asarray(E_eV, dtype=float)
    if np.any(~np.isfinite(energy)) or np.any(energy < 0.0):
        raise ValueError("E_eV must be finite and non-negative")
    return energy


def _validate_population(n_cm3, T_eV, *, kappa=None):
    if not np.isfinite(n_cm3) or n_cm3 < 0.0:
        raise ValueError("density must be finite and non-negative")
    if n_cm3 == 0.0:
        return
    if not np.isfinite(T_eV) or T_eV <= 0.0:
        raise ValueError("temperature must be finite and positive")
    if kappa is not None and (not np.isfinite(kappa) or kappa <= 1.5):
        raise ValueError("Kappa populations require kappa > 1.5")


def _validate_speed_kms(speed_kms):
    speed = float(speed_kms)
    if not np.isfinite(speed):
        raise ValueError("drift speed must be finite")
    return abs(speed)


def j_maxwell(E_eV, n_cm3, T_eV):
    """Return rest-frame Maxwellian particle flux per eV.

    Parameters use eV and cm^-3. The returned differential particle flux has
    units cm^-2 s^-1 sr^-1 eV^-1.
    """

    energy = _validate_energy_eV(E_eV)
    _validate_population(n_cm3, T_eV)
    if n_cm3 == 0.0:
        return np.zeros_like(energy)
    return (
        _PARTICLE_FLUX_PREFACTOR
        * n_cm3
        * energy
        / T_eV**1.5
        * np.exp(-energy / T_eV)
    )


def j_kappa(E_eV, n_cm3, T_eV, kappa):
    """Return rest-frame Kappa particle flux per eV.

    ``T_eV`` is the kinetic temperature of the Kappa population, so the
    distribution is defined only for ``kappa > 1.5``.
    """

    energy = _validate_energy_eV(E_eV)
    _validate_population(n_cm3, T_eV, kappa=kappa)
    if n_cm3 == 0.0:
        return np.zeros_like(energy)
    normalization = np.exp(
        gammaln(kappa + 1.0)
        - gammaln(kappa - 0.5)
        - 1.5 * np.log(kappa - 1.5)
    )
    return (
        _PARTICLE_FLUX_PREFACTOR
        * n_cm3
        * energy
        / T_eV**1.5
        * normalization
        * (1.0 + energy / ((kappa - 1.5) * T_eV)) ** (-(kappa + 1.0))
    )


def _omnidirectional_drift(E_eV, rest_flux, drift_kms):
    """Pitch-angle average an isotropic rest-frame particle flux."""

    energy = _validate_energy_eV(E_eV)
    speed = _validate_speed_kms(drift_kms)
    if speed == 0.0:
        return rest_flux(energy)

    drift_energy = 0.5 * mass_H_kg * (speed * 1.0e3) ** 2 / e_charge
    shape = (1,) * energy.ndim + (-1,)
    mu = _DRIFT_MU.reshape(shape)
    weight = _DRIFT_WEIGHT.reshape(shape)
    observed_energy = energy[..., None]
    rest_energy = (
        observed_energy
        + drift_energy
        - 2.0 * np.sqrt(observed_energy * drift_energy) * mu
    )
    rest_energy = np.maximum(rest_energy, np.finfo(float).tiny)
    transformed_flux = observed_energy / rest_energy * rest_flux(rest_energy)
    return 0.5 * np.sum(weight * transformed_flux, axis=-1)


def j_dual_model(
    E_eV,
    n_core_cm3,
    T_core_eV,
    n_kappa_cm3=0.0,
    T_kappa_eV=None,
    kappa=None,
    drift_kms=0.0,
):
    """Return an omnidirectional Maxwellian-core plus Kappa-tail flux per eV.

    The Maxwellian core and Kappa tail share one bulk drift speed. This models
    a single convecting proton population with a suprathermal Kappa tail.
    """

    energy = _validate_energy_eV(E_eV)
    core = _omnidirectional_drift(
        energy,
        lambda rest_energy: j_maxwell(rest_energy, n_core_cm3, T_core_eV),
        drift_kms,
    )
    tail = _omnidirectional_drift(
        energy,
        lambda rest_energy: j_kappa(
            rest_energy, n_kappa_cm3, T_kappa_eV, kappa
        ),
        drift_kms,
    )
    return core + tail


def fit_dual_moment_constrained(
    E_eV,
    observed_flux,
    obs_n_cm3,
    obs_T_eV,
    *,
    en=False,
    drift_kms=None,
    drift_bounds=(0.0, 3000.0),
    R=None,
    tail_fraction_bounds=(1.0e-5, 0.8),
    temperature_ratio_bounds=(1.001, 100.0),
    kappa_bounds=(1.6, 20.0),
    log_flux_sigma_dex=0.15,
):
    """Fit a shared-drift Maxwellian core plus Kappa tail with exact moments.

    ``obs_n_cm3`` and ``obs_T_eV`` are imposed exactly through

    ``n_maxwell + n_kappa = obs_n_cm3`` and
    ``n_maxwell*T_maxwell + n_kappa*T_kappa = obs_n_cm3*T_esa_eV``.

    ``R`` fixes the Maxwellian density fraction ``n_maxwell / obs_n_cm3`` when
    supplied and must lie in ``[0, 1]``. ``R=1`` selects a pure Maxwellian,
    ``R=0`` selects a pure Kappa distribution, and ``R=None`` retains the
    original behavior of fitting the component fraction. For a mixed
    population, the Kappa term is constrained to be hotter than the Maxwellian
    core. Both components share one drift-speed magnitude. ``drift_kms=None``
    fits that speed inside ``drift_bounds``; supplying a numerical
    ``drift_kms`` keeps it fixed for a control fit. ``observed_flux`` is
    differential particle flux when ``en=False`` and differential energy flux
    when ``en=True``. Energy and all returned temperatures use eV.
    """

    energy = _validate_energy_eV(E_eV)
    measured = np.asarray(observed_flux, dtype=float)
    if energy.ndim != 1 or measured.ndim != 1 or energy.shape != measured.shape:
        raise ValueError("E_eV and observed_flux must be one-dimensional arrays of equal length")
    if not isinstance(en, (bool, np.bool_)):
        raise TypeError("en must be bool")
    if not np.isfinite(obs_n_cm3) or obs_n_cm3 <= 0.0:
        raise ValueError("obs_n_cm3 must be finite and positive")
    if not np.isfinite(obs_T_eV) or obs_T_eV <= 0.0:
        raise ValueError("obs_T_eV must be finite and positive")
    if not np.isfinite(log_flux_sigma_dex) or log_flux_sigma_dex <= 0.0:
        raise ValueError("log_flux_sigma_dex must be finite and positive")

    if R is not None:
        R = float(R)
        if not np.isfinite(R) or not 0.0 <= R <= 1.0:
            raise ValueError("R must lie in [0, 1]")

    fit_drift = drift_kms is None
    if fit_drift:
        drift_limits = np.asarray(drift_bounds, dtype=float)
        if (
            drift_limits.shape != (2,)
            or np.any(~np.isfinite(drift_limits))
            or not (0.0 <= drift_limits[0] < drift_limits[1])
        ):
            raise ValueError("drift_bounds must be increasing and non-negative")
        fixed_drift_kms = None
    else:
        fixed_drift_kms = _validate_speed_kms(drift_kms)

    tail_bounds = np.asarray(tail_fraction_bounds, dtype=float)
    ratio_bounds = np.asarray(temperature_ratio_bounds, dtype=float)
    kappa_limits = np.asarray(kappa_bounds, dtype=float)
    if R is None and (
        tail_bounds.shape != (2,)
        or not (0.0 < tail_bounds[0] < tail_bounds[1] < 1.0)
    ):
        raise ValueError("tail_fraction_bounds must lie strictly inside (0, 1)")
    if R != 0.0 and R != 1.0 and (
        ratio_bounds.shape != (2,)
        or not (1.0 < ratio_bounds[0] < ratio_bounds[1])
    ):
        raise ValueError("temperature_ratio_bounds must be increasing and greater than one")
    if R != 1.0 and (
        kappa_limits.shape != (2,)
        or not (1.5 < kappa_limits[0] < kappa_limits[1])
    ):
        raise ValueError("kappa_bounds must be increasing and greater than 1.5")

    valid = np.isfinite(measured) & (energy > 0.0) & (measured > 0.0)
    if np.count_nonzero(valid) < 7:
        raise ValueError("at least seven finite positive energy channels are required")
    energy_fit = energy[valid]
    particle_flux_fit = measured[valid]
    if en:
        particle_flux_fit = particle_flux_fit / energy_fit
    if np.any(~np.isfinite(particle_flux_fit)) or np.any(particle_flux_fit <= 0.0):
        raise ValueError("observed_flux cannot be converted to positive particle flux")
    log_measured = np.log10(particle_flux_fit)

    def parameters_from_free(params):
        if R == 1.0:
            return obs_n_cm3, 0.0, obs_T_eV, None, None, None
        if R == 0.0:
            return 0.0, obs_n_cm3, None, obs_T_eV, None, params[0]
        if R is None:
            tail_fraction, log10_temperature_ratio, kappa = params
            maxwell_fraction = 1.0 - tail_fraction
        else:
            maxwell_fraction = R
            tail_fraction = 1.0 - R
            log10_temperature_ratio, kappa = params
        temperature_ratio = 10.0**log10_temperature_ratio
        n_kappa = tail_fraction * obs_n_cm3
        n_maxwell = maxwell_fraction * obs_n_cm3
        T_maxwell = obs_T_eV / (
            maxwell_fraction + tail_fraction * temperature_ratio
        )
        T_kappa = temperature_ratio * T_maxwell
        return n_maxwell, n_kappa, T_maxwell, T_kappa, temperature_ratio, kappa

    def split_free(params):
        if fit_drift:
            return params[:-1], float(params[-1])
        return params, fixed_drift_kms

    def model_from_free(params, energy_eval):
        distribution_params, fitted_drift_kms = split_free(params)
        n_maxwell, n_kappa, T_maxwell, T_kappa, _, kappa = (
            parameters_from_free(distribution_params)
        )
        return j_dual_model(
            energy_eval,
            n_maxwell,
            T_maxwell,
            n_kappa,
            T_kappa,
            kappa,
            drift_kms=fitted_drift_kms,
        )

    def residual(params):
        model = np.maximum(model_from_free(params, energy_fit), np.finfo(float).tiny)
        return (np.log10(model) - log_measured) / log_flux_sigma_dex

    if R is None:
        distribution_lower = np.array([
            tail_bounds[0],
            np.log10(ratio_bounds[0]),
            kappa_limits[0],
        ])
        distribution_upper = np.array([
            tail_bounds[1],
            np.log10(ratio_bounds[1]),
            kappa_limits[1],
        ])
        distribution_starts = [
            np.array([fraction, np.log10(ratio), kappa])
            for fraction in (0.01, 0.05, 0.20)
            for ratio in (2.0, 6.0, 20.0)
            for kappa in (2.5, 4.0, 7.0)
        ]
    elif R == 0.0:
        distribution_lower = np.array([kappa_limits[0]])
        distribution_upper = np.array([kappa_limits[1]])
        distribution_starts = [np.array([kappa]) for kappa in (2.5, 4.0, 7.0)]
    elif R == 1.0:
        distribution_lower = distribution_upper = np.empty(0)
        distribution_starts = [np.empty(0)]
    else:
        distribution_lower = np.array([
            np.log10(ratio_bounds[0]),
            kappa_limits[0],
        ])
        distribution_upper = np.array([
            np.log10(ratio_bounds[1]),
            kappa_limits[1],
        ])
        distribution_starts = [
            np.array([np.log10(ratio), kappa])
            for ratio in (2.0, 6.0, 20.0)
            for kappa in (2.5, 4.0, 7.0)
        ]

    if fit_drift:
        lower = np.append(distribution_lower, drift_limits[0])
        upper = np.append(distribution_upper, drift_limits[1])
        trial_drifts = np.linspace(drift_limits[0], drift_limits[1], 4)
        starts = [
            np.append(distribution_start, trial_drift)
            for distribution_start in distribution_starts
            for trial_drift in trial_drifts
        ]
    else:
        lower = distribution_lower
        upper = distribution_upper
        starts = distribution_starts

    def interior(values):
        return np.clip(values, lower + 1.0e-8, upper - 1.0e-8)

    starts = [interior(start) for start in starts]

    if lower.size == 0:
        fitted_parameters = np.empty(0)
        fit_residual = residual(fitted_parameters)
        fit_cost = float(np.sum(np.sqrt(1.0 + fit_residual**2) - 1.0))
        fit_success = True
        fit_message = "all model parameters were fixed; no fit was required"
    else:
        candidates = [
            least_squares(
                residual,
                start,
                bounds=(lower, upper),
                x_scale="jac",
                loss="soft_l1",
                f_scale=1.0,
                max_nfev=20_000,
            )
            for start in starts
        ]
        successful = [
            result
            for result in candidates
            if result.success and np.isfinite(result.cost)
        ]
        if not successful:
            raise RuntimeError("all constrained dual-population fits failed")
        result = min(successful, key=lambda candidate: candidate.cost)
        fitted_parameters = result.x
        fit_cost = float(result.cost)
        fit_success = bool(result.success)
        fit_message = result.message

    distribution_parameters, fitted_drift_kms = split_free(fitted_parameters)
    n_maxwell, n_kappa, T_maxwell, T_kappa, temperature_ratio, kappa = (
        parameters_from_free(distribution_parameters)
    )
    model_particle_flux = model_from_free(fitted_parameters, energy)
    density_sum = n_maxwell + n_kappa
    temperature_sum = (
        n_maxwell * (T_maxwell or 0.0) + n_kappa * (T_kappa or 0.0)
    ) / density_sum
    log_rms_dex = float(
        np.sqrt(np.mean(residual(fitted_parameters) ** 2)) * log_flux_sigma_dex
    )
    maxwell_fraction = n_maxwell / density_sum

    return {
        "n_maxwell_cm3": float(n_maxwell),
        "n_kappa_cm3": float(n_kappa),
        "T_maxwell_eV": None if T_maxwell is None else float(T_maxwell),
        "T_kappa_eV": None if T_kappa is None else float(T_kappa),
        "kappa": None if kappa is None else float(kappa),
        "R": float(maxwell_fraction),
        "tail_fraction": float(1.0 - maxwell_fraction),
        "temperature_ratio": (
            None if temperature_ratio is None else float(temperature_ratio)
        ),
        "drift_kms": float(fitted_drift_kms),
        "density_sum_cm3": float(density_sum),
        "temperature_sum_eV": float(temperature_sum),
        "model_particle_flux": model_particle_flux,
        "model_energy_flux": energy * model_particle_flux,
        "log_rms_dex": log_rms_dex,
        "cost": fit_cost,
        "success": fit_success,
        "message": fit_message,
    }


"""Lowest-order estimate of ENA counts from a local ion spectrum.

The default unit convention is:

* ion energy: eV
* differential ion flux: cm^-2 s^-1 sr^-1 eV^-1
* neutral density: cm^-3
* charge-exchange cross section: cm^2
* detector area: cm^2
* path length: Earth radii
* accumulation time: s

The calculation assumes that the supplied ion spectrum is an omnidirectional
average expressed per steradian and that it is isotropic within the detector
field of view.  ``directional_enhancement`` can be used for a later, simple
anisotropy sensitivity test; its default value is one.
"""



def rectangular_solid_angle(width_deg: float, height_deg: float) -> float:
    """Return the exact solid angle of a centered rectangular field of view."""

    width = np.deg2rad(width_deg)
    height = np.deg2rad(height_deg)
    return float(
        4.0 * np.arcsin(np.sin(width / 2.0) * np.sin(height / 2.0))
    )


#return callable
#use scipy interp1d to do the interpolation
def interpolate_flux(
    sample_energy_eV: np.ndarray,
    measured_energy_eV: np.ndarray,
    measured_flux_per_eV: np.ndarray,
    method: str,
) -> np.ndarray:
    if method == "loglog":
        if np.any(measured_flux_per_eV <= 0.0):
            raise ValueError("Log-log interpolation requires strictly positive flux.")
        return 10.0 ** np.interp(
            np.log10(sample_energy_eV),
            np.log10(measured_energy_eV),
            np.log10(measured_flux_per_eV),
        )

    if method == "linear":
        return np.interp(
            sample_energy_eV,
            measured_energy_eV,
            measured_flux_per_eV,
        )

    raise ValueError("interpolation must be 'loglog' or 'linear'.")

from scipy.interpolate import interp1d
def fun_interpolate_flux(
    measured_energy_eV: np.ndarray,
    measured_flux_per_eV: np.ndarray,
    method: str,
) -> np.ndarray:
    if method == "loglog":
        if np.any(measured_flux_per_eV <= 0.0):
            raise ValueError("Log-log interpolation requires strictly positive flux.")
        return lambda E: 10**interp1d(np.log10(measured_energy_eV), np.log10(measured_flux_per_eV))(np.log10(E))

    if method == "linear":
        return lambda E: np.interp(measured_energy_eV, measured_flux_per_eV)(E)

    raise ValueError("interpolation must be 'loglog' or 'linear'.")


def qk_estimate_ena_counts(
    energy_eV: Sequence[float],
    ion_flux_per_eV: Sequence[float],
    energy_windows_eV: Sequence[tuple[float, float]] = ((1.0e3, 6.0e4),),
    *,
    neutral_density_cm3: float = 0.2,
    path_length_RE: float = 10.0,
    efficiency: float = 0.1,
    detector_area_cm2: float = 3.75,
    fov_deg: tuple[float, float] = (20.0, 10.0),
    accumulation_time_s: float = 60.0,
    directional_enhancement: float = 1.0,
    interpolation: str = "loglog",
    integration_points_per_decade: int = 400,
) -> dict[str, float | np.ndarray]:
    """Estimate ENA counts in one or more energy windows.

    Parameters
    ----------
    energy_eV, ion_flux_per_eV
        Discrete local ion spectrum.  Flux units must be
        cm^-2 s^-1 sr^-1 eV^-1.
    energy_windows_eV
        Sequence of ``(lower, upper)`` energy windows in eV.
    neutral_density_cm3
        Constant neutral-H density used by the quick estimate.
    path_length_RE
        Uniform source-region length in Earth radii.
    efficiency
        Overall dimensionless detection efficiency.
    detector_area_cm2
        Effective detector-plane area.
    fov_deg
        Full rectangular field-of-view widths ``(width, height)`` in degrees.
    accumulation_time_s
        Count accumulation time.
    directional_enhancement
        Multiplicative enhancement relative to the isotropic use of the
        omnidirectional-average spectrum.  Keep at 1 for the baseline model.
    interpolation
        ``"loglog"`` (default) or ``"linear"`` one-dimensional interpolation.
    integration_points_per_decade
        Density of the logarithmic numerical integration grid.

    Returns
    -------
    dict
        Includes ``counts_per_window`` and ``total_counts``, together with the
        main intermediate quantities used in the estimate.

    Notes
    -----
    A requested window must lie inside the supplied energy range.  The
    function intentionally does not extrapolate the measured spectrum.
    """

    energy, flux = energy_eV, ion_flux_per_eV
    windows = np.asarray(energy_windows_eV, dtype=float)

    if integration_points_per_decade < 20:
        raise ValueError("integration_points_per_decade must be at least 20.")
    solid_angle_sr = rectangular_solid_angle(*fov_deg)
    maximum_directional_enhancement = 4.0 * np.pi / solid_angle_sr
    if directional_enhancement > maximum_directional_enhancement * (1.0 + 1.0e-12):
        raise ValueError(
            "directional_enhancement exceeds the 4*pi/FOV limit implied by "
            "an omnidirectional-average input spectrum."
        )

    path_length_cm = path_length_RE * RE * 100
    column_density_cm2 = neutral_density_cm3 * path_length_cm
    exposure_cm2_s_sr = (
        accumulation_time_s
        * detector_area_cm2
        * solid_angle_sr
        * efficiency
        * directional_enhancement
    )

    counts = []
    sigma_j_integrals = []
    for lower_eV, upper_eV in windows:
        decades = np.log10(upper_eV / lower_eV)
        grid_size = max(2, int(np.ceil(integration_points_per_decade * decades)) + 1)
        integration_energy = np.geomspace(lower_eV, upper_eV, grid_size)
        interpolated_flux = interpolate_flux(
            integration_energy,
            energy,
            flux,
            interpolation,
        )
        sigma = np.asarray(sigma_HH(integration_energy / 1.0e3), dtype=float)
        if sigma.shape != integration_energy.shape:
            sigma = np.broadcast_to(sigma, integration_energy.shape)
        if np.any(~np.isfinite(sigma)) or np.any(sigma < 0.0):
            raise ValueError("The cross-section function returned invalid values.")

        integrand = sigma * interpolated_flux
        sigma_j = float(np.trapezoid(integrand, integration_energy))
        sigma_j_integrals.append(sigma_j)
        counts.append(exposure_cm2_s_sr * column_density_cm2 * sigma_j)

    counts_array = np.asarray(counts)
    return {
        "energy_windows_eV": windows,
        "counts_per_window": counts_array,
        "total_counts": float(counts_array.sum()),
        "sigma_j_integrals_s-1_sr-1": np.asarray(sigma_j_integrals),
        "solid_angle_sr": solid_angle_sr,
        "column_density_cm-2": column_density_cm2,
        "exposure_cm2_s_sr": exposure_cm2_s_sr,
        "maximum_directional_enhancement": maximum_directional_enhancement,
    }

#get directon vector from azimuth and elevation angles in degrees
def get_los_dir(az_deg, el_deg):
    az_rad, el_rad = np.radians(az_deg), np.radians(el_deg)
    return np.array(
        [
            -np.cos(el_rad) * np.cos(az_rad),
            -np.cos(el_rad) * np.sin(az_rad),
            np.sin(el_rad),
        ]
    )


def line_of_sight(S, U , dl=0.1, Lend=50):
    """
    Calculate coordinate points along the line of sight
    
    Parameters:
        Sx, Sy, Sz: Observer position (units: RE)
        Ax (A_long): Angle along Z-axis (-22.5 to 22.5 degrees)
        Ay (A_short): Angle along Y-axis (-5 to 5 degrees)
        dl: Length of each line segment (units: RE)
        Lend: Total length of the line of sight (units: RE)
    
    Returns:
        x, y, z: Arrays of coordinates along the line of sight (units: RE)
    """
    N = int(Lend / dl)
    x = np.zeros(N)
    y = np.zeros(N)
    z = np.zeros(N)
    
    # Generate coordinates along the line of sight
    distances = dl * np.arange(N)  # Distances from observer
    x = S[0] + U[0] * distances
    y = S[1] + U[1] * distances
    z = S[2] + U[2] * distances
    
    return x, y, z



import math
from scipy.integrate import quad
def sigma_HpH(E_ev):
    E_kev = E_ev / 1000
    return sigma_HH(E_kev)

#assuming J = k E^(-gamma)
def J(E_ev, k,gamma):
    return k * (E_ev ** (-gamma))

from scipy.interpolate import interp1d
def J_int(E_ev,Flux):
    # Create an interpolation function for the flux
    interp_func = interp1d(Flux[:, 0], Flux[:, 1], kind='linear', fill_value="extrapolate")
    return interp_func(E_ev)


# 被积函数
def integrand(E_ev, k,gamma):
    return J(E_ev,k,gamma)*sigma_HpH(E_ev)


def integrand_real(E_ev, J):
    return J(E_ev)*sigma_HpH(E_ev)

if __name__ == "__main__":
    e1,e2 = 1000, 60000
    j1,j2= 100, 0.01
    gamma = -(math.log10(j2) - math.log10(j1)) / (math.log10(e2) - math.log10(e1))
    k = j1 * e1**(gamma)
    Phi, err = quad(integrand, e1, e2, args=(k,gamma))
    print(f"积分结果 Phi = {Phi:.3e}, 误差估计 {err:.3e}")

from typing import Callable
def integral_ena_counts(energy_eV: np.ndarray, flux_per_eV: np.ndarray, density: float):
    """Unfinished draft; density alone does not define detector counts.

    Use qk_estimate_ena_counts with explicit path length, geometry, efficiency
    and accumulation time until this draft's physical inputs are specified.
    """
    raise NotImplementedError("Use qk_estimate_ena_counts with explicit geometry and exposure.")

__all__ = [
    "qk_estimate_ena_counts",
    "fit_dual_moment_constrained",
    "j_dual_model",
    "j_kappa",
    "j_maxwell",
]
