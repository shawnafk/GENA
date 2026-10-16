"""Joint ENA and in-situ estimates of magnetotail energy transport.

The routines in this module keep three physically different operations
separate:

1. particle energy flux from a measured velocity distribution;
2. electromagnetic energy flux from electric and magnetic fields; and
3. the optically thin charge-exchange mapping between parent-ion and ENA
   directional intensities.

All public functions use SI units unless a parameter name explicitly says
otherwise.  ENA energy is expressed in keV to match GENA products.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import lsq_linear

from libs.plasma import MU_0, sigma_HH


KEV_TO_J = 1e3 * 1.602176634e-19


@dataclass(frozen=True)
class ENAInversionResult:
    """Result of an energy-by-energy non-negative ENA inversion."""

    parent_intensity: np.ndarray
    modeled_intensity: np.ndarray
    normalized_residual: np.ndarray
    reduced_chi2: np.ndarray
    success: np.ndarray
    messages: tuple[str, ...]


@dataclass(frozen=True)
class TransportBudget:
    """Power crossing a surface, separated into particle and EM terms."""

    particle_power_w: np.ndarray
    electromagnetic_power_w: np.ndarray
    total_power_w: np.ndarray


def kinetic_energy_flux(
    distribution_si,
    velocity_mps,
    velocity_cell_volume,
    *,
    particle_mass_kg,
):
    """Calculate ``Q = integral(0.5 m v^2 v f d^3v)`` in W m^-2.

    Parameters
    ----------
    distribution_si : array_like, shape (...)
        Phase-space density in s^3 m^-6.
    velocity_mps : array_like, shape (..., 3)
        Velocity at each distribution cell in m s^-1.
    velocity_cell_volume : float or array_like, shape (...)
        The corresponding ``d^3v`` cell volume in m^3 s^-3.
    particle_mass_kg : float
        Particle mass in kg.
    """
    f = np.asarray(distribution_si, dtype=float)
    velocity = np.asarray(velocity_mps, dtype=float)
    dv3 = np.asarray(velocity_cell_volume, dtype=float)
    if velocity.shape != f.shape + (3,):
        raise ValueError("velocity_mps must have shape distribution.shape + (3,)")
    try:
        dv3 = np.broadcast_to(dv3, f.shape)
    except ValueError as exc:
        raise ValueError("velocity_cell_volume must broadcast to distribution") from exc
    speed2 = np.einsum("...i,...i->...", velocity, velocity)
    cell_flux = (
        0.5
        * float(particle_mass_kg)
        * f[..., None]
        * speed2[..., None]
        * velocity
        * dv3[..., None]
    )
    axes = tuple(range(f.ndim))
    return np.sum(cell_flux, axis=axes)


def mhd_particle_energy_flux(
    density_m3,
    bulk_velocity_mps,
    pressure_pa,
    *,
    particle_mass_kg,
    heat_flux_w_m2=None,
):
    """Return the bulk/pressure particle-energy flux in W m^-2.

    A scalar ``pressure_pa`` is treated as isotropic and gives the familiar
    Angelopoulos/Cao MHD expression

    ``0.5 rho U^2 U + 2.5 p U``.

    A pressure tensor uses the exact second-moment decomposition

    ``0.5 rho U^2 U + P dot U + 0.5 trace(P) U``.

    Supplying the third central moment ``heat_flux_w_m2`` makes the expression
    equal to the full kinetic third moment (up to measurement coverage).
    """
    density = np.asarray(density_m3, dtype=float)
    velocity = np.asarray(bulk_velocity_mps, dtype=float)
    pressure = np.asarray(pressure_pa, dtype=float)
    if velocity.shape[-1:] != (3,):
        raise ValueError("bulk_velocity_mps must end with three components")

    rho = density * float(particle_mass_kg)
    speed2 = np.einsum("...i,...i->...", velocity, velocity)
    bulk = 0.5 * rho[..., None] * speed2[..., None] * velocity

    if pressure.ndim >= 2 and pressure.shape[-2:] == (3, 3):
        pressure_work = np.einsum("...ij,...j->...i", pressure, velocity)
        internal = 0.5 * np.trace(pressure, axis1=-2, axis2=-1)[..., None] * velocity
        thermal = pressure_work + internal
    else:
        thermal = 2.5 * pressure[..., None] * velocity

    result = bulk + thermal
    if heat_flux_w_m2 is not None:
        heat_flux = np.asarray(heat_flux_w_m2, dtype=float)
        if heat_flux.shape != result.shape:
            raise ValueError("heat_flux_w_m2 must match the output vector shape")
        result = result + heat_flux
    return result


def cao_two_component_correction(fast_ion_fraction):
    """Cao et al. (2013) two-component ratio ``Q_K / Q_MHD``.

    The model ratio is ``1 / (3 r - 2 r^2)``, where ``r = N_fast/N``.
    It is a diagnostic closure, not a replacement for integrating a measured
    velocity distribution.  A zero fraction returns infinity.
    """
    fraction = np.asarray(fast_ion_fraction, dtype=float)
    if np.any((fraction < 0) | (fraction > 1)):
        raise ValueError("fast_ion_fraction must lie between zero and one")
    with np.errstate(divide="ignore", invalid="ignore"):
        return 1.0 / (3.0 * fraction - 2.0 * fraction**2)


def poynting_flux(electric_field_v_m, magnetic_field_t):
    """Return the SI Poynting vector ``E cross B / mu0`` in W m^-2."""
    electric = np.asarray(electric_field_v_m, dtype=float)
    magnetic = np.asarray(magnetic_field_t, dtype=float)
    if electric.shape[-1:] != (3,) or magnetic.shape[-1:] != (3,):
        raise ValueError("electric and magnetic fields must end in 3 components")
    return np.cross(electric, magnetic) / MU_0


def convective_electric_field(bulk_velocity_mps, magnetic_field_t):
    """Return the ideal-MHD estimate ``E = -U cross B`` in V m^-1."""
    velocity = np.asarray(bulk_velocity_mps, dtype=float)
    magnetic = np.asarray(magnetic_field_t, dtype=float)
    if velocity.shape[-1:] != (3,) or magnetic.shape[-1:] != (3,):
        raise ValueError("velocity and magnetic field must end in 3 components")
    return -np.cross(velocity, magnetic)


def build_los_response(
    path_length_m,
    neutral_density_m3,
    energy_keV,
    *,
    survival_probability=1.0,
    cross_section_cm2=None,
):
    """Build an optically thin ENA response matrix.

    Parameters
    ----------
    path_length_m : array_like, shape (n_rays, n_voxels)
        Distance travelled by each ray in each source voxel.
    neutral_density_m3 : array_like
        Neutral density with shape ``(n_voxels,)``, ``(n_rays, n_voxels)``,
        or ``(n_energy, n_rays, n_voxels)``.
    energy_keV : array_like, shape (n_energy,)
        ENA energy.  Charge exchange is approximately energy preserving.
    survival_probability : float or array_like
        Optional loss/attenuation factor broadcastable to the response.
    cross_section_cm2 : array_like, optional
        Charge-exchange cross-section for each energy in cm^2.  The default is
        :func:`libs.plasma.sigma_HH` for H+ + H.  Pass a species-appropriate
        array for oxygen or other parent ions.

    Returns
    -------
    ndarray, shape (n_energy, n_rays, n_voxels)
        Dimensionless response multiplying parent directional intensity.
    """
    path = np.asarray(path_length_m, dtype=float)
    energy = np.atleast_1d(np.asarray(energy_keV, dtype=float))
    neutral = np.asarray(neutral_density_m3, dtype=float)
    if path.ndim != 2:
        raise ValueError("path_length_m must have shape (n_rays, n_voxels)")
    n_rays, n_voxels = path.shape
    if neutral.shape == (n_voxels,):
        neutral = neutral[None, None, :]
    elif neutral.shape == (n_rays, n_voxels):
        neutral = neutral[None, :, :]
    elif neutral.shape != (len(energy), n_rays, n_voxels):
        raise ValueError("neutral_density_m3 has an incompatible shape")

    sigma_cm2 = sigma_HH(energy) if cross_section_cm2 is None else cross_section_cm2
    sigma_cm2 = np.asarray(sigma_cm2, dtype=float)
    try:
        sigma_cm2 = np.broadcast_to(sigma_cm2, energy.shape)
    except ValueError as exc:
        raise ValueError("cross_section_cm2 must broadcast to energy_keV") from exc
    sigma_m2 = sigma_cm2 * 1e-4
    response = sigma_m2[:, None, None] * path[None, :, :] * neutral
    try:
        response = response * np.broadcast_to(
            np.asarray(survival_probability, dtype=float), response.shape
        )
    except ValueError as exc:
        raise ValueError("survival_probability cannot broadcast to the response") from exc
    return response


def forward_ena_intensity(parent_intensity, response_matrix):
    """Apply a LOS response to parent directional intensity.

    ``parent_intensity`` has shape ``(n_energy, n_voxels)`` and the response
    has shape ``(n_energy, n_rays, n_voxels)``.  The output is ENA intensity
    with shape ``(n_energy, n_rays)`` in the same differential-intensity units.
    """
    parent = np.asarray(parent_intensity, dtype=float)
    response = np.asarray(response_matrix, dtype=float)
    if response.ndim != 3 or parent.shape != (response.shape[0], response.shape[2]):
        raise ValueError("parent_intensity and response_matrix shapes are incompatible")
    return np.einsum("erv,ev->er", response, parent)


def column_averaged_parent_intensity(ena_intensity, optical_depth):
    """Recover the charge-exchange-weighted parent-ion intensity.

    For an optically thin line of sight,

    ``j_ENA = integral(j_i n_H sigma P_surv ds)``.

    Dividing by ``tau = integral(n_H sigma P_surv ds)`` therefore returns the
    kernel-weighted column mean of ``j_i``.  Both input and output use the same
    differential-intensity units, normally m^-2 s^-1 sr^-1 keV^-1.  Spatially
    varying parent intensity requires :func:`invert_ena_intensity` instead.
    """
    ena = np.asarray(ena_intensity, dtype=float)
    tau = np.asarray(optical_depth, dtype=float)
    try:
        tau = np.broadcast_to(tau, ena.shape)
    except ValueError as exc:
        raise ValueError("optical_depth must broadcast to ena_intensity") from exc
    if np.any(~np.isfinite(tau)) or np.any(tau <= 0):
        raise ValueError("optical_depth must be finite and strictly positive")
    return ena / tau


def invert_ena_intensity(
    observed_intensity,
    response_matrix,
    *,
    uncertainty=None,
    regularization=0.0,
    regularization_matrix=None,
    nonnegative=True,
):
    """Invert ENA intensity independently in each energy channel.

    This solves a weighted, Tikhonov-regularized least-squares problem.  The
    response should include charge exchange, neutral density, path length,
    survival probability, and any known instrumental mixing.  Smoothness can
    be imposed by passing a finite-difference ``regularization_matrix``.
    """
    observed = np.asarray(observed_intensity, dtype=float)
    response = np.asarray(response_matrix, dtype=float)
    if response.ndim != 3 or observed.shape != response.shape[:2]:
        raise ValueError("observed_intensity must match response energy/ray axes")
    if regularization < 0:
        raise ValueError("regularization must be non-negative")

    sigma = np.ones_like(observed) if uncertainty is None else np.asarray(uncertainty, dtype=float)
    try:
        sigma = np.broadcast_to(sigma, observed.shape)
    except ValueError as exc:
        raise ValueError("uncertainty must broadcast to observed_intensity") from exc

    n_energy, _, n_voxels = response.shape
    if regularization_matrix is None:
        reg = np.eye(n_voxels)
    else:
        reg = np.asarray(regularization_matrix, dtype=float)
        if reg.ndim != 2 or reg.shape[1] != n_voxels:
            raise ValueError("regularization_matrix must have n_voxels columns")

    parent = np.full((n_energy, n_voxels), np.nan)
    modeled = np.full_like(observed, np.nan)
    residual = np.full_like(observed, np.nan)
    chi2 = np.full(n_energy, np.nan)
    success = np.zeros(n_energy, dtype=bool)
    messages: list[str] = []
    bounds = (0.0, np.inf) if nonnegative else (-np.inf, np.inf)

    for index in range(n_energy):
        valid = (
            np.isfinite(observed[index])
            & np.isfinite(sigma[index])
            & (sigma[index] > 0)
            & np.all(np.isfinite(response[index]), axis=1)
        )
        if not np.any(valid):
            messages.append("no finite rays")
            continue
        matrix = response[index, valid] / sigma[index, valid, None]
        target = observed[index, valid] / sigma[index, valid]
        if regularization > 0:
            matrix = np.vstack((matrix, np.sqrt(regularization) * reg))
            target = np.concatenate((target, np.zeros(reg.shape[0])))
        fit = lsq_linear(matrix, target, bounds=bounds, method="trf")
        parent[index] = fit.x
        modeled[index] = response[index] @ fit.x
        residual[index, valid] = (
            modeled[index, valid] - observed[index, valid]
        ) / sigma[index, valid]
        dof = max(int(np.count_nonzero(valid)) - n_voxels, 1)
        chi2[index] = np.nansum(residual[index, valid] ** 2) / dof
        success[index] = fit.success
        messages.append(fit.message)

    return ENAInversionResult(
        parent_intensity=parent,
        modeled_intensity=modeled,
        normalized_residual=residual,
        reduced_chi2=chi2,
        success=success,
        messages=tuple(messages),
    )


def velocity_distribution_from_intensity(
    differential_intensity,
    energy_keV,
    *,
    particle_mass_kg,
):
    """Convert directional number intensity ``j_E`` to phase-space density.

    The input energy axis is the first axis and intensity units are
    m^-2 s^-1 sr^-1 keV^-1.  The non-relativistic relation is

    ``j_E[J^-1] = v^2 f / m``  and hence  ``f = m j_E[J^-1] / v^2``.

    The returned SI phase-space density has units s^3 m^-6.  In particular,
    ``f`` is not ``v^2 j_E``; conventions using energy flux rather than number
    intensity require another factor of energy.
    """
    intensity = np.asarray(differential_intensity, dtype=float)
    energy = np.atleast_1d(np.asarray(energy_keV, dtype=float))
    if intensity.ndim < 1 or intensity.shape[0] != energy.size:
        raise ValueError("the first intensity axis must match energy_keV")
    if np.any(energy <= 0):
        raise ValueError("energy_keV must be strictly positive")
    mass = float(particle_mass_kg)
    energy_j = energy * KEV_TO_J
    speed2 = 2.0 * energy_j / mass
    reshape = (energy.size,) + (1,) * (intensity.ndim - 1)
    intensity_per_joule = intensity / KEV_TO_J
    return mass * intensity_per_joule / speed2.reshape(reshape)


def _validate_directional_intensity(
    differential_intensity, energy_edges_keV, direction_vectors, solid_angle_sr
):
    intensity = np.asarray(differential_intensity, dtype=float)
    edges = np.asarray(energy_edges_keV, dtype=float)
    directions = np.asarray(direction_vectors, dtype=float)
    omega = np.asarray(solid_angle_sr, dtype=float)
    if intensity.ndim < 2:
        raise ValueError("differential_intensity must end in (energy, direction)")
    if edges.shape != (intensity.shape[-2] + 1,):
        raise ValueError("energy_edges_keV must have n_energy + 1 values")
    if np.any(np.diff(edges) <= 0) or np.any(edges < 0):
        raise ValueError("energy bin edges must be non-negative and increasing")
    if directions.shape != (intensity.shape[-1], 3):
        raise ValueError("direction_vectors shape is incompatible with intensity")
    omega = np.broadcast_to(omega, (intensity.shape[-1],))
    if np.any(omega < 0):
        raise ValueError("solid_angle_sr must be non-negative")
    norms = np.linalg.norm(directions, axis=1)
    if np.any(~np.isfinite(norms)) or np.any(norms == 0):
        raise ValueError("direction_vectors must be finite and non-zero")
    unit = directions / norms[:, None]
    widths_keV = np.diff(edges)
    centers_j = 0.5 * (edges[:-1] + edges[1:]) * KEV_TO_J
    return intensity, widths_keV, centers_j, unit, omega


def directional_number_density(
    differential_intensity,
    energy_edges_keV,
    direction_vectors,
    solid_angle_sr,
    *,
    particle_mass_kg,
):
    """Integrate directional number intensity into number density in m^-3.

    The evaluated angular/energy domain contributes
    ``n = integral(j_E / v dE dOmega)``.  Missing directions or energies are
    not silently extrapolated; supply closure weights explicitly if desired.
    """
    intensity, widths, centers_j, _, omega = _validate_directional_intensity(
        differential_intensity, energy_edges_keV, direction_vectors, solid_angle_sr
    )
    speed = np.sqrt(2.0 * centers_j / float(particle_mass_kg))
    weights = widths / speed
    return np.einsum("...ed,e,d->...", intensity, weights, omega)


def directional_energy_density(
    differential_intensity,
    energy_edges_keV,
    direction_vectors,
    solid_angle_sr,
    *,
    particle_mass_kg,
):
    """Integrate directional number intensity into kinetic energy density.

    Returns J m^-3 from ``u = integral(E j_E / v dE dOmega)`` over the supplied
    energy and angular coverage.
    """
    intensity, widths, centers_j, _, omega = _validate_directional_intensity(
        differential_intensity, energy_edges_keV, direction_vectors, solid_angle_sr
    )
    speed = np.sqrt(2.0 * centers_j / float(particle_mass_kg))
    weights = widths * centers_j / speed
    return np.einsum("...ed,e,d->...", intensity, weights, omega)


def directional_energy_flux(
    differential_intensity,
    energy_edges_keV,
    direction_vectors,
    solid_angle_sr,
):
    """Integrate directional parent-ion intensity into energy flux.

    Parameters
    ----------
    differential_intensity : array_like, shape (..., n_energy, n_direction)
        Parent-ion intensity in m^-2 s^-1 sr^-1 keV^-1.
    energy_edges_keV : array_like, shape (n_energy + 1,)
        Energy-bin edges.
    direction_vectors : array_like, shape (n_direction, 3)
        Unit particle-velocity directions in the desired coordinate frame.
    solid_angle_sr : float or array_like, shape (n_direction,)
        Angular quadrature weights.  For incomplete viewing these must include
        an angular closure inferred from an in-situ VDF; detector pixel solid
        angles alone do not recover the unobserved distribution.
    """
    intensity, widths_keV, centers_j, unit, omega = _validate_directional_intensity(
        differential_intensity, energy_edges_keV, direction_vectors, solid_angle_sr
    )
    energy_weights = widths_keV * centers_j
    return np.einsum("...ed,e,d,di->...i", intensity, energy_weights, omega, unit)


def scale_spatial_template(relative_template, anchor_flux, template_at_anchor):
    """Scale an ENA-derived spatial template to a kinetic in-situ anchor.

    The last axis is interpreted as spatial cell.  ``anchor_flux`` and
    ``template_at_anchor`` must broadcast to the leading axes.  This operation
    uses ENA for relative extent/structure and the in-situ VDF for absolute
    particle-energy flux, avoiding an unjustified absolute ENA normalization.
    """
    template = np.asarray(relative_template, dtype=float)
    anchor = np.asarray(anchor_flux, dtype=float)
    template_anchor = np.asarray(template_at_anchor, dtype=float)
    if np.any(~np.isfinite(template_anchor)) or np.any(template_anchor == 0):
        raise ValueError("template_at_anchor must be finite and non-zero")
    scale = anchor / template_anchor
    return template * scale[..., None]


def integrate_cross_section(
    particle_flux_w_m2,
    electromagnetic_flux_w_m2,
    cell_area_m2,
    *,
    surface_normal=(1.0, 0.0, 0.0),
    activity_fraction=1.0,
):
    """Integrate particle and EM flux maps over a Y-Z-like surface.

    Flux arrays must end in ``(n_cells, 3)``.  Area and activity fraction are
    applied cell by cell; the leading axes (commonly time) are preserved.
    """
    particle = np.asarray(particle_flux_w_m2, dtype=float)
    electromagnetic = np.asarray(electromagnetic_flux_w_m2, dtype=float)
    area = np.asarray(cell_area_m2, dtype=float)
    normal = np.asarray(surface_normal, dtype=float)
    fraction = np.asarray(activity_fraction, dtype=float)
    if particle.shape != electromagnetic.shape or particle.shape[-1] != 3:
        raise ValueError("particle and electromagnetic flux maps must share (..., cells, 3)")
    if area.shape != (particle.shape[-2],):
        raise ValueError("cell_area_m2 must contain one area per spatial cell")
    if normal.shape != (3,) or not np.isfinite(normal).all() or np.linalg.norm(normal) == 0:
        raise ValueError("surface_normal must be a finite three-vector")
    normal = normal / np.linalg.norm(normal)
    try:
        fraction = np.broadcast_to(fraction, particle.shape[:-1])
    except ValueError as exc:
        raise ValueError("activity_fraction must broadcast to (..., cells)") from exc
    weights = area * fraction
    particle_normal = np.einsum("...ci,i->...c", particle, normal)
    em_normal = np.einsum("...ci,i->...c", electromagnetic, normal)
    particle_power = np.sum(particle_normal * weights, axis=-1)
    em_power = np.sum(em_normal * weights, axis=-1)
    return TransportBudget(
        particle_power_w=particle_power,
        electromagnetic_power_w=em_power,
        total_power_w=particle_power + em_power,
    )


def integrate_event_energy(time_s, power_w):
    """Integrate a power time series with the trapezoid rule and return joules."""
    time = np.asarray(time_s, dtype=float)
    power = np.asarray(power_w, dtype=float)
    if time.ndim != 1 or power.shape[0] != time.size:
        raise ValueError("power_w must use time as its first axis")
    if np.any(np.diff(time) <= 0):
        raise ValueError("time_s must be strictly increasing")
    # ``np.trapezoid`` is unavailable in the NumPy 1.x runtime used by some
    # existing GENA environments.
    return np.trapezoid(power, x=time, axis=0)
