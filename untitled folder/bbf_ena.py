"""Data-driven adiabatic forward model for BBF energetic neutral atoms.

The model maps a time series of rest-frame ion differential intensity and
bulk velocity onto a moving, finite-width BBF source.  It then integrates the
optically thin charge-exchange emissivity along detector lines of sight.  The
ENA time of flight and the sign reversal between a detector-to-source line of
sight and the arriving particle velocity are both included explicitly.

All vectors use a single caller-selected Cartesian frame, normally GSM.
Physical inputs and outputs use SI units except for the explicitly named RE,
keV, km/s, and cm^2 quantities.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np
from scipy.interpolate import RegularGridInterpolator

from common.constants import RE, e_charge, mass_H_kg
from libs.plasma import density_atom, sigma_HH


@dataclass(frozen=True)
class FOVRays:
    """Quadrature rays spanning a rectangular tangent-plane field of view."""

    direction: np.ndarray
    solid_angle_sr: np.ndarray


@dataclass(frozen=True)
class ENAForwardResult:
    """Physical output of :func:`forward_bbf_ena`."""

    observation_time_s: np.ndarray
    energy_keV: np.ndarray
    ena_intensity_m2_s_sr_keV: np.ndarray
    source_time_s: np.ndarray
    source_center_re: np.ndarray
    los_unit_vectors: np.ndarray
    ray_solid_angle_sr: np.ndarray


@dataclass(frozen=True)
class ENACountResult:
    """Detector-folded ENA count rates."""

    energy_keV: np.ndarray
    count_rate_by_energy_s: np.ndarray
    total_count_rate_s: np.ndarray
    pixel_geometric_factor_cm2_sr: np.ndarray


@dataclass(frozen=True)
class ENASpectrumCountResult:
    """Expected detector counts folded from an ENA energy-flux spectrum."""

    energy_eV: np.ndarray
    energy_edges_eV: np.ndarray
    particle_intensity_cm2_s_sr_eV: np.ndarray
    count_rate_by_energy_s: np.ndarray
    expected_counts_by_energy: np.ndarray
    total_count_rate_s: np.ndarray
    total_expected_counts: np.ndarray


def integrate_bulk_trajectory(time_s, bulk_velocity_kms, initial_position_re):
    """Integrate a measured bulk-velocity series into a source-center track."""

    time = np.asarray(time_s, dtype=float)
    velocity = np.asarray(bulk_velocity_kms, dtype=float)
    initial = np.asarray(initial_position_re, dtype=float)
    if time.ndim != 1 or time.size < 2 or np.any(np.diff(time) <= 0):
        raise ValueError("time_s must be a strictly increasing one-dimensional array")
    if velocity.shape != (time.size, 3) or not np.all(np.isfinite(velocity)):
        raise ValueError("bulk_velocity_kms must have shape (n_time, 3) and be finite")
    if initial.shape != (3,) or not np.all(np.isfinite(initial)):
        raise ValueError("initial_position_re must be a finite three-vector")

    increments_re = (
        0.5
        * (velocity[:-1] + velocity[1:])
        * np.diff(time)[:, None]
        * 1e3
        / RE
    )
    track = np.empty_like(velocity, dtype=float)
    track[0] = initial
    track[1:] = initial + np.cumsum(increments_re, axis=0)
    return track


def adiabatic_rescale_intensity(
    reference_energy_keV,
    reference_intensity,
    output_energy_keV,
    density_ratio,
    *,
    gamma=5.0 / 3.0,
):
    """Rescale an isotropic spectrum under polytropic adiabatic evolution.

    For ``c = n/n0`` and ``a = T/T0 = c**(gamma - 1)``, conservation of
    phase-space density gives

    ``j(E) = c / sqrt(a) * j0(E/a)``.

    Intensities may use any consistent differential-intensity area and energy
    units.  Values queried outside the reference energy range are returned as
    zero rather than extrapolated.
    """

    energy0 = np.asarray(reference_energy_keV, dtype=float)
    intensity0 = np.asarray(reference_intensity, dtype=float)
    energy = np.asarray(output_energy_keV, dtype=float)
    ratio = np.asarray(density_ratio, dtype=float)
    if energy0.ndim != 1 or energy0.size < 2 or np.any(np.diff(energy0) <= 0):
        raise ValueError("reference_energy_keV must be strictly increasing")
    if intensity0.shape != energy0.shape or np.any(intensity0 < 0):
        raise ValueError("reference_intensity must be non-negative and match energy")
    if np.any(energy <= 0) or np.any(ratio < 0) or gamma <= 1:
        raise ValueError("energies must be positive, density_ratio non-negative, gamma > 1")

    scale = np.power(ratio, gamma - 1.0)
    query = np.divide(energy, scale, out=np.full(np.broadcast_shapes(energy.shape, scale.shape), np.inf), where=scale > 0)
    positive = intensity0 > 0
    if np.count_nonzero(positive) < 2:
        return np.zeros_like(query, dtype=float)
    log_value = np.interp(
        np.log(query),
        np.log(energy0[positive]),
        np.log(intensity0[positive]),
        left=np.nan,
        right=np.nan,
    )
    result = np.exp(log_value)
    result *= np.divide(ratio, np.sqrt(scale), out=np.zeros_like(query), where=scale > 0)
    return np.where(np.isfinite(result), result, 0.0)


def rectangular_fov_rays(
    boresight,
    horizontal_axis,
    *,
    width_deg,
    height_deg,
    n_horizontal,
    n_vertical,
):
    """Return normalized center rays and exact pixel solid angles.

    The rectangle is uniform in tangent-plane coordinates.  ``boresight`` and
    ``horizontal_axis`` use the desired output coordinate frame.
    """

    bore = _unit_vector(boresight, "boresight")
    horizontal = np.asarray(horizontal_axis, dtype=float)
    if horizontal.shape != (3,):
        raise ValueError("horizontal_axis must be a three-vector")
    horizontal = horizontal - np.dot(horizontal, bore) * bore
    horizontal = _unit_vector(horizontal, "horizontal_axis")
    vertical = _unit_vector(np.cross(bore, horizontal), "vertical axis")
    if width_deg <= 0 or height_deg <= 0 or width_deg >= 180 or height_deg >= 180:
        raise ValueError("FOV width and height must lie between 0 and 180 degrees")
    if n_horizontal < 1 or n_vertical < 1:
        raise ValueError("FOV grid dimensions must be positive")

    u_edges = np.linspace(
        -np.tan(np.deg2rad(width_deg) / 2.0),
        np.tan(np.deg2rad(width_deg) / 2.0),
        n_horizontal + 1,
    )
    v_edges = np.linspace(
        -np.tan(np.deg2rad(height_deg) / 2.0),
        np.tan(np.deg2rad(height_deg) / 2.0),
        n_vertical + 1,
    )
    u_center = 0.5 * (u_edges[:-1] + u_edges[1:])
    v_center = 0.5 * (v_edges[:-1] + v_edges[1:])
    uu, vv = np.meshgrid(u_center, v_center, indexing="xy")
    direction = (
        bore[None, None, :]
        + uu[..., None] * horizontal[None, None, :]
        + vv[..., None] * vertical[None, None, :]
    )
    direction /= np.linalg.norm(direction, axis=-1, keepdims=True)

    u0, v0 = np.meshgrid(u_edges[:-1], v_edges[:-1], indexing="xy")
    u1, v1 = np.meshgrid(u_edges[1:], v_edges[1:], indexing="xy")
    solid_angle = (
        _tangent_solid_angle_primitive(u1, v1)
        - _tangent_solid_angle_primitive(u0, v1)
        - _tangent_solid_angle_primitive(u1, v0)
        + _tangent_solid_angle_primitive(u0, v0)
    )
    return FOVRays(
        direction=direction.reshape(-1, 3),
        solid_angle_sr=solid_angle.ravel(),
    )


def forward_bbf_ena(
    source_time_s,
    source_energy_keV,
    source_intensity_m2_s_sr_keV,
    bulk_velocity_kms,
    observation_time_s,
    ena_energy_keV,
    observer_position_re,
    los_unit_vectors,
    *,
    initial_center_re,
    source_center_re=None,
    fwhm_re=(8.0, 3.0, 2.0),
    ray_solid_angle_sr=None,
    max_distance_re=80.0,
    step_re=0.25,
    gamma=5.0 / 3.0,
    flow_aligned=True,
    include_time_of_flight=True,
    neutral_density_model: Callable = density_atom,
    cross_section_cm2: Callable = sigma_HH,
    particle_mass_kg=mass_H_kg,
    ray_chunk_size=16,
):
    """Forward model ENA differential intensity from a moving BBF source.

    ``source_intensity_m2_s_sr_keV`` must describe the isotropic rest-frame
    intensity at the BBF center.  The finite source uses a Gaussian FWHM in
    the parallel/cross-flow/normal directions.  Its spectrum is cooled and
    diluted toward the edge with an isotropic polytropic closure.  Supplying
    spacecraft-frame omnidirectional spectra is only an approximation unless
    they have first been transformed or fitted in the bulk frame.

    ``source_center_re`` may provide an explicit ``(n_source_time, 3)`` source
    envelope track.  When it is omitted, the historical Lagrangian behavior
    is retained and the measured bulk velocity is integrated from
    ``initial_center_re``.  Supplying an explicit track is useful for an
    Eulerian BBF channel, because the plasma bulk velocity need not equal the
    propagation speed of the channel boundary or dipolarization front.

    ``los_unit_vectors`` point from detector to source.  The ion/ENA velocity
    evaluated at the source is therefore ``-v * los_unit_vectors``.
    """

    source_time = np.asarray(source_time_s, dtype=float)
    source_energy = np.asarray(source_energy_keV, dtype=float)
    source_intensity = np.asarray(source_intensity_m2_s_sr_keV, dtype=float)
    bulk_velocity = np.asarray(bulk_velocity_kms, dtype=float)
    observation_time = np.asarray(observation_time_s, dtype=float)
    ena_energy = np.asarray(ena_energy_keV, dtype=float)
    los = np.asarray(los_unit_vectors, dtype=float)
    fwhm = np.asarray(fwhm_re, dtype=float)
    _validate_forward_inputs(
        source_time,
        source_energy,
        source_intensity,
        bulk_velocity,
        observation_time,
        ena_energy,
        los,
        fwhm,
        max_distance_re,
        step_re,
        gamma,
        particle_mass_kg,
        ray_chunk_size,
    )
    los = los / np.linalg.norm(los, axis=1, keepdims=True)
    observer = _broadcast_observer(observer_position_re, observation_time.size)
    ray_solid_angle = _validate_ray_solid_angle(ray_solid_angle_sr, los.shape[0])
    if source_center_re is None:
        center_track = integrate_bulk_trajectory(
            source_time, bulk_velocity, initial_center_re
        )
    else:
        center_track = np.asarray(source_center_re, dtype=float)
        if (
            center_track.shape != (source_time.size, 3)
            or np.any(~np.isfinite(center_track))
        ):
            raise ValueError(
                "source_center_re must be finite with shape (n_source_time, 3)"
            )
    spectrum_interpolator = _log_spectrum_interpolator(
        source_time, source_energy, source_intensity
    )

    distance_re = np.arange(step_re / 2.0, max_distance_re, step_re)
    distance_m = distance_re * RE
    ds_m = step_re * RE
    n_time = observation_time.size
    n_energy = ena_energy.size
    n_ray = los.shape[0]
    result = np.zeros((n_time, n_energy, n_ray), dtype=float)

    for energy_index, energy_keV in enumerate(ena_energy):
        speed_mps = np.sqrt(2.0 * energy_keV * 1e3 * e_charge / particle_mass_kg)
        if include_time_of_flight:
            production_time = observation_time[:, None] - distance_m[None, :] / speed_mps
        else:
            production_time = np.broadcast_to(
                observation_time[:, None], (n_time, distance_re.size)
            )
        center = _interp_series(source_time, center_track, production_time)
        bulk_ms = _interp_series(source_time, bulk_velocity * 1e3, production_time)
        valid_time = np.isfinite(center[..., 0]) & np.isfinite(bulk_ms[..., 0])
        center = np.where(valid_time[..., None], center, 0.0)
        bulk_ms = np.where(valid_time[..., None], bulk_ms, 0.0)

        for ray_start in range(0, n_ray, ray_chunk_size):
            ray_stop = min(ray_start + ray_chunk_size, n_ray)
            ray = los[ray_start:ray_stop]
            position = (
                observer[:, None, None, :]
                + ray[None, :, None, :] * distance_re[None, None, :, None]
            )
            delta = position - center[:, None, :, :]
            local_coordinate = _source_coordinates(
                delta, bulk_ms[:, None, :, :], flow_aligned=flow_aligned
            )
            radius2 = np.sum((local_coordinate / fwhm) ** 2, axis=-1)
            envelope = np.where(
                valid_time[:, None, :],
                np.exp(-4.0 * np.log(2.0) * radius2),
                0.0,
            )

            velocity_lab = -speed_mps * ray[None, :, None, :]
            velocity_rest = velocity_lab - bulk_ms[:, None, :, :]
            energy_rest_keV = (
                0.5
                * particle_mass_kg
                * np.sum(velocity_rest**2, axis=-1)
                / (1e3 * e_charge)
            )
            thermal_scale = np.power(envelope, gamma - 1.0)
            reference_energy = np.divide(
                energy_rest_keV,
                thermal_scale,
                out=np.full_like(energy_rest_keV, np.inf),
                where=thermal_scale > 1e-12,
            )
            query_time = np.broadcast_to(
                production_time[:, None, :], reference_energy.shape
            )
            center_intensity = _evaluate_log_spectrum(
                spectrum_interpolator, query_time, reference_energy
            )
            local_rest_intensity = center_intensity * np.divide(
                envelope,
                np.sqrt(thermal_scale),
                out=np.zeros_like(envelope),
                where=thermal_scale > 1e-12,
            )
            lab_intensity = local_rest_intensity * np.divide(
                energy_keV,
                energy_rest_keV,
                out=np.zeros_like(energy_rest_keV),
                where=energy_rest_keV > 1e-12,
            )

            neutral_m3 = np.asarray(
                neutral_density_model(
                    position[..., 0], position[..., 1], position[..., 2]
                ),
                dtype=float,
            )
            if neutral_m3.shape != lab_intensity.shape:
                neutral_m3 = np.broadcast_to(neutral_m3, lab_intensity.shape)
            sigma_m2 = float(np.asarray(cross_section_cm2(energy_keV))) * 1e-4
            result[:, energy_index, ray_start:ray_stop] = np.sum(
                lab_intensity * neutral_m3 * sigma_m2 * ds_m, axis=-1
            )

    return ENAForwardResult(
        observation_time_s=observation_time,
        energy_keV=ena_energy,
        ena_intensity_m2_s_sr_keV=result,
        source_time_s=source_time,
        source_center_re=center_track,
        los_unit_vectors=los,
        ray_solid_angle_sr=ray_solid_angle,
    )


def fold_count_rate(
    ena_intensity_m2_s_sr_keV,
    energy_edges_keV,
    ray_solid_angle_sr,
    *,
    full_fov_geometric_factor_cm2_sr,
    relative_efficiency=1.0,
):
    """Fold physical ENA intensity through an effective detector response.

    The input intensity shape is ``(..., n_energy, n_ray)``.  The geometric
    factor is distributed between rays in proportion to their solid angle,
    which assumes uniform effective area across the modeled field of view.
    """

    intensity = np.asarray(ena_intensity_m2_s_sr_keV, dtype=float)
    edges = np.asarray(energy_edges_keV, dtype=float)
    solid_angle = np.asarray(ray_solid_angle_sr, dtype=float)
    if intensity.ndim < 2:
        raise ValueError("ENA intensity must end in (energy, ray)")
    if (
        edges.shape != (intensity.shape[-2] + 1,)
        or np.any(edges <= 0)
        or np.any(np.diff(edges) <= 0)
    ):
        raise ValueError("energy_edges_keV must contain n_energy + 1 increasing edges")
    if solid_angle.shape != (intensity.shape[-1],) or np.any(solid_angle <= 0):
        raise ValueError("ray_solid_angle_sr must contain one positive value per ray")
    if full_fov_geometric_factor_cm2_sr <= 0:
        raise ValueError("full_fov_geometric_factor_cm2_sr must be positive")
    efficiency = np.asarray(relative_efficiency, dtype=float)
    try:
        efficiency = np.broadcast_to(efficiency, (intensity.shape[-2],))
    except ValueError as exc:
        raise ValueError("relative_efficiency must broadcast to the energy axis") from exc
    if np.any(~np.isfinite(efficiency)) or np.any(efficiency < 0):
        raise ValueError("relative_efficiency must be finite and non-negative")

    pixel_gf = (
        full_fov_geometric_factor_cm2_sr * solid_angle / np.sum(solid_angle)
    )
    widths_keV = np.diff(edges)
    intensity_cm2 = intensity / 1e4
    count_by_energy = np.einsum(
        "...er,r,e,e->...e",
        intensity_cm2,
        pixel_gf,
        widths_keV,
        efficiency,
    )
    centers = np.sqrt(edges[:-1] * edges[1:])
    return ENACountResult(
        energy_keV=centers,
        count_rate_by_energy_s=count_by_energy,
        total_count_rate_s=np.sum(count_by_energy, axis=-1),
        pixel_geometric_factor_cm2_sr=pixel_gf,
    )


def energy_flux_to_counts(
    energy_eV,
    energy_flux_eV_cm2_s_sr_eV,
    exposure_s,
    *,
    full_fov_geometric_factor_cm2_sr,
    relative_efficiency=1.0,
    energy_edges_eV=None,
):
    """Fold an ENA differential energy-flux spectrum into expected counts.

    ``energy_flux_eV_cm2_s_sr_eV`` must be the directional (or explicitly
    angle-averaged) ENA energy intensity in
    ``eV cm^-2 s^-1 sr^-1 eV^-1``.  In particular, an omnidirectional
    spectrum integrated over solid angle must be converted to an intensity
    before it is supplied here.  The final axis is energy; leading axes may
    represent time or other independent spectra.

    The conversion is ``j(E) = Eflux(E) / E`` followed by
    ``N(E) = j(E) * dE * G * efficiency(E) * exposure``.  The returned
    counts are expectation values, not Poisson-realized integer samples.
    When bin edges are omitted, they are inferred from positive, monotonic
    energy centres with logarithmic midpoints.  Descending energy channels are
    accepted and are returned in increasing-energy order.
    """

    energy = np.asarray(energy_eV, dtype=float)
    energy_flux = np.asarray(energy_flux_eV_cm2_s_sr_eV, dtype=float)
    if energy.ndim != 1 or energy.size == 0 or np.any(~np.isfinite(energy)) or np.any(energy <= 0):
        raise ValueError("energy_eV must be a non-empty one-dimensional array of positive finite values")
    if energy_flux.ndim < 1 or energy_flux.shape[-1] != energy.size:
        raise ValueError("energy_flux_eV_cm2_s_sr_eV must end with the energy axis")
    if np.any(np.isfinite(energy_flux) & (energy_flux < 0)):
        raise ValueError("energy_flux_eV_cm2_s_sr_eV cannot contain negative finite values")

    descending = energy.size > 1 and np.all(np.diff(energy) < 0)
    if energy.size > 1 and not (np.all(np.diff(energy) > 0) or descending):
        raise ValueError("energy_eV must be strictly monotonic")
    if descending:
        energy = energy[::-1]
        energy_flux = energy_flux[..., ::-1]

    if energy_edges_eV is None:
        if energy.size < 2:
            raise ValueError("energy_edges_eV is required when only one energy centre is supplied")
        edges = np.empty(energy.size + 1, dtype=float)
        edges[1:-1] = np.sqrt(energy[:-1] * energy[1:])
        edges[0] = energy[0] ** 2 / edges[1]
        edges[-1] = energy[-1] ** 2 / edges[-2]
    else:
        edges = np.asarray(energy_edges_eV, dtype=float)
        if edges.shape != (energy.size + 1,):
            raise ValueError("energy_edges_eV must contain one more value than energy_eV")
        if descending:
            edges = edges[::-1]
    if np.any(~np.isfinite(edges)) or np.any(edges <= 0) or np.any(np.diff(edges) <= 0):
        raise ValueError("energy_edges_eV must be positive, finite, and strictly increasing")

    geometric_factor = np.asarray(full_fov_geometric_factor_cm2_sr, dtype=float)
    if geometric_factor.ndim != 0 or not np.isfinite(geometric_factor) or geometric_factor <= 0:
        raise ValueError("full_fov_geometric_factor_cm2_sr must be a positive finite scalar")
    efficiency = np.asarray(relative_efficiency, dtype=float)
    try:
        efficiency = np.broadcast_to(efficiency, energy.shape)
    except ValueError as exc:
        raise ValueError("relative_efficiency must broadcast to the energy axis") from exc
    if descending:
        efficiency = efficiency[::-1]
    if np.any(~np.isfinite(efficiency)) or np.any(efficiency < 0):
        raise ValueError("relative_efficiency must be finite and non-negative")

    try:
        exposure = np.broadcast_to(np.asarray(exposure_s, dtype=float), energy_flux.shape[:-1])
    except ValueError as exc:
        raise ValueError("exposure_s must be scalar or broadcast to the leading spectrum axes") from exc
    if np.any(~np.isfinite(exposure)) or np.any(exposure <= 0):
        raise ValueError("exposure_s must be positive and finite")

    particle_intensity = energy_flux / energy
    count_rate = particle_intensity * np.diff(edges) * geometric_factor * efficiency
    expected_counts = count_rate * exposure[..., np.newaxis]
    return ENASpectrumCountResult(
        energy_eV=energy,
        energy_edges_eV=edges,
        particle_intensity_cm2_s_sr_eV=particle_intensity,
        count_rate_by_energy_s=count_rate,
        expected_counts_by_energy=expected_counts,
        total_count_rate_s=np.sum(count_rate, axis=-1),
        total_expected_counts=np.sum(expected_counts, axis=-1),
    )


def _validate_forward_inputs(
    source_time,
    source_energy,
    source_intensity,
    bulk_velocity,
    observation_time,
    ena_energy,
    los,
    fwhm,
    max_distance_re,
    step_re,
    gamma,
    particle_mass_kg,
    ray_chunk_size,
):
    if source_time.ndim != 1 or source_time.size < 2 or np.any(np.diff(source_time) <= 0):
        raise ValueError("source_time_s must be strictly increasing")
    if source_energy.ndim != 1 or source_energy.size < 2 or np.any(np.diff(source_energy) <= 0):
        raise ValueError("source_energy_keV must be strictly increasing")
    if np.any(source_energy <= 0):
        raise ValueError("source energies must be positive")
    if source_intensity.shape != (source_time.size, source_energy.size):
        raise ValueError("source intensity must have shape (n_time, n_energy)")
    if np.any(~np.isfinite(source_intensity)) or np.any(source_intensity < 0):
        raise ValueError("source intensity must be finite and non-negative")
    if bulk_velocity.shape != (source_time.size, 3) or np.any(~np.isfinite(bulk_velocity)):
        raise ValueError("bulk velocity must have shape (n_time, 3) and be finite")
    if observation_time.ndim != 1 or np.any(~np.isfinite(observation_time)):
        raise ValueError("observation_time_s must be a finite one-dimensional array")
    if ena_energy.ndim != 1 or np.any(ena_energy <= 0) or np.any(~np.isfinite(ena_energy)):
        raise ValueError("ena_energy_keV must be a positive one-dimensional array")
    if los.ndim != 2 or los.shape[1] != 3:
        raise ValueError("los_unit_vectors must have shape (n_ray, 3)")
    norm = np.linalg.norm(los, axis=1)
    if np.any(~np.isfinite(norm)) or np.any(norm == 0):
        raise ValueError("LOS vectors must be finite and non-zero")
    if fwhm.shape != (3,) or np.any(fwhm <= 0) or np.any(~np.isfinite(fwhm)):
        raise ValueError("fwhm_re must contain three positive finite values")
    if max_distance_re <= 0 or step_re <= 0 or step_re > max_distance_re:
        raise ValueError("LOS distance and step must be positive with step <= distance")
    if gamma <= 1 or particle_mass_kg <= 0 or ray_chunk_size < 1:
        raise ValueError("gamma > 1, positive mass, and positive ray_chunk_size are required")


def _broadcast_observer(observer_position_re, n_time):
    observer = np.asarray(observer_position_re, dtype=float)
    if observer.shape == (3,):
        observer = np.broadcast_to(observer, (n_time, 3))
    if observer.shape != (n_time, 3) or np.any(~np.isfinite(observer)):
        raise ValueError("observer_position_re must be (3,) or (n_time, 3) and finite")
    return observer


def _validate_ray_solid_angle(ray_solid_angle_sr, n_ray):
    if ray_solid_angle_sr is None:
        return np.full(n_ray, 1.0 / n_ray)
    solid_angle = np.asarray(ray_solid_angle_sr, dtype=float)
    if solid_angle.shape != (n_ray,) or np.any(solid_angle <= 0) or np.any(~np.isfinite(solid_angle)):
        raise ValueError("ray_solid_angle_sr must contain one positive finite value per ray")
    return solid_angle


def _log_spectrum_interpolator(time, energy, intensity):
    positive = intensity[intensity > 0]
    if positive.size == 0:
        log_intensity = np.full_like(intensity, -np.inf)
    else:
        floor = max(float(np.min(positive)) * 1e-12, np.finfo(float).tiny)
        log_intensity = np.log(np.maximum(intensity, floor))
    return RegularGridInterpolator(
        (time, np.log(energy)),
        log_intensity,
        method="linear",
        bounds_error=False,
        fill_value=np.nan,
    )


def _evaluate_log_spectrum(interpolator, time, energy):
    query = np.column_stack((time.ravel(), np.log(energy.ravel())))
    log_value = interpolator(query).reshape(energy.shape)
    return np.exp(np.where(np.isfinite(log_value), log_value, -np.inf))


def _interp_series(time, values, query_time):
    output = np.empty(query_time.shape + (values.shape[1],), dtype=float)
    for component in range(values.shape[1]):
        output[..., component] = np.interp(
            query_time,
            time,
            values[:, component],
            left=np.nan,
            right=np.nan,
        )
    return output


def _source_coordinates(delta, velocity_ms, *, flow_aligned):
    if not flow_aligned:
        return delta
    speed = np.linalg.norm(velocity_ms, axis=-1, keepdims=True)
    parallel = np.divide(
        velocity_ms,
        speed,
        out=np.zeros_like(velocity_ms),
        where=speed > 1.0,
    )
    stationary = speed[..., 0] <= 1.0
    parallel[stationary] = np.array([1.0, 0.0, 0.0])
    reference = np.zeros_like(parallel)
    reference[..., 2] = 1.0
    near_vertical = np.abs(np.sum(parallel * reference, axis=-1)) > 0.9
    reference[near_vertical] = np.array([0.0, 1.0, 0.0])
    cross_flow = np.cross(reference, parallel)
    cross_flow /= np.linalg.norm(cross_flow, axis=-1, keepdims=True)
    normal = np.cross(parallel, cross_flow)
    return np.stack(
        (
            np.sum(delta * parallel, axis=-1),
            np.sum(delta * cross_flow, axis=-1),
            np.sum(delta * normal, axis=-1),
        ),
        axis=-1,
    )


def _unit_vector(vector, name):
    vector = np.asarray(vector, dtype=float)
    if vector.shape != (3,) or not np.all(np.isfinite(vector)):
        raise ValueError(f"{name} must be a finite three-vector")
    norm = np.linalg.norm(vector)
    if norm == 0:
        raise ValueError(f"{name} must be non-zero")
    return vector / norm


def _tangent_solid_angle_primitive(u, v):
    return np.arctan2(u * v, np.sqrt(1.0 + u**2 + v**2))
