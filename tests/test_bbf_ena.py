import numpy as np
import pytest

from common.constants import e_charge, mass_H_kg
from libs.bbf_ena import (
    adiabatic_rescale_intensity,
    energy_flux_to_counts,
    fold_count_rate,
    forward_bbf_ena,
    integrate_bulk_trajectory,
    rectangular_fov_rays,
)


def test_bulk_trajectory_uses_measured_velocity_series():
    time_s = np.array([0.0, 10.0, 20.0])
    velocity_kms = np.array([[100.0, 0.0, 0.0]] * 3)
    track = integrate_bulk_trajectory(time_s, velocity_kms, [-60.0, 1.0, 2.0])
    expected_dx_re = 100.0 * 20.0 / 6378.1
    np.testing.assert_allclose(track[-1], [-60.0 + expected_dx_re, 1.0, 2.0])


def test_adiabatic_rescaling_preserves_requested_density_ratio():
    energy = np.logspace(-4, 3, 20000)
    reference = energy * np.exp(-energy / 2.0)
    density_ratio = 0.4
    scaled = adiabatic_rescale_intensity(
        energy, reference, energy, density_ratio, gamma=5.0 / 3.0
    )
    speed = np.sqrt(2.0 * energy * 1e3 * e_charge / mass_H_kg)
    reference_density_proxy = np.trapezoid(reference / speed, energy)
    scaled_density_proxy = np.trapezoid(scaled / speed, energy)
    np.testing.assert_allclose(
        scaled_density_proxy / reference_density_proxy,
        density_ratio,
        rtol=2e-3,
    )


def test_rectangular_fov_solid_angle_matches_analytic_pyramid():
    width_deg = 45.0
    height_deg = 20.0
    rays = rectangular_fov_rays(
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        width_deg=width_deg,
        height_deg=height_deg,
        n_horizontal=9,
        n_vertical=5,
    )
    a = np.deg2rad(width_deg) / 2.0
    b = np.deg2rad(height_deg) / 2.0
    expected = 4.0 * np.arctan2(
        np.tan(a) * np.tan(b),
        np.sqrt(1.0 + np.tan(a) ** 2 + np.tan(b) ** 2),
    )
    np.testing.assert_allclose(np.linalg.norm(rays.direction, axis=1), 1.0)
    np.testing.assert_allclose(rays.solid_angle_sr.sum(), expected)


def test_time_of_flight_makes_high_energy_pulse_arrive_first():
    source_time_s = np.arange(-300.0, 301.0, 10.0)
    source_energy_keV = np.logspace(-2, 2, 80)
    pulse = np.exp(-0.5 * (source_time_s / 20.0) ** 2)
    spectrum = np.exp(-source_energy_keV / 10.0)
    source_intensity = 1e12 * pulse[:, None] * spectrum[None, :]
    observation_time_s = np.arange(0.0, 251.0, 5.0)

    result = forward_bbf_ena(
        source_time_s,
        source_energy_keV,
        source_intensity,
        np.zeros((source_time_s.size, 3)),
        observation_time_s,
        [1.0, 4.0],
        [0.0, 0.0, 0.0],
        [[1.0, 0.0, 0.0]],
        initial_center_re=[10.0, 0.0, 0.0],
        fwhm_re=[0.5, 0.5, 0.5],
        max_distance_re=12.0,
        step_re=0.1,
        neutral_density_model=lambda x, y, z: np.ones_like(x) * 1e6,
        cross_section_cm2=lambda energy: 1e-15,
    )
    low_peak_s = observation_time_s[np.argmax(result.ena_intensity_m2_s_sr_keV[:, 0, 0])]
    high_peak_s = observation_time_s[np.argmax(result.ena_intensity_m2_s_sr_keV[:, 1, 0])]
    assert high_peak_s < low_peak_s
    np.testing.assert_allclose(low_peak_s, 145.0, atol=15.0)
    np.testing.assert_allclose(high_peak_s, 75.0, atol=10.0)


def test_detector_fold_uses_full_fov_geometric_factor_once():
    intensity = np.ones((3, 2, 2)) * 1e4
    folded = fold_count_rate(
        intensity,
        [1.0, 2.0, 4.0],
        [1.0, 3.0],
        full_fov_geometric_factor_cm2_sr=20.0,
        relative_efficiency=[1.0, 0.5],
    )
    np.testing.assert_allclose(folded.pixel_geometric_factor_cm2_sr, [5.0, 15.0])
    np.testing.assert_allclose(folded.count_rate_by_energy_s, 20.0)
    np.testing.assert_allclose(folded.total_count_rate_s, 40.0)


def test_energy_flux_to_counts_converts_energy_flux_and_exposure():
    result = energy_flux_to_counts(
        [1000.0, 2000.0],
        [[1000.0, 4000.0], [2000.0, 2000.0]],
        [2.0, 3.0],
        energy_edges_eV=[500.0, 1500.0, 2500.0],
        full_fov_geometric_factor_cm2_sr=10.0,
        relative_efficiency=[1.0, 0.5],
    )

    np.testing.assert_allclose(
        result.particle_intensity_cm2_s_sr_eV,
        [[1.0, 2.0], [2.0, 1.0]],
    )
    np.testing.assert_allclose(
        result.count_rate_by_energy_s,
        [[10000.0, 10000.0], [20000.0, 5000.0]],
    )
    np.testing.assert_allclose(
        result.expected_counts_by_energy,
        [[20000.0, 20000.0], [60000.0, 15000.0]],
    )
    np.testing.assert_allclose(result.total_expected_counts, [40000.0, 75000.0])


def test_energy_flux_to_counts_accepts_descending_energy_channels():
    result = energy_flux_to_counts(
        [2000.0, 1000.0],
        [4000.0, 1000.0],
        2.0,
        energy_edges_eV=[2500.0, 1500.0, 500.0],
        full_fov_geometric_factor_cm2_sr=10.0,
        relative_efficiency=[0.5, 1.0],
    )

    np.testing.assert_allclose(result.energy_eV, [1000.0, 2000.0])
    np.testing.assert_allclose(result.expected_counts_by_energy, [20000.0, 20000.0])


def test_forward_model_accepts_explicit_eulerian_source_track():
    source_time = np.array([0.0, 10.0, 20.0])
    source_center = np.array([
        [-12.0, 0.0, 0.0],
        [-12.0, 0.5, 0.0],
        [-12.0, 1.0, 0.0],
    ])
    result = forward_bbf_ena(
        source_time,
        np.array([5.0, 10.0]),
        np.ones((3, 2)) * 1e8,
        np.ones((3, 3)) * np.array([400.0, 0.0, 0.0]),
        np.array([5.0, 15.0]),
        np.array([7.0]),
        np.array([-10.0, 0.0, 0.0]),
        np.array([[-1.0, 0.0, 0.0]]),
        initial_center_re=np.array([-12.0, 0.0, 0.0]),
        source_center_re=source_center,
        fwhm_re=(4.0, 2.0, 2.0),
        max_distance_re=8.0,
        step_re=1.0,
        include_time_of_flight=False,
        neutral_density_model=lambda x, y, z: np.ones_like(x) * 1e5,
        cross_section_cm2=lambda energy: 1e-15,
    )

    assert np.allclose(result.source_center_re, source_center)
    assert np.all(np.isfinite(result.ena_intensity_m2_s_sr_keV))
    assert np.any(result.ena_intensity_m2_s_sr_keV > 0.0)


def test_forward_model_rejects_invalid_explicit_source_track():
    with pytest.raises(ValueError, match="source_center_re"):
        forward_bbf_ena(
            np.array([0.0, 10.0]),
            np.array([5.0, 10.0]),
            np.ones((2, 2)),
            np.zeros((2, 3)),
            np.array([5.0]),
            np.array([7.0]),
            np.zeros(3),
            np.array([[1.0, 0.0, 0.0]]),
            initial_center_re=np.zeros(3),
            source_center_re=np.zeros((3, 3)),
        )
