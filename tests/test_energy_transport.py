import numpy as np

from libs.energy_transport import (
    build_los_response,
    cao_two_component_correction,
    column_averaged_parent_intensity,
    convective_electric_field,
    directional_energy_density,
    directional_energy_flux,
    directional_number_density,
    forward_ena_intensity,
    integrate_cross_section,
    integrate_event_energy,
    invert_ena_intensity,
    kinetic_energy_flux,
    mhd_particle_energy_flux,
    poynting_flux,
    scale_spatial_template,
    velocity_distribution_from_intensity,
)
from libs.plasma import M_P


def test_kinetic_third_moment_for_one_velocity_cell():
    density = 2.0e6
    velocity = np.array([[3.0e5, 0.0, 0.0]])
    result = kinetic_energy_flux(
        np.array([density]), velocity, 1.0, particle_mass_kg=M_P
    )
    expected = 0.5 * M_P * density * 3.0e5**3
    np.testing.assert_allclose(result, [expected, 0.0, 0.0])


def test_isotropic_mhd_expression_and_cao_closure():
    density = 1.0e6
    velocity = np.array([4.0e5, 0.0, 0.0])
    pressure = 2.0e-10
    result = mhd_particle_energy_flux(
        density, velocity, pressure, particle_mass_kg=M_P
    )
    expected_x = 0.5 * density * M_P * velocity[0] ** 3 + 2.5 * pressure * velocity[0]
    np.testing.assert_allclose(result, [expected_x, 0.0, 0.0])
    np.testing.assert_allclose(cao_two_component_correction([1.0, 0.2]), [1.0, 1.923076923])


def test_poynting_flux_from_direct_and_convective_electric_field():
    velocity = np.array([1.0e5, 0.0, 0.0])
    magnetic = np.array([0.0, 0.0, 10.0e-9])
    electric = convective_electric_field(velocity, magnetic)
    np.testing.assert_allclose(electric, [0.0, 1.0e-3, 0.0])
    result = poynting_flux(electric, magnetic)
    assert result[0] > 0
    np.testing.assert_allclose(result[1:], 0.0)


def test_los_forward_inverse_recovers_nonnegative_source():
    energy = np.array([5.0, 20.0])
    paths = 6.371e6 * np.array(
        [[1.0, 0.1], [0.2, 1.0], [0.7, 0.4], [0.4, 0.8]]
    )
    response = build_los_response(paths, [2.0e7, 1.0e7], energy)
    truth = np.array([[2.0e10, 1.0e10], [0.8e10, 2.5e10]])
    observed = forward_ena_intensity(truth, response)
    fit = invert_ena_intensity(
        observed,
        response,
        uncertainty=np.maximum(observed * 1e-3, 1.0),
        regularization=0.0,
    )
    assert fit.success.all()
    np.testing.assert_allclose(fit.parent_intensity, truth, rtol=1e-8)
    np.testing.assert_allclose(fit.modeled_intensity, observed, rtol=1e-10)


def test_direct_ena_deconvolution_and_intensity_moments():
    ena = np.array([[2.0, 6.0]])
    parent = column_averaged_parent_intensity(ena, [[0.1, 0.2]])
    np.testing.assert_allclose(parent, [[20.0, 30.0]])

    energy_center_keV = np.array([1.5])
    intensity = np.array([[1.0e10]])
    psd = velocity_distribution_from_intensity(
        intensity, energy_center_keV, particle_mass_kg=M_P
    )
    energy_j = energy_center_keV[0] * 1e3 * 1.602176634e-19
    speed2 = 2.0 * energy_j / M_P
    intensity_per_joule = intensity[0, 0] / (1e3 * 1.602176634e-19)
    np.testing.assert_allclose(psd[0, 0], M_P * intensity_per_joule / speed2)

    directions = [[1.0, 0.0, 0.0]]
    density = directional_number_density(
        intensity, [1.0, 2.0], directions, [1.0], particle_mass_kg=M_P
    )
    energy_density = directional_energy_density(
        intensity, [1.0, 2.0], directions, [1.0], particle_mass_kg=M_P
    )
    energy_flux = directional_energy_flux(intensity, [1.0, 2.0], directions, [1.0])
    speed = np.sqrt(2.0 * energy_j / M_P)
    np.testing.assert_allclose(energy_density, density * energy_j)
    np.testing.assert_allclose(energy_flux[0], energy_density * speed)
    np.testing.assert_allclose(energy_flux[1:], 0.0)


def test_directional_flux_spatial_scaling_and_surface_integration():
    intensity = np.array([[1.0e10, 0.0], [2.0e10, 0.0]])
    q = directional_energy_flux(
        intensity,
        [1.0, 2.0, 4.0],
        [[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]],
        [1.0, 1.0],
    )
    assert q[0] > 0
    np.testing.assert_allclose(q[1:], 0.0)

    template = np.array([[1.0, 2.0]])
    scaled = scale_spatial_template(template, np.array([0.4]), np.array([2.0]))
    np.testing.assert_allclose(scaled, [[0.2, 0.4]])
    particle = np.stack((scaled, np.zeros_like(scaled), np.zeros_like(scaled)), axis=-1)
    electromagnetic = np.zeros_like(particle)
    electromagnetic[..., 0] = 0.1
    budget = integrate_cross_section(particle, electromagnetic, [2.0, 3.0])
    np.testing.assert_allclose(budget.particle_power_w, [1.6])
    np.testing.assert_allclose(budget.electromagnetic_power_w, [0.5])
    np.testing.assert_allclose(budget.total_power_w, [2.1])
    np.testing.assert_allclose(integrate_event_energy([0.0, 10.0], [2.0, 4.0]), 30.0)
