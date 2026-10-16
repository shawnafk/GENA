import numpy as np

from common.constants import e_charge, mass_H_kg
from libs.distribution import rest_2_lab, rest_flux, vb_from_T_lab


def test_rest_flux_pure_maxwellian_has_expected_normalization_and_peak():
    energy_keV = np.array([0.0, 2.0, 4.0, 8.0])
    density_cm3 = 0.5
    temperature_keV = 4.0

    flux = rest_flux(energy_keV, density_cm3, temperature_keV)

    energy_eV = 1.0e3 * energy_keV
    temperature_eV = 1.0e3 * temperature_keV
    prefactor = 1.0e2 * np.sqrt(e_charge / (2.0 * mass_H_kg)) / np.pi**1.5
    expected = (
        1.0e3
        * prefactor
        * density_cm3
        * energy_eV
        / temperature_eV**1.5
        * np.exp(-energy_eV / temperature_eV)
    )
    np.testing.assert_allclose(flux, expected)
    assert np.argmax(flux) == 2


def test_rest_flux_adds_maxwellian_and_kappa_components():
    energy_keV = np.array([2.0, 5.0, 10.0, 20.0])

    combined = rest_flux(energy_keV, 0.45, 4.5, 0.12, 10.0, 4.0)
    maxwellian = rest_flux(energy_keV, 0.45, 4.5)
    kappa = rest_flux(energy_keV, 0.0, 1.0, 0.12, 10.0, 4.0)

    np.testing.assert_allclose(combined, maxwellian + kappa)


def test_rest_flux_requires_physical_kappa_tail_parameters():
    with np.testing.assert_raises_regex(ValueError, "Kappa tail"):
        rest_flux([5.0], 0.4, 5.0, 0.1, 8.0, 1.5)


def test_rest_2_lab_zero_bulk_is_the_rest_flux():
    energy_keV = np.geomspace(0.1, 30.0, 20)

    expected = rest_flux(energy_keV, 0.06, 2.6)
    actual = rest_2_lab(energy_keV, 0.0, 0.06, 2.6)

    np.testing.assert_allclose(actual, expected)


def test_rest_2_lab_omnidirectional_flux_depends_only_on_bulk_speed():
    energy_keV = np.geomspace(0.1, 30.0, 20)

    positive = rest_2_lab(energy_keV, 197.0, 0.06, 2.6)
    negative = rest_2_lab(energy_keV, -197.0, 0.06, 2.6)

    np.testing.assert_allclose(positive, negative)
    assert np.all(np.isfinite(positive))
    assert np.all(positive > 0.0)


def test_vb_from_T_lab_recovers_synthetic_bulk_speed_independent_of_scale():
    energy_keV = np.geomspace(0.5, 25.0, 80)
    temperature_keV = 2.6
    expected_speed_kms = 798.0
    energy_flux = (
        37.0
        * energy_keV
        * rest_2_lab(energy_keV, expected_speed_kms, 1.0, temperature_keV)
    )

    fitted_speed_kms = vb_from_T_lab(
        energy_keV,
        energy_flux,
        temperature_keV,
        [0.7, 20.0],
    )

    assert abs(fitted_speed_kms - expected_speed_kms) < 1.0


def test_vb_from_T_lab_requires_three_valid_interval_samples():
    with np.testing.assert_raises_regex(ValueError, "at least three"):
        vb_from_T_lab(
            [1.0, 2.0, 3.0],
            [10.0, np.nan, 30.0],
            2.0,
            [1.0, 3.0],
        )
