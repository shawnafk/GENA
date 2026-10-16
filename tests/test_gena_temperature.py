import numpy as np
import pandas as pd

from common.constants import e_charge, mass_H_kg
from gena.post.statistic import (
    SCIME_ENERGY_BINS_KEV,
    counts_nd,
    counts_to_flux,
    flux_to_temperature,
)


def unit_cross_section(values):
    return np.ones_like(values, dtype=float)


def test_counts_to_flux_uses_weighted_counts_and_energy_axis_position():
    events = pd.DataFrame(
        {
            "PX": [0.5, 0.5, 0.5, 0.5, 1.5, 1.5, 1.5],
            "Energy": [1.25, 1.75, 3.0, 6.0, 1.5, 3.0, 6.0],
        }
    )
    weights = np.array([1.0, 2.0, 6.0, 16.0, 4.0, 10.0, 24.0])
    result = counts_nd(
        events,
        {"PX": [0.0, 1.0, 2.0], "Energy": [1.0, 2.0, 4.0, 8.0]},
        weights=weights,
    )

    energy, flux, flux_errors = counts_to_flux(
        result,
        exposure=2.0,
        geometric_factor=0.5,
    )

    assert result.dims == ("PX", "Energy")
    np.testing.assert_allclose(energy, np.sqrt([2.0, 8.0, 32.0]))
    np.testing.assert_allclose(flux, [[3.0, 3.0, 4.0], [4.0, 5.0, 6.0]])
    np.testing.assert_allclose(flux_errors, [[np.sqrt(5.0), 3.0, 4.0], [4.0, 5.0, 6.0]])


def test_flux_to_temperature_removes_an_arbitrary_energy_axis():
    energy = np.sqrt(SCIME_ENERGY_BINS_KEV[:-1] * SCIME_ENERGY_BINS_KEV[1:])
    expected = np.array([[6.0, 8.0, 10.0], [12.0, 14.0, 16.0]])
    amplitude = np.array([[2.0, 3.0, 4.0], [5.0, 6.0, 7.0]])
    flux = amplitude[:, None, :] * energy[None, :, None] * np.exp(
        -energy[None, :, None] / expected[:, None, :]
    )
    flux_errors = 0.05 * flux

    temperature, temperature_error, slope, slope_error, r_squared = flux_to_temperature(
        (energy, flux, flux_errors),
        energy_axis=1,
        sigma_func=unit_cross_section,
    )

    assert temperature.shape == expected.shape
    np.testing.assert_allclose(temperature, expected)
    np.testing.assert_allclose(slope, -1.0 / expected)
    np.testing.assert_allclose(temperature_error, slope_error * expected**2)
    np.testing.assert_allclose(r_squared, 1.0)
    assert np.all(temperature_error > 0)
    assert np.all(slope_error > 0)


def test_flux_to_temperature_corrects_known_bulk_drift_per_spectrum():
    energy = np.sqrt(SCIME_ENERGY_BINS_KEV[:-1] * SCIME_ENERGY_BINS_KEV[1:])
    expected = np.array([[6.0, 9.0], [12.0, 15.0]])
    drift_velocity_kms = np.array([[400.0, -300.0], [0.0, 500.0]])
    direction_cosine = np.array([[1.0, 0.6], [-0.5, -0.8]])
    drift_energy = (
        0.5 * mass_H_kg * (np.abs(drift_velocity_kms) * 1.0e3) ** 2
        / e_charge
        / 1.0e3
    )
    plasma_energy = (
        energy[:, None, None]
        + drift_energy[None, :, :]
        - 2.0
        * np.sqrt(energy[:, None, None] * drift_energy[None, :, :])
        * direction_cosine[None, :, :]
        * np.sign(drift_velocity_kms)[None, :, :]
    )
    flux = 5.0 * energy[:, None, None] * np.exp(
        -plasma_energy / expected[None, :, :]
    )
    flux_errors = 0.05 * flux

    temperature, temperature_error, slope, slope_error, r_squared = flux_to_temperature(
        (energy, flux, flux_errors),
        drift_velocity_kms=drift_velocity_kms,
        look_direction_cosine=direction_cosine,
        sigma_func=unit_cross_section,
    )

    np.testing.assert_allclose(temperature, expected)
    np.testing.assert_allclose(slope, -1.0 / expected)
    np.testing.assert_allclose(temperature_error, slope_error * expected**2)
    np.testing.assert_allclose(r_squared, 1.0)

    uncorrected, _, _, _, _ = flux_to_temperature(
        (energy, flux[:, 0, 0], flux_errors[:, 0, 0]),
        sigma_func=unit_cross_section,
    )
    assert not np.isclose(uncorrected.item(), expected[0, 0], rtol=0.05)


def test_flux_to_temperature_validates_drift_geometry():
    energy = np.array([8.0, 10.0, 12.0])
    flux = np.ones((3, 2))
    flux_errors = np.ones_like(flux)

    with np.testing.assert_raises_regex(ValueError, "not broadcastable"):
        flux_to_temperature(
            (energy, flux, flux_errors),
            drift_velocity_kms=np.ones(3),
            sigma_func=unit_cross_section,
        )

    with np.testing.assert_raises_regex(ValueError, "between -1 and 1"):
        flux_to_temperature(
            (energy, flux, flux_errors),
            drift_velocity_kms=300.0,
            look_direction_cosine=1.1,
            sigma_func=unit_cross_section,
        )


def test_counts_flux_temperature_pipeline_handles_any_retained_dimensions():
    edges = SCIME_ENERGY_BINS_KEV
    energy = np.sqrt(edges[:-1] * edges[1:])
    widths = np.diff(edges)
    px = np.array([-0.5, 0.5])
    py = np.array([-0.5, 0.5])
    expected = np.array([[6.0, 8.0], [10.0, 12.0]])
    rows = []
    weights = []

    for px_index, px_value in enumerate(px):
        for py_index, py_value in enumerate(py):
            temperature = expected[px_index, py_index]
            for energy_value, width in zip(energy, widths):
                rows.append((energy_value, px_value, py_value))
                flux_value = 5.0 * energy_value * np.exp(-energy_value / temperature)
                weights.append(flux_value * width)

    events = pd.DataFrame(rows, columns=["Energy", "PX", "PY"])
    result = counts_nd(
        events,
        {
            "Energy": edges,
            "PX": [-1.0, 0.0, 1.0],
            "PY": [-1.0, 0.0, 1.0],
        },
        weights=np.asarray(weights),
    )
    combo_flux = counts_to_flux(result)
    energy_out, flux, flux_errors = combo_flux
    temperature, temperature_error, slope, slope_error, r_squared = flux_to_temperature(
        combo_flux,
        energy_axis=result.dims.index("Energy"),
        sigma_func=unit_cross_section,
    )

    assert result.dims == ("Energy", "PX", "PY")
    np.testing.assert_allclose(energy_out, energy)
    assert flux.shape == (len(energy), 2, 2)
    assert flux_errors.shape == flux.shape
    assert temperature.shape == (2, 2)
    np.testing.assert_allclose(temperature, expected)
    np.testing.assert_allclose(slope, -1.0 / expected)
    np.testing.assert_allclose(temperature_error, slope_error * expected**2)
    np.testing.assert_allclose(r_squared, 1.0)


def test_flux_to_temperature_returns_nan_for_insufficient_or_positive_slope():
    energy = np.array([8.0, 10.0, 12.0])
    flux = np.column_stack(
        [
            np.array([1.0, 0.0, 1.0]),
            energy * np.exp(energy),
        ]
    )
    flux_errors = np.ones_like(flux)

    temperature, temperature_error, slope, slope_error, r_squared = flux_to_temperature(
        (energy, flux, flux_errors),
        sigma_func=unit_cross_section,
    )

    np.testing.assert_array_equal(np.isnan(temperature), [True, True])
    np.testing.assert_array_equal(np.isnan(temperature_error), [True, True])
    assert np.isnan(slope[0])
    assert np.isnan(slope_error[0])
    assert np.isnan(r_squared[0])
    np.testing.assert_allclose(slope[1], 1.0)
    assert slope_error[1] > 0
    np.testing.assert_allclose(r_squared[1], 1.0)


def test_flux_to_temperature_rejects_misaligned_flux_errors():
    energy = np.array([8.0, 10.0, 12.0])
    flux = np.ones((3, 2))

    with np.testing.assert_raises_regex(ValueError, "same shape"):
        flux_to_temperature(
            (energy, flux, np.ones(3)),
            sigma_func=unit_cross_section,
        )


def test_counts_to_flux_rejects_nonbroadcastable_exposure():
    result = counts_nd(
        pd.DataFrame({"Energy": [1.5, 3.0]}),
        {"Energy": [1.0, 2.0, 4.0]},
    )

    with np.testing.assert_raises_regex(ValueError, "not broadcastable"):
        counts_to_flux(result, exposure=np.ones(3))
