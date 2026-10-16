from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import numpy as np
import pandas as pd


SCRIPT = Path(__file__).resolve().parents[1] / "1.Paper_BBF" / "ena_cg_fov_peak" / "fit_cg_alpha.py"
SPEC = spec_from_file_location("fit_cg_alpha", SCRIPT)
MODULE = module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def synthetic_cube(*, gamma=2.0, speed_kms=300.0, seed=20260803, sample_poisson=True):
    energy = np.geomspace(5.0, 40.0, 6)
    ratio = energy[1] / energy[0]
    edges = np.concatenate([[energy[0] / np.sqrt(ratio)], np.sqrt(energy[:-1] * energy[1:]), [energy[-1] * np.sqrt(ratio)]])
    widths = np.diff(edges)
    centers = np.arange(-5.0, 6.0)
    px, py = np.meshgrid(centers, centers, indexing="ij")
    vectors = MODULE.fov_2_vec(px.ravel(), py.ravel()).reshape(*px.shape, 3)
    direction = MODULE.fov_2_vec(0.0, 0.0)[0]
    mu = vectors @ direction
    geometric_factor = np.full(px.shape, 0.01)
    expected = []
    for energy_value in energy:
        log_factor, _ = MODULE.cg_log_angular_factor(energy_value, gamma, speed_kms, mu.ravel())
        weight = geometric_factor.ravel() * np.exp(log_factor - np.max(log_factor))
        expected.append((60000.0 * weight / weight.sum()).reshape(px.shape))
    expected = np.asarray(expected)
    rng = np.random.default_rng(seed)
    raw_counts = rng.poisson(expected).astype(float) if sample_poisson else expected
    corrected = np.empty_like(raw_counts)
    for index, energy_value in enumerate(energy):
        corrected[index] = energy_value ** (-gamma) * widths[index] * geometric_factor
    corrected_errors = np.sqrt(np.maximum(corrected, 1e-12))
    return {
        "energy_centers_keV": energy,
        "energy_widths_keV": widths,
        "theta1_deg": px,
        "theta2_deg": py,
        "raw_counts": raw_counts,
        "efficiency_corrected_counts": corrected,
        "efficiency_corrected_count_errors": corrected_errors,
        "pixel_geometric_factor_cm2_sr": geometric_factor,
        "peak_px_deg": np.asarray(0.0),
        "peak_py_deg": np.asarray(0.0),
        "exposure_seconds": np.asarray(1.0),
    }


def test_proton_speed_and_alpha_have_expected_energy_scaling():
    speed = MODULE.proton_speed_kms([4.0, 16.0])

    assert np.isclose(speed[1], 2.0 * speed[0])
    assert np.isclose(speed[0], 875.4, rtol=2e-3)
    _, alpha_low = MODULE.cg_log_angular_factor(4.0, 2.0, 400.0, [1.0])
    _, alpha_high = MODULE.cg_log_angular_factor(16.0, 2.0, 400.0, [1.0])
    assert np.isclose(alpha_low, 2.0 * alpha_high)


def test_local_log_slope_recovers_power_law_spectrum():
    cube = synthetic_cube(gamma=2.3)
    spectrum = MODULE.estimate_spectral_index(cube)

    assert np.allclose(spectrum["gamma"], 2.3, atol=1e-10)


def test_profiled_prediction_preserves_counts_in_every_selected_energy_bin():
    cube = synthetic_cube()
    spectrum = pd.DataFrame({"gamma": np.full(6, 2.0)})
    arrays = MODULE._prepare_fit_arrays(cube, spectrum, 5.0, (4.0, 60.0))
    direction = MODULE.fov_2_vec(0.0, 0.0)[0]
    mu = arrays["vectors"] @ direction
    fitted = MODULE.profiled_poisson_prediction(
        arrays["counts"],
        arrays["pixel_gf"],
        mu,
        arrays["energy"],
        arrays["gamma"],
        300.0,
        arrays["energy_mask"],
    )

    assert np.allclose(
        np.nansum(fitted["prediction"], axis=1),
        np.sum(arrays["counts"], axis=1),
    )


def test_fixed_direction_poisson_fit_recovers_injected_speed():
    injected_speed = 300.0
    cube = synthetic_cube(speed_kms=injected_speed, sample_poisson=False)
    spectrum = pd.DataFrame({"gamma": np.full(6, 2.0)})
    fitted = MODULE.fit_cg_speed_fixed_direction(
        cube,
        spectrum,
        radius_deg=5.0,
        energy_range_keV=(4.0, 60.0),
        speed_bounds_kms=(0.0, 500.0),
    )

    assert abs(fitted["speed_kms"] - injected_speed) < 25.0
    assert not fitted["upper_boundary_hit"]


def test_likelihood_interval_marks_open_grid_boundary():
    profile = pd.DataFrame(
        {
            "speed_kms": [0.0, 100.0, 200.0, 300.0],
            "delta_negative_log_likelihood": [8.0, 4.0, 1.0, 0.0],
        }
    )
    interval = MODULE.likelihood_interval(profile)

    assert not interval["lower_open"]
    assert interval["upper_open"]
    assert interval["upper_kms"] == 300.0
