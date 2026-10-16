from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import numpy as np
import pandas as pd


SCRIPT = Path(__file__).resolve().parents[1] / "1.Paper_BBF" / "ena_cg_fov_peak" / "fit_conditional_cg_short_windows.py"
SPEC = spec_from_file_location("fit_conditional_cg_short_windows", SCRIPT)
MODULE = module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def synthetic_conditional_cube(*, speed_kms=300.0, pixel_amplitudes=None):
    energy = np.sqrt(MODULE.ENERGY_EDGES_KEV[:-1] * MODULE.ENERGY_EDGES_KEV[1:])
    gamma = np.array([1.5, 2.0, 2.5, 2.0])
    coordinates = np.arange(-6.0, 7.0, 2.0)
    px, py = np.meshgrid(coordinates, coordinates, indexing="ij")
    radius = np.hypot(px, py)
    selected = radius <= 7.5
    active_px = px[selected]
    active_py = py[selected]
    vectors = MODULE.fov_2_vec(active_px, active_py)
    direction = MODULE.fov_2_vec(0.0, 0.0)[0]
    mu = vectors @ direction
    if pixel_amplitudes is None:
        pixel_amplitudes = np.exp(0.12 * active_px - 0.08 * active_py + 0.02 * active_px * active_py)
    pixel_amplitudes = np.asarray(pixel_amplitudes, dtype=float)
    energy_amplitudes = np.array([8.0, 12.0, 5.0, 2.0])
    log_factor = np.empty((len(energy), len(mu)))
    for index in range(len(energy)):
        log_factor[index], _ = MODULE.cg_log_angular_factor(energy[index], gamma[index], speed_kms, mu)
    counts = pixel_amplitudes[None, :] * energy_amplitudes[:, None] * np.exp(log_factor)
    return {
        "active_counts": counts,
        "active_vectors": vectors,
        "active_px_deg": active_px,
        "active_py_deg": active_py,
        "active_theta1_deg": active_px,
        "active_theta2_deg": active_py,
        "active_radius_deg": np.hypot(active_px, active_py),
        "energy_centers_keV": energy,
        "gamma": gamma,
        "centroid": {"px_deg": 0.0, "py_deg": 0.0},
    }


def test_strongest_window_is_full_duration_and_maximizes_count():
    start = pd.Timestamp("2025-01-01T00:00:00")
    timestamps = np.concatenate(
        [
            start.timestamp() + np.arange(5) * 60.0,
            start.timestamp() + 1800.0 + np.arange(20) * 20.0,
        ]
    )
    events = pd.DataFrame({"Timestamp": timestamps})
    window = MODULE.find_strongest_window(
        events,
        interval_start=start,
        interval_end=start + pd.Timedelta(hours=1),
    )

    assert window["n_events"] == 20
    assert np.isclose(window["end_unix"] - window["start_unix"], 600.0)


def test_conditional_fit_preserves_every_pixel_total():
    cube = synthetic_conditional_cube()
    zero_log_factor = np.zeros_like(cube["active_counts"])
    fitted = MODULE.conditional_multinomial_fit(cube["active_counts"], zero_log_factor)

    assert np.allclose(fitted["expected_counts"].sum(axis=0), cube["active_counts"].sum(axis=0))
    assert np.allclose(fitted["probability"].sum(axis=0), 1.0)


def test_conditional_speed_fit_recovers_injected_speed_with_arbitrary_zero_order_peak():
    injected_speed = 300.0
    cube = synthetic_conditional_cube(speed_kms=injected_speed)
    fitted = MODULE.fit_conditional_speed_fixed_direction(
        cube,
        speed_bounds_kms=(0.0, 500.0),
    )

    assert abs(fitted["speed_kms"] - injected_speed) < 1.0
    assert not fitted["upper_boundary_hit"]


def test_changing_zero_order_pixel_brightness_does_not_change_recovered_speed():
    first = synthetic_conditional_cube(speed_kms=350.0)
    n_pixels = first["active_counts"].shape[1]
    second_amplitude = np.geomspace(0.1, 20.0, n_pixels)
    second = synthetic_conditional_cube(speed_kms=350.0, pixel_amplitudes=second_amplitude)
    first_fit = MODULE.fit_conditional_speed_fixed_direction(first)
    second_fit = MODULE.fit_conditional_speed_fixed_direction(second)

    assert abs(first_fit["speed_kms"] - 350.0) < 1.0
    assert abs(second_fit["speed_kms"] - 350.0) < 1.0


def test_two_degree_grid_covers_the_fixed_fov_exactly():
    px_edges, py_edges, _, _, px_grid, py_grid = MODULE.make_coarse_fov_grid()

    assert px_grid.shape == (30, 10)
    assert px_edges[0] == -30.0 and px_edges[-1] == 30.0
    assert py_edges[0] == -10.0 and py_edges[-1] == 10.0
