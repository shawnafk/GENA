from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import numpy as np

from gena.frame import ParticleDataFrame


SCRIPT = Path(__file__).resolve().parents[1] / "1.Paper_BBF" / "ena_cg_fov_peak" / "analyze_local_fov_peak.py"
SPEC = spec_from_file_location("analyze_local_fov_peak", SCRIPT)
MODULE = module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def synthetic_flux_map(flux):
    flux = np.asarray(flux, dtype=float)
    ones = np.ones_like(flux)
    return {
        "differential_flux": flux,
        "differential_flux_error": np.sqrt(np.clip(flux, 0.0, None)),
        "raw_counts": flux,
        "raw_count_errors": np.sqrt(np.clip(flux, 0.0, None)),
        "efficiency_corrected_counts": flux,
        "efficiency_corrected_count_errors": np.sqrt(np.clip(flux, 0.0, None)),
        "pixel_solid_angle_sr": ones,
        "pixel_geometric_factor_cm2_sr": ones,
        "flux_unit": MODULE.FLUX_UNIT,
    }


def test_load_fov_events_applies_each_agreed_selection(monkeypatch):
    raw = ParticleDataFrame(
        {
            "Energy": [10.0, 10.0, np.nan, 3.0, 60.0, 10.0, 10.0, 10.0],
            "Spec": [1, 2, 1, 1, 1, 1, 1, 1],
            "PX": [0.0, 0.0, 0.0, 0.0, 0.0, 30.0, 0.0, 0.0],
            "PY": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -10.0, 10.0],
            "X1": np.zeros(8),
            "Y1": np.zeros(8),
        }
    )
    requested = {}

    def fake_load(trange, product):
        requested["trange"] = trange
        requested["product"] = product
        return raw

    monkeypatch.setattr(MODULE, "load_gena", fake_load)
    monkeypatch.setattr(MODULE.flt, "block_eff", lambda distance: np.ones(len(raw), dtype=bool))
    selected = MODULE.load_fov_events(["2025-01-01", "2025-01-02"])

    assert requested["product"] == "l2g_a"
    assert selected["Energy"].tolist() == [10.0, 10.0]
    assert selected["PY"].tolist() == [0.0, -10.0]
    assert selected.attrs["selection_counts"] == {
        "l2g_raw": 8,
        "after_block_eff": 8,
        "after_excluding_spec_2": 7,
        "after_finite_energy_px_py": 6,
        "after_energy_range": 4,
        "after_fixed_fov": 2,
    }


def test_peak_pixel_and_five_degree_centroid_recover_gaussian_center():
    rng = np.random.default_rng(20260803)
    true_center = np.array([4.2, -1.7])
    points = rng.normal(true_center, [1.25, 1.1], size=(30000, 2))
    events = ParticleDataFrame({"PX": points[:, 0], "PY": points[:, 1], "Energy": np.full(len(points), 10.0)})
    px_edges, py_edges, px_centers, py_centers, _, _ = MODULE.make_fov_grid()
    flux_map = MODULE.build_fov_differential_flux_map(events, px_edges, py_edges, exposure_seconds=1.0)
    peak = MODULE.find_peak_pixel(flux_map["differential_flux"], px_centers, py_centers)
    weights = MODULE.event_flux_weights(events, px_edges, py_edges, flux_map["pixel_geometric_factor_cm2_sr"])
    centroid = MODULE.refine_peak_centroid(events, peak["px_deg"], peak["py_deg"], weights=weights)

    assert np.hypot(peak["px_deg"] - true_center[0], peak["py_deg"] - true_center[1]) < 1.0
    assert np.hypot(centroid["px_deg"] - true_center[0], centroid["py_deg"] - true_center[1]) < 0.08


def test_constant_map_has_flat_mean_profile_when_rings_are_clipped():
    _, _, _, _, px_grid, py_grid = MODULE.make_fov_grid()
    counts = np.full(px_grid.shape, 10.0)
    distance = MODULE.distance_from_peak(px_grid, py_grid, 27.4, 7.1)
    cells = MODULE.build_local_cell_table(synthetic_flux_map(counts), px_grid, py_grid, distance)
    profile = MODULE.build_radial_profile(cells)

    assert np.allclose(profile["mean_differential_flux"], 10.0)
    assert (profile["coverage_fraction"] < MODULE.MIN_RING_COVERAGE).any()
    assert profile["boundary_truncated"].any()


def test_monotonic_synthetic_peak_has_monotonic_radial_decay():
    _, _, _, _, px_grid, py_grid = MODULE.make_fov_grid()
    distance = MODULE.distance_from_peak(px_grid, py_grid, 0.2, -0.3)
    counts = 100.0 * np.exp(-distance / 2.5)
    cells = MODULE.build_local_cell_table(synthetic_flux_map(counts), px_grid, py_grid, distance)
    profile = MODULE.build_radial_profile(cells)

    assert np.all(np.diff(profile["mean_differential_flux"]) < 0)
    assert profile["usable_coverage"].all()


def test_differential_flux_uses_energy_efficiency_solid_angle_and_exposure():
    events = ParticleDataFrame(
        {
            "PX": [0.1, 0.2, 0.3],
            "PY": [0.1, 0.2, 0.3],
            "Energy": [5.0, 10.0, 20.0],
        }
    )
    px_edges, py_edges, _, _, _, _ = MODULE.make_fov_grid()
    first = MODULE.build_fov_differential_flux_map(events, px_edges, py_edges, exposure_seconds=10.0)
    second = MODULE.build_fov_differential_flux_map(events, px_edges, py_edges, exposure_seconds=20.0)

    assert np.all(first["pixel_solid_angle_sr"] > 0)
    assert np.isclose(first["pixel_geometric_factor_cm2_sr"].sum(), MODULE.GENA_GEOMETRIC_FACTOR_CM2_SR)
    assert np.allclose(second["differential_flux"], 0.5 * first["differential_flux"])
    assert np.allclose(second["differential_flux_error"], 0.5 * first["differential_flux_error"])


def test_energy_theta_cube_has_expected_shape_counts_and_exposure_scaling():
    events = ParticleDataFrame(
        {
            "PX": [0.1, 0.2, -2.4, 12.7, 29.9],
            "PY": [0.1, 0.2, -1.8, 4.6, 9.9],
            "Energy": [5.0, 7.0, 10.0, 20.0, 59.0],
        }
    )
    px_edges, py_edges, _, _, _, _ = MODULE.make_fov_grid()
    energy_edges = np.array([4.0, 8.0, 16.0, 60.0])
    first = MODULE.build_energy_fov_differential_flux_cube(
        events,
        energy_edges,
        px_edges,
        py_edges,
        exposure_seconds=10.0,
    )
    second = MODULE.build_energy_fov_differential_flux_cube(
        events,
        energy_edges,
        px_edges,
        py_edges,
        exposure_seconds=20.0,
    )

    assert first["differential_flux"].shape == (3, 60, 20)
    assert np.array_equal(first["raw_counts"].sum(axis=(1, 2)), [2.0, 1.0, 2.0])
    assert np.allclose(first["raw_count_errors"], np.sqrt(first["raw_counts"]))
    assert np.isclose(first["raw_counts"].sum(), len(events))
    assert np.isclose(first["pixel_geometric_factor_cm2_sr"].sum(), MODULE.GENA_GEOMETRIC_FACTOR_CM2_SR)
    assert np.allclose(second["differential_flux"], 0.5 * first["differential_flux"])
    assert np.allclose(second["differential_flux_error"], 0.5 * first["differential_flux_error"])


def test_theta_offset_grids_are_centered_on_refined_peak():
    _, _, _, _, px_grid, py_grid = MODULE.make_fov_grid()
    peak_px = 4.25
    peak_py = -1.75
    theta1, theta2 = MODULE.theta_offset_grids(px_grid, py_grid, peak_px, peak_py)

    assert np.allclose(theta1, px_grid - peak_px)
    assert np.allclose(theta2, py_grid - peak_py)
    assert np.allclose(np.hypot(theta1, theta2), MODULE.distance_from_peak(px_grid, py_grid, peak_px, peak_py))


def test_energy_theta_cube_table_preserves_every_voxel_and_local_mask():
    events = ParticleDataFrame(
        {
            "PX": [0.1, 1.1, -2.4],
            "PY": [0.1, 1.1, -1.8],
            "Energy": [5.0, 10.0, 20.0],
        }
    )
    px_edges, py_edges, _, _, px_grid, py_grid = MODULE.make_fov_grid()
    energy_edges = np.array([4.0, 8.0, 16.0, 60.0])
    cube = MODULE.build_energy_fov_differential_flux_cube(
        events,
        energy_edges,
        px_edges,
        py_edges,
        exposure_seconds=10.0,
    )
    theta1, theta2 = MODULE.theta_offset_grids(px_grid, py_grid, 0.25, -0.25)
    table = MODULE.energy_theta_cube_table("synthetic", cube, px_grid, py_grid, theta1, theta2)

    assert len(table) == 3 * 60 * 20
    assert table["energy_index"].nunique() == 3
    assert np.isclose(table["raw_counts"].sum(), len(events))
    assert np.allclose(table["raw_count_error"], np.sqrt(table["raw_counts"]))
    assert np.array_equal(table["inside_7p5deg"].to_numpy(), (table["theta_radius_deg"] <= 7.5).to_numpy())
    assert (table["flux_unit"] == MODULE.FLUX_UNIT).all()


def test_interval_comparison_uses_normalized_maps_and_detector_peak_distances():
    first = np.zeros((3, 3))
    first[1, 1] = 10
    second = first * 5
    peaks = {
        "first": {
            "peak_pixel": {"px_deg": 0.0, "py_deg": 0.0},
            "centroid": {"px_deg": 0.1, "py_deg": 0.1},
        },
        "second": {
            "peak_pixel": {"px_deg": 1.0, "py_deg": 0.0},
            "centroid": {"px_deg": 1.1, "py_deg": 0.1},
        },
    }
    comparison = MODULE.compare_interval_maps({"first": first, "second": second}, peaks)

    assert np.isclose(comparison.loc[0, "normalized_map_correlation"], 1.0)
    assert np.isclose(comparison.loc[0, "peak_pixel_distance_deg"], 1.0)
    assert bool(comparison.loc[0, "fixed_hotspot_suspect"])


def test_first_stage_source_has_no_external_coordinate_or_pose_calls():
    source = SCRIPT.read_text(encoding="utf-8")
    forbidden = ["common.coord", "get_spacecraft_pos", "get_pos(", "get_rot(", ".rot(", "GSM"]
    assert not any(term in source for term in forbidden)
