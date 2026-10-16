import importlib.util
from pathlib import Path

import numpy as np


SCRIPT = Path(__file__).parents[1] / "1.Paper_BBF" / "ENA_pitch_angle_estimation.py"
SPEC = importlib.util.spec_from_file_location("ena_pitch_angle_estimation", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_isotropic_count_coefficients():
    energy = np.array([10_000.0])
    denergy = np.array([1_000.0])
    j_pa = np.full((1, len(MODULE.PA_CENTERS)), 2.0)
    result = MODULE.estimate_count_rates(energy, denergy, j_pa, np.array([2.0]), 60.0)
    omega = result["domega_sr"].sum()
    expected_ion = MODULE.AREA_CM2 * MODULE.EXPOSURE_S * MODULE.ION_LEAK * omega * 2.0 * 1_000.0 * MODULE.EFFICIENCY_AT_10KEV
    expected_ena = MODULE.AREA_CM2 * MODULE.EXPOSURE_S * MODULE.NEUTRAL_H_CM3 * MODULE.LOS_CM * omega * 2.0 * 1_000.0 * MODULE.EFFICIENCY_AT_10KEV * MODULE.sigma_HH(10.0)
    np.testing.assert_allclose(result["ion_counts_min"], expected_ion, rtol=1e-12)
    np.testing.assert_allclose(result["ena_counts_min"], expected_ena, rtol=1e-12)
    np.testing.assert_allclose(result["ion_counts_min_iso"], expected_ion, rtol=1e-12)


def test_synthetic_gyrotropic_distribution_reduces_to_pa_energy():
    pitch = np.deg2rad(np.repeat(np.arange(7.5, 180.0, 15.0), 24))
    gyro = np.tile(np.linspace(0, 2 * np.pi, 24, endpoint=False), 12)
    look = np.column_stack((np.sin(pitch) * np.cos(gyro), np.sin(pitch) * np.sin(gyro), np.cos(pitch)))
    energy = np.array([6_000.0, 12_000.0, 20_000.0])
    particle_flux = (energy[:, None] / energy[0]) ** -2 * (1.0 + 0.8 * np.cos(pitch)[None, :] ** 2)
    record = {
        "energy_ev": energy,
        "denergy_ev": np.array([2_000.0, 4_000.0, 6_000.0]),
        "eflux": particle_flux * energy[:, None],
        "bins": np.ones_like(particle_flux, dtype=bool),
        "look": np.broadcast_to(look, particle_flux.shape + (3,)),
        "domega": np.ones_like(particle_flux),
        "b_dsl": np.array([0.0, 0.0, 1.0]),
    }
    reduced = MODULE.reduce_to_pa_energy(record)
    metrics = MODULE.gyrotropy_metrics(record, reduced)
    assert metrics["gyro_log_rms_dex"] < 0.03
    assert metrics["gyro_first_harmonic"] < 0.03
    assert metrics["pa_explained_fraction"] > 0.95
