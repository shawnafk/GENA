import numpy as np

from libs import flux


def test_j_dual_model_is_the_sum_of_maxwellian_and_kappa_terms():
    energy_eV = np.geomspace(100.0, 100_000.0, 50)

    combined = flux.j_dual_model(energy_eV, 0.04, 1.2e3, 0.002, 12.0e3, 3.5)
    components = (
        flux.j_maxwell(energy_eV, 0.04, 1.2e3)
        + flux.j_kappa(energy_eV, 0.002, 12.0e3, 3.5)
    )

    np.testing.assert_allclose(combined, components)


def test_shared_drift_transforms_both_components():
    energy_eV = np.geomspace(100.0, 100_000.0, 50)

    combined = flux.j_dual_model(
        energy_eV,
        0.04,
        1.2e3,
        0.002,
        12.0e3,
        3.5,
        drift_kms=450.0,
    )
    components = flux._omnidirectional_drift(
        energy_eV,
        lambda energy: flux.j_maxwell(energy, 0.04, 1.2e3),
        450.0,
    ) + flux._omnidirectional_drift(
        energy_eV,
        lambda energy: flux.j_kappa(energy, 0.002, 12.0e3, 3.5),
        450.0,
    )

    np.testing.assert_allclose(combined, components)


def test_constrained_dual_fit_recovers_moments_and_synthetic_spectrum():
    energy_eV = np.geomspace(30.0, 100_000.0, 80)
    n_total_cm3 = 0.05
    T_total_eV = 1_000.0
    tail_fraction = 0.05
    temperature_ratio = 10.0
    kappa = 3.5
    T_maxwell_eV = T_total_eV / (
        1.0 - tail_fraction + tail_fraction * temperature_ratio
    )
    T_kappa_eV = temperature_ratio * T_maxwell_eV
    synthetic_flux = flux.j_dual_model(
        energy_eV,
        (1.0 - tail_fraction) * n_total_cm3,
        T_maxwell_eV,
        tail_fraction * n_total_cm3,
        T_kappa_eV,
        kappa,
    )

    fitted = flux.fit_dual_moment_constrained(
        energy_eV,
        synthetic_flux,
        n_total_cm3,
        T_total_eV,
    )

    np.testing.assert_allclose(fitted["density_sum_cm3"], n_total_cm3)
    np.testing.assert_allclose(fitted["temperature_sum_eV"], T_total_eV)
    np.testing.assert_allclose(fitted["model_particle_flux"], synthetic_flux, rtol=1e-5)


def test_constrained_dual_fit_respects_fixed_maxwell_fraction():
    energy_eV = np.geomspace(30.0, 100_000.0, 80)
    n_total_cm3 = 0.05
    T_total_eV = 1_000.0
    R = 0.9
    temperature_ratio = 10.0
    expected_kappa = 3.5
    expected_drift_kms = 450.0
    T_maxwell_eV = T_total_eV / (
        R + (1.0 - R) * temperature_ratio
    )
    T_kappa_eV = temperature_ratio * T_maxwell_eV
    synthetic_flux = flux.j_dual_model(
        energy_eV,
        R * n_total_cm3,
        T_maxwell_eV,
        (1.0 - R) * n_total_cm3,
        T_kappa_eV,
        expected_kappa,
        drift_kms=expected_drift_kms,
    )
    synthetic_energy_flux = energy_eV * synthetic_flux

    fitted = flux.fit_dual_moment_constrained(
        energy_eV,
        synthetic_energy_flux,
        n_total_cm3,
        T_total_eV,
        en=True,
        R=R,
    )

    np.testing.assert_allclose(fitted["R"], R)
    np.testing.assert_allclose(fitted["n_maxwell_cm3"], R * n_total_cm3)
    np.testing.assert_allclose(fitted["n_kappa_cm3"], (1.0 - R) * n_total_cm3)
    np.testing.assert_allclose(fitted["temperature_sum_eV"], T_total_eV)
    np.testing.assert_allclose(fitted["drift_kms"], expected_drift_kms, rtol=1e-5)
    np.testing.assert_allclose(
        fitted["model_energy_flux"], synthetic_energy_flux, rtol=1e-5
    )


def test_constrained_dual_fit_supports_pure_maxwellian_R_one():
    energy_eV = np.geomspace(30.0, 50_000.0, 70)
    expected_drift_kms = 450.0
    synthetic_flux = flux.j_dual_model(
        energy_eV,
        0.05,
        1_000.0,
        drift_kms=expected_drift_kms,
    )

    fitted = flux.fit_dual_moment_constrained(
        energy_eV,
        synthetic_flux,
        0.05,
        1_000.0,
        R=1.0,
    )

    assert fitted["R"] == 1.0
    assert fitted["n_kappa_cm3"] == 0.0
    assert fitted["T_kappa_eV"] is None
    assert fitted["kappa"] is None
    np.testing.assert_allclose(fitted["drift_kms"], expected_drift_kms, rtol=1e-5)
    np.testing.assert_allclose(fitted["model_particle_flux"], synthetic_flux)


def test_constrained_dual_fit_supports_pure_kappa_R_zero():
    energy_eV = np.geomspace(30.0, 100_000.0, 80)
    expected_kappa = 3.5
    expected_drift_kms = 350.0
    synthetic_flux = flux.j_dual_model(
        energy_eV,
        0.0,
        None,
        0.05,
        1_000.0,
        expected_kappa,
        drift_kms=expected_drift_kms,
    )

    fitted = flux.fit_dual_moment_constrained(
        energy_eV,
        synthetic_flux,
        0.05,
        1_000.0,
        R=0.0,
    )

    assert fitted["R"] == 0.0
    assert fitted["n_maxwell_cm3"] == 0.0
    assert fitted["T_maxwell_eV"] is None
    np.testing.assert_allclose(fitted["kappa"], expected_kappa, rtol=1e-5)
    np.testing.assert_allclose(fitted["drift_kms"], expected_drift_kms, rtol=1e-5)
    np.testing.assert_allclose(fitted["model_particle_flux"], synthetic_flux, rtol=1e-5)


def test_constrained_dual_fit_rejects_R_outside_unit_interval():
    energy_eV = np.geomspace(30.0, 50_000.0, 20)
    measured = flux.j_maxwell(energy_eV, 0.05, 1_000.0)

    with np.testing.assert_raises_regex(ValueError, "R must lie"):
        flux.fit_dual_moment_constrained(
            energy_eV,
            measured,
            0.05,
            1_000.0,
            R=1.1,
        )


def test_manual_partition_preserves_total_moments():
    partition = flux.dual_moment_partition(
        0.05, 1_000.0, maxwell_fraction=0.9, temperature_ratio=6.0
    )

    np.testing.assert_allclose(partition["density_sum_cm3"], 0.05)
    np.testing.assert_allclose(partition["temperature_sum_eV"], 1_000.0)
    np.testing.assert_allclose(partition["n_maxwell_cm3"], 0.045)
    np.testing.assert_allclose(partition["n_kappa_cm3"], 0.005)


def test_fixed_partition_jointly_recovers_kappa_and_drift_from_energy_flux():
    energy_eV = np.geomspace(30.0, 100_000.0, 80)
    partition = flux.dual_moment_partition(
        0.05, 1_000.0, maxwell_fraction=0.9, temperature_ratio=6.0
    )
    expected_kappa = 3.7
    expected_speed_kms = 450.0
    energy_flux = energy_eV * flux.j_dual_model(
        energy_eV,
        partition["n_maxwell_cm3"],
        partition["T_maxwell_eV"],
        partition["n_kappa_cm3"],
        partition["T_kappa_eV"],
        expected_kappa,
        drift_kms=expected_speed_kms,
    )

    fitted = flux.fit_dual_fixed_partition(
        energy_eV,
        energy_flux,
        0.05,
        1_000.0,
        maxwell_fraction=0.9,
        temperature_ratio=6.0,
        en=True,
        drift_bounds=(0.0, 1_000.0),
    )

    np.testing.assert_allclose(fitted["kappa"], expected_kappa, rtol=1.0e-4)
    np.testing.assert_allclose(fitted["drift_kms"], expected_speed_kms, rtol=1.0e-4)
    np.testing.assert_allclose(fitted["model_energy_flux"], energy_flux, rtol=1.0e-4)


def test_fixed_partition_handles_a_pure_maxwellian_without_a_kappa_fit():
    energy_eV = np.geomspace(30.0, 50_000.0, 70)
    expected_speed_kms = 450.0
    energy_flux = energy_eV * flux.j_dual_model(
        energy_eV, 0.05, 1_000.0, drift_kms=expected_speed_kms
    )

    fitted = flux.fit_dual_fixed_partition(
        energy_eV,
        energy_flux,
        0.05,
        1_000.0,
        maxwell_fraction=1.0,
        en=True,
        drift_bounds=(0.0, 1_000.0),
    )

    assert fitted["kappa"] is None
    np.testing.assert_allclose(fitted["drift_kms"], expected_speed_kms, rtol=1.0e-4)


def test_fixed_partition_handles_a_pure_kappa_population():
    energy_eV = np.geomspace(30.0, 100_000.0, 80)
    expected_kappa = 3.5
    expected_speed_kms = 350.0
    energy_flux = energy_eV * flux.j_dual_model(
        energy_eV,
        0.0,
        1.0,
        0.05,
        1_000.0,
        expected_kappa,
        drift_kms=expected_speed_kms,
    )

    fitted = flux.fit_dual_fixed_partition(
        energy_eV,
        energy_flux,
        0.05,
        1_000.0,
        maxwell_fraction=0.0,
        en=True,
        drift_bounds=(0.0, 1_000.0),
    )

    np.testing.assert_allclose(fitted["kappa"], expected_kappa, rtol=1.0e-4)
    np.testing.assert_allclose(fitted["drift_kms"], expected_speed_kms, rtol=1.0e-4)
