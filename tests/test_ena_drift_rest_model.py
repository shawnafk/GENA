import numpy as np

from gena.post.ena_drift import _rest_model_parameters
from libs import flux


def test_ena_drift_accepts_fixed_partition_fit_as_its_rest_model():
    partition = flux.dual_moment_partition(
        0.05, 1_000.0, maxwell_fraction=0.9, temperature_ratio=6.0
    )
    rest_fit = {**partition, "kappa": 3.5}

    actual = _rest_model_parameters(None, rest_fit)

    expected = (
        partition["n_maxwell_cm3"],
        partition["T_maxwell_eV"],
        partition["n_kappa_cm3"],
        partition["T_kappa_eV"],
        3.5,
    )
    np.testing.assert_allclose(actual, expected)


def test_ena_drift_keeps_the_legacy_pure_maxwellian_fallback():
    actual = _rest_model_parameters(2.5, None)

    np.testing.assert_allclose(actual[:3], (1.0, 2_500.0, 0.0))
    assert actual[3:] == (None, None)
