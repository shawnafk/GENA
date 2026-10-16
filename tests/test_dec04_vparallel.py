import importlib.util
from pathlib import Path

import numpy as np


SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "1.Paper_BBF"
    / "ena_Vbulk"
    / "fit_dec04_secondary_vparallel.py"
)
SPEC = importlib.util.spec_from_file_location("fit_dec04_secondary_vparallel", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_kappa_drift_ratio_is_zero_without_drift():
    energy = np.array([4.0, 8.0, 20.0])
    mu = np.array([-0.5, 0.0, 0.8])
    result = MODULE.kappa_drift_log_ratio(energy, mu, 0.0, 3.5, 3.0)
    np.testing.assert_allclose(result, 0.0)


def test_profile_interval_recovers_quadratic_limits():
    speed = np.arange(0.0, 1001.0, 10.0)
    log_likelihood = -0.5 * ((speed - 500.0) / 100.0) ** 2
    result = MODULE.profile_interval(speed, log_likelihood)
    assert result["best_speed_kms"] == 500.0
    assert 300.0 < result["lower_95_kms"] < 310.0
    assert 690.0 < result["upper_95_kms"] < 700.0
    assert not result["zero_in_95"]
