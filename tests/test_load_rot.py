"""Exercise PDF attitude loading across legacy and flagged files."""

import numpy as np
import pandas as pd
import pytest
from scipy.spatial.transform import Rotation

from load import load_rot
from load.loaders import gena


def test_load_rot_interpolates_valid_pose_packets_across_pdf_formats(tmp_path, monkeypatch):
    monkeypatch.setattr(gena, "get_configured_paths", lambda: {"pqt_base": tmp_path})
    midnight = pd.Timestamp("2024-10-10T00:00:00Z").timestamp()
    angle_90 = np.sqrt(0.5)
    legacy = pd.DataFrame({
        "Timestamp": [midnight - 2, midnight - 1],
        "MJ2000_q1": [angle_90, 0.0],
        "MJ2000_q2": [0.0, 0.0],
        "MJ2000_q3": [0.0, 0.0],
    })
    current = pd.DataFrame({
        "Timestamp": midnight + np.arange(4),
        "MJ2000_q1": [0.0, angle_90, 1.01, angle_90],
        "MJ2000_q2": [0.0, 0.0, 0.0, 0.0],
        "MJ2000_q3": [0.0, 0.0, 0.0, 0.0],
        "PosePresent": [True, False, True, True],
    })
    legacy.to_parquet(tmp_path / "20241009_pdf.pqt")
    current.to_parquet(tmp_path / "20241010_pdf.pqt")

    times = midnight + np.array([1.5, -1, 0, 1.5])
    quaternions = load_rot("queqiao", times)
    assert quaternions.shape == (4, 4)
    np.testing.assert_allclose(
        Rotation.from_quat(quaternions).as_euler("xyz", degrees=True)[:, 0],
        [45, 45, 0, 45], atol=1e-8,
    )
    np.testing.assert_allclose(load_rot("gena", midnight), [0, 0, 0, 1], atol=1e-12)
    with pytest.raises(ValueError, match="outside PDF coverage"):
        load_rot("qq", midnight + 4)
    partial = load_rot("qq", [midnight + 4, midnight], bounds_error=False)
    assert np.isnan(partial[0]).all()
    np.testing.assert_allclose(partial[1], [0, 0, 0, 1], atol=1e-12)
