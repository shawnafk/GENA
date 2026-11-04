import numpy as np
import pandas as pd

from gena.post.project import get_projection, get_ray
from gena.post import project


def test_projection_preserves_rows_with_missing_pose():
    events = pd.DataFrame(
        {
            "Timestamp": 1720483200.0,
            "X1": 0.0,
            "Y1": 0.0,
            "X2": 1.0,
            "Y2": 2.0,
            "MJ2000_q1": np.array([0.0, np.nan, 0.6, 1.1, 0.0], dtype=np.float32),
            "MJ2000_q2": np.array([0.0, np.nan, 0.8, 0.0, 0.0], dtype=np.float32),
            "MJ2000_q3": 0.0,
            "MJ2000_Star_X": [1000.0, 1000.0, 1000.0, 1000.0, np.nan],
            "MJ2000_Star_Y": 2000.0,
            "MJ2000_Star_Z": 3000.0,
        },
        index=[3, 8, 3, 1, 0],
    )

    xyz = get_projection(events)

    assert xyz.shape == (5, 3)
    assert np.isfinite(xyz[[0, 2]]).all()
    np.testing.assert_allclose(xyz[[0, 2], 1], 0.0, atol=1e-12)
    assert np.isnan(xyz[[1, 3, 4]]).all()
    assert np.isnan(get_ray(events)["GSM"][[1, 3]]).all()
    for row in (0, 2):
        np.testing.assert_allclose(get_projection(events.iloc[[row]]), xyz[[row]])
    assert np.isnan(get_projection(events.iloc[[1]])).all()
    assert get_projection(events.iloc[:0]).shape == (0, 3)


def test_instrument_projection_uses_the_apparent_earth_limb_and_mlt(monkeypatch):
    from scipy.spatial.transform import Rotation

    transform = project.gcd.transform_vectors
    monkeypatch.setattr(project.gcd, "transform_vectors", lambda points, **kwargs:
                        points if kwargs["source"] == "SM" else transform(points, **kwargs))
    def moon(time, *, geocentric):
        assert geocentric
        return np.zeros(3)

    monkeypatch.setattr(project.epms, "moon", moon)
    monkeypatch.setattr(project, "load_pos", lambda mission, time: np.array([0, 0, 10]) * project.RE)
    monkeypatch.setattr(project, "load_rot", lambda mission, time:
                        Rotation.from_euler("x", 5, degrees=True).as_quat())
    result = project.get_instrument_projection("2025-12-04T09:00:00Z", l_shells=(4,))

    np.testing.assert_allclose(result["earth_center"], [0, 0], atol=1e-12)
    angular_radius = np.linalg.norm(np.tan(np.radians(result["earth_limb"])), axis=1)
    np.testing.assert_allclose(angular_radius, np.tan(np.arcsin(0.1)))
    np.testing.assert_allclose(
        [line["px_py"][90] for line in result["field_lines"]],
        np.degrees(np.arctan(np.array([[-.4, 0], [0, -.4], [.4, 0], [0, .4]]))),
        atol=1e-12,
    )
    monkeypatch.setattr(project, "load_rot", lambda mission, time: np.array([0., 0., 0., 1.]))
    for probe, expected_y in (("A", 5), ("B", -5)):
        result = project.get_instrument_projection("2025-12-04T09:00:00Z", prob=probe, l_shells=())
        np.testing.assert_allclose(result["earth_center"], [0, expected_y], atol=1e-12)


def test_fov_geometry_overlay_preserves_the_original_histogram_and_extent():
    from matplotlib import pyplot as plt
    from plots.gena import fov

    events = pd.DataFrame({"PX": [-1, 0, 1, 1], "PY": [0, 0, 1, 1]})
    projection = {
        "earth_limb": np.array([[0., 0.], [0., 1.]]),
        "field_lines": [{"px_py": np.array([[0., 0.], [70., 40.]])}],
    }
    base, base_ax = fov(events, vmin=1, vmax=4)
    overlay, overlay_ax = fov(events, projection=projection, vmin=1, vmax=4)
    try:
        np.testing.assert_array_equal(overlay_ax.collections[0].get_array(), base_ax.collections[0].get_array())
        assert overlay_ax.collections[0].get_clim() == base_ax.collections[0].get_clim()
        assert overlay_ax.get_xlim() == base_ax.get_xlim()
        assert overlay_ax.get_ylim() == base_ax.get_ylim()
        assert len(overlay_ax.lines) == 2
    finally:
        plt.close(base)
        plt.close(overlay)
