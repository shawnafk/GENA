from datetime import timedelta
import io
import json
import zipfile
import numpy as np
import pytest
from scipy.spatial.transform import Rotation
from pydantic import ValidationError
from digitalQQ.backend.models import Scene, Orbit, OrbitSample, ExportRequest
from digitalQQ.backend.physics import (
    attitude,
    solve,
    rotation,
    environment,
    environment_cached,
    convert,
    pose,
)
from digitalQQ.backend.service import export_report


def test_middle_20261021_regression_and_quaternion():
    s = Scene()
    a = attitude(s, True)
    assert a["solution"]["roll"] == pytest.approx(-30.932073968, abs=1e-7)
    assert a["solution"]["pitch"] == pytest.approx(-1.643111289, abs=1e-7)
    assert a["pointing_error"] < 1e-10
    assert a["ray_miss"] < 1e-10
    assert a["target_forward"]
    assert a["sun_angle"] == pytest.approx(96.487978, abs=1e-5)
    q = Rotation.from_quat(a["quaternion"])
    np.testing.assert_allclose(
        convert(q.apply([0, 0, -1]), "ECI2000", "GSM", s.time),
        -np.array(a["axes"])[2],
        atol=1e-12,
    )
    for corners in [a["fov"][:4], a["fov"][4:]]:
        np.testing.assert_allclose(np.linalg.norm(corners, axis=1), 1, atol=1e-12)


def test_rotation_order_and_unreachable():
    v = np.array([0.3, 0.4, 0.5])
    np.testing.assert_allclose(
        rotation(26, -12).apply(v),
        Rotation.from_euler("x", -12, degrees=True).apply(
            Rotation.from_euler("z", 26, degrees=True).apply(v)
        ),
    )
    r, p, e = solve([0.2, 0.1, 1], [-1, 0, 0], True)
    assert p == 0 and e > 70
    assert abs(rotation(r, p).apply(np.array([0.2, 0.1, 1]))[1]) < 1e-12


def test_random_reachable_targets():
    rng = np.random.default_rng(42)
    for _ in range(30):
        a = rng.normal(size=3)
        a /= np.linalg.norm(a)
        b = rotation(*rng.uniform(-180, 180, 2)).apply(a)
        r, p, e = solve(a, b)
        assert e < 1e-10
        np.testing.assert_allclose(rotation(r, p).apply(a), b, atol=1e-12)


def imported():
    s = Scene()
    return s.model_copy(
        update={
            "orbit": Orbit(
                source="imported",
                unit="RE",
                samples=[
                    OrbitSample(
                        time=s.time - timedelta(hours=1),
                        position=(50, 10, 0),
                        quaternion=(0, 0, 0, 1),
                    ),
                    OrbitSample(
                        time=s.time + timedelta(hours=1),
                        position=(52, 12, 0),
                        quaternion=(0, 0, 1, 0),
                    ),
                ],
            )
        }
    )


def test_orbit_interpolation_coverage_and_measured():
    s = imported()
    p, _, _, q = pose(s, s.time)
    np.testing.assert_allclose(p, [51, 11, 0])
    np.testing.assert_allclose(q.apply([1, 0, 0]), [0, 1, 0], atol=1e-12)
    with pytest.raises(ValueError, match="覆盖"):
        attitude(s.model_copy(update={"time": s.time + timedelta(hours=2)}))
    assert attitude(s)["measured_axes"] is not None


def test_hold_modes():
    s = Scene()
    t = s.time + timedelta(hours=6)
    ref = attitude(s)
    inertial = attitude(s.model_copy(update={"time": t, "hold": "inertial"}))
    np.testing.assert_allclose(
        Rotation.from_quat(ref["quaternion"]).as_matrix(),
        Rotation.from_quat(inertial["quaternion"]).as_matrix(),
        atol=1e-12,
    )
    tracked = attitude(s.model_copy(update={"time": t, "hold": "tracking"}))
    assert tracked["pointing_error"] < 1e-10
    earth = attitude(s.model_copy(update={"time": t}))
    assert earth["pointing_error"] > 0.01


def test_validation_and_scene_roundtrip():
    s = Scene()
    assert Scene.model_validate_json(s.model_dump_json()) == s
    with pytest.raises(ValidationError):
        Scene.model_validate({"time": "2026-10-21T00:00"})
    with pytest.raises(ValidationError):
        Scene.model_validate({"weather": {"pressure": -1}})
    with pytest.raises(ValidationError):
        Orbit(source="imported", samples=[])
    with pytest.raises(ValidationError):
        Orbit(
            source="imported",
            samples=[
                OrbitSample(time=s.time, position=(10, 0, 0)),
                OrbitSample(time=s.time, position=(11, 0, 0)),
            ],
        )
    with pytest.raises(ValidationError):
        OrbitSample(time=s.time, position=(10, 0, 0), quaternion=(0, 0, 0, 2))


def test_environment_finite_and_cache_parameter_changes():
    s = Scene(magnetic=False)
    environment_cached.cache_clear()
    e = environment(s)
    for key in ["bowshock", "magnetopause"]:
        grid = np.array(e[key])
        assert np.isfinite(grid).any()
        valid = np.isfinite(grid).all(axis=0)
        assert np.linalg.norm(grid[:, valid], axis=0).max() <= 100
    environment(s.model_copy(update={"roll": 0}))
    assert environment_cached.cache_info().hits == 1
    w = s.weather.model_copy(update={"pressure": 4})
    assert (
        environment(s.model_copy(update={"weather": w}))["magnetopause"]
        != e["magnetopause"]
    )


def test_export_archive_and_exact_midpoint():
    s = Scene(magnetic=False)
    req = ExportRequest(
        scene=s, start=s.time, stop=s.time + timedelta(minutes=25), cadence_minutes=10
    )
    data = export_report(req.model_dump(mode="json"))
    with zipfile.ZipFile(io.BytesIO(data)) as z:
        assert {
            "scene.json",
            "parameters_10min.csv",
            "parameters_midpoint.csv",
            "middle_3d.html",
            "middle_2d.png",
            "middle_3d.png",
            "parameter_history.png",
            "sun_avoidance.png",
            "calculation_metadata.json",
            "观测姿态计算说明.md",
        } <= set(z.namelist())
        assert Scene.model_validate_json(z.read("scene.json")) == s
        assert "16:12:30" in z.read("parameters_midpoint.csv").decode()
        assert z.read("middle_3d.png").startswith(b"\x89PNG")


def test_local_t96_sheet_roots_without_network():
    from digitalQQ.backend.physics import space_models

    s = Scene()
    sheet = space_models().T96NeutralSheet(
        s.time.timestamp(),
        2,
        -20,
        0,
        -3,
        nx=3,
        ny=3,
        nz=7,
        x_range=(-35, -15),
        y_range=(-5, 5),
    )
    x, y, z = sheet.surface_xyz()
    valid = np.isfinite(z)
    assert valid.any()
    for xx, yy, zz in zip(x[valid], y[valid], z[valid]):
        assert abs(sheet.bx(xx, yy, zz)) < 1e-7
    lower, upper = sheet.thickness_xyz()
    np.testing.assert_allclose(upper[2][valid] - lower[2][valid], 6)
