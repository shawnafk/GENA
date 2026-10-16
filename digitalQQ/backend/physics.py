"""Pure scene calculations. Called sequentially inside isolated worker processes."""

from contextlib import contextmanager
from functools import lru_cache
from pathlib import Path
import json
import numpy as np
from astropy.time import Time as AstroTime
from scipy.spatial.transform import Rotation, Slerp
from common.constants import RE
from common.coord import transform_vectors
from gena.constants import va1, va2, va3, va4, vb1, vb2, vb3, vb4
from libs.ephemeris import kernel
from .models import Scene, Orbit

DATA = Path(__file__).resolve().parents[1] / "data"
FOV = np.array([va1, va2, va3, va4, vb1, vb2, vb3, vb4])


@lru_cache(maxsize=1)
def space_models():
    """Disable geopack's unbounded import-time coefficient update in this worker."""
    import importlib.util
    from unittest.mock import patch

    spec = importlib.util.find_spec("geopack")
    if spec is None or not spec.origin:
        raise ValueError("请安装 geopack，并配置本地 IGRF14 系数")
    coefficient = Path(spec.origin).parent / "igrf_coeffs" / "igrf14coeffs.txt"
    if not coefficient.exists():
        raise ValueError(
            f"缺少本地 IGRF14：请将 digitalQQ/data/igrf14coeffs.txt 复制到 {coefficient}"
        )

    def offline(*args, **kwargs):
        raise OSError("digitalQQ uses local IGRF14 coefficients")

    with patch("urllib.request.urlopen", offline):
        from geopack import geopack
        from libs import space
    original = geopack.urlopen
    try:
        geopack.urlopen = offline
        geopack.init_igrf("14")
    finally:
        geopack.urlopen = original
    return space


@lru_cache(maxsize=1)
def coefficients():
    import spacepy.igrf as igrf

    return igrf.IGRFCoefficients(str(DATA / "igrf14coeffs.txt"))


@contextmanager
def igrf14():
    import spacepy.igrf as igrf

    old = igrf.igrfcoeffs
    igrf.igrfcoeffs = coefficients()
    try:
        yield
    finally:
        igrf.igrfcoeffs = old


def convert(v, source, target, time):
    with igrf14():
        return transform_vectors(v, source=source, target=target, times=time)


def ephemeris(time):
    jd = AstroTime(time).tdb.jd
    moon = (kernel[3, 301].compute(jd) - kernel[3, 399].compute(jd)) * 1000 / RE
    sun = (
        (
            kernel[0, 10].compute(jd)
            - kernel[0, 3].compute(jd)
            - kernel[3, 399].compute(jd)
        )
        * 1000
        / RE
    )
    return moon, sun


def initial_rotation(look):
    z = -np.asarray(look) / np.linalg.norm(look)
    x = np.array([z[1], -z[0], 0.0])
    if np.linalg.norm(x) < 1e-14:
        x = np.array([1.0, 0.0, 0.0])
    x /= np.linalg.norm(x)
    return Rotation.from_matrix(np.column_stack([x, np.cross(z, x), z]))


def rotation(roll, pitch):
    return Rotation.from_euler("x", pitch, degrees=True) * Rotation.from_euler(
        "z", roll, degrees=True
    )


def angle(a, b):
    return float(np.degrees(np.arctan2(np.linalg.norm(np.cross(a, b)), np.dot(a, b))))


def solve(a, b, far=False):
    a, b = np.array(a) / np.linalg.norm(a), np.array(b) / np.linalg.norm(b)
    rho, phi = np.hypot(a[0], a[1]), np.arctan2(a[1], a[0])
    wrap = lambda x: (x + 180) % 360 - 180
    if rho < 1e-14:
        rolls = [0.0]
    elif far:
        rolls = [wrap(np.degrees(np.pi - phi))]
    else:
        delta = np.arccos(np.clip(b[0] / rho, -1, 1))
        rolls = [wrap(np.degrees(sign * delta - phi)) for sign in (1, -1)]
    candidates = []
    for roll in rolls:
        u = rotation(roll, 0).apply(a)
        pitch = (
            0.0
            if far or np.hypot(*b[1:]) < 1e-14
            else wrap(np.degrees(np.arctan2(b[2], b[1]) - np.arctan2(u[2], u[1])))
        )
        candidates.append(
            (float(roll), float(pitch), angle(rotation(roll, pitch).apply(a), b))
        )
    err = min(v[2] for v in candidates)
    return min(
        (v for v in candidates if v[2] <= err + 1e-9),
        key=lambda v: (abs(v[1]), abs(v[0])),
    )


@lru_cache(maxsize=4)
def imported_interpolators(orbit_json):
    """Convert imported sample frames once, instead of on every playback tick."""
    orbit = Orbit.model_validate_json(orbit_json)
    samples = orbit.samples
    ts = np.array([s.time.timestamp() for s in samples])
    factor = {"km": 1000 / RE, "m": 1 / RE, "RE": 1}[orbit.unit]
    positions = convert(
        np.array([s.position for s in samples]) * factor,
        orbit.frame,
        "ECI2000",
        [s.time for s in samples],
    )
    rotations = (
        None
        if samples[0].quaternion is None
        else Slerp(ts, Rotation.from_quat([s.quaternion for s in samples]))
    )
    return ts, positions, rotations


def pose(scene, time):
    moon, sun = ephemeris(time)
    measured = None
    if scene.orbit.source == "moon":
        position = moon
    else:
        samples = scene.orbit.samples
        ts, positions, rotations = imported_interpolators(scene.orbit.model_dump_json())
        t = time.timestamp()
        if t < ts[0] or t > ts[-1]:
            raise ValueError(
                f"超出导入轨道覆盖范围：{samples[0].time.isoformat()} 至 {samples[-1].time.isoformat()}"
            )
        position = np.array([np.interp(t, ts, positions[:, i]) for i in range(3)])
        if rotations is not None:
            measured = rotations([t])[0]
    if np.linalg.norm(position) < 1:
        raise ValueError("卫星位置不能位于地球内部")
    return position, sun, moon, measured


def target_direction(scene, pos_gsm, time):
    if scene.target_mode == "far":
        return np.array([-1.0, 0, 0]), None
    target = np.array(scene.target, dtype=float)
    if scene.target_mode == "sheet":
        T96NeutralSheet = space_models().T96NeutralSheet
        w = scene.weather
        sheet = T96NeutralSheet(time.timestamp(), w.pressure, w.dst, w.by, w.bz)
        target[2] = sheet.z_at(target[0], target[1])
        if not np.isfinite(target[2]):
            raise ValueError(
                "该目标 XY 位置没有可用的 T96 中性片根，请选择 X=-70…-5、Y=-30…30 R_E 内的位置"
            )
    vec = target - pos_gsm
    if np.linalg.norm(vec) < 1e-10:
        raise ValueError("目标不能与卫星位置重合")
    return vec / np.linalg.norm(vec), target


def attitude(scene: Scene, force_solve=False):
    p, sun, moon, measured = pose(scene, scene.time)
    pg = convert(p, "ECI2000", "GSM", scene.time)
    basis = convert(initial_rotation(-p).apply(np.eye(3)), "ECI2000", "GSM", scene.time)
    target, point = target_direction(scene, pg, scene.time)
    roll, pitch, error = solve(-basis[2], target, scene.target_mode == "far")
    used_roll, used_pitch = (
        (roll, pitch)
        if force_solve or scene.hold == "tracking"
        else (scene.roll, scene.pitch)
    )
    final = rotation(used_roll, used_pitch).apply(basis)
    if scene.hold == "inertial" and not force_solve:
        p0, _, _, _ = pose(scene, scene.reference_time)
        b0 = convert(
            initial_rotation(-p0).apply(np.eye(3)),
            "ECI2000",
            "GSM",
            scene.reference_time,
        )
        final_eci = convert(
            rotation(scene.roll, scene.pitch).apply(b0),
            "GSM",
            "ECI2000",
            scene.reference_time,
        )
        final = convert(final_eci, "ECI2000", "GSM", scene.time)
    look = -final[2]
    ray = point - pg if point is not None else target
    distance = float(np.dot(ray, look))
    miss = (
        float(np.linalg.norm(ray - max(distance, 0) * look))
        if point is not None
        else None
    )
    sg = convert(sun - p, "ECI2000", "GSM", scene.time)
    sun_angle = angle(look, sg)
    q = Rotation.from_matrix(convert(final, "GSM", "ECI2000", scene.time).T).as_quat()
    # Time-local moon track or imported orbit track, with no extrapolation.
    from datetime import timedelta

    times = [scene.time + timedelta(hours=h) for h in np.linspace(-24, 24, 49)]
    if scene.orbit.source == "imported":
        times = [
            t
            for t in times
            if scene.orbit.samples[0].time <= t <= scene.orbit.samples[-1].time
        ]
        times = sorted(set([scene.time] + times))
    track = [convert(pose(scene, t)[0], "ECI2000", "GSM", t).tolist() for t in times]
    return dict(
        time=scene.time.isoformat(),
        position=pg.tolist(),
        moon=convert(moon, "ECI2000", "GSM", scene.time).tolist(),
        axes=final.tolist(),
        reference_axes=basis.tolist(),
        fov=(FOV @ final).tolist(),
        measured_axes=(
            None
            if measured is None
            else convert(
                measured.apply(np.eye(3)), "ECI2000", "GSM", scene.time
            ).tolist()
        ),
        sun_direction=(sg / np.linalg.norm(sg)).tolist(),
        target_direction=target.tolist(),
        target=None if point is None else point.tolist(),
        roll=used_roll,
        pitch=used_pitch,
        solution={
            "roll": roll,
            "pitch": pitch,
            "error": error,
            "reachable": error < 1e-7,
        },
        pointing_error=angle(look, target),
        sun_angle=sun_angle,
        sun_ok=sun_angle > 45,
        ray_miss=miss,
        target_forward=distance > 0,
        quaternion=q.tolist(),
        track=track,
        source="月心近似" if scene.orbit.source == "moon" else "导入轨道",
        measured_error=(
            None
            if measured is None
            else angle(
                look, convert(measured.apply([0, 0, -1]), "ECI2000", "GSM", scene.time)
            )
        ),
    )


def finite_grid(grid):
    return [np.where(np.isfinite(v), v, np.nan).tolist() for v in grid]


@lru_cache(maxsize=24)
def environment_cached(key: str):
    scene = Scene.model_validate_json(key)
    space = space_models()
    Bowshock, Magnetopause, T96_fieldlines, T96NeutralSheet = (
        space.Bowshock,
        space.Magnetopause,
        space.T96_fieldlines,
        space.T96NeutralSheet,
    )
    w = scene.weather
    errors = []
    result = {
        "time": scene.time.isoformat(),
        "weather": w.model_dump(),
        "models": {
            "magnetopause": "Repository empirical magnetopause",
            "bowshock": "Repository empirical bow shock",
            "sheet": "T96 + IGRF Bx=0; thickness is assumed",
        },
        "errors": errors,
    }
    for name, model in [
        (
            "magnetopause",
            Magnetopause(np.nextafter(0.0, 1.0) if w.bz == 0 else w.bz, w.pressure),
        ),
        ("bowshock", Bowshock(w.bz, w.beta, w.mach, w.pressure)),
    ]:
        # Avoid the bow-shock asymptote; clip beyond finite scene domain.
        model.lat = np.linspace(0, np.radians(130 if name == "bowshock" else 165), 48)
        model.long = np.linspace(0, 2 * np.pi, 48)
        with np.errstate(invalid="ignore", over="ignore"):
            grid = np.array(model.surface_xyz())
            radius = model.r_func(np.meshgrid(model.lat, model.long)[0])
        valid = (
            (radius > 0)
            & np.isfinite(grid).all(axis=0)
            & (np.linalg.norm(grid, axis=0) <= 100)
            & (grid[0] >= -80)
        )
        if not valid.any():
            errors.append(f"{name}: 参数组合没有有效曲面")
        grid[:, ~valid] = np.nan
        result[name] = finite_grid(grid)
    result["field_lines"] = []
    result["sheet"] = None
    if scene.magnetic:
        try:
            result["field_lines"] = [
                np.asarray(l).T.tolist()
                for l in T96_fieldlines(
                    scene.time.timestamp(),
                    w.pressure,
                    w.dst,
                    w.by,
                    w.bz,
                    Nlat=9 if scene.quality == "preview" else 21,
                ).trace_lines_y0()
            ]
        except Exception as exc:
            errors.append(f"T96 磁力线：{exc}")
        try:
            n = 11 if scene.quality == "preview" else 25
            sheet = T96NeutralSheet(
                scene.time.timestamp(),
                w.pressure,
                w.dst,
                w.by,
                w.bz,
                nx=n,
                ny=n,
                nz=13,
                half_thickness=scene.half_thickness,
            )
            result["sheet"] = finite_grid(sheet.surface_xyz())
            if not np.isfinite(np.asarray(result["sheet"][2])).any():
                errors.append("T96 中性片：当前网格未找到 Bx=0 根")
        except Exception as exc:
            errors.append(f"T96 中性片：{exc}")
    return result


def environment(scene, magnetic=None):
    if magnetic is not None:
        scene = scene.model_copy(update={"magnetic": magnetic})
    # Attitude, target, orbit and thickness do not affect the central magnetic grid.
    key = json.dumps(
        {
            "time": scene.time.isoformat(),
            "weather": scene.weather.model_dump(),
            "magnetic": scene.magnetic,
            "quality": scene.quality,
        },
        sort_keys=True,
    )
    return environment_cached(key)


def json_safe(value):
    if isinstance(value, dict):
        return {k: json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(v) for v in value]
    if isinstance(value, (float, np.floating)) and not np.isfinite(value):
        return None
    return value
