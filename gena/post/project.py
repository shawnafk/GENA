import common.coord as gcd
from common.coord import fov_to_vector, vector_to_fov
from common.time import Time
from load import load, load_pos, load_rot
from gena.frame import ParameterDataFrame, _row_quaternions
from gena import utils as gena_utils
from common.constants import RE, RM
from common.frame import split_time_intervals
from libs import ephemeris as epms

import numpy as np
import pandas  as pd

RM_RE = RM / RE
FOV_AXES_MINUS_Z = np.array([[1, 0, 0], [0, 1, 0], [0, 0, -1]])

__all__ = [
    "event_vectors",
    "axial_mean_angle",
    "projected_fov_axis_angle",
    "get_spacecraft_pos",
    "get_ray",
    "get_frame",
    "get_projection",
    "get_instrument_projection",
    "ray_sphere_intersections",
    "moon_blocking_mask",
    "moon_blocking_intervals",
]


def _event_attitude(events):
    quaternions = _row_quaternions(events)
    valid = np.isfinite(quaternions[:, 3])
    # The coordinate API requires finite rotations. Mask vectors on these rows
    # so the placeholder rotations cannot produce a physical direction.
    quaternions[~valid] = (0.0, 0.0, 0.0, 1.0)
    return quaternions, valid


def _event_coord(dates, vectors, probe, coord, attitude):
    quaternions, valid = attitude
    ray = gcd.coord(dates, np.where(valid[:, None], vectors, np.nan), probe, quaternions)
    ray.convert(coord)
    return ray


def axial_mean_angle(angle_deg, weights=None) -> float:
    """Return the circular mean of unoriented line angles in degrees.

    A projected FOV axis is an unoriented line, so angles separated by 180
    degrees are equivalent.  The mean is therefore evaluated after doubling
    the angles and is returned in the interval ``[-90, 90)``.

    Parameters
    ----------
    angle_deg : array-like
        Line angles in degrees.
    weights : array-like, optional
        Non-negative weights.  Passing one angle per event with no explicit
        weights gives an event-count-weighted mean automatically.
    """
    angle = np.asarray(angle_deg, dtype=float)
    if weights is None:
        weight = np.ones(angle.shape, dtype=float)
    else:
        weight = np.broadcast_to(np.asarray(weights, dtype=float), angle.shape)
        if np.any(np.isfinite(weight) & (weight < 0)):
            raise ValueError("weights must be non-negative.")

    valid = np.isfinite(angle) & np.isfinite(weight) & (weight > 0)
    if not np.any(valid):
        return np.nan

    doubled = np.deg2rad(2.0 * angle[valid])
    sine = np.sum(weight[valid] * np.sin(doubled))
    cosine = np.sum(weight[valid] * np.cos(doubled))
    if np.isclose(sine, 0.0) and np.isclose(cosine, 0.0):
        return np.nan
    return float(0.5 * np.rad2deg(np.arctan2(sine, cosine)))


def projected_fov_axis_angle(
    l2df,
    *,
    prob="A",
    coord="AGSM",
    norm=(1.0, 0.0, 0.0),
    point=(-65.0, 0.0, 0.0),
    fov_axis="PX",
    offset_deg=1.0,
    plot_axes=(1, 2),
    weights=None,
):
    """Return per-row and mean projected FOV-axis angles.

    Two rays at ``-offset_deg`` and ``+offset_deg`` on the requested angular
    axis are transformed with each row's attitude and intersected with the
    requested plane.  The line angle is measured from ``plot_axes[0]`` toward
    ``plot_axes[1]``; for the usual AGSM Y-Z map use the default ``(1, 2)``.

    ``PX`` and ``PY`` in L2 describe the source angular coordinates, whereas
    :func:`get_ray` follows particle flight.  Consequently a source direction
    ``(px, py)`` is represented here by the flight vector
    ``(-tan(px), -tan(py), +1)``.  Negating the complete ray would give the
    same line-plane intersection.

    Parameters
    ----------
    l2df : ParticleDataFrame
        L2 rows containing timestamps, row-wise attitude, and position.
    prob : {'A', 'B'}, optional
        Probe coordinate frame.
    coord : str, optional
        Output coordinate frame.
    norm, point : array-like, optional
        Projection-plane normal and one point on the plane.
    fov_axis : {'PX', 'PY'}, optional
        Instrument angular axis whose projected orientation is measured.
    offset_deg : float, optional
        Symmetric angular offset used to construct the two rays.
    plot_axes : tuple[int, int], optional
        Coordinate components used as horizontal and vertical plot axes.
    weights : array-like, optional
        Weights for the axial circular mean.  With one input row per event,
        omitting weights gives an event-count-weighted mean.

    Returns
    -------
    angle_deg : np.ndarray
        Projected line angle for every input row, normalized to ``[-90, 90)``.
    mean_angle_deg : float
        Axial circular mean of ``angle_deg``.
    """
    if len(l2df) == 0:
        return np.array([], dtype=float), np.nan
    if fov_axis.upper() not in {"PX", "PY"}:
        raise ValueError("fov_axis must be 'PX' or 'PY'.")
    if not np.isfinite(offset_deg) or offset_deg <= 0:
        raise ValueError("offset_deg must be a positive finite scalar.")
    if len(plot_axes) != 2 or any(axis not in (0, 1, 2) for axis in plot_axes):
        raise ValueError("plot_axes must contain two coordinate indices from (0, 1, 2).")

    n_rows = len(l2df)
    px = np.zeros((2, n_rows), dtype=float)
    py = np.zeros((2, n_rows), dtype=float)
    if fov_axis.upper() == "PX":
        px[0] = -offset_deg
        px[1] = offset_deg
    else:
        py[0] = -offset_deg
        py[1] = offset_deg

    dates = pd.to_datetime(l2df["Timestamp"].to_numpy(dtype=float), unit="s").to_list()
    attitude = _event_attitude(l2df)
    start = np.asarray(get_spacecraft_pos(l2df, coord)[coord], dtype=float)
    intersections = []
    for px_side, py_side in zip(px, py):
        flight = np.column_stack((
            -np.tan(np.deg2rad(px_side)),
            -np.tan(np.deg2rad(py_side)),
            np.ones(n_rows, dtype=float),
        ))
        ray = _event_coord(dates, flight, prob, coord, attitude)
        intersections.append(project(start, ray[coord], norm, point))

    delta = intersections[1] - intersections[0]
    angle = np.rad2deg(np.arctan2(delta[:, plot_axes[1]], delta[:, plot_axes[0]]))
    angle = (angle + 90.0) % 180.0 - 90.0
    return angle, axial_mean_angle(angle, weights=weights)


def event_vectors(events, *, coord="GSM", probe="A"):
    """Return a copy with unit particle-flight direction and velocity in ``coord``.

    The added columns are always ``DX/DY/DZ`` and ``VX/VY/VZ``.  The input is
    never mutated, and no implicit velocity scaling is supported.
    """
    required = ["X1", "Y1", "X2", "Y2", "Tof"]
    missing = [name for name in required if name not in events]
    if missing:
        raise KeyError(f"Missing required event columns: {missing}")
    target = events.copy()
    if len(target) == 0:
        for name in ("DX", "DY", "DZ", "VX", "VY", "VZ"):
            target[name] = np.array([], dtype=float)
        return target
    dates = target["Datetime"].to_list()
    attitude = _event_attitude(target)
    ray_probe = np.column_stack((
        target["X2"].to_numpy(dtype=float) - target["X1"].to_numpy(dtype=float),
        target["Y2"].to_numpy(dtype=float) - target["Y1"].to_numpy(dtype=float),
        np.full(len(target), 50.0),
    ))
    ray = _event_coord(dates, ray_probe, probe, coord, attitude)
    direction = np.asarray(ray[coord], dtype=float)
    norms = np.linalg.norm(direction, axis=1)
    direction = np.divide(
        direction,
        norms[:, np.newaxis],
        out=np.full_like(direction, np.nan),
        where=norms[:, np.newaxis] > 0,
    )
    velocity = _event_coord(dates, gena_utils.get_velocity(target), probe, coord, attitude)
    target.loc[:, ["DX", "DY", "DZ"]] = direction
    target.loc[:, ["VX", "VY", "VZ"]] = np.asarray(velocity[coord], dtype=float)
    target.attrs = dict(events.attrs)
    target.attrs.update(vector_coordinate=str(coord).upper(), probe=str(probe).upper())
    return target

#* for verification
#* project A B and QQ center to GSM

'''
#! deprecated
def project(start, dir, plane_normal=np.array([0,1,0]), plane_point=np.array([0,0,0])):
    """
    计算由起点和方向矢量定义的直线与Y=0平面的交点
    
    参数:
        direction_vec: N * 3 numpy数组，表示直线的方向矢量
        start_point: N * 3 numpy数组，表示直线的起点（默认原点）
        
    返回:
        intersection: 3x1 numpy数组，表示交点坐标（x, 0, z）
        若直线与Y=0平面平行或在平面上，则返回None
    """
    # 提取方向矢量的Y分量和起点的Y分量
    dir_y = dir[:,1]
    start_y = start[:, 1]
    
    t = -start_y / dir_y
    
    # 计算交点（此时Y分量必为0）
    intersection = start + t[:,np.newaxis] * dir
    #return intersection[0,0],intersection[2,0]
    return intersection[:,0],intersection[:,2]
'''

def project(start, dir, plane_normal=np.array([0,1,0]), plane_point=np.array([0,0,0])):
    """
    计算由起点和方向矢量定义的射线（或直线）与任意平面的交点
    
    参数:
        start: N * 3 numpy数组，表示N条射线的起点坐标（x,y,z）
        dir: N * 3 numpy数组，表示N条射线的方向矢量（x,y,z）
        plane_normal: 1 * 3 numpy数组，表示平面的法向量（默认Y=0平面法向量[0,1,0]）
        plane_point: 1 * 3 numpy数组，表示平面上的任意一点（默认Y=0平面上的原点[0,0,0]）
        
    返回:
        x: N维numpy数组，交点的x坐标（仅当存在有效交点时）
        z: N维numpy数组，交点的z坐标（仅当存在有效交点时）
        若射线与平面平行/在平面上/交点在射线反向延长线，对应位置返回np.nan
    """
    # 确保输入数组维度正确
    start = np.atleast_2d(start)
    dir = np.atleast_2d(dir)
    plane_normal = np.atleast_1d(plane_normal)
    plane_point = np.atleast_1d(plane_point)
    
    # 1. 计算法向量与方向向量的点积（n·d）
    n_dot_d = np.dot(dir, plane_normal.T).reshape(-1)  # 形状：(N,)
    
    # 2. 计算分子部分：n·(r0 - s)
    r0_minus_s = plane_point - start  # 形状：(N, 3)
    n_dot_r0s = np.dot(r0_minus_s, plane_normal.T).reshape(-1)  # 形状：(N,)
    
    # 3. 初始化t数组，默认设为nan（无效值）
    t = np.full_like(n_dot_d, np.nan, dtype=np.float64)
    
    # 4. 筛选出非平行的射线（n·d ≠ 0），计算t
    valid_mask = np.abs(n_dot_d) > 1e-9  # 避免浮点精度问题，不直接用==0
    t[valid_mask] = n_dot_r0s[valid_mask] / n_dot_d[valid_mask]
    
    # 5. （可选）筛选出射线有效交点（t≥0），若需直线交点可注释此行
    #ray_valid_mask = valid_mask & (t >= 0)
    #t[~ray_valid_mask] = np.nan
    
    # 6. 计算交点坐标
    intersection = start + t[:, np.newaxis] * dir  # 形状：(N, 3)
    
    return intersection




#verify the focus is aiming the right place
def verify_point(ts_list,coord='GSM'):
    #& datetime
    timestamps = np.asarray(ts_list, dtype=float)
    if timestamps.ndim != 1 or len(timestamps) == 0 or not np.isfinite(timestamps).all():
        raise ValueError("ts_list must contain finite one-dimensional Unix timestamps")
    pd_date = pd.to_datetime(timestamps,unit='s').to_list()
    start = pd.to_datetime(timestamps.min(), unit="s", utc=True) - pd.Timedelta("1h")
    stop = pd.to_datetime(timestamps.max(), unit="s", utc=True) + pd.Timedelta("1h")
    pdf = load("gena", 
        trange=[start, stop],
        product="pdf",
        columns=[
            "Timestamp",
            "MJ2000_Star_X",
            "MJ2000_Star_Y",
            "MJ2000_Star_Z",
            "MJ2000_q1",
            "MJ2000_q2",
            "MJ2000_q3",
        ],
    )
    qq_loc = pdf.pos_at(timestamps)/RE
    qq_rot = pdf.rot_at(timestamps)
    moon_loc = epms.moon(pd_date)/RE*1000
    
    start  = gcd.coord(pd_date,qq_loc+moon_loc,'ECI2000',0)
    ray_qq = gcd.coord(pd_date,fov_to_vector(np.zeros(len(qq_rot)),np.zeros(len(qq_rot))),'QQ',qq_rot)
    ray_a = gcd.coord(pd_date,fov_to_vector(np.zeros(len(qq_rot)),np.zeros(len(qq_rot))),'A',qq_rot)
    ray_b = gcd.coord(pd_date,fov_to_vector(np.zeros(len(qq_rot)),np.zeros(len(qq_rot))),'B',qq_rot)
    
    start.convert(coord)
    ray_qq.convert(coord)
    ray_a.convert(coord)
    ray_b.convert(coord)
    
    c_qq = project(start[coord],ray_qq[coord])[:]
    c_a = project(start[coord],ray_a[coord])[:]
    c_b = project(start[coord],ray_b[coord])[:]
    return c_qq,c_a,c_b

#may be use read_ to infer directly read, and derive_ to infer calculated
def get_spacecraft_pos(l2df, coord='GSM'):
    """Convert spacecraft positions stored on individual L2 event rows."""
    ts_list = l2df['Timestamp'].to_list()
    pd_date = pd.to_datetime(ts_list,unit='s').to_list()
    position_columns = ['MJ2000_Star_X', 'MJ2000_Star_Y', 'MJ2000_Star_Z']
    missing = [name for name in position_columns if name not in l2df.columns]
    if missing:
        raise KeyError(f"position columns are missing: {missing}")
    qq_loc = l2df[position_columns].to_numpy(dtype=float)/RE
    moon_loc = epms.moon(pd_date)/RE*1000
    start  = gcd.coord(pd_date,qq_loc+moon_loc,'ECI2000',0)
    start.convert(coord)
    return start

def get_ray(l2df,prob='A',coord='GSM'):
    """Convert event flight vectors, retaining NaN rows for invalid attitude."""
    ts_list = l2df['Timestamp'].to_list()
    pd_date = pd.to_datetime(ts_list,unit='s').to_list()
    attitude = _event_attitude(l2df)
    # Particle flight is layer 1 -> layer 2.  The source-looking vector defined
    # in the GENA geometry note is its exact negative,
    # (X1-X2, Y1-Y2, -50).  A line-plane projection is invariant when the full
    # direction vector is negated, so this convention cannot mirror Y or Z.
    vec_in_qq = np.column_stack((l2df['X2']-l2df['X1'],l2df['Y2']-l2df['Y1'],np.ones(len(l2df))*50))
    return _event_coord(pd_date, vec_in_qq, prob, coord, attitude)

def get_frame(l2df,prob='A',coord='GSM'):
    """Convert each event's probe axes, retaining NaN rows for invalid attitude."""
    ts_list = l2df['Timestamp'].to_list()
    pd_date = pd.to_datetime(ts_list,unit='s').to_list()
    attitude = _event_attitude(l2df)
    #end - start
    qqx = np.column_stack((np.ones(len(l2df)),np.ones(len(l2df))*0,np.ones(len(l2df))*0))
    qqy = np.column_stack((np.ones(len(l2df))*0,np.ones(len(l2df)),np.ones(len(l2df))*0))
    qqz = np.column_stack((np.ones(len(l2df))*0,np.ones(len(l2df))*0,np.ones(len(l2df))))
    rayx = _event_coord(pd_date, qqx, prob, coord, attitude)
    rayy = _event_coord(pd_date, qqy, prob, coord, attitude)
    rayz = _event_coord(pd_date, qqz, prob, coord, attitude)
    return rayx,rayy,rayz

def get_projection(l2df,prob='A',coord='GSM',norm=[0,1,0],point=[0,0,0]):
    """Return one plane intersection per event, with NaN for invalid geometry."""
    ray = get_ray(l2df,prob,coord)
    start = get_spacecraft_pos(l2df, coord)
    proj = project(start[coord],ray[coord],norm,point)[:]
    return proj


def _earth_limb(observer, distance, samples):
    """Sample the tangent circle on a spherical Earth of radius 1 R_E."""
    axis = observer / distance
    reference = np.eye(3)[np.argmin(np.abs(axis))]
    u = np.cross(axis, reference)
    u /= np.linalg.norm(u)
    v = np.cross(axis, u)
    phase = np.linspace(0, 2 * np.pi, samples)
    circle = np.cos(phase)[:, None] * u + np.sin(phase)[:, None] * v
    return axis / distance + np.sqrt(1 - distance**-2) * circle


def get_instrument_projection(time, *, prob="A", l_shells=(4, 8),
                              mlt_hours=(0, 6, 12, 18), samples=181,
                              magnetotail=False):
    """Project Earth and dipole reference lines into a GENA angular image.

    Evaluate one scalar UTC time (or Unix timestamp), usually the center of
    the ENA accumulation interval. Use interpolated Queqiao position and
    attitude and the probe's source-looking -Z convention. Magnetic local
    times are specified in SM, with noon along +X and dusk along +Y.

    Return `time_utc`, `probe`, `earth_distance_re`, `earth_center` (shape
    (2,)), `earth_limb` (samples by 2), and `field_lines` (a list of records
    with `l_shell`, `mlt_hours` and `px_py`). All pixel coordinates are
    PX/PY angles in degrees, ready to overlay on the event histogram. Points
    behind the probe are NaN. Detector aperture clipping and occultation of
    the reference curves by Earth are left to the caller.

    With magnetotail=True, trace nightside T96 field lines using local
    space-weather data (with the usual load_sw defaults for missing inputs).
    Return them as a list of (N, 2) angular tracks in magnetotail_lines.
    """
    instant = Time(time)
    if not instant.isscalar:
        raise ValueError("time must be one UTC instant")
    probe = str(prob).upper()
    if probe not in {"A", "B"}:
        raise ValueError("prob must be 'A' or 'B'")
    if not isinstance(samples, (int, np.integer)) or samples < 3:
        raise ValueError("samples must be an integer >= 3")
    shells = np.atleast_1d(np.asarray(l_shells, dtype=float))
    local_times = np.atleast_1d(np.asarray(mlt_hours, dtype=float))
    if shells.ndim != 1 or not np.isfinite(shells).all() or np.any(shells <= 1):
        raise ValueError("l_shells must be a one-dimensional sequence of finite values > 1")
    if local_times.ndim != 1 or not np.isfinite(local_times).all() or np.any((local_times < 0) | (local_times >= 24)):
        raise ValueError("mlt_hours must be a one-dimensional sequence in [0, 24)")

    stamp = instant.datetime
    observer = (load_pos("queqiao", stamp) + 1000 * epms.moon(stamp, geocentric=True)) / RE
    attitude = load_rot("queqiao", stamp)
    distance = np.linalg.norm(observer)
    if not np.isfinite(distance) or distance <= 1 or not np.isfinite(attitude).all():
        raise ValueError("projection requires valid attitude and a spacecraft position outside Earth")

    def angular(points):
        sight = gcd.transform_vectors(points - observer, source="ECI2000",
                                      target=probe, attitude=attitude)
        front = np.isfinite(sight).all(axis=1) & (sight[:, 2] < 0)
        pixels = np.full((len(points), 2), np.nan)
        if front.any():
            pixels[front] = np.column_stack(vector_to_fov(sight[front], direction="look"))
        return pixels

    field_lines = []
    curves = []
    for shell in shells:
        footpoint = np.arccos(1 / np.sqrt(shell))
        latitude = np.linspace(-footpoint, footpoint, samples)
        radius = shell * np.cos(latitude)**2
        rho = radius * np.cos(latitude)
        z = radius * np.sin(latitude)
        for hour in local_times:
            phi = (hour - 12) * np.pi / 12
            curves.append(np.column_stack((rho * np.cos(phi), rho * np.sin(phi), z)))
            field_lines.append({"l_shell": float(shell), "mlt_hours": float(hour)})
    if curves:
        points = gcd.transform_vectors(np.concatenate(curves), source="SM", target="ECI2000", times=stamp)
        for line, pixels in zip(field_lines, angular(points).reshape(len(curves), samples, 2)):
            line["px_py"] = pixels

    magnetotail_lines = []
    if magnetotail:
        from libs.space import T96_fieldlines
        from load import load_sw

        weather = load_sw(stamp, download=False)
        parameters = [float(np.asarray(weather[name]).item())
                      for name in ("Pressure", "dst", "BY_GSM", "BZ_GSM")]
        if not np.isfinite(parameters).all() or parameters[0] <= 0:
            raise ValueError("T96 requires positive solar-wind pressure and finite Dst, By and Bz")
        model = T96_fieldlines(stamp.tz_localize("UTC").timestamp(), *parameters, Nlat=11)
        tail = [np.asarray(line).T for line in model.trace_lines_y0()]
        if tail:
            points = gcd.transform_vectors(np.concatenate(tail), source="GSM",
                                          target="ECI2000", times=stamp)
            magnetotail_lines = np.split(angular(points), np.cumsum([len(line) for line in tail])[:-1])

    return {
        "time_utc": stamp.tz_localize("UTC"),
        "probe": probe,
        "earth_distance_re": float(distance),
        "earth_center": angular(np.zeros((1, 3)))[0],
        "earth_limb": angular(_earth_limb(observer, distance, samples)),
        "field_lines": field_lines,
        "magnetotail_lines": magnetotail_lines,
    }


def moon_surface_points(
    center=(0.0, 0.0, 0.0),
    *,
    radius: float = RM_RE,
    n_theta: int = 32,
    n_phi: int = 15,
) -> np.ndarray:
    """Return sampled points on a lunar sphere in Earth-radius units."""
    theta = np.linspace(0, 2 * np.pi, n_theta)
    phi = np.linspace(-np.pi / 2, np.pi / 2, n_phi)
    th, ph = np.meshgrid(theta, phi, indexing="ij")
    center = np.asarray(center, dtype=float)
    return np.column_stack([
        center[0] + radius * np.sin(th.ravel()) * np.cos(ph.ravel()),
        center[1] + radius * np.sin(th.ravel()) * np.sin(ph.ravel()),
        center[2] + radius * np.cos(th.ravel()),
    ])


def is_in_fov(vec, *, ax_deg: float = 22.5, ay_deg: float = 10.0) -> np.ndarray:
    """Return whether vector(s) fall inside the rectangular GENA FOV."""
    vec = np.atleast_2d(np.asarray(vec, dtype=float))
    px, py = vector_to_fov(vec, FOV_AXES_MINUS_Z)
    return (np.abs(px) < ax_deg) & (np.abs(py) < ay_deg)


def ray_sphere_intersections(
    origins,
    directions,
    centers=(0.0, 0.0, 0.0),
    *,
    radius: float = RM_RE,
) -> np.ndarray:
    """Return nearest forward ray-sphere intersection for each ray."""
    origins = np.atleast_2d(np.asarray(origins, dtype=float))
    directions = np.atleast_2d(np.asarray(directions, dtype=float))
    centers = np.atleast_2d(np.asarray(centers, dtype=float))
    n = max(len(origins), len(directions), len(centers))

    def _broadcast_rows(arr, name):
        if len(arr) == n:
            return arr
        if len(arr) == 1:
            return np.repeat(arr, n, axis=0)
        raise ValueError(f"{name} must have either 1 row or {n} rows.")

    origins = _broadcast_rows(origins, "origins")
    directions = _broadcast_rows(directions, "directions")
    centers = _broadcast_rows(centers, "centers")

    a = np.einsum("ij,ij->i", directions, directions)
    if np.any(np.isclose(a, 0)):
        raise ValueError("directions must not contain zero vectors.")
    delta = origins - centers
    b = 2 * np.einsum("ij,ij->i", directions, delta)
    c = np.einsum("ij,ij->i", delta, delta) - radius**2
    discriminant = b**2 - 4 * a * c

    out = np.full((n, 3), np.nan, dtype=float)
    hit = discriminant >= 0
    if not np.any(hit):
        return out

    sqrt_disc = np.sqrt(np.clip(discriminant[hit], 0, None))
    a_hit = a[hit]
    b_hit = b[hit]
    t0 = (-b_hit - sqrt_disc) / (2 * a_hit)
    t1 = (-b_hit + sqrt_disc) / (2 * a_hit)
    ts = np.column_stack([t0, t1])
    ts[ts < 0] = np.nan
    has_forward = np.any(np.isfinite(ts), axis=1)
    nearest = np.full(len(ts), np.nan, dtype=float)
    nearest[has_forward] = np.nanmin(ts[has_forward], axis=1)
    hit_indices = np.flatnonzero(hit)
    forward = np.isfinite(nearest)
    out[hit_indices[forward]] = (
        origins[hit_indices[forward]]
        + nearest[forward, np.newaxis] * directions[hit_indices[forward]]
    )
    return out


def moon_blocking_mask(
    pdf: ParameterDataFrame,
    *,
    gap: int = 600,
    ax_deg: float = 22.5,
    ay_deg: float = 10.0,
    n_theta: int = 32,
    n_phi: int = 15,
) -> pd.DataFrame:
    """Return per-sample lunar partial/full FOV blocking flags."""
    if not isinstance(pdf, ParameterDataFrame):
        raise TypeError("pdf must be a ParameterDataFrame")
    st_time, st_pos = pdf._loc_from_pdf(gap=gap)
    rot_time, st_rot = pdf._rot_from_pdf(gap=gap)
    n = min(len(st_time), len(rot_time), len(st_pos), len(st_rot))
    if n == 0:
        return pd.DataFrame(columns=["Datetime", "partial_block", "full_block"])

    st_time = np.asarray(st_time, dtype=float)[:n]
    st_pos = np.asarray(st_pos, dtype=float)[:n]
    dates = pd.to_datetime(st_time, unit="s")
    sat_pos_re = st_pos / RE

    partial = np.zeros(n, dtype=bool)
    full = np.zeros(n, dtype=bool)
    for i in range(n):
        surface_vec_eci = moon_surface_points(
            -sat_pos_re[i],
            radius=RM_RE,
            n_theta=n_theta,
            n_phi=n_phi,
        )
        surface_vec_qq = st_rot[i].inv().apply(surface_vec_eci)
        inside = is_in_fov(surface_vec_qq, ax_deg=ax_deg, ay_deg=ay_deg)
        partial[i] = bool(np.any(inside))
        full[i] = bool(np.all(inside))
    return pd.DataFrame({"Datetime": dates, "partial_block": partial, "full_block": full})


def moon_blocking_intervals(
    pdf: ParameterDataFrame,
    *,
    gap: int = 600,
    max_gap="60min",
    full_disk: bool = False,
    **mask_kwargs,
) -> list[tuple[pd.Timestamp, pd.Timestamp]]:
    """Return contiguous moon-blocking intervals from a pose PDF frame."""
    mask = moon_blocking_mask(pdf, gap=gap, **mask_kwargs)
    key = "full_block" if full_disk else "partial_block"
    return split_time_intervals(mask.loc[mask[key], "Datetime"], max_gap=max_gap)
