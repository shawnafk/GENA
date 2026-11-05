"""观测几何：数组几何 → 事件几何 → 目标定义 → 通用视场投影 → 月球遮挡。

本模块不加载工程数据或太阳风、不插值、不构造环境模型。事件表自带逐行
位置和姿态；位置偏移和仪器快照输入均由调用方准备，不读取星历。

空间位置统一使用 R_E；工程表 MJ2000_Star_X/Y/Z 的输入仍为月心米。
姿态遵循 QQ → ECI2000 的 xyzw 四元数约定。粒子飞行方向为层 1 → 层 2，
仪器观测方向为其反向（探头 -Z）。事件坐标结果保留 coord 对象，
求交结果返回 (N, 3) 数组。
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from common.coord import coord as Coordinate, fov_to_vector, transform_vectors, vector_to_fov
from common.constants import RE, RM
from common.frame import split_time_intervals
from common.time import Time

__all__ = [
    "plane_intersections", "ray_sphere_intersections", "sphere_limb", "axial_mean_angle",
    "position_as_coord", "events_rays", "event_axes",
    "project_events", "projected_fov_axis_angle", "dipole_fieldlines", "obj_in_view",
    "moon_blocking_mask", "moon_blocking_intervals",
]

_ROTATION_COLUMNS = ["MJ2000_q1", "MJ2000_q2", "MJ2000_q3"]
_MOON_RADIUS_RE = RM / RE


# 1. 基础几何：只处理数组，不依赖 GENA 表格。

def plane_intersections(origins, directions, *, plane_normal=(0, 1, 0),
                        plane_point=(0, 0, 0), forward=False):
    """Return line-plane intersections as (N, 3), with NaN for parallel rows.

    Inputs share a coordinate frame and length unit. Single rows broadcast.
    The default extends lines in both directions; forward=True keeps t >= 0.
    """
    origins, directions = np.broadcast_arrays(
        np.atleast_2d(np.asarray(origins, dtype=float)),
        np.atleast_2d(np.asarray(directions, dtype=float)),
    )
    denominator = np.sum(directions * np.asarray(plane_normal), axis=1)
    numerator = np.sum((np.asarray(plane_point) - origins) * plane_normal, axis=1)
    distance = np.divide(numerator, denominator, out=np.full(len(origins), np.nan),
                         where=np.abs(denominator) > 1e-9)
    if forward:
        distance[distance < 0] = np.nan
    return origins + distance[:, None] * directions


def ray_sphere_intersections(origins, directions, centers=(0, 0, 0), *, radius=_MOON_RADIUS_RE):
    """Return the nearest forward sphere hit for each ray, shape (N, 3).

    Single rows broadcast. Misses, backward hits and zero directions yield NaN.
    All positions and radius must use the same coordinate frame and length unit.
    """
    origins, directions, centers = np.broadcast_arrays(*[
        np.atleast_2d(np.asarray(values, dtype=float))
        for values in (origins, directions, centers)
    ])
    delta = origins - centers
    a = np.sum(directions**2, axis=1)
    b = 2 * np.sum(directions * delta, axis=1)
    c = np.sum(delta**2, axis=1) - radius**2
    discriminant = b**2 - 4 * a * c
    valid = (a > 0) & (discriminant >= 0)
    root = np.sqrt(np.clip(discriminant, 0, None))
    distances = np.stack([
        np.divide(-b - root, 2 * a, out=np.full(len(a), np.inf), where=valid),
        np.divide(-b + root, 2 * a, out=np.full(len(a), np.inf), where=valid),
    ])
    distances[distances < 0] = np.inf
    nearest = np.min(distances, axis=0)
    nearest[~np.isfinite(nearest)] = np.nan
    return origins + nearest[:, None] * directions


def sphere_limb(center, radius, observer=(0, 0, 0), *, samples=181):
    """Return the exact tangent circle of a sphere as (samples, 3) positions.

    Center, observer and radius use the same frame and length unit. The
    observer must be outside the sphere. The returned curve is closed.
    """
    center, observer = np.asarray(center, dtype=float), np.asarray(observer, dtype=float)
    delta = observer - center
    distance = np.linalg.norm(delta)
    if radius <= 0 or distance <= radius:
        raise ValueError("radius must be positive and observer must be outside the sphere")
    axis = delta / distance
    u = np.cross(axis, np.eye(3)[np.argmin(np.abs(axis))])
    u /= np.linalg.norm(u)
    v = np.cross(axis, u)
    phase = np.linspace(0, 2 * np.pi, samples)
    circle = np.cos(phase)[:, None] * u + np.sin(phase)[:, None] * v
    return center + axis * (radius**2 / distance) + radius * np.sqrt(1 - (radius / distance)**2) * circle


def axial_mean_angle(angle_deg, weights=None):
    """Mean unoriented line angles in degrees, returned in [-90, 90).

    Uses doubled angles because a line is unchanged by a 180-degree reversal.
    Missing angles and zero weights are ignored; an undefined mean is NaN.
    """
    angles = np.asarray(angle_deg, dtype=float)
    weights = np.ones_like(angles) if weights is None else np.broadcast_to(weights, angles.shape)
    if np.any(weights < 0):
        raise ValueError("weights must be non-negative")
    valid = np.isfinite(angles) & np.isfinite(weights) & (weights > 0)
    doubled = np.deg2rad(2 * angles[valid])
    sine = np.sum(weights[valid] * np.sin(doubled))
    cosine = np.sum(weights[valid] * np.cos(doubled))
    return np.nan if np.isclose(sine, 0) and np.isclose(cosine, 0) else float(
        np.rad2deg(np.arctan2(sine, cosine)) / 2
    )


# 2. 事件几何：一行对应一个粒子，始终保留行顺序和无效行。
#general vector to coordinate class.
def _vectors_as_coord(dates, vectors, quaternions, *, ref_coord, res_coord='GSM'):
    """Rotate event vectors and preserve invalid-attitude rows as NaN."""
    valid = np.isfinite(quaternions).all(axis=1)
    rotations = np.where(valid[:, None], quaternions, (0, 0, 0, 1))
    result = Coordinate(dates, np.where(valid[:, None], vectors, np.nan), ref_coord, rotations)
    result.convert(res_coord)
    return result


def position_as_coord(dates, locations, *, pos_offset=None,
                      ref_coord="ECI2000", res_coord="GSM"):
    """Add an optional position offset, then convert to a coord object.

    locations and pos_offset share ref_coord and units of R_E. pos_offset
    can be one (3,) vector or an aligned (N, 3) array. No data are loaded.
    """
    position = np.asarray(locations, dtype=float)
    if pos_offset is not None:
        position = position + pos_offset
    result = Coordinate(dates, position, ref_coord, unit="Re")
    result.convert(res_coord)
    return result


def events_rays(L2, *, probe="A", coord="GSM"):
    """Return unnormalised layer-1 → layer-2 flight vectors as a coord object.

    Native components are (X2-X1, Y2-Y1, 50) in millimetres. Attitude comes
    directly from each event's MJ2000_q1/q2/q3 columns, without interpolation.
    """
    vectors = np.column_stack((L2["X2"] - L2["X1"],
                               L2["Y2"] - L2["Y1"], np.full(len(L2), 50.)))
    return _vectors_as_coord(L2["Datetime"], vectors, L2.rot(),
                                    ref_coord=probe, res_coord=coord)


def event_axes(L2, *, probe="A", coord="GSM"):
    """Return the event-wise probe X/Y/Z unit axes as three coord objects."""
    quaternions = L2.rot()
    return tuple(
        _vectors_as_coord(L2["Datetime"], np.broadcast_to(axis, (len(L2), 3)),
                                 quaternions, ref_coord=probe, res_coord=coord)
        for axis in np.eye(3)
    )


def project_events(events, *, pos_offset, probe="A", coord="GSM",
                   plane_normal=(0, 1, 0), plane_point=(0, 0, 0), forward=False):
    """Intersect event flight lines with a plane; return (N, 3) points in R_E.

    pos_offset is the prepared geocentric Moon position in ECI2000 and R_E.
    The plane uses coord and R_E.
    forward=False preserves bidirectional line projection; True uses flight rays.
    """
    position = position_as_coord(events["Datetime"], events.pos() / RE,
                                 pos_offset=pos_offset, ref_coord="ECI2000", res_coord=coord)[coord]
    direction = events_rays(events, probe=probe, coord=coord)[coord]
    return plane_intersections(position, direction, plane_normal=plane_normal,
                               plane_point=plane_point, forward=forward)


def projected_fov_axis_angle(events, *, pos_offset, probe="A", coord="AGSM",
                             plane_normal=(1, 0, 0), plane_point=(-65, 0, 0),
                             fov_axis="PX", offset_deg=1., plot_axes=(1, 2), weights=None):
    """Return each event's projected PX/PY line angle and its axial mean.

    Two flight rays at ±offset_deg on the chosen source-angle axis intersect
    the plane. Angles are measured from plot_axes[0] toward plot_axes[1], in
    degrees within [-90, 90). pos_offset is the prepared geocentric Moon
    position in ECI2000 and R_E; the plane also uses R_E.
    """
    if fov_axis not in {"PX", "PY"} or offset_deg <= 0:
        raise ValueError("choose fov_axis='PX'/'PY' and a positive offset_deg")
    position = position_as_coord(events["Datetime"], events.pos() / RE,
                                 pos_offset=pos_offset, ref_coord="ECI2000", res_coord=coord)[coord]
    quaternions = events.rot()
    intersections = []
    for offset in (-offset_deg, offset_deg):
        px, py = (offset, 0.) if fov_axis == "PX" else (0., offset)
        flight = -fov_to_vector(px, py, direction="look")
        ray = _vectors_as_coord(events["Datetime"], np.broadcast_to(flight, (len(events), 3)),
                                       quaternions, ref_coord=probe, res_coord=coord)
        intersections.append(plane_intersections(position, ray[coord], plane_normal=plane_normal,
                                                 plane_point=plane_point))
    delta = intersections[1] - intersections[0]
    angle = np.rad2deg(np.arctan2(delta[:, plot_axes[1]], delta[:, plot_axes[0]]))
    angle = (angle + 90.) % 180. - 90.
    return angle, axial_mean_angle(angle, weights)


# 3. 目标定义：只产生空间几何，不依赖观察者、姿态或视场。
# 球体用 {"kind": "sphere", "center": ..., "radius": ..., "coord": ...}。
# 点/曲线用 {"kind": "points"/"curve", "points": ..., "coord": ...}。

def dipole_fieldlines(l_shells=(4, 8), mlt_hours=(0, 6, 12, 18), *, samples=181):
    """Define named reference dipole curves in SM/R_E, independently of a view."""
    objects = {}
    for shell in l_shells:
        latitude = np.linspace(-np.arccos(1 / np.sqrt(shell)), np.arccos(1 / np.sqrt(shell)), samples)
        radius = shell * np.cos(latitude)**2
        for hour in mlt_hours:
            phi = (hour - 12) * np.pi / 12
            points = np.column_stack((radius * np.cos(latitude) * np.cos(phi),
                                      radius * np.cos(latitude) * np.sin(phi), radius * np.sin(latitude)))
            objects[f"dipole_L{shell}_MLT{hour}"] = {
                "kind": "curve", "points": points, "coord": "SM",
                "l_shell": shell, "mlt_hours": hour,
            }
    return objects


# 4. 通用视场投影：不识别天体名称，不生成物理模型，不加载或插值数据。

def _project_to_fov(points, normals):
    """Project camera-frame points; retain NaN gaps and reject the rear hemisphere."""
    points = np.atleast_2d(points)
    front = np.isfinite(points).all(axis=1) & (points[:, 2] < 0)
    px_py = np.full((len(points), 2), np.nan)
    if front.any():
        px_py[front] = np.column_stack(vector_to_fov(points[front], direction="look"))
    forward = points * (1, 1, -1)
    inside = front & (forward @ normals.T >= -1e-12).all(axis=1)
    return px_py, inside


def _clip_view_segments(points, normals):
    """Clip straight curve segments to the four planes of the viewing cone."""
    starts, ends = points[:-1], points[1:]
    forward = starts * (1, 1, -1)
    delta = (ends - starts) * (1, 1, -1)
    value, slope = forward @ normals.T, delta @ normals.T
    crossing = np.divide(-value, slope, out=np.zeros_like(value), where=slope != 0)
    enter = np.maximum(0., np.max(np.where(slope > 0, crossing, -np.inf), axis=1))
    leave = np.minimum(1., np.min(np.where(slope < 0, crossing, np.inf), axis=1))
    valid = np.isfinite(starts).all(axis=1) & np.isfinite(ends).all(axis=1)
    valid &= ~((slope == 0) & (value < 0)).any(axis=1) & (enter <= leave)
    offsets = np.stack((enter[valid], leave[valid]), axis=1)
    clipped = starts[valid, None, :] + offsets[:, :, None] * (ends - starts)[valid, None, :]
    angles, _ = _project_to_fov(clipped.reshape(-1, 3), normals)
    segments = angles.reshape(-1, 2, 2)
    return segments[np.isfinite(segments).all(axis=(1, 2))]


def _sphere_view_flags(center, radius, normals, corners):
    """Test spherical-disk/cone overlap analytically, independently of limb sampling."""
    distance = np.linalg.norm(center)
    direction = center * (1, 1, -1) / distance
    margin = normals @ direction
    sine = radius / distance
    cosine = np.sqrt(1 - sine**2)
    fully_in_view = bool((margin >= sine - 1e-12).all())
    covers_view = bool((corners @ direction >= cosine - 1e-12).all())
    # The closest direction in the cone is the center direction, an edge
    # projection, or a corner. This also detects disks enclosing the FOV.
    edges = direction - margin[:, None] * normals
    norms = np.linalg.norm(edges, axis=1)
    edges = np.divide(edges, norms[:, None], out=np.zeros_like(edges), where=norms[:, None] > 0)
    valid_edges = (edges @ normals.T >= -1e-12).all(axis=1) & (norms > 0)
    best = np.max(corners @ direction)
    if valid_edges.any():
        best = max(best, np.max(edges[valid_edges] @ direction))
    if (margin >= 0).all():
        best = 1.
    return bool(best >= cosine - 1e-12), fully_in_view, covers_view


def obj_in_view(time, position, attitude, objects, *, probe="A", ref_coord="ECI2000",
                ax_deg=22.5, ay_deg=10., samples=181):
    """Project named sphere, curve and point objects into one observer's FOV.

    position is the observer's (3,) location in ref_coord, in R_E. attitude
    is QQ→ECI2000 xyzw or Rotation. The view looks along probe -Z; ax_deg
    and ay_deg are rectangular half angles. All object positions share the
    observer's origin and R_E unit; each object may specify its own coord.

    objects maps names to definitions: sphere uses kind='sphere', center and
    radius; curve/points use kind='curve'/'points' and an (N, 3) points array.
    A (3,) point is accepted. coord defaults to ref_coord. Extra metadata is
    retained, including optional plotting style. No objects are mutated.

    Each result contains px_py (front-hemisphere angles), in_view (per-sample
    mask), visible_px_py (outside samples NaN), any_in_view and all_in_view.
    Curves/sphere limbs also return clipped segments_px_py, (M, 2, 2), so a
    segment crossing the view is retained even if both endpoints are outside.
    Spheres also return center_px_py and covers_view. Their overlap, whole-disk
    containment and whole-FOV coverage flags are analytic, not sampled.
    This function tests each object independently; it does not hide one object
    behind another. time is one UTC instant; no ephemeris or models are loaded.
    """
    if not (0 < ax_deg < 90 and 0 < ay_deg < 90):
        raise ValueError("FOV half angles must be between 0 and 90 degrees")
    stamp = Time(time).datetime
    position = np.asarray(position, dtype=float)
    ax, ay = np.deg2rad([ax_deg, ay_deg])
    normals = np.array([[np.cos(ax), 0, np.sin(ax)], [-np.cos(ax), 0, np.sin(ax)],
                        [0, np.cos(ay), np.sin(ay)], [0, -np.cos(ay), np.sin(ay)]])
    corners = fov_to_vector(np.array([-ax_deg, -ax_deg, ax_deg, ax_deg]),
                            np.array([-ay_deg, ay_deg, -ay_deg, ay_deg]))
    results = {}
    for name, obj in objects.items():
        kind = obj["kind"]
        if kind not in {"sphere", "curve", "points"}:
            raise ValueError("object kind must be 'sphere', 'curve' or 'points'")
        points = np.atleast_2d(np.asarray(obj["center"] if kind == "sphere" else obj["points"], dtype=float))
        points = transform_vectors(points, source=obj.get("coord", ref_coord), target=ref_coord,
                                   times=stamp, attitude=attitude, unit="Re")
        camera = transform_vectors(points - position, source=ref_coord, target=probe,
                                   times=stamp, attitude=attitude, unit="Re")
        result = dict(obj)
        if kind == "sphere":
            center = camera[0]
            camera = sphere_limb(center, obj["radius"], samples=samples)
            result["center_px_py"] = _project_to_fov(center, normals)[0][0]
        px_py, inside = _project_to_fov(camera, normals)
        result.update(px_py=px_py, in_view=inside, visible_px_py=np.where(inside[:, None], px_py, np.nan),
                      any_in_view=bool(inside.any()), all_in_view=bool(len(inside) and inside.all()))
        if kind in {"sphere", "curve"}:
            result["segments_px_py"] = _clip_view_segments(camera, normals)
            result["any_in_view"] |= bool(len(result["segments_px_py"]))
        if kind == "sphere":
            overlap, contained, covers = _sphere_view_flags(center, obj["radius"], normals, corners)
            result.update(any_in_view=overlap, all_in_view=contained, covers_view=covers)
        results[name] = result
    return results


# 5. 月球遮挡：只准备月球目标和观测姿态，复用 obj_in_view。

def moon_blocking_mask(pose, *, probe="A", gap=600, ax_deg=22.5, ay_deg=10., samples=181):
    """Return Datetime/partial_block/full_block from preloaded pose rows.

    pose is a GENA frame supplying Datetime, pos() in metres, and rot().
    gap is a stride through valid pose rows, not a time step.
    A valid PosePresent flag takes precedence; legacy rows without the flag
    require nonzero stored quaternion components. No interpolation is performed.

    partial_block means the lunar disk intersects the FOV; full_block means
    the whole disk is inside it. covers_view means the Moon covers the entire
    FOV. These disk/FOV tests are analytic; samples only controls limb drawing.
    """
    position = pose.pos()
    quaternions = pose.rot()
    present = pose[_ROTATION_COLUMNS].ne(0).any(axis=1)
    if "PosePresent" in pose:
        present = pose["PosePresent"].fillna(present).astype(bool)
    valid = np.isfinite(position).all(axis=1) & np.isfinite(quaternions).all(axis=1)
    valid &= (np.linalg.norm(position, axis=1) > 0) & np.asarray(present)
    selected = np.flatnonzero(valid)[::gap]
    dates = pose["Datetime"].iloc[selected].to_numpy()

    objects = {"moon": {"kind": "sphere", "center": (0, 0, 0), "radius": _MOON_RADIUS_RE}}
    partial, full, covers = [], [], []
    for index in selected:
        moon = obj_in_view(pose["Datetime"].iloc[index], position[index] / RE, quaternions[index],
                           objects, probe=probe, ax_deg=ax_deg, ay_deg=ay_deg, samples=samples)["moon"]
        partial.append(moon["any_in_view"])
        full.append(moon["all_in_view"])
        covers.append(moon["covers_view"])
    return pd.DataFrame({"Datetime": dates, "partial_block": np.asarray(partial, dtype=bool),
                         "full_block": np.asarray(full, dtype=bool), "covers_view": np.asarray(covers, dtype=bool)})


def moon_blocking_intervals(pose, *, probe="A", gap=600, max_gap="60min", full_disk=False, **kwargs):
    """Return contiguous blocking intervals; full_disk selects the full_block flag."""
    mask = moon_blocking_mask(pose, probe=probe, gap=gap, **kwargs)
    selected = mask["full_block" if full_disk else "partial_block"]
    return split_time_intervals(mask.loc[selected, "Datetime"], max_gap=max_gap)
