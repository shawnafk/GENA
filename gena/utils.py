"""Reusable GENA event, environment, and observation-geometry helpers."""

from __future__ import annotations

import matplotlib.dates as mdates
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd

from common.constants import RE
from common.coord import coord as Coord, quaternion_array
from gena.frame import _row_rotations
from plots.tools import axes_for_columns as _axes_for_columns, column_groups as _column_groups

__all__ = [
    "get_velocity",
    "orbit",
    "pointing_geometry",
    "velocity",
]


_ORBIT_LABELS = {
    "spacecraft_x": r"$X_{QQ}$",
    "spacecraft_y": r"$Y_{QQ}$",
    "spacecraft_z": r"$Z_{QQ}$",
    "xy_x": r"$X_{FOV\cap XY}$",
    "xy_y": r"$Y_{FOV\cap XY}$",
    "xz_x": r"$X_{FOV\cap XZ}$",
    "xz_z": r"$Z_{FOV\cap XZ}$",
    "yz_y": r"$Y_{FOV\cap YZ}$",
    "yz_z": r"$Z_{FOV\cap YZ}$",
}






def _ray_plane_intersection(origins, directions, normal, point, *, forward=True):
    origins = np.atleast_2d(np.asarray(origins, dtype=float))
    directions = np.atleast_2d(np.asarray(directions, dtype=float))
    normal = np.asarray(normal, dtype=float)
    point = np.asarray(point, dtype=float)
    denominator = directions @ normal
    numerator = (point - origins) @ normal
    distance = np.full(len(origins), np.nan, dtype=float)
    valid = np.abs(denominator) > 1e-9
    distance[valid] = numerator[valid] / denominator[valid]
    if forward:
        distance[distance < 0] = np.nan
    return origins + distance[:, np.newaxis] * directions


def pointing_geometry(pdf, *, probe="A", coord="GSM", forward=True, attitude=None):
    """Return QQ orbit and nominal source-looking FOV intersections.

    PDF positions are Moon-relative ECI2000 metres. They are converted to
    Earth-centred ``coord`` positions in Earth radii. The nominal source-looking
    direction is probe-frame ``-Z`` and is intersected with the XY, XZ, and YZ
    coordinate planes. Invalid pose rows are omitted. Pass unit xyzw
    quaternions in attitude to use rotations aligned to the input rows.
    """
    probe = str(probe).strip().upper()
    if probe not in {"A", "B"}:
        raise ValueError("probe must be 'A' or 'B'")
    from libs import ephemeris as epms

    columns = [
        "Timestamp", "Datetime",
        "spacecraft_x", "spacecraft_y", "spacecraft_z",
        "look_x", "look_y", "look_z",
        "xy_x", "xy_y", "xy_z",
        "xz_x", "xz_y", "xz_z",
        "yz_x", "yz_y", "yz_z",
    ]
    if len(pdf) == 0:
        return pd.DataFrame(columns=columns)

    position_columns = ["MJ2000_Star_X", "MJ2000_Star_Y", "MJ2000_Star_Z"]
    vector_columns = ["MJ2000_q1", "MJ2000_q2", "MJ2000_q3"]
    required = ["Timestamp", *position_columns]
    if attitude is None:
        required.extend(vector_columns)
    missing = [name for name in required if name not in pdf.columns]
    if missing:
        raise KeyError(f"missing pose columns: {missing}")
    positions = pdf[position_columns].to_numpy(dtype=float)
    valid = np.isfinite(pdf["Timestamp"]) & np.isfinite(positions).all(axis=1)
    valid &= (positions != 0).any(axis=1)
    if attitude is None:
        vectors = pdf[vector_columns].to_numpy(dtype=float)
        valid &= np.isfinite(vectors).all(axis=1) & (vectors != 0).any(axis=1)
        valid &= np.sum(vectors ** 2, axis=1) <= 1.0 + 1e-8
    pose = pdf.loc[valid].reset_index(drop=True)
    if len(pose) == 0:
        return pd.DataFrame(columns=columns)

    dates = pose["Datetime"].to_list()
    if attitude is None:
        rotation = _row_rotations(pose)
    else:
        rotation = quaternion_array(attitude, count=len(pdf))
        if len(rotation) != 1:
            rotation = rotation[valid]
    moon_eci2000_m = np.asarray(epms.moon(dates), dtype=float) * 1000.0
    earth_centered_eci2000_re = (
        positions[valid] + moon_eci2000_m
    ) / RE
    spacecraft = Coord(dates, earth_centered_eci2000_re, "ECI2000", 0)
    spacecraft.convert(coord)
    spacecraft_xyz = np.asarray(spacecraft[coord], dtype=float)

    source_look = Coord(
        dates,
        np.tile([0.0, 0.0, -1.0], (len(pose), 1)),
        probe,
        rotation,
    )
    source_look.convert(coord)
    look_xyz = np.asarray(source_look[coord], dtype=float)
    look_norm = np.linalg.norm(look_xyz, axis=1)
    look_xyz = np.divide(
        look_xyz,
        look_norm[:, np.newaxis],
        out=np.full_like(look_xyz, np.nan),
        where=look_norm[:, np.newaxis] > 0,
    )

    intersections = {
        "xy": _ray_plane_intersection(
            spacecraft_xyz, look_xyz, [0, 0, 1], [0, 0, 0], forward=forward
        ),
        "xz": _ray_plane_intersection(
            spacecraft_xyz, look_xyz, [0, 1, 0], [0, 0, 0], forward=forward
        ),
        "yz": _ray_plane_intersection(
            spacecraft_xyz, look_xyz, [1, 0, 0], [0, 0, 0], forward=forward
        ),
    }
    result = pd.DataFrame(
        {
            "Timestamp": pose["Timestamp"].to_numpy(dtype=float),
            "Datetime": pose["Datetime"],
            "spacecraft_x": spacecraft_xyz[:, 0],
            "spacecraft_y": spacecraft_xyz[:, 1],
            "spacecraft_z": spacecraft_xyz[:, 2],
            "look_x": look_xyz[:, 0],
            "look_y": look_xyz[:, 1],
            "look_z": look_xyz[:, 2],
        }
    )
    for plane, values in intersections.items():
        for index, component in enumerate("xyz"):
            result[f"{plane}_{component}"] = values[:, index]
    result.attrs.update(coordinate=str(coord).upper(), probe=probe, units="R_E")
    return result.reindex(columns=columns)


def orbit(
    pdf=None,
    *,
    geometry=None,
    probe="A",
    coord="GSM",
    cols=None,
    f=None,
    ax=None,
    figsize=None,
    xlim=None,
    grid=True,
):
    """Plot the QQ orbit and nominal FOV intersections with coordinate planes.

    ``cols`` follows the dynamic-panel convention of ``plots.joint.env``. The
    default layout groups the three spacecraft coordinates and the two varying
    coordinates on each of the XY, XZ, and YZ intersection planes.
    """
    if geometry is None:
        if pdf is None:
            raise ValueError("pdf is required when geometry is not supplied")
        geometry = pointing_geometry(pdf, probe=probe, coord=coord)
    if cols is None:
        cols = [
            ("spacecraft_x", "spacecraft_y", "spacecraft_z"),
            ("xy_x", "xy_y"),
            ("xz_x", "xz_z"),
            ("yz_y", "yz_z"),
        ]
    groups = _column_groups(cols)
    f, axes = _axes_for_columns(
        len(groups), f=f, ax=ax, figsize=figsize, sharex=True
    )
    times = geometry["Datetime"]
    palette = ("black", "tab:blue", "tab:red")
    for axis, group in zip(axes, groups):
        plotted = 0
        for index, name in enumerate(group):
            if name not in geometry.columns:
                continue
            axis.plot(
                times,
                geometry[name].to_numpy(dtype=float),
                color=palette[index % len(palette)],
                label=_ORBIT_LABELS.get(name, name),
            )
            plotted += 1
        if plotted:
            axis.legend(loc="upper right", ncol=min(3, plotted))
        else:
            axis.text(
                0.5, 0.5, f"No data: {', '.join(group)}",
                ha="center", va="center", transform=axis.transAxes,
            )
        axis.set_ylabel(f"{str(coord).upper()} ($R_E$)")
        if xlim is not None:
            axis.set_xlim(xlim)
        if grid:
            axis.grid(alpha=0.2, linewidth=0.5)
    locator = mdates.AutoDateLocator(minticks=4, maxticks=10)
    axes[-1].xaxis.set_major_locator(locator)
    axes[-1].xaxis.set_major_formatter(mdates.ConciseDateFormatter(locator))
    axes[-1].set_xlabel("UTC")
    f.autofmt_xdate()
    return f, axes[0] if len(axes) == 1 else axes

def get_velocity(df):
    dx = df['X2'] - df['X1']
    dy = df['Y2'] - df['Y1']
    dz = 50
    dt = df['Tof']
    #mm to km, and ns to s
    vx = (dx*1e-3*1e-3)/(dt*1e-9)
    vy = (dy*1e-3*1e-3)/(dt*1e-9)
    vz = (dz*1e-3*1e-3)/(dt*1e-9)
    return np.column_stack([vx,vy,vz])

def velocity(
    df: pd.DataFrame,
    *,
    coord="GSM",
    probe="A",
    vscale=1.0,
    inplace=False,
    direction_cols=("DX", "DY", "DZ"),
    velocity_cols=("VX", "VY", "VZ"),
) -> pd.DataFrame:
    """Calculate event look-direction and velocity columns in a target frame."""
    target = df if inplace else df.copy()
    required = ["X1", "Y1", "X2", "Y2", "Tof"]
    for column in required:
        if column not in target.columns:
            raise KeyError(f"Missing required column {column!r}.")
    dates = target["Datetime"].to_list()
    rotation = _row_rotations(target)
    ray_in_probe = np.column_stack(
        (
            target["X2"].to_numpy(dtype=float)
            - target["X1"].to_numpy(dtype=float),
            target["Y2"].to_numpy(dtype=float)
            - target["Y1"].to_numpy(dtype=float),
            np.full(len(target), 50.0),
        )
    )
    ray = Coord(dates, ray_in_probe, probe, rotation)
    ray.convert(coord)
    direction = np.asarray(ray[coord], dtype=float)
    direction_norm = np.linalg.norm(direction, axis=1)
    valid_direction = direction_norm > 0
    direction_unit = np.full_like(direction, np.nan, dtype=float)
    direction_unit[valid_direction] = direction[valid_direction] / direction_norm[
        valid_direction, np.newaxis
    ]

    event_velocity = Coord(dates, get_velocity(target), probe, rotation)
    event_velocity.convert(coord)
    velocity_values = np.asarray(event_velocity[coord], dtype=float)

    scale = vscale(target) if callable(vscale) else vscale
    scale = np.asarray(scale, dtype=float)
    if scale.ndim == 0:
        velocity_values = velocity_values * float(scale)
    else:
        if len(scale) != len(target):
            raise ValueError(
                "vscale must be scalar or have the same length as the DataFrame."
            )
        velocity_values = velocity_values * scale[:, np.newaxis]

    target.loc[:, list(direction_cols)] = direction_unit
    target.loc[:, list(velocity_cols)] = velocity_values
    return target



def deflect_volt(E):
    E_thres = 0.0 * E
    return E_thres
