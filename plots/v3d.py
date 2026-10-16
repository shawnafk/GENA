"""Interactive GENA 3-D geometry scenes and HTML export."""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path
import warnings

import pandas as pd
from plotly import graph_objects as go
import numpy as np
from scipy.interpolate import interp1d
from scipy.spatial.transform import Slerp

from libs import utils as common_utils
from common.config import get_configured_paths, resolve_output_path
from common.constants import RE
from common.coord import coord
from gena.constants import va1, va2, va3, va4, vb1, vb2, vb3, vb4
from gena.frame import _row_rotations
from load import load_gena
from libs import ephemeris as epms
from libs.environment import sw_keys
from libs.space import Bowshock, Magnetopause, T96NeutralSheet, T96_fieldlines

__all__ = [
    "SceneContext",
    "arrow",
    "dumpGeo",
    "line",
    "make_scene",
    "plane",
    "point",
    "prepare_scene",
    "pyramid",
    "scene",
    "sphere",
    "surface",
    "traces",
]

def _format_time(value):
    return pd.Timestamp(value).strftime('%Y-%m-%d %H:%M:%S')


def point(p, color='black', size=2, symbol='diamond', **kwargs):
    x, y, z = p
    return [go.Scatter3d(
        x=[x], y=[y], z=[z],
        mode='markers',
        marker=dict(symbol=symbol, size=size, color=color, **kwargs),
        showlegend=False,
    )]


def line(l, color='black', width=2, **kwargs):
    x, y, z = l
    return [go.Scatter3d(
        x=x, y=y, z=z,
        mode='lines',
        line=dict(width=width, color=color, **kwargs),
        showlegend=False,
    )]


def arrow(start, direction, c='blue'):
    end = np.array(start) + np.array(direction)
    trace = []
    trace.append(go.Scatter3d(
        x=[start[0], end[0]], y=[start[1], end[1]], z=[start[2], end[2]],
        mode='lines', line=dict(color=c, width=3), showlegend=False,
    ))
    trace.append(go.Cone(
        x=[end[0]], y=[end[1]], z=[end[2]],
        u=[direction[0]], v=[direction[1]], w=[direction[2]],
        sizemode='scaled', sizeref=0.2,
        colorscale=[[0, c], [1, c]],
        showscale=False, anchor='tip', showlegend=False,
    ))
    return trace


def pyramid(apex, base, facecolor='lightblue'):
    for i in range(3):
        base[0][i] = apex[i] + base[0][i]
        base[1][i] = apex[i] + base[1][i]
        base[2][i] = apex[i] + base[2][i]
        base[3][i] = apex[i] + base[3][i]
    vertices = np.vstack([apex, base])
    faces = [
        [0, 1, 2],
        [0, 2, 3],
        [0, 3, 4],
        [0, 4, 1],
        [1, 4, 3, 2],
    ]
    trace = []
    trace.append(go.Mesh3d(
        x=vertices[:, 0], y=vertices[:, 1], z=vertices[:, 2],
        i=[f[0] for f in faces[:-1]],
        j=[f[1] for f in faces[:-1]],
        k=[f[2] for f in faces[:-1]],
        color=facecolor, opacity=0.8,
    ))
    trace.append(go.Scatter3d(
        x=[vertices[i, 0] for i in faces[-1]],
        y=[vertices[i, 1] for i in faces[-1]],
        z=[vertices[i, 2] for i in faces[-1]],
        mode='lines', line=dict(color='red', width=1), showlegend=False,
    ))
    trace.append(go.Scatter3d(
        x=[apex[0]], y=[apex[1]], z=[apex[2]],
        mode='markers', marker=dict(color='red', size=1), name='apex',
    ))
    trace.append(go.Scatter3d(
        x=[p[0] for p in base], y=[p[1] for p in base], z=[p[2] for p in base],
        mode='markers', marker=dict(color='red', size=1), name='base',
    ))
    for bp in base:
        trace.append(go.Scatter3d(
            x=[apex[0], bp[0]], y=[apex[1], bp[1]], z=[apex[2], bp[2]],
            mode='lines', line=dict(color='black', width=3), showlegend=False,
        ))
    return trace


def sphere(center, radius, opacity=1):
    phi = np.linspace(0, np.pi, 15)
    theta = np.linspace(0, 2 * np.pi, 15)
    phi, theta = np.meshgrid(phi, theta)
    x = center[0] + radius * np.sin(phi) * np.cos(theta)
    y = center[1] + radius * np.sin(phi) * np.sin(theta)
    z = center[2] + radius * np.cos(phi)
    relative_x = x - center[0]
    return [go.Surface(
        x=x, y=y, z=z,
        surfacecolor=relative_x,
        colorscale=[
            [0, 'black'],
            [0.5, 'black'],
            [0.5, 'red'],
            [1, 'red'],
        ],
        cmin=-radius, cmax=radius,
        colorbar=dict(
            title='Relative X-Coordinate',
            tickvals=[-radius, 0, radius],
            ticktext=['Negative', 'Zero', 'Positive'],
        ),
        opacity=opacity,
        showscale=False,
    )]


def plane(
    coef,
    x_range=(-80, 20),
    y_range=(-30, 30),
    z_range=(-30, 30),
    color='lightblue',
    opacity=0.1,
    name=None,
):
    A, B, C, D = coef
    if C != 0:
        x = np.linspace(x_range[0], x_range[1], 5)
        y = np.linspace(y_range[0], y_range[1], 5)
        X, Y = np.meshgrid(x, y)
        Z = (-(A * X + B * Y + D)) / C
        mask = (Z >= z_range[0]) & (Z <= z_range[1])
        X = np.where(mask, X, np.nan)
        Y = np.where(mask, Y, np.nan)
        Z = np.where(mask, Z, np.nan)
    elif B != 0:
        x = np.linspace(x_range[0], x_range[1], 5)
        z = np.linspace(z_range[0], z_range[1], 5)
        X, Z = np.meshgrid(x, z)
        Y = (-(A * X + C * Z + D)) / B
        mask = (Y >= y_range[0]) & (Y <= y_range[1])
        X = np.where(mask, X, np.nan)
        Y = np.where(mask, Y, np.nan)
        Z = np.where(mask, Z, np.nan)
    elif A != 0:
        y = np.linspace(y_range[0], y_range[1], 5)
        z = np.linspace(z_range[0], z_range[1], 5)
        Y, Z = np.meshgrid(y, z)
        X = (-(B * Y + C * Z + D)) / A
        mask = (X >= x_range[0]) & (X <= x_range[1])
        X = np.where(mask, X, np.nan)
        Y = np.where(mask, Y, np.nan)
        Z = np.where(mask, Z, np.nan)
    else:
        raise ValueError("plane coefficients A, B, C cannot all be zero")
    return [go.Surface(
        x=X, y=Y, z=Z,
        colorscale=[[0, color], [1, color]],
        opacity=opacity,
        showscale=False,
        name=name,
        legendgroup=name,
        showlegend=name is not None,
        hoverinfo='skip',
    )]


def surface(xyz_data, color=0, line_width=1, opacity=0.5, **kwargs):
    xx, yy, zz = xyz_data
    color_array = np.full_like(zz, color)
    return [go.Surface(
        x=xx, y=yy, z=zz,
        contours=dict(
            x=dict(show=True, width=line_width),
            y=dict(show=True, width=line_width),
            z=dict(show=True, width=line_width),
        ),
        opacity=opacity,
        surfacecolor=color_array,
        colorscale=[
            [0, 'gray'],
            [1, 'lightblue'],
        ],
        showscale=False,
        **kwargs,
    )]


def traces(qq, qax, qvxa, qvxb, vscale, gap, moon=None, sats=None):
    trace_qq = point([qq[0], qq[1], qq[2]], color='black', symbol='circle')

    qx, qy, qz = qax
    trace_qq_x = arrow([qq[0], qq[1], qq[2]], qx * vscale, c='red')
    trace_qq_y = arrow([qq[0], qq[1], qq[2]], qy * vscale, c='green')
    trace_qq_z = arrow([qq[0], qq[1], qq[2]], qz * vscale, c='blue')

    da = -gap
    db = gap
    v1a, v2a, v3a, v4a = qvxa * vscale
    v1b, v2b, v3b, v4b = qvxb * vscale

    trace_fov_a = pyramid(
        [qq[0] + da * qy[0], qq[1] + da * qy[1], qq[2] + da * qy[2]],
        [
            [v1a[0], v1a[1], v1a[2]],
            [v2a[0], v2a[1], v2a[2]],
            [v3a[0], v3a[1], v3a[2]],
            [v4a[0], v4a[1], v4a[2]],
        ],
        facecolor='green',
    )
    trace_fov_b = pyramid(
        [qq[0] + db * qy[0], qq[1] + db * qy[1], qq[2] + db * qy[2]],
        [
            [v1b[0], v1b[1], v1b[2]],
            [v2b[0], v2b[1], v2b[2]],
            [v3b[0], v3b[1], v3b[2]],
            [v4b[0], v4b[1], v4b[2]],
        ],
        facecolor='blue',
    )

    all_traces = [
        *trace_qq,
        *trace_qq_x,
        *trace_qq_y,
        *trace_qq_z,
        *trace_fov_a,
        *trace_fov_b,
    ]
    if moon is not None:
        all_traces.extend(sphere(moon, 0.27))
    if sats is not None:
        for sat in sats:
            trace_sat = point(sat)
            all_traces.append(*trace_sat)

    fig = go.Figure()
    fig.add_traces(all_traces)
    return fig


def _traces_raw(
    qq,
    qax,
    qvxa,
    qvxb,
    vscale,
    gap,
    moon=None,
    sats=None,
    sat_colors=None,
    satellite_names=None,
    satellite_marker_size=6,
):
    trace_qq = point([qq[0], qq[1], qq[2]], color='black', symbol='circle')
    qx, qy, qz = qax
    trace_qq_x = arrow([qq[0], qq[1], qq[2]], qx * vscale, c='red')
    trace_qq_y = arrow([qq[0], qq[1], qq[2]], qy * vscale, c='green')
    trace_qq_z = arrow([qq[0], qq[1], qq[2]], qz * vscale, c='blue')
    da = -gap
    db = gap
    v1a, v2a, v3a, v4a = qvxa * vscale
    v1b, v2b, v3b, v4b = qvxb * vscale
    trace_fov_a = pyramid(
        [qq[0] + da * qy[0], qq[1] + da * qy[1], qq[2] + da * qy[2]],
        [
            [v1a[0], v1a[1], v1a[2]],
            [v2a[0], v2a[1], v2a[2]],
            [v3a[0], v3a[1], v3a[2]],
            [v4a[0], v4a[1], v4a[2]],
        ],
        facecolor='green',
    )
    trace_fov_b = pyramid(
        [qq[0] + db * qy[0], qq[1] + db * qy[1], qq[2] + db * qy[2]],
        [
            [v1b[0], v1b[1], v1b[2]],
            [v2b[0], v2b[1], v2b[2]],
            [v3b[0], v3b[1], v3b[2]],
            [v4b[0], v4b[1], v4b[2]],
        ],
        facecolor='blue',
    )
    all_traces = [
        *trace_qq,
        *trace_qq_x,
        *trace_qq_y,
        *trace_qq_z,
        *trace_fov_a,
        *trace_fov_b,
    ]
    if moon is not None:
        all_traces.extend(sphere(moon, 0.27))
    if sats is not None:
        sat_colors = sat_colors or ['black'] * len(sats)
        satellite_names = satellite_names or [
            f'Satellite {index + 1}' for index in range(len(sats))
        ]
        for sat, color, name in zip(sats, sat_colors, satellite_names):
            all_traces.append(go.Scatter3d(
                x=[sat[0]],
                y=[sat[1]],
                z=[sat[2]],
                mode='markers',
                marker={
                    'symbol': 'diamond',
                    'size': satellite_marker_size,
                    'color': color,
                    'line': {'color': 'white', 'width': 1},
                },
                name=f'{name} current position',
                legendgroup=f'satellite-{name}',
                showlegend=False,
                hovertemplate=(
                    f'{name}<br>X=%{{x:.2f}} R_E<br>Y=%{{y:.2f}} R_E'
                    '<br>Z=%{z:.2f} R_E<extra></extra>'
                ),
            ))
    return all_traces


def _satellite_trajectory_traces(trajectory, *, name, color):
    trajectory = np.asarray(trajectory, dtype=float)
    traces_out = [go.Scatter3d(
        x=trajectory[0],
        y=trajectory[1],
        z=trajectory[2],
        mode='lines',
        line={'width': 4, 'color': color},
        opacity=0.8,
        name=f'{name} trajectory',
        legendgroup=f'satellite-{name}',
        showlegend=True,
        hovertemplate=(
            f'{name} trajectory<br>X=%{{x:.2f}} R_E<br>Y=%{{y:.2f}} R_E'
            '<br>Z=%{z:.2f} R_E<extra></extra>'
        ),
    )]
    if trajectory.shape[1] >= 2:
        direction = trajectory[:, -1] - trajectory[:, -2]
        norm = np.linalg.norm(direction)
        if np.isfinite(norm) and norm > 0:
            direction = direction / norm
            traces_out.append(go.Cone(
                x=[trajectory[0, -1]],
                y=[trajectory[1, -1]],
                z=[trajectory[2, -1]],
                u=[direction[0]],
                v=[direction[1]],
                w=[direction[2]],
                anchor='tip',
                sizemode='absolute',
                sizeref=0.8,
                colorscale=[[0, color], [1, color]],
                showscale=False,
                name=f'{name} direction',
                legendgroup=f'satellite-{name}',
                showlegend=False,
                hoverinfo='skip',
            ))
    return traces_out


def _field_line_traces(field_lines, color='#E3A018'):
    traces_out = []
    for index, coordinates in enumerate(field_lines):
        coordinates = np.asarray(coordinates, dtype=float)
        if coordinates.ndim != 2 or coordinates.shape[0] != 3:
            raise ValueError('each T96 field line must have shape (3, N)')
        traces_out.append(go.Scatter3d(
            x=coordinates[0],
            y=coordinates[1],
            z=coordinates[2],
            mode='lines',
            line={'width': 2, 'color': color},
            opacity=0.75,
            name='T96 field lines',
            legendgroup='t96-field-lines',
            showlegend=index == 0,
        ))
    return traces_out


def _neutral_sheet_traces(sheet):
    if sheet is None:
        return []
    X, Y, Z = sheet.surface_xyz()
    lower, upper = sheet.thickness_xyz()
    central_color = np.zeros_like(Z)
    lower_color = np.zeros_like(lower[2])
    upper_color = np.zeros_like(upper[2])
    common = {
        'showscale': False,
        'hovertemplate': (
            'T96+IGRF Bx=0 neutral sheet'
            '<br>X=%{x:.2f} R_E<br>Y=%{y:.2f} R_E'
            '<br>Z=%{z:.2f} R_E<extra></extra>'
        ),
    }
    return [
        go.Surface(
            x=X,
            y=Y,
            z=Z,
            surfacecolor=central_color,
            colorscale=[[0, '#C0392B'], [1, '#C0392B']],
            opacity=0.32,
            name='T96 Bx=0 neutral sheet',
            legendgroup='t96-neutral-sheet',
            showlegend=True,
            **common,
        ),
        go.Surface(
            x=lower[0],
            y=lower[1],
            z=lower[2],
            surfacecolor=lower_color,
            colorscale=[[0, '#F5B7B1'], [1, '#F5B7B1']],
            opacity=0.08,
            name='T96 sheet lower visual bound',
            legendgroup='t96-neutral-sheet',
            showlegend=False,
            **common,
        ),
        go.Surface(
            x=upper[0],
            y=upper[1],
            z=upper[2],
            surfacecolor=upper_color,
            colorscale=[[0, '#F5B7B1'], [1, '#F5B7B1']],
            opacity=0.08,
            name='T96 sheet upper visual bound',
            legendgroup='t96-neutral-sheet',
            showlegend=False,
            **common,
        ),
    ]


def _fov_z0_intersection(spacecraft_position, fov_corners, *, vscale=60):
    direction = np.mean(np.asarray(fov_corners, dtype=float), axis=0) * vscale
    if not np.isfinite(direction).all() or abs(direction[2]) < 1e-12:
        return None
    scale = -spacecraft_position[2] / direction[2]
    if not 0 <= scale <= 1:
        return None
    return spacecraft_position + scale * direction


def _scene_annotation(time, sheet, spacecraft_position, fov_a_corners):
    lines = [
        '<b>T96 geometry</b>',
        f'UTC: {_format_time(time)}',
    ]
    if sheet is None:
        lines.append('Neutral sheet: disabled')
    else:
        pressure, dst, by, bz = sheet.parameters
        lines.extend([
            f'Dipole tilt ψ: {sheet.dipole_tilt_degrees:+.2f}°',
            f'Pdyn: {pressure:.2f} nPa · Dst: {dst:.1f} nT',
            f'IMF By/Bz: {by:+.2f} / {bz:+.2f} nT',
            'Sheet: total T96+IGRF Bx=0',
            f'Visual half-thickness: ±{sheet.half_thickness:.1f} R_E',
        ])
        intersection = _fov_z0_intersection(
            np.asarray(spacecraft_position, dtype=float),
            fov_a_corners,
        )
        if intersection is None:
            lines.append('A FOV has no Z=0 intersection in displayed cone')
        else:
            sheet_z = sheet.z_at(intersection[0], intersection[1])
            lines.append(
                'A LOS at Z=0: '
                f'X={intersection[0]:.1f}, Y={intersection[1]:.1f} R_E'
            )
            if np.isfinite(sheet_z):
                offset = intersection[2] - sheet_z
                state = (
                    'inside visual band'
                    if abs(offset) <= sheet.half_thickness
                    else 'outside visual band'
                )
                lines.append(
                    f'Bx=0 there: Z={sheet_z:+.2f} R_E '
                    f'(ΔZ={offset:+.2f}; {state})'
                )
            else:
                lines.append('Bx=0 there: outside rendered sheet domain')
    return {
        'xref': 'paper',
        'yref': 'paper',
        'x': 0.99,
        'y': 0.99,
        'xanchor': 'right',
        'yanchor': 'top',
        'align': 'left',
        'showarrow': False,
        'text': '<br>'.join(lines),
        'font': {'size': 11, 'color': '#202020'},
        'bgcolor': 'rgba(255, 255, 255, 0.84)',
        'bordercolor': 'rgba(80, 80, 80, 0.55)',
        'borderwidth': 1,
        'borderpad': 6,
    }


def scene(
    st_datetime_list,
    st_pose,
    sat_orb,
    bs_ave,
    mp_ave,
    field_lines=None,
    neutral_sheets=None,
    *,
    sat_trajectories=None,
    satellite_names=None,
    include_reference_planes=True,
    satellite_marker_size=6,
):
    """Build the original interactive GENA geometry scene.

    ``st_pose`` positions are Moon-relative MJ2000 metres. ``sat_orb`` entries
    contain the current satellite positions with shape ``(3, N)``.
    ``sat_trajectories`` contains the wider orbit context around those scene
    times. Both and every rendered coordinate are in Earth radii (``R_E``).

    ``field_lines`` and ``neutral_sheets`` contain time-dependent T96 magnetic
    geometry. Field lines are preserved alongside the total T96+IGRF ``Bx=0``
    sheet; GSM ``Z=0`` and ``Y=0`` remain fixed coordinate references.
    """
    Nt = len(st_datetime_list)
    if Nt < 1:
        raise ValueError('scene requires at least one pose sample')
    st_pos, st_rot = st_pose
    st_timestamp = [t.timestamp() for t in st_datetime_list]
    checked_satellite_orbits = []
    for orbit in sat_orb:
        orbit = np.asarray(orbit, dtype=float)
        if orbit.shape != (3, len(st_datetime_list)):
            raise ValueError(
                'each satellite orbit must have shape (3, N) in GSM R_E units'
            )
        orbit_norm_re = np.linalg.norm(orbit, axis=0)
        if not np.isfinite(orbit_norm_re).all() or np.max(orbit_norm_re) > 100:
            raise ValueError(
                'satellite orbit contains nonphysical coordinates; scene '
                'expects GSM positions in R_E and does not permit distances '
                'greater than 100 R_E'
            )
        checked_satellite_orbits.append(orbit)
    sat_orb = checked_satellite_orbits
    if sat_trajectories is None:
        sat_trajectories = sat_orb
    if len(sat_trajectories) != len(sat_orb):
        raise ValueError(
            'sat_trajectories and current satellite positions must have the '
            'same length'
        )
    checked_trajectories = []
    for trajectory in sat_trajectories:
        trajectory = np.asarray(trajectory, dtype=float)
        if trajectory.ndim != 2 or trajectory.shape[0] != 3:
            raise ValueError('each satellite trajectory must have shape (3, M)')
        trajectory_norm_re = np.linalg.norm(trajectory, axis=0)
        if (
            trajectory.shape[1] < 1
            or not np.isfinite(trajectory_norm_re).all()
            or np.max(trajectory_norm_re) > 100
        ):
            raise ValueError(
                'satellite trajectory contains nonphysical coordinates; scene '
                'expects GSM positions in R_E and does not permit distances '
                'greater than 100 R_E'
            )
        checked_trajectories.append(trajectory)
    sat_trajectories = checked_trajectories
    if satellite_names is None:
        satellite_names = [
            f'Satellite {index + 1}' for index in range(len(sat_orb))
        ]
    if len(satellite_names) != len(sat_orb):
        raise ValueError('satellite_names must match the satellite orbit count')
    if field_lines is None:
        field_lines = [[] for _ in range(Nt)]
    if len(field_lines) != Nt:
        raise ValueError('field_lines must contain one line list per scene time')
    field_line_counts = [len(lines_at_time) for lines_at_time in field_lines]
    if len(set(field_line_counts)) > 1:
        raise ValueError('every scene time must contain the same number of field lines')
    if neutral_sheets is None:
        neutral_sheets = [None] * Nt
    if len(neutral_sheets) != Nt:
        raise ValueError('neutral_sheets must contain one sheet per scene time')
    moon_pos_re = epms.moon(pd.to_datetime(st_timestamp, unit='s')) * 1000 / RE
    spacecraft_offset_re = np.asarray(st_pos, dtype=float) / RE
    spacecraft_pos_ej2000_re = moon_pos_re + spacecraft_offset_re

    t_reference = st_datetime_list[0]
    date_series = pd.date_range(
        start=t_reference,
        periods=28 * 24,
        freq='h',
    ).to_list()
    full_moon_pos = coord(date_series, epms.moon(date_series) * 1000 / RE, 'ECI2000', 0)
    full_moon_pos.convert('GSM')

    coord_moon = coord(st_datetime_list, moon_pos_re, 'ECI2000', 0)
    coord_moon.convert('GSM')

    coord_qq = coord(st_datetime_list, spacecraft_pos_ej2000_re, 'ECI2000', st_rot)
    coord_qq.convert('GSM')

    coord_qq_minus_z = coord(st_datetime_list, np.tile([0, 0, -1], (Nt, 1)), 'QQ', st_rot)
    coord_qq_minus_z.convert('GSM')

    coord_xyz = [coord(st_datetime_list, np.tile(v, (Nt, 1)), 'QQ', st_rot) for v in [[1, 0, 0], [0, 1, 0], [0, 0, 1]]]
    [_.convert('GSM') for _ in coord_xyz]

    va = [coord(st_datetime_list, np.tile(v, (Nt, 1)), 'QQ', st_rot) for v in [va1, va2, va3, va4]]
    [_.convert('GSM') for _ in va]

    vb = [coord(st_datetime_list, np.tile(v, (Nt, 1)), 'QQ', st_rot) for v in [vb1, vb2, vb3, vb4]]
    [_.convert('GSM') for _ in vb]

    fig = go.Figure()
    tid_list = range(Nt)
    tid_start = 0
    sat_palette = ['red', 'blue', 'green', 'orange', 'purple', 'black']
    sat_colors = [
        sat_palette[index % len(sat_palette)] for index in range(len(sat_orb))
    ]
    Trace_satellites = []
    for trajectory, name, color in zip(
        sat_trajectories,
        satellite_names,
        sat_colors,
    ):
        Trace_satellites.extend(
            _satellite_trajectory_traces(
                trajectory,
                name=name,
                color=color,
            )
        )

    Trace_moon_in_month = line(full_moon_pos['GSM'].transpose(), color='black')
    Trace_earth = sphere([0, 0, 0], 1)
    trace_bs = surface(bs_ave.surface_xyz(), color=0, opacity=0.1)
    trace_mp = surface(mp_ave.surface_xyz(), color=1, opacity=0.1)
    trace_reference_planes = []
    if include_reference_planes:
        trace_reference_planes.extend(plane(
            [0, 0, 1, 0],
            color='#B8C4CE',
            opacity=0.10,
            name='GSM Z=0 reference',
        ))
        trace_reference_planes.extend(plane(
            [0, 1, 0, 0],
            color='#D7B5E5',
            opacity=0.08,
            name='GSM Y=0 reference',
        ))

    fig.add_traces(
        Trace_earth
        + trace_bs
        + trace_mp
        + trace_reference_planes
        + Trace_moon_in_month
        + Trace_satellites
    )

    initial_dynamic_traces = _traces_raw(
        coord_qq['GSM'][tid_start],
        np.array([_['GSM'][tid_start] for _ in coord_xyz]),
        np.array([_['GSM'][tid_start] for _ in va]),
        np.array([_['GSM'][tid_start] for _ in vb]),
        vscale=60,
        gap=0,
        moon=coord_moon['GSM'][tid_start],
        sats=[_[:, tid_start] for _ in sat_orb],
        sat_colors=sat_colors,
        satellite_names=satellite_names,
        satellite_marker_size=satellite_marker_size,
    )
    initial_dynamic_traces.extend(_field_line_traces(field_lines[tid_start]))
    initial_dynamic_traces.extend(_neutral_sheet_traces(neutral_sheets[tid_start]))
    dynamic_start = len(fig.data)
    fig.add_traces(initial_dynamic_traces)
    initial_dynamic_count = len(initial_dynamic_traces)
    dynamic_indices = list(range(dynamic_start, dynamic_start + initial_dynamic_count))

    frames = []
    for tid in tid_list:
        current_time = st_datetime_list[tid]
        formatted_time = _format_time(current_time)

        dynamic_traces = _traces_raw(
            coord_qq['GSM'][tid],
            np.array([_['GSM'][tid] for _ in coord_xyz]),
            np.array([_['GSM'][tid] for _ in va]),
            np.array([_['GSM'][tid] for _ in vb]),
            vscale=60,
            gap=0,
            moon=coord_moon['GSM'][tid],
            sats=[_[:, tid] for _ in sat_orb],
            sat_colors=sat_colors,
            satellite_names=satellite_names,
            satellite_marker_size=satellite_marker_size,
        )
        dynamic_traces.extend(_field_line_traces(field_lines[tid]))
        dynamic_traces.extend(_neutral_sheet_traces(neutral_sheets[tid]))

        frames.append(go.Frame(
            data=dynamic_traces,
            traces=dynamic_indices,
            name=formatted_time,
            layout=go.Layout(annotations=[_scene_annotation(
                current_time,
                neutral_sheets[tid],
                coord_qq['GSM'][tid],
                np.array([_['GSM'][tid] for _ in va]),
            )]),
        ))

    if Nt > 1:
        fig.frames = frames
        controls = {
            'sliders': [{
            'active': 0,
            'currentvalue': {
                'prefix': 'time: ',
                'visible': True,
                'font': {'size': 14},
            },
            'steps': [{
                'method': 'animate',
                'args': [[f.name], {'frame': {'duration': 300, 'redraw': True}, 'mode': 'immediate'}],
                'label': f.name,
            } for f in fig.frames],
            }],
            'updatemenus': [{
            'type': 'buttons',
            'direction': 'left',
            'x': 0.02,
            'y': 0.02,
            'buttons': [
                {
                    'label': 'Play',
                    'method': 'animate',
                    'args': [None, {
                        'fromcurrent': True,
                        'frame': {'duration': 180, 'redraw': True},
                        'transition': {'duration': 0},
                    }],
                },
                {
                    'label': 'Pause',
                    'method': 'animate',
                    'args': [[None], {
                        'mode': 'immediate',
                        'frame': {'duration': 0, 'redraw': False},
                        'transition': {'duration': 0},
                    }],
                },
            ],
            }],
        }
        title = (
            f'timerange: {_format_time(st_datetime_list[0])} to '
            f'{_format_time(st_datetime_list[-1])}'
        )
    else:
        controls = {'sliders': [], 'updatemenus': []}
        title = f'time: {_format_time(st_datetime_list[0])}'
    fig.update_layout(
        **controls,
        scene=dict(
            aspectmode='data',
            xaxis_title='X GSM (R_E)',
            yaxis_title='Y GSM (R_E)',
            zaxis_title='Z GSM (R_E)',
        ),
        title={
            'text': title,
            'font': {'size': 16},
        },
        annotations=[_scene_annotation(
            st_datetime_list[tid_start],
            neutral_sheets[tid_start],
            coord_qq['GSM'][tid_start],
            np.array([_['GSM'][tid_start] for _ in va]),
        )],
    )
    return fig


@dataclass
class SceneContext:
    """Preloaded GENA pose, satellite, and space-weather scene data.

    Construct this object once for a day or time range, then call
    :meth:`make_scene` repeatedly. Parquet and OMNI data are read only during
    construction; subsequent calls use cached interpolators in memory.
    """

    pose_timestamps: np.ndarray
    pose_positions: np.ndarray
    pose_rotations: object
    cache_timestamps: np.ndarray
    orbit_timestamps: np.ndarray
    satellite_orbits: dict[str, np.ndarray]
    space_weather: dict[str, np.ndarray]
    include_field_lines: bool = True
    include_neutral_sheet: bool = True
    include_reference_planes: bool = True
    field_line_nlat: int = 21
    field_line_footprint_latitudes: tuple[float, float] = (70.0, 82.0)
    field_line_rlim: float = 80.0
    field_line_cadence_seconds: float = 3600.0
    neutral_sheet_half_thickness: float = 3.0
    satellite_trajectory_padding_seconds: float = 3600.0
    satellite_marker_size: float = 6.0
    omitted_satellites: tuple[str, ...] = ()
    _position_interpolators: tuple = field(init=False, repr=False)
    _rotation_interpolator: object = field(init=False, repr=False)
    _field_line_cache: dict[float, list[np.ndarray]] = field(
        default_factory=dict,
        init=False,
        repr=False,
    )
    _neutral_sheet_cache: dict[float, T96NeutralSheet] = field(
        default_factory=dict,
        init=False,
        repr=False,
    )

    def __post_init__(self):
        self._position_interpolators = tuple(
            interp1d(
                self.pose_timestamps,
                self.pose_positions[:, axis],
                bounds_error=True,
            )
            for axis in range(3)
        )
        self._rotation_interpolator = Slerp(
            self.pose_timestamps,
            self.pose_rotations,
        )

    @classmethod
    def load(
        cls,
        trange,
        *,
        satellites=('THA', 'ATMS1', 'ATMS2', 'THD', 'THE', 'MMS1'),
        cache_cadence='5min',
        use_space_weather=True,
        include_field_lines=True,
        include_neutral_sheet=True,
        include_reference_planes=True,
        field_line_nlat=21,
        field_line_footprint_latitudes=(70.0, 82.0),
        field_line_rlim=80,
        field_line_cadence='1h',
        neutral_sheet_half_thickness=3.0,
        satellite_trajectory_padding='1h',
        satellite_marker_size=6,
        update_ephemeris=False,
    ):
        """Preload every data source needed for scenes in ``trange``."""
        position_columns = [
            'MJ2000_Star_X', 'MJ2000_Star_Y', 'MJ2000_Star_Z'
        ]
        rotation_columns = ['MJ2000_q1', 'MJ2000_q2', 'MJ2000_q3']
        pose_columns = ['Timestamp', *position_columns, *rotation_columns]
        pose = load_gena(
            trange=trange,
            product="pdf",
            columns=pose_columns,
        )
        values = pose[pose_columns].to_numpy(dtype=float)
        valid = np.isfinite(values).all(axis=1)
        valid &= (values[:, 1:4] != 0).all(axis=1)
        valid &= (values[:, 4:7] != 0).all(axis=1)
        valid &= np.sum(values[:, 4:7] ** 2, axis=1) <= 1.0 + 1e-8
        pose = pose.loc[valid].reset_index(drop=True)
        if len(pose) < 2:
            raise ValueError('fewer than two valid position/attitude samples')

        pose_timestamps = pose['Timestamp'].to_numpy(dtype=float)
        pose_positions = pose[position_columns].to_numpy(dtype=float)
        pose_rotations = _row_rotations(pose)
        cache_start = pd.to_datetime(pose_timestamps[0], unit='s', utc=True)
        cache_stop = pd.to_datetime(pose_timestamps[-1], unit='s', utc=True)
        cadence = pd.Timedelta(cache_cadence)
        if cadence <= pd.Timedelta(0):
            raise ValueError('cache_cadence must be positive')
        if not isinstance(field_line_nlat, int) or field_line_nlat < 1:
            raise ValueError('field_line_nlat must be a positive integer')
        field_line_footprint_latitudes = tuple(
            float(value) for value in field_line_footprint_latitudes
        )
        if (
            len(field_line_footprint_latitudes) != 2
            or not 0
            < field_line_footprint_latitudes[0]
            < field_line_footprint_latitudes[1]
            < 90
        ):
            raise ValueError(
                'field_line_footprint_latitudes must be two increasing values '
                'between 0 and 90 degrees'
            )
        if not np.isfinite(field_line_rlim) or field_line_rlim <= 3.5:
            raise ValueError('field_line_rlim must be greater than 3.5 R_E')
        if (
            not np.isfinite(neutral_sheet_half_thickness)
            or neutral_sheet_half_thickness <= 0
        ):
            raise ValueError('neutral_sheet_half_thickness must be positive')
        field_line_cadence = pd.Timedelta(field_line_cadence)
        if field_line_cadence <= pd.Timedelta(0):
            raise ValueError('field_line_cadence must be positive')
        satellite_trajectory_padding = pd.Timedelta(
            satellite_trajectory_padding
        )
        if satellite_trajectory_padding < pd.Timedelta(0):
            raise ValueError('satellite_trajectory_padding must not be negative')
        if (
            not np.isfinite(satellite_marker_size)
            or satellite_marker_size <= 0
        ):
            raise ValueError('satellite_marker_size must be positive')
        cache_dates = pd.date_range(cache_start, cache_stop, freq=cadence)
        cache_dates = cache_dates.union(pd.DatetimeIndex([cache_start, cache_stop]))
        cache_timestamps = common_utils.Time(cache_dates).to_timestamp()
        orbit_start = cache_start - satellite_trajectory_padding
        orbit_stop = cache_stop + satellite_trajectory_padding
        orbit_dates = pd.date_range(orbit_start, orbit_stop, freq=cadence)
        orbit_dates = orbit_dates.union(
            pd.DatetimeIndex([orbit_start, orbit_stop])
        )
        orbit_timestamps = common_utils.Time(orbit_dates).to_timestamp()

        satellite_orbits = {}
        omitted = []
        for satellite in dict.fromkeys(str(name).upper() for name in satellites):
            try:
                orbit_gse_re = epms.satellites(
                    orbit_timestamps,
                    satellite,
                    update=update_ephemeris,
                )
                orbit = coord(orbit_dates.to_list(), orbit_gse_re, 'GSE', 0)
                orbit.convert('GSM')
                orbit_gsm_re = np.asarray(orbit['GSM'], dtype=float)
                norms = np.linalg.norm(orbit_gsm_re, axis=1)
                if not np.isfinite(norms).all() or np.max(norms) > 100:
                    raise ValueError('nonphysical orbit coordinates')
                satellite_orbits[satellite] = orbit_gsm_re
            except (KeyError, OSError, RuntimeError, TypeError, ValueError):
                omitted.append(satellite)

        if omitted:
            warnings.warn(
                'omitted unavailable satellite orbit(s): ' + ', '.join(omitted),
                stacklevel=2,
            )

        if use_space_weather:
            raw_weather = sw_keys(cache_timestamps, download=False)
            space_weather = {
                name: np.asarray(raw_weather.get(name, []), dtype=float)
                for name in (
                    'BY_GSM',
                    'BZ_GSM',
                    'Beta',
                    'Mgs_mach_num',
                    'Pressure',
                    'dst',
                )
            }
        else:
            space_weather = {}

        return cls(
            pose_timestamps=pose_timestamps,
            pose_positions=pose_positions,
            pose_rotations=pose_rotations,
            cache_timestamps=np.asarray(cache_timestamps, dtype=float),
            orbit_timestamps=np.asarray(orbit_timestamps, dtype=float),
            satellite_orbits=satellite_orbits,
            space_weather=space_weather,
            include_field_lines=include_field_lines,
            include_neutral_sheet=include_neutral_sheet,
            include_reference_planes=include_reference_planes,
            field_line_nlat=field_line_nlat,
            field_line_footprint_latitudes=field_line_footprint_latitudes,
            field_line_rlim=float(field_line_rlim),
            field_line_cadence_seconds=field_line_cadence.total_seconds(),
            neutral_sheet_half_thickness=float(neutral_sheet_half_thickness),
            satellite_trajectory_padding_seconds=(
                satellite_trajectory_padding.total_seconds()
            ),
            satellite_marker_size=float(satellite_marker_size),
            omitted_satellites=tuple(omitted),
        )

    def _normalize_times(self, times):
        normalized = common_utils.Time(times).datetime
        if isinstance(normalized, pd.Timestamp):
            dates = pd.DatetimeIndex([normalized])
        else:
            dates = pd.DatetimeIndex(normalized)
        dates = dates.tz_localize('UTC') if dates.tz is None else dates.tz_convert('UTC')
        timestamps = np.asarray(common_utils.Time(dates).to_timestamp(), dtype=float)
        if len(timestamps) == 0:
            raise ValueError('at least one scene time is required')
        if timestamps.min() < self.pose_timestamps[0] or timestamps.max() > self.pose_timestamps[-1]:
            coverage = (
                pd.to_datetime(self.pose_timestamps[0], unit='s', utc=True),
                pd.to_datetime(self.pose_timestamps[-1], unit='s', utc=True),
            )
            raise ValueError(
                f'scene time is outside loaded pose coverage {coverage[0]} to '
                f'{coverage[1]}'
            )
        return dates, timestamps

    def _evaluate_weather(self, timestamps):
        defaults = {
            'BY_GSM': 0.0,
            'BZ_GSM': 0.0,
            'Beta': 5.0,
            'Mgs_mach_num': 8.0,
            'Pressure': 5.0,
            'dst': -40.0,
        }
        values = {}
        for name, default in defaults.items():
            source = self.space_weather.get(name)
            if source is None or len(source) != len(self.cache_timestamps):
                evaluated = np.full(len(timestamps), default, dtype=float)
            else:
                valid = np.isfinite(source)
                if valid.sum() < 2:
                    evaluated = np.full(len(timestamps), default, dtype=float)
                else:
                    evaluated = np.interp(
                        timestamps,
                        self.cache_timestamps[valid],
                        source[valid],
                    )
            values[name] = evaluated
        return values

    def _boundaries(self, timestamps, *, average):
        evaluated = self._evaluate_weather(timestamps)
        values = {
            name: float(np.mean(value)) if average else float(value[0])
            for name, value in evaluated.items()
        }
        return (
            Bowshock(
                values['BZ_GSM'],
                values['Beta'],
                values['Mgs_mach_num'],
                values['Pressure'],
            ),
            Magnetopause(values['BZ_GSM'], values['Pressure']),
        )

    def _field_lines(self, timestamps, *, exact):
        if not self.include_field_lines:
            return None
        if exact:
            cache_keys = np.asarray(timestamps, dtype=float)
        else:
            cache_keys = (
                np.floor(
                    np.asarray(timestamps, dtype=float)
                    / self.field_line_cadence_seconds
                )
                * self.field_line_cadence_seconds
            )
        for cache_key in np.unique(cache_keys):
            cache_key = float(cache_key)
            if cache_key in self._field_line_cache:
                continue
            weather = self._evaluate_weather(np.array([cache_key]))
            model = T96_fieldlines(
                cache_key,
                float(weather['Pressure'][0]),
                float(weather['dst'][0]),
                float(weather['BY_GSM'][0]),
                float(weather['BZ_GSM'][0]),
                Nlat=self.field_line_nlat,
            )
            self._field_line_cache[cache_key] = model.trace_lines_y0(
                footprint_latitude_range=self.field_line_footprint_latitudes,
                rlim=self.field_line_rlim,
            )
        return [self._field_line_cache[float(key)] for key in cache_keys]

    def _neutral_sheets(self, timestamps, *, exact):
        if not self.include_neutral_sheet:
            return [None] * len(timestamps)
        if exact:
            cache_keys = np.asarray(timestamps, dtype=float)
        else:
            cache_keys = (
                np.floor(
                    np.asarray(timestamps, dtype=float)
                    / self.field_line_cadence_seconds
                )
                * self.field_line_cadence_seconds
            )
        for cache_key in np.unique(cache_keys):
            cache_key = float(cache_key)
            if cache_key in self._neutral_sheet_cache:
                continue
            weather = self._evaluate_weather(np.array([cache_key]))
            self._neutral_sheet_cache[cache_key] = T96NeutralSheet(
                cache_key,
                float(weather['Pressure'][0]),
                float(weather['dst'][0]),
                float(weather['BY_GSM'][0]),
                float(weather['BZ_GSM'][0]),
                half_thickness=self.neutral_sheet_half_thickness,
            )
        return [self._neutral_sheet_cache[float(key)] for key in cache_keys]

    def _build(self, times, *, average_space_weather):
        dates, timestamps = self._normalize_times(times)
        positions = np.column_stack(
            [interpolator(timestamps) for interpolator in self._position_interpolators]
        )
        rotations = self._rotation_interpolator(timestamps)
        satellite_positions = []
        satellite_trajectories = []
        trajectory_start = (
            timestamps.min() - self.satellite_trajectory_padding_seconds
        )
        trajectory_stop = (
            timestamps.max() + self.satellite_trajectory_padding_seconds
        )
        trajectory_mask = (
            (self.orbit_timestamps >= trajectory_start)
            & (self.orbit_timestamps <= trajectory_stop)
        )
        for orbit in self.satellite_orbits.values():
            evaluated = np.column_stack(
                [
                    np.interp(timestamps, self.orbit_timestamps, orbit[:, axis])
                    for axis in range(3)
                ]
            )
            satellite_positions.append(evaluated.transpose(1, 0))
            satellite_trajectories.append(
                orbit[trajectory_mask].transpose(1, 0)
            )
        bowshock, magnetopause = self._boundaries(
            timestamps,
            average=average_space_weather,
        )
        field_lines = self._field_lines(
            timestamps,
            exact=len(timestamps) == 1,
        )
        neutral_sheets = self._neutral_sheets(
            timestamps,
            exact=len(timestamps) == 1,
        )
        return scene(
            dates.to_list(),
            (positions, rotations),
            satellite_positions,
            bowshock,
            magnetopause,
            field_lines=field_lines,
            neutral_sheets=neutral_sheets,
            sat_trajectories=satellite_trajectories,
            satellite_names=list(self.satellite_orbits),
            include_reference_planes=self.include_reference_planes,
            satellite_marker_size=self.satellite_marker_size,
        )

    def make_scene(self, time):
        """Build one instantaneous scene entirely from preloaded data."""
        dates, _ = self._normalize_times(time)
        if len(dates) != 1:
            raise ValueError('make_scene accepts exactly one time')
        return self._build(dates, average_space_weather=False)

    def make_timerange(self, trange=None, *, max_frames=48):
        """Build an animated scene using vectorized evaluation of this cache."""
        if not isinstance(max_frames, int) or max_frames < 2:
            raise ValueError('max_frames must be an integer of at least 2')
        timestamps = self.pose_timestamps
        if trange is not None:
            times = common_utils.Time(trange)
            if times.isscalar or len(times.datetime) != 2:
                raise ValueError("trange must contain exactly two time values")
            start, stop = times.to_timestamp()
            timestamps = timestamps[(timestamps >= start) & (timestamps <= stop)]
        if len(timestamps) < 2:
            raise ValueError('fewer than two pose samples in requested timerange')
        if len(timestamps) > max_frames:
            indices = np.unique(
                np.linspace(0, len(timestamps) - 1, max_frames, dtype=int)
            )
            timestamps = timestamps[indices]
        dates = pd.to_datetime(timestamps, unit='s', utc=True)
        return self._build(dates, average_space_weather=True)


def prepare_scene(trange, **kwargs):
    """Preload a reusable :class:`SceneContext` for a time range."""
    return SceneContext.load(trange, **kwargs)


def make_scene(
    time,
    *,
    context=None,
    satellites=('THA', 'ATMS1', 'ATMS2', 'THD', 'THE', 'MMS1'),
    use_space_weather=True,
    include_field_lines=True,
    include_neutral_sheet=True,
    include_reference_planes=True,
    field_line_nlat=21,
    field_line_footprint_latitudes=(70.0, 82.0),
    field_line_rlim=80,
    neutral_sheet_half_thickness=3.0,
    satellite_trajectory_padding='1h',
    satellite_marker_size=6,
    update_ephemeris=False,
):
    """Build one scene, automatically loading its UTC day when necessary.

    Repeated calls should pass a preloaded ``context`` to avoid repeated I/O.
    """
    if context is None:
        normalized = common_utils.Time(time).datetime
        if not isinstance(normalized, pd.Timestamp):
            raise ValueError('make_scene accepts exactly one time')
        normalized = (
            normalized.tz_localize('UTC')
            if normalized.tzinfo is None
            else normalized.tz_convert('UTC')
        )
        start = normalized.floor('D')
        stop = start + pd.Timedelta(days=1) - pd.Timedelta(nanoseconds=1)
        context = prepare_scene(
            [start, stop],
            satellites=satellites,
            use_space_weather=use_space_weather,
            include_field_lines=include_field_lines,
            include_neutral_sheet=include_neutral_sheet,
            include_reference_planes=include_reference_planes,
            field_line_nlat=field_line_nlat,
            field_line_footprint_latitudes=field_line_footprint_latitudes,
            field_line_rlim=field_line_rlim,
            neutral_sheet_half_thickness=neutral_sheet_half_thickness,
            satellite_trajectory_padding=satellite_trajectory_padding,
            satellite_marker_size=satellite_marker_size,
            update_ephemeris=update_ephemeris,
        )
    return context.make_scene(time)


def _output_one_day(
    day,
    *,
    output_dir,
    max_frames,
    use_average_sw_params,
    include_field_lines,
    include_neutral_sheet,
    include_reference_planes,
    field_line_nlat,
    field_line_footprint_latitudes,
    field_line_rlim,
    neutral_sheet_half_thickness,
    field_line_cadence,
    satellite_trajectory_padding,
    satellite_marker_size,
    update_ephemeris,
    include_plotlyjs,
    overwrite,
):
    destination = Path(output_dir) / day / '3d.html'
    if destination.exists() and not overwrite:
        return destination
    start = pd.Timestamp(day, tz='UTC')
    stop = start + pd.Timedelta(days=1) - pd.Timedelta(nanoseconds=1)
    context = prepare_scene(
        [start, stop],
        use_space_weather=use_average_sw_params,
        include_field_lines=include_field_lines,
        include_neutral_sheet=include_neutral_sheet,
        include_reference_planes=include_reference_planes,
        field_line_nlat=field_line_nlat,
        field_line_footprint_latitudes=field_line_footprint_latitudes,
        field_line_rlim=field_line_rlim,
        neutral_sheet_half_thickness=neutral_sheet_half_thickness,
        field_line_cadence=field_line_cadence,
        satellite_trajectory_padding=satellite_trajectory_padding,
        satellite_marker_size=satellite_marker_size,
        update_ephemeris=update_ephemeris,
    )
    figure = context.make_timerange(max_frames=max_frames)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name('.3d.html.tmp')
    try:
        figure.write_html(
            temporary,
            include_plotlyjs=include_plotlyjs,
            full_html=True,
            auto_play=False,
        )
        temporary.replace(destination)
    finally:
        temporary.unlink(missing_ok=True)
    return destination


def dumpGeo(
    *,
    daterange=None,
    dates=None,
    output_dir=None,
    max_frames=48,
    use_average_sw_params=True,
    include_field_lines=True,
    include_neutral_sheet=True,
    include_reference_planes=True,
    field_line_nlat=21,
    field_line_footprint_latitudes=(70.0, 82.0),
    field_line_rlim=80,
    neutral_sheet_half_thickness=3.0,
    field_line_cadence='1h',
    satellite_trajectory_padding='1h',
    satellite_marker_size=6,
    update_ephemeris=False,
    include_plotlyjs=True,
    overwrite=True,
    workers=1,
    verbose=False,
):
    """Generate standalone daily Plotly 3-D observation geometry files.

    With neither ``daterange`` nor ``dates``, every local ``*_pdf.pqt`` day is
    processed. Outputs are written as ``<output_dir>/YYYYMMDD/3d.html``. The
    default output directory is the configured ``GENA_RES`` directory.
    """
    pqt_dir = Path(get_configured_paths()['pqt_base'])
    output_dir = resolve_output_path(output_dir if output_dir is not None else "")
    if not isinstance(max_frames, int) or max_frames < 2:
        raise ValueError('max_frames must be an integer of at least 2')
    if not isinstance(workers, int) or workers < 1:
        raise ValueError('workers must be a positive integer')
    if daterange is None and dates is None:
        days = sorted({
            pd.Timestamp(path.name[:8]).strftime('%Y%m%d')
            for path in pqt_dir.glob('*_pdf.pqt')
        })
    else:
        if daterange is not None and dates is not None:
            raise ValueError("provide exactly one of daterange or dates")
        time = common_utils.Time(daterange if daterange is not None else dates)
        if daterange is not None:
            if time.isscalar or len(time.datetime) != 2:
                raise ValueError("day range must contain exactly two time values")
            start, stop = time.datetime
            if start > stop:
                raise ValueError("day range start must not be later than end")
            day_times = pd.date_range(start.floor("D"), stop.floor("D"), freq="D")
        else:
            native = pd.DatetimeIndex([time.datetime]) if time.isscalar else time.datetime
            if len(native) == 0:
                raise ValueError("dates must contain at least one date")
            day_times = native.floor("D").drop_duplicates().sort_values()
        days = [day.strftime('%Y%m%d') for day in day_times]
    kwargs = {
        'output_dir': output_dir,
        'max_frames': max_frames,
        'use_average_sw_params': use_average_sw_params,
        'include_field_lines': include_field_lines,
        'include_neutral_sheet': include_neutral_sheet,
        'include_reference_planes': include_reference_planes,
        'field_line_nlat': field_line_nlat,
        'field_line_footprint_latitudes': field_line_footprint_latitudes,
        'field_line_rlim': field_line_rlim,
        'neutral_sheet_half_thickness': neutral_sheet_half_thickness,
        'field_line_cadence': field_line_cadence,
        'satellite_trajectory_padding': satellite_trajectory_padding,
        'satellite_marker_size': satellite_marker_size,
        'update_ephemeris': update_ephemeris,
        'include_plotlyjs': include_plotlyjs,
        'overwrite': overwrite,
    }
    generated = {}
    failures = {}
    if workers == 1:
        for index, day in enumerate(days, start=1):
            try:
                generated[day] = _output_one_day(day, **kwargs)
                if verbose:
                    print(f'[{index}/{len(days)}] generated {generated[day]}', flush=True)
            except Exception as exc:
                failures[day] = exc
                warnings.warn(f'failed to generate {day}/3d.html: {exc}')
    else:
        with ProcessPoolExecutor(max_workers=workers) as executor:
            futures = {
                executor.submit(_output_one_day, day, **kwargs): day
                for day in days
            }
            for index, future in enumerate(as_completed(futures), start=1):
                day = futures[future]
                try:
                    generated[day] = future.result()
                    if verbose:
                        print(f'[{index}/{len(days)}] generated {generated[day]}', flush=True)
                except Exception as exc:
                    failures[day] = exc
                    warnings.warn(f'failed to generate {day}/3d.html: {exc}')
    if failures:
        failed_days = ', '.join(sorted(failures))
        raise RuntimeError(f'3D output failed for {len(failures)} day(s): {failed_days}')
    return {day: generated[day] for day in sorted(generated)}


# Backward-compatible name for callers that have not migrated to ``dumpGeo``.
output_3d = dumpGeo
