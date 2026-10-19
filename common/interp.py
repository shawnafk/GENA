"""Time-series, vector, and attitude interpolation helpers."""

from __future__ import annotations

import numpy as np
from scipy import interpolate
from scipy.spatial.transform import Slerp
from scipy.spatial.transform import Rotation as Rot

from common.time import Time


i1d = interpolate.interp1d


def mid_ary(ary):
    """Return the midpoint of consecutive array elements."""
    return (ary[:-1] + ary[1:]) / 2


def f_interp_vars(dat):
    """Build scalar or three-component interpolators for a backend variable."""
    t_intp = np.asarray(Time(dat["x"]).to_timestamp(), dtype=float)
    interp_kw = dict(bounds_error=False, fill_value=np.nan)
    values = np.asarray(dat["y"])
    if values.ndim == 1:
        idx = ~np.isnan(values)
        order = np.argsort(t_intp[idx])
        return i1d(t_intp[idx][order], values[idx][order], **interp_kw)
    if values.ndim == 2 and values.shape[-1] == 3:
        idx = np.all(np.isfinite(values), axis=1)
        order = np.argsort(t_intp[idx])
        times, valid = t_intp[idx][order], values[idx][order]
        return tuple(i1d(times, valid[:, axis], **interp_kw) for axis in range(3))
    return None


def i_rot(t, rot, t_eval):
    """Interpolate rotations with SLERP and constant edge extension."""
    t_array = np.asarray(t)
    quat_array = rot.as_quat()
    sorted_idx = np.argsort(t_array)
    t_sorted, quat_sorted = t_array[sorted_idx], quat_array[sorted_idx]
    _, unique_idx = np.unique(t_sorted, return_index=True)
    t_clean, quat_clean = t_sorted[unique_idx], quat_sorted[unique_idx]
    if len(t_clean) < 2:
        return Rot.concatenate([Rot.from_quat([quat_clean[0]]) for _ in t_eval])
    rot_clean = Rot.from_quat(quat_clean)
    if t_eval[0] < t_clean[0]:
        t_clean = np.insert(t_clean, 0, t_eval[0])
        rot_clean = Rot.concatenate([Rot.from_quat([quat_clean[0]]), rot_clean])
    if t_eval[-1] > t_clean[-1]:
        t_clean = np.append(t_clean, t_eval[-1])
        rot_clean = Rot.concatenate([rot_clean, Rot.from_quat([quat_clean[-1]])])
    return Slerp(t_clean, rot_clean)(t_eval)


def i_pos(t, pos, teval, extrapolate=True):
    """Interpolate three-component positions at evaluation timestamps."""
    fill_value = "extrapolate" if extrapolate else np.nan
    functions = [i1d(t, pos[:, axis], bounds_error=not extrapolate, fill_value=fill_value) for axis in range(3)]
    return np.column_stack([function(teval) for function in functions])


def clean_nan(x, y):
    """Remove paired non-finite samples from two arrays."""
    valid_mask = ~np.isnan(x) & ~np.isnan(y)
    return x[valid_mask], y[valid_mask]
