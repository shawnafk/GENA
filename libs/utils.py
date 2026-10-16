from common.time import Time
import numpy as np
import pandas as pd
from scipy import interpolate
from scipy.spatial.transform import Slerp
from scipy.spatial.transform import Rotation as Rot
from functools import wraps

"""
Shared utility decorators and formatters for the GENA codebase.
Includes:
    - Decorators for timing, logging, and repeating function execution.
"""

import functools
import time
import datetime


# ── Decorators ────────────────────────────────────────────────────────
def timer(func):
    """Print function execution time in seconds."""
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        start_time = time.perf_counter()
        result = func(*args, **kwargs)
        end_time = time.perf_counter()
        run_time = end_time - start_time
        print(f"Function '{func.__name__}' executed,耗时: {run_time:.4f} s")
        return result
    return wrapper


def logger(func):
    """Log function entry and exit."""
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        print(f"Preparing to execute: {func.__name__}")
        result = func(*args, **kwargs)
        print(f"Function {func.__name__} executed.")
        return result
    return wrapper

def repeat(num_times):
    """一个带参数的装饰器，用于重复执行函数指定的次数。"""
    def decorator_repeat(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            for _ in range(num_times):
                result = func(*args, **kwargs)
            return result
        return wrapper
    return decorator_repeat




# ── Array Utilities ───────────────────────────────────────────────────
def mid_ary(ary):
    """Return the midpoint of consecutive array elements."""
    return (ary[:-1] + ary[1:]) / 2




# Alias for scipy 1-D interpolation (used throughout the codebase)
i1d = interpolate.interp1d

def f_interp_vars(dat):
    """
    1D or 3D interpolation of data variables.

    Parameters
    ----------
    dat : dict
        Dictionary with ``'x'`` (timestamps) and ``'y'`` (values) keys.

    Returns
    -------
    callable or tuple of callable
        - 1-D data -> single ``interp1d`` callable
        - 3-D (N, 3) data -> ``(fx, fy, fz)`` tuple of ``interp1d`` callables
    """
    t_intp = np.asarray(Time(dat['x']).to_timestamp(), dtype=float)
    interp_kw = dict(bounds_error=False, fill_value=np.nan)
    if 1 == dat['y'].ndim:
        idx = ~np.isnan(dat['y'])
        y_valid = dat['y'][idx]
        # ns to s
        t_valid = t_intp[idx]
        order = np.argsort(t_valid)
        f = i1d(t_valid[order], y_valid[order], **interp_kw)
        return f
    if 2 == dat['y'].ndim and 3 == dat['y'].shape[-1]:
        idx_x = ~np.isnan(dat['y'][:, 0])
        idx_y = ~np.isnan(dat['y'][:, 1])
        idx_z = ~np.isnan(dat['y'][:, 2])
        idx = idx_x & idx_y & idx_z
        y_valid = dat['y'][idx]
        # ns to s
        t_valid = t_intp[idx]
        order = np.argsort(t_valid)
        t_valid, y_valid = t_valid[order], y_valid[order]
        fx = i1d(t_valid, y_valid[:, 0], **interp_kw)
        fy = i1d(t_valid, y_valid[:, 1], **interp_kw)
        fz = i1d(t_valid, y_valid[:, 2], **interp_kw)
        return fx, fy, fz


# ── Rotation / Quaternion Interpolation ──────────────────────────────
def i_rot(t, rot, t_eval):
    """Interpolate rotation via SLERP with boundary extension.

    Input timestamps are sorted and deduplicated before interpolation.
    If *t_eval* extends beyond the input range the edge rotation value
    is held constant.

    Parameters
    ----------
    t : array-like
        Source timestamps.
    rot : scipy.spatial.transform.Rotation
        Source rotation object.
    t_eval : array-like
        Timestamps at which to evaluate the interpolated rotation.

    Returns
    -------
    scipy.spatial.transform.Rotation
        Interpolated rotation at *t_eval*.
    """
    t_array = np.array(t)
    quat_array = rot.as_quat()

    # 1. Sort by time
    sorted_idx = np.argsort(t_array)
    t_sorted = t_array[sorted_idx]
    quat_sorted = quat_array[sorted_idx]

    # 2. Deduplicate (keep first occurrence)
    _, unique_idx = np.unique(t_sorted, return_index=True)
    t_clean = t_sorted[unique_idx]
    quat_clean = quat_sorted[unique_idx]

    # 3. Guard: insufficient data
    if len(t_clean) < 2:
        return Rot.concatenate([
            Rot.from_quat([quat_clean[0]]) for _ in t_eval
        ])

    rot_clean = Rot.from_quat(quat_clean)

    # 4. Extend left boundary
    if t_eval[0] < t_clean[0]:
        t_clean = np.insert(t_clean, 0, t_eval[0])
        rot_clean = Rot.concatenate([
            Rot.from_quat([quat_clean[0]]), rot_clean
        ])

    # 5. Extend right boundary
    if t_eval[-1] > t_clean[-1]:
        t_clean = np.append(t_clean, t_eval[-1])
        rot_clean = Rot.concatenate([
            rot_clean, Rot.from_quat([quat_clean[-1]])
        ])

    # 6. SLERP interpolation
    slerp = Slerp(t_clean, rot_clean)
    return slerp(t_eval)

#interpolation always take timestamp
def i_pos(t,pos,teval,extrapolate=True):
    """Interpolate position array columns at new evaluation times.

    Parameters
    ----------
    t : array-like, shape ``(M,)``
        Source timestamps.
    pos : np.ndarray, shape ``(M, D)``
        Source positions.
    t_eval : array-like, shape ``(N,)``
        Target timestamps.

    Returns
    -------
    np.ndarray, shape ``(N, D)``
        Interpolated positions at *t_eval*.
    """
    if extrapolate:
        f_x = i1d(t, pos[:,0],bounds_error=False, fill_value="extrapolate")
        f_y = i1d(t, pos[:,1],bounds_error=False, fill_value="extrapolate")
        f_z = i1d(t, pos[:,2],bounds_error=False, fill_value="extrapolate")
    else:
        f_x = i1d(t, pos[:,0])
        f_y = i1d(t, pos[:,1])
        f_z = i1d(t, pos[:,2])
    return np.column_stack([f_x(teval), f_y(teval), f_z(teval)])

def clean_nan(x,y):
    valid_mask = ~np.isnan(x) & ~np.isnan(y)
    x_clean = x[valid_mask]
    y_clean = y[valid_mask]
    return x_clean, y_clean

from functools import wraps

def is_multi_layer(trange):
    return isinstance(trange, (list, tuple)) and len(trange) and isinstance(trange[0], (list, tuple))

def multi_trange_loader(mode: str):
    """Add transparent multi-interval support to a loader.

    A normal ``[start, end]`` range is passed through unchanged.  For input of
    the form ``[[start1, end1], [start2, end2], ...]``, each interval is loaded
    separately and then merged:

    ``mode='df'``
        Concatenate DataFrames with ``ignore_index=True``.
    ``mode='dict'``
        Merge ``{'x', 'y'}`` variable dictionaries; output ``x`` is a
        ``DatetimeIndex`` and ``y`` (plus optional ``v``) is an ndarray.
    """
    if mode not in {"df", "dict"}:
        raise ValueError(f"mode only supports 'df' or 'dict'; got {mode!r}")

    def decorator(func):
        @wraps(func)
        def wrapper(trange, *args, **kwargs):
            if not is_multi_layer(trange):
                return func(trange, *args, **kwargs)

            chunk_list = []
            for sub_trange in trange:
                chunk = func(sub_trange, *args, **kwargs)
                chunk_list.append(chunk)

            if mode == "df":
                concat_df = pd.concat(chunk_list, ignore_index=True)
                return concat_df
            elif mode == "dict":
                merged_data = {}
                for chunk in chunk_list:
                    for var_name, var_dict in chunk.items():
                        if var_name not in merged_data:
                            merged_data[var_name] = {"x": [], "y": []}
                        
                            if "v" in var_dict:
                                merged_data[var_name]["v"] = []
                        merged_data[var_name]["x"].extend(var_dict["x"])
                        merged_data[var_name]["y"].extend(var_dict["y"])
                        if "v" in var_dict:
                            merged_data[var_name]["v"].extend(var_dict["v"])
                for var_name, var_dict in merged_data.items():
                    merged_data[var_name]['x'] =pd.DatetimeIndex(merged_data[var_name]['x'])
                    merged_data[var_name]['y'] = np.array(merged_data[var_name]['y'])
                    if "v" in var_dict:
                        merged_data[var_name]['v'] = np.array(merged_data[var_name]['v'])
                return merged_data
        return wrapper
    return decorator




# ---- useful functions 
def normalize_probe(prob):
    prob = str(prob).strip().lower()
    if prob not in {"a", "b"}:
        raise ValueError("prob must be 'A' or 'B'")
    return prob


