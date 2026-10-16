"""Higher-level space-environment data products."""

from __future__ import annotations

import numpy as np
import pandas as pd

from libs import utils as _common_utils
from load import load_indices, load_spedas

__all__ = ["sw_keys"]


def sw_keys(t, *, trange_pad="1h", verbose=False, **kwargs):
    """Interpolate common solar-wind/context keys at time(s).

    This is a collective application-level loader: it combines OMNI and Dst
    products, normalizes flexible time inputs, and returns arrays evaluated at
    the requested times.  Lower-level loading uses ``load_spedas`` and ``load_indices``.
    """
    times = _common_utils.Time(t)
    scalar = times.isscalar
    datetimes = pd.DatetimeIndex([times.datetime]) if scalar else times.datetime
    if datetimes.empty:
        return {}

    float_t = np.atleast_1d(times.to_timestamp())
    pad = pd.Timedelta(trange_pad)
    trange = [
        (datetimes.min() - pad).isoformat(),
        (datetimes.max() + pad).isoformat(),
    ]

    namelist = ["Pressure", "BY_GSM", "BZ_GSM", "Beta", "Mgs_mach_num", "SYM_H"]
    values = {}
    try:
        omni_kwargs = dict(kwargs)
        omni_kwargs.pop("variables", None)
        omni_kwargs.pop("return_type", None)
        omni_data = load_spedas(
            "omni", product="hro",
            trange=trange, variables=namelist, return_type="dict", **omni_kwargs
        )
        for name in namelist:
            values[name] = _common_utils.f_interp_vars(omni_data[name])(float_t)
        if verbose:
            print("load solar wind parameters successfully")
    except Exception as exc:
        if verbose:
            print("Warning: cannot load omni data, set all default.", exc)
        values["Pressure"] = np.full_like(float_t, 1, dtype=np.float32)
        values["BY_GSM"] = np.full_like(float_t, 0, dtype=np.float32)
        values["BZ_GSM"] = np.full_like(float_t, 0, dtype=np.float32)
        values["Beta"] = np.full_like(float_t, 1, dtype=np.float32)
        values["Mgs_mach_num"] = np.full_like(float_t, 5, dtype=np.float32)
        values["SYM_H"] = np.full_like(float_t, 0, dtype=np.float32)

    try:
        index_kwargs = dict(kwargs)
        index_kwargs.pop("return_type", None)
        dst_vars = load_indices(
            trange=trange, product="dst", return_type="dict", **index_kwargs
        )
        values["dst"] = _common_utils.f_interp_vars(dst_vars)(float_t)
    except Exception as exc:
        if verbose:
            print("Warning: cannot load dst data, set all default.", exc)
        values["dst"] = np.full_like(float_t, 0, dtype=np.float32)

    if scalar:
        return {k: np.asarray(v) for k, v in values.items()}
    return values
