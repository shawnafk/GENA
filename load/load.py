"""Public loading facade for GENA and external mission products."""

from __future__ import annotations

import numpy as np
import pandas as pd
from functools import wraps
from collections.abc import Mapping

from common.time import Time
from common.interp import f_interp_vars
from common.time import Time

from load.loaders.gena import load_gena as _load_gena
from load.loaders.indices import load_indices as _load_indices
from load.loaders.spedas import load_spedas as _load_spedas
from load.api.spedas import read_spedas
from load.api.indices import read_indices
from common.frame import concat
from load.api.public import _merge_variable_dicts


def is_multi_interval(trange):
    """Return whether *trange* is a nested sequence of intervals."""
    return isinstance(trange, (list, tuple)) and bool(trange) and isinstance(trange[0], (list, tuple))


def normalize_interval(trange):
    """Normalize one inclusive interval to two ordered UTC timestamps."""
    dates = Time(trange).datetime
    if isinstance(dates, pd.Timestamp) or len(dates) != 2:
        raise ValueError("trange must contain exactly two time values")
    if dates[0] > dates[1]:
        raise ValueError("trange start must not be later than end")
    return dates


def iter_intervals(trange):
    """Normalize one interval or a nested sequence of intervals."""
    ranges = trange if is_multi_interval(trange) else [trange]
    return [normalize_interval(interval) for interval in ranges]


def _merge_interval_results(results):
    if not results:
        return None
    first = results[0]
    if isinstance(first, Mapping):
        return _merge_variable_dicts(results)
    if isinstance(first, np.ndarray):
        return np.concatenate(results, axis=0)
    if isinstance(first, pd.DataFrame):
        return concat(results)
    return results


def _multi_interval(func):
    """Apply a loader or location function once per requested interval."""
    @wraps(func)
    def wrapped(*args, **kwargs):
        trange = kwargs.get("trange")
        position = None
        if trange is None and len(args) > 1:
            position = 1
            trange = args[position]
        if trange is None or not is_multi_interval(trange):
            return func(*args, **kwargs)
        results = []
        for interval in iter_intervals(trange):
            if position is None:
                call_kwargs = dict(kwargs)
                call_kwargs["trange"] = interval
                results.append(func(*args, **call_kwargs))
            else:
                call_args = list(args)
                call_args[position] = interval
                results.append(func(*call_args, **kwargs))
        return _merge_interval_results(results)
    return wrapped


load_gena = _multi_interval(_load_gena)
load_indices = _multi_interval(_load_indices)
load_spedas = _multi_interval(_load_spedas)


def load_sw(t, *, trange_pad="1h", verbose=False, **kwargs):
    """Load and interpolate common solar-wind and geomagnetic keys."""
    times = Time(t)
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
            "omni", product="hro", trange=trange,
            variables=namelist, return_type="dict", **omni_kwargs,
        )
        for name in namelist:
            values[name] = f_interp_vars(omni_data[name])(float_t)
        if verbose:
            print("load solar wind parameters successfully")
    except Exception as exc:
        if verbose:
            print("Warning: cannot load omni data, set all default.", exc)
        values.update({
            "Pressure": np.full_like(float_t, 1, dtype=np.float32),
            "BY_GSM": np.full_like(float_t, 0, dtype=np.float32),
            "BZ_GSM": np.full_like(float_t, 0, dtype=np.float32),
            "Beta": np.full_like(float_t, 1, dtype=np.float32),
            "Mgs_mach_num": np.full_like(float_t, 5, dtype=np.float32),
            "SYM_H": np.full_like(float_t, 0, dtype=np.float32),
        })

    try:
        index_kwargs = dict(kwargs)
        index_kwargs.pop("return_type", None)
        dst_vars = load_indices(
            trange=trange, product="dst", return_type="dict", **index_kwargs,
        )
        values["dst"] = f_interp_vars(dst_vars)(float_t)
    except Exception as exc:
        if verbose:
            print("Warning: cannot load dst data, set all default.", exc)
        values["dst"] = np.full_like(float_t, 0, dtype=np.float32)

    if scalar:
        return {key: np.asarray(value) for key, value in values.items()}
    return values


@_multi_interval
def load(mission, trange, product=None, *, probe=None, key=None, **kwargs):
    """Load a selected mission product through its source loader."""
    mission = str(mission).lower()
    if mission == "gena":
        return load_gena(trange=trange, product=product, **kwargs)
    if mission in {"themis", "mms", "omni"}:
        return load_spedas(mission, trange, product, probe=probe, **kwargs)
    if mission in {"indices", "index", "kyoto"}:
        return load_indices(trange, product=key or product, **kwargs)
    raise ValueError("mission must be one of 'gena', 'themis', 'mms', 'omni', or 'indices'")


__all__ = [
    "get_goes_loc",
    "get_mms_loc",
    "get_queqiao_loc",
    "get_themis_loc",
    "load",
    "load_sw",
    "load_pos",
    "read_spedas",
    "read_indices",
]


# Explicit convenience operations using the public loaders.
def _position_times(value):
    """Normalize position evaluation times and return dates plus Unix seconds."""
    times = Time(value)
    dates = pd.DatetimeIndex([times.datetime]) if times.isscalar else times.datetime
    if dates.empty:
        raise ValueError("at least one probe-position time is required")
    return times, dates, np.asarray(times.to_timestamp(), dtype=float)


def _padded_position_range(dates):
    """Return a one-hour-padded ISO range for a position backend request."""
    pad = pd.Timedelta("1h")
    return [(dates.min() - pad).isoformat(), (dates.max() + pad).isoformat()]


def _interpolate_position(variable, evaluation):
    """Interpolate one backend vector variable at Unix-second timestamps."""
    values = np.asarray(variable["y"])
    if values.ndim != 2 or values.shape[-1] != 3:
        raise ValueError(f"position variable must have shape (N, 3), got {values.shape}")
    normalized = dict(variable)
    normalized["y"] = values
    interpolators = f_interp_vars(normalized)
    if interpolators is None:
        raise ValueError("position variable has an unsupported sample layout")
    return np.column_stack([interpolator(evaluation) for interpolator in interpolators])


def _backend_position(value, variable_name, reader):
    """Load and interpolate a named backend position variable."""
    times, dates, evaluation = _position_times(value)
    raw = reader(_padded_position_range(dates), variable_name)
    if variable_name not in raw:
        available = ", ".join(sorted(raw)[:8])
        raise KeyError(f"position variable {variable_name!r} was not returned; available: {available}")
    result = _interpolate_position(raw[variable_name], evaluation)
    return result[0] if times.isscalar else result


def get_themis_loc(t, probe="b", coord="gsm"):
    """Interpolate a THEMIS probe position in the common ``(N, 3)`` layout.

    The source is the THEMIS state product. ``coord`` is the state-variable
    coordinate suffix, such as ``gsm`` or ``gse``; values retain the source
    product's native position unit. A scalar time returns one ``(3,)`` vector.
    """
    probe = str(probe).lower()
    if probe not in {"a", "b", "c", "d", "e"}:
        raise ValueError("THEMIS probe must be one of 'a', 'b', 'c', 'd', or 'e'")
    variable_name = f"th{probe}_pos_{str(coord).lower()}"

    def reader(trange, name):
        return load_spedas(
            "themis", trange, "state", variables=[name], probe=probe,
            level="l1", return_type="dict",
        )

    return _backend_position(t, variable_name, reader)


def get_mms_loc(t, probe="1", coord="gsm", *, download=None, **kwargs):
    """Interpolate an MMS MEC position at arbitrary UTC times.

    MMS MEC position variables are named ``mmsN_mec_r_COORD``. The returned
    values retain the MEC product's native unit, normally kilometres. Nested
    ``[[start, end], ...]`` input is evaluated interval by interval and
    concatenated in the input order.
    """
    probe = str(probe)
    if probe not in {"1", "2", "3", "4"}:
        raise ValueError("MMS probe must be one of '1', '2', '3', or '4'")
    coord = str(coord).lower()
    if coord not in {"gsm", "gse", "sm", "geo", "gei", "dsc", "dbcs"}:
        raise ValueError("unsupported MMS position coordinate")
    variable_name = f"mms{probe}_mec_r_{coord}"

    def reader(trange, name):
        from pyspedas.projects import mms

        extra = dict(kwargs)
        call = dict(
            trange=trange,
            probe=probe,
            data_rate=extra.pop("data_rate", "srvy"),
            level=extra.pop("level", "l2"),
            datatype=extra.pop("datatype", "epht89q"),
            varnames=[name],
            notplot=True,
            time_clip=True,
            no_update=download in (False, "offline", "local"),
            spdf=extra.pop("spdf", False),
        )
        call.update(extra)
        return mms.mms_load_mec(**call) or {}

    return _backend_position(t, variable_name, reader)


def get_goes_loc(t, satellite, coord="gsm", *, download=None, **kwargs):
    """Interpolate a GOES orbit position at arbitrary UTC times.

    ``satellite`` is the GOES spacecraft number. Orbit vectors are loaded
    from pySPEDAS and retain the ephemeris product's native unit, normally
    kilometres. Nested ``[[start, end], ...]`` input is evaluated interval by
    interval and concatenated in the input order.
    """
    try:
        satellite = str(int(satellite))
    except (TypeError, ValueError) as exc:
        raise ValueError("GOES satellite id must be an integer") from exc
    if int(satellite) < 1:
        raise ValueError("GOES satellite id must be positive")
    coord = str(coord).upper()
    if coord not in {"GSM", "GSE", "SM", "GEO"}:
        raise ValueError("GOES orbit coordinate must be GSM, GSE, SM, or GEO")
    variable_name = f"g{satellite}_orbit_XYZ_{coord}"

    def reader(trange, name):
        from pyspedas.projects import goes

        call = dict(
            trange=trange,
            probe=satellite,
            notplot=True,
            time_clip=True,
            no_update=download in (False, "offline", "local"),
            force_download=download == "refresh",
        )
        call.update(kwargs)
        return goes.load_orbit(**call) or {}

    return _backend_position(t, variable_name, reader)


def get_queqiao_loc(t, satellite=None):
    """Interpolate the Queqiao position recorded in GENA PDF engineering data.

    The PDF position columns are Moon-relative MJ2000 coordinates in metres.
    ``satellite`` is accepted for symmetry with the other mission helpers and
    is currently ignored because the configured PDF contains one Queqiao
    spacecraft trajectory. Nested ``[[start, end], ...]`` input is evaluated
    interval by interval and concatenated in the input order.
    """
    times, dates, _ = _position_times(t)
    frame = load_gena(
        trange=_padded_position_range(dates),
        product="pdf",
        columns=["MJ2000_Star_X", "MJ2000_Star_Y", "MJ2000_Star_Z"],
    )
    result = frame.pos_at(times)
    return result[0] if times.isscalar else result


@_multi_interval
def load_pos(mission, trange, id=None, *, coord="gsm", **kwargs):
    """Return interpolated spacecraft positions for a supported mission.

    ``trange`` may be a scalar, a one-dimensional sequence of UTC times, or a
    nested sequence of two-value intervals. Nested intervals are evaluated and
    concatenated in input order. The result is one ``(x, y, z)`` vector for a
    scalar and an ``(N, 3)`` array otherwise. THEMIS ids are ``a``-``e``, MMS
    ids are ``1``-``4``, and GOES ids are spacecraft numbers. Queqiao uses the
    local GENA PDF orbit.
    """
    key = str(mission).strip().lower()
    if key == "themis":
        return get_themis_loc(trange, probe="a" if id is None else id, coord=coord)
    if key == "mms":
        return get_mms_loc(trange, probe="1" if id is None else id, coord=coord, **kwargs)
    if key in {"goes", "goes-r"}:
        if id is None:
            raise ValueError("id is required for GOES")
        return get_goes_loc(trange, id, coord=coord, **kwargs)
    if key in {"queqiao", "qq", "gena"}:
        return get_queqiao_loc(trange, id)
    raise ValueError("mission must be 'themis', 'mms', 'goes', or 'queqiao'")
