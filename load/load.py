"""Public loading facade for GENA and external mission products."""

from __future__ import annotations

import numpy as np
import pandas as pd

from common.time import Time
from libs import utils as _common_utils

from load.loaders.gena import load_gena
from load.loaders.indices import load_indices
from load.loaders.spedas import load_spedas
from load.loaders.common import is_multi_interval
from load.api.spedas import read_spedas
from load.api.indices import read_indices


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
    interpolators = _common_utils.f_interp_vars(normalized)
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
    if is_multi_interval(t):
        return np.concatenate(
            [get_themis_loc(interval, probe=probe, coord=coord) for interval in t],
            axis=0,
        )
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
    if is_multi_interval(t):
        return np.concatenate(
            [get_mms_loc(interval, probe=probe, coord=coord, download=download, **kwargs) for interval in t],
            axis=0,
        )
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
    if is_multi_interval(t):
        return np.concatenate(
            [get_goes_loc(interval, satellite, coord=coord, download=download, **kwargs) for interval in t],
            axis=0,
        )
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
    if is_multi_interval(t):
        return np.concatenate(
            [get_queqiao_loc(interval, satellite=satellite) for interval in t],
            axis=0,
        )
    times, dates, _ = _position_times(t)
    frame = load_gena(
        trange=_padded_position_range(dates),
        product="pdf",
        columns=["MJ2000_Star_X", "MJ2000_Star_Y", "MJ2000_Star_Z"],
    )
    result = frame.pos_at(times)
    return result[0] if times.isscalar else result


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
