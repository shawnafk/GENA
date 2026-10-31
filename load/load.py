"""Public loading facade for GENA and external mission products."""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd

from common.time import Time
from common.interp import f_interp_vars

from load.loaders.gena import load_gena
from load.loaders.indices import load_indices
from load.loaders.spedas import load_spedas
from load.api.spedas import read_spedas
from load.api.indices import read_indices


logger = logging.getLogger(__name__)


def load_sw(t, *, trange_pad="1h", **kwargs):
    """Load and interpolate common solar-wind and geomagnetic keys.

    Each OMNI key falls back independently: an unusable variable keeps its
    documented constant while the others retain real values. Every fallback is
    logged, so synthetic inputs are never silently mistaken for measurements.
    """
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
    defaults = {
        "Pressure": 1.0,
        "BY_GSM": 0.0,
        "BZ_GSM": 0.0,
        "Beta": 1.0,
        "Mgs_mach_num": 5.0,
        "SYM_H": 0.0,
    }
    values = {
        name: np.full_like(float_t, value, dtype=np.float32)
        for name, value in defaults.items()
    }
    omni_error = None
    try:
        omni_kwargs = dict(kwargs)
        omni_kwargs.pop("variables", None)
        omni_kwargs.pop("return_type", None)
        omni_data = load_spedas(
            "omni", product="hro", trange=trange,
            variables=namelist, return_type="dict", **omni_kwargs,
        )
    except Exception as exc:  # loader boundary: fall back per variable below
        omni_data = {}
        omni_error = exc
    resolved = []
    for name in namelist:
        record = omni_data.get(name)
        if record is None:
            continue
        try:
            interpolated = f_interp_vars(record)(float_t)
        except Exception as exc:  # one bad variable must not hide the others
            omni_error = omni_error or exc
            continue
        if interpolated is None or not np.any(np.isfinite(interpolated)):
            continue
        values[name] = np.asarray(interpolated, dtype=float)
        resolved.append(name)
    missing = [name for name in namelist if name not in resolved]
    if missing:
        logger.warning(
            "load_sw: fallback constants used for %s (last error: %r)",
            ", ".join(missing), omni_error,
        )

    try:
        index_kwargs = dict(kwargs)
        index_kwargs.pop("return_type", None)
        dst_vars = load_indices(
            trange=trange, product="dst", return_type="dict", **index_kwargs,
        )
        interpolator = f_interp_vars(dst_vars)
        if interpolator is None:
            raise ValueError("dst variable has an unsupported sample layout")
        values["dst"] = interpolator(float_t)
    except Exception as exc:
        logger.warning("load_sw: fallback dst=0 used (%r)", exc)
        values["dst"] = np.full_like(float_t, 0, dtype=np.float32)

    if scalar:
        return {key: np.asarray(value) for key, value in values.items()}
    return values


#: What ``product``, ``probe`` and ``key`` mean for each mission. ``key`` is
#: the "what part of this product do I want" selector: a variable name for
#: SPEDAS missions, a column for GENA, and the index name for indices.
_MISSION_ROLES = {
    "gena": {"probe": False, "key": "columns"},
    "themis": {"probe": True, "key": "variables"},
    "mms": {"probe": True, "key": "variables"},
    "omni": {"probe": False, "key": "variables"},
    "indices": {"probe": False, "key": "product"},
}
_MISSION_ALIASES = {"index": "indices", "kyoto": "indices"}


def load(mission, trange, product=None, *, probe=None, key=None, **kwargs):
    """Load one mission product through its source loader.

    ``mission`` selects the role of the other arguments:

    - ``gena``: ``product`` is a Parquet suffix such as ``l2s_a``; there is no
      separate probe (it is encoded in the suffix); ``key`` selects columns
      (``key="Energy"``) and is equivalent to ``columns=[...]``.
    - ``themis`` / ``mms``: ``product`` is the instrument product (``esa``,
      ``fpi``, ...), ``probe`` is the spacecraft id (``b``, ``1``, ...), and
      ``key`` selects variables (``key="peir_flux"``) instead of the product's
      default list; it is equivalent to ``variables=[...]``.
    - ``omni``: ``product`` defaults to ``hro``, there is no probe, and
      ``key`` selects variables.
    - ``indices``: ``product`` (or ``key``) is the index name (``dst``,
      ``ae``, ``kp``, ...); there is no probe or instrument level.

    Use ``key`` for a short selection and ``variables``/``columns`` when the
    explicit loader-style name is clearer; passing both is rejected.
    """
    name = _MISSION_ALIASES.get(str(mission).lower(), str(mission).lower())
    role = _MISSION_ROLES.get(name)
    if role is None:
        choices = ", ".join(sorted(_MISSION_ROLES))
        raise ValueError(f"mission must be one of {choices}")
    if not role["probe"] and probe is not None:
        raise ValueError(f"{name} does not accept probe")
    if key is not None:
        if role["key"] == "variables":
            if "variables" in kwargs:
                raise ValueError("pass either key or variables, not both")
            kwargs["variables"] = [key] if isinstance(key, str) else list(key)
        elif role["key"] == "columns":
            if "columns" in kwargs:
                raise ValueError("pass either key or columns, not both")
            kwargs["columns"] = [key] if isinstance(key, str) else list(key)
        elif product is None:
            product = key
    if name == "gena":
        return load_gena(trange=trange, product=product, **kwargs)
    if name in {"themis", "mms", "omni"}:
        return load_spedas(name, trange, product, probe=probe, **kwargs)
    return load_indices(trange, product=product, **kwargs)


__all__ = [
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


def _get_themis_loc(t, probe="b", coord="gsm"):
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


def _get_mms_loc(t, probe="1", coord="gsm", *, download=None, **kwargs):
    """Interpolate an MMS MEC position at arbitrary UTC times.

    MMS MEC position variables are named ``mmsN_mec_r_COORD``. The returned
    values retain the MEC product's native unit, normally kilometres.
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


def _get_goes_loc(t, satellite, coord="gsm", *, download=None, **kwargs):
    """Interpolate a GOES orbit position at arbitrary UTC times.

    ``satellite`` is the GOES spacecraft number. Orbit vectors are loaded
    from pySPEDAS and retain the ephemeris product's native unit, normally
    kilometres.
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


def _get_queqiao_loc(t, satellite=None):
    """Interpolate the Queqiao position recorded in GENA PDF engineering data.

    The PDF position columns are Moon-relative MJ2000 coordinates in metres.
    ``satellite`` is accepted for symmetry with the other mission helpers and
    is currently ignored because the configured PDF contains one Queqiao
    spacecraft trajectory.
    """
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

    ``trange`` is a scalar time or a one-dimensional sequence of UTC times.
    The result is one ``(x, y, z)`` vector for a scalar and an ``(N, 3)``
    array otherwise. THEMIS ids are ``a``-``e``, MMS
    ids are ``1``-``4``, and GOES ids are spacecraft numbers. Queqiao uses the
    local GENA PDF orbit.
    """
    key = str(mission).strip().lower()
    if key == "themis":
        return _get_themis_loc(trange, probe="a" if id is None else id, coord=coord)
    if key == "mms":
        return _get_mms_loc(trange, probe="1" if id is None else id, coord=coord, **kwargs)
    if key in {"goes", "goes-r"}:
        if id is None:
            raise ValueError("id is required for GOES")
        return _get_goes_loc(trange, id, coord=coord, **kwargs)
    if key in {"queqiao", "qq", "gena"}:
        return _get_queqiao_loc(trange, id)
    raise ValueError("mission must be 'themis', 'mms', 'goes', or 'queqiao'")
