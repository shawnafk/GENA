"""Public loading facade for GENA and external mission products."""

from __future__ import annotations

import numpy as np
import pandas as pd

from common.time import Time
from common.frame import TimeSeriesDataFrame

from load.loaders.gena import load_gena
from load.loaders.indices import load_indices
from load.loaders.spedas import load_spedas
from load.api.spedas import SPEDAS_OPTIONS, normalize_probe

__all__ = ["load_data", "load_pos", "load_rot", "load_sw"]


#: What ``product``, ``probe`` and ``key`` mean for each mission. ``key`` is
#: the "what part of this product do I want" selector: a variable name for
#: SPEDAS missions, a column for GENA, and the index name for indices.
_MISSION_ROLES = {
    "gena": {"probe": False, "key": "columns"},
    "themis": {"probe": True, "key": "variables"},
    "mms": {"probe": True, "key": "variables"},
    "goes": {"probe": True, "key": "variables"},
    "omni": {"probe": False, "key": "variables"},
    "indices": {"probe": False, "key": "product"},
}
_MISSION_ALIASES = {"index": "indices", "kyoto": "indices", "goes-r": "goes"}


def load_data(mission, trange, product=None, *, probe=None, key=None, **kwargs):
    """Load one mission product through its source loader.

    ``mission`` selects the role of the other arguments:

    - ``gena``: ``product`` is a Parquet suffix such as ``l2m_a``; there is no
      separate probe (it is encoded in the suffix); ``key`` selects columns
      (``key="Energy"``) and is equivalent to ``columns=[...]``.
    - ``themis`` / ``mms`` / ``goes``: ``product`` is an instrument or orbit product (``esa``,
      ``fpi``, ``mec``, ``orbit``, ...), ``probe`` is the spacecraft id (``b``, ``1``, ...), and
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
    if name in SPEDAS_OPTIONS:
        return load_spedas(name, trange, product, probe=probe, **kwargs)
    return load_indices(trange, product=product, **kwargs)


def load_sw(trange, **kwargs):
    """Load native OMNI solar-wind parameters and Dst over [start, end].

    Return a TimeSeriesDataFrame on the union of the source timestamps, with
    Pressure, BY_GSM, BZ_GSM, Beta, Mgs_mach_num, SYM_H, and dst columns.
    Exactly two ordered UTC endpoints are required. No time padding,
    interpolation, extrapolation, or constant fallback is performed; missing
    measurements remain NaN and source-loading errors propagate.

    Common options such as prefix and download are forwarded to both load_data()
    calls. Omit download to retain each source loader's default policy.
    To evaluate parameters at specific times, call frame.resample_at(times)
    explicitly after loading a covering interval.
    """
    times = Time(trange)
    if times.isscalar or len(times.datetime) != 2:
        raise ValueError("trange must contain exactly two UTC times: [start, end]")
    if times.datetime[0] > times.datetime[1]:
        raise ValueError("trange start must be <= end")
    dates = [value.isoformat() for value in times.datetime]
    if kwargs.pop("return_type", "dataframe") != "dataframe":
        raise ValueError("load_sw returns a TimeSeriesDataFrame")
    namelist = ["Pressure", "BY_GSM", "BZ_GSM", "Beta", "Mgs_mach_num", "SYM_H"]
    omni = load_data("omni", dates, "hro", variables=namelist, **kwargs)
    dst = load_data("indices", dates, "dst", **kwargs)

    tables, variables, sources = [], {}, {}
    for source, frame in (("omni", omni), ("indices", dst)):
        table = pd.DataFrame(frame).drop(columns="Timestamp")
        table.attrs = {}
        tables.append(table.set_index("Datetime"))
        variables.update(frame.attrs.get("variables", {}))
        sources[source] = {key: value for key, value in frame.attrs.items() if key != "variables"}
    joined = pd.concat(tables, axis=1, join="outer").sort_index()
    joined = joined.reindex(columns=[*namelist, "dst"]).reset_index()
    joined.insert(1, "Timestamp", Time(joined["Datetime"]).to_timestamp())
    result = TimeSeriesDataFrame(joined)
    result.attrs = {"source": ["omni", "indices"], "product": "solar_wind",
                    "variables": variables, "sources": sources}
    return result


# Explicit convenience operations using the public loaders.
def _evaluation_times(value):
    """Normalize evaluation times and return dates plus Unix seconds."""
    times = Time(value)
    dates = pd.DatetimeIndex([times.datetime]) if times.isscalar else times.datetime
    if dates.empty:
        raise ValueError("at least one evaluation time is required")
    return times, dates, np.asarray(times.to_timestamp(), dtype=float)


def _padded_evaluation_range(dates):
    """Return a one-hour-padded ISO range for a backend request."""
    pad = pd.Timedelta("1h")
    return [(dates.min() - pad).isoformat(), (dates.max() + pad).isoformat()]


def _get_themis_loc(trange, probe="b", coord="gsm", **kwargs):
    """Read native THEMIS STATE positions over an inclusive UTC interval."""
    probe = normalize_probe("themis", probe)
    variable_name = f"th{probe}_pos_{str(coord).lower()}"

    return load_spedas(
        "themis", trange, "state", variables=[variable_name], probe=probe,
        **kwargs,
    )


def _get_mms_loc(trange, probe="1", coord="gsm", **kwargs):
    """Read native MMS MEC positions through the shared SPEDAS loader."""
    probe = normalize_probe("mms", probe)
    coord = str(coord).lower()
    if coord not in {"gsm", "gse", "sm", "geo", "gei", "dsc", "dbcs"}:
        raise ValueError("unsupported MMS position coordinate")
    return load_spedas(
        "mms", trange, "mec", probe=probe,
        variables=[f"mms{probe}_mec_r_{coord}"], **kwargs,
    )


def _get_goes_loc(trange, satellite, coord="gsm", **kwargs):
    """Read native GOES orbit positions through the shared SPEDAS loader."""
    satellite = normalize_probe("goes", satellite)
    coord = str(coord).upper()
    if coord not in {"GSM", "GSE", "SM", "GEO"}:
        raise ValueError("GOES orbit coordinate must be GSM, GSE, SM, or GEO")
    return load_spedas(
        "goes", trange, "orbit", probe=satellite,
        variables=[f"g{satellite}_orbit_XYZ_{coord}"], **kwargs,
    )


def _get_queqiao_loc(trange, satellite=None):
    """Read native Queqiao positions recorded in GENA PDF engineering data.

    The PDF position columns are Moon-relative MJ2000 coordinates in metres.
    ``satellite`` is accepted for symmetry with the other mission helpers and
    is currently ignored because the configured PDF contains one Queqiao
    spacecraft trajectory.
    """
    frame = load_gena(
        trange=trange,
        product="pdf",
        columns=["MJ2000_Star_X", "MJ2000_Star_Y", "MJ2000_Star_Z"],
    )
    result = TimeSeriesDataFrame(frame)
    result.attrs = dict(frame.attrs)
    result.attrs.update(coord="MJ2000", origin="Moon", units="m")
    return result


def load_pos(mission, trange, id=None, *, coord="gsm", **kwargs):
    """Load all native position samples in the inclusive interval [start, end].

    ``trange`` must contain exactly two ordered UTC times. No interpolation,
    padding, or extrapolation is performed. All missions return a
    TimeSeriesDataFrame with source variable component columns (GENA/Queqiao
    uses MJ2000_Star_X/Y/Z). Native units
    are preserved (CDF units for external missions, Moon-relative MJ2000 metres
    for Queqiao; ``coord`` only selects external-mission variables).

    To interpolate explicitly for any mission, call frame.resample_at(times).
    THEMIS ids are a--e, MMS ids
    are 1--4, and GOES requires its spacecraft number. All SPEDAS missions
    accept the same download/prefix options as load_data(mission, ...).
    """
    times = Time(trange)
    if times.isscalar or len(times.datetime) != 2:
        raise ValueError("trange must contain exactly two UTC times: [start, end]")
    if times.datetime[0] > times.datetime[1]:
        raise ValueError("trange start must be <= end")
    dates = [value.isoformat() for value in times.datetime]
    if kwargs.pop("return_type", "dataframe") != "dataframe":
        raise ValueError("load_pos returns a DataFrame; use load_data(..., return_type='dict') for raw variables")
    key = str(mission).strip().lower()
    if key == "themis":
        return _get_themis_loc(dates, probe="a" if id is None else id, coord=coord, **kwargs)
    if key == "mms":
        return _get_mms_loc(dates, probe="1" if id is None else id, coord=coord, **kwargs)
    if key in {"goes", "goes-r"}:
        if id is None:
            raise ValueError("id is required for GOES")
        return _get_goes_loc(dates, id, coord=coord, **kwargs)
    if key in {"queqiao", "qq", "gena"}:
        if kwargs:
            raise TypeError(f"unsupported Queqiao position options: {', '.join(kwargs)}")
        return _get_queqiao_loc(dates, id)
    raise ValueError("mission must be 'themis', 'mms', 'goes', or 'queqiao'")


def load_rot(mission, trange, *, bounds_error=True):
    """Return interpolated Queqiao MJ2000 attitude as unit xyzw quaternions.

    Accepts the mission names queqiao, qq, and gena. A scalar evaluation
    time returns shape (4,), and a sequence returns shape (N, 4) in the
    requested order. Interpolation uses SLERP within valid PDF coverage.
    Set bounds_error=False to return NaN rows outside that coverage.
    """
    if str(mission).strip().lower() not in {"queqiao", "qq", "gena"}:
        raise ValueError("load_rot supports the Queqiao attitude recorded in GENA PDF")
    times, dates, evaluation = _evaluation_times(trange)
    frame = load_gena(
        trange=_padded_evaluation_range(dates),
        product="pdf",
        columns=["MJ2000_q1", "MJ2000_q2", "MJ2000_q3"],
        optional_columns=["PosePresent"],
    )
    if bounds_error:
        result = frame.rot_at(times).as_quat()
    else:
        from common.coord import interpolate_quaternions

        source, rotations = frame._rot_from_pdf()
        covered = (evaluation >= source[0]) & (evaluation <= source[-1])
        result = np.full((len(evaluation), 4), np.nan)
        result[covered] = interpolate_quaternions(source, rotations, evaluation[covered])
    return result[0] if times.isscalar else result
