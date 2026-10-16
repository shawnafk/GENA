"""Pandas-backed UTC times with Astropy scientific-format conversion."""

import datetime
from collections.abc import Iterable

import numpy as np
import pandas as pd
from astropy.time import Time as AstroTime

__all__ = ["Time"]


class Time:
    """Convert scalar or 1-D time input to native pandas UTC nanoseconds.

    datetime is a Timestamp or DatetimeIndex without a timezone marker.
    Bare numbers mean Unix seconds; other numeric formats must be explicit.
    Indexing, slicing and arithmetic belong to the native pandas object.
    """

    def __init__(self, value, *, format=None, scale=None):
        self.datetime = self._convert(value, format=format, scale=scale)

    def __getattr__(self, name):
        """Delegate missing attributes and methods to the native pandas time.

        pandas methods return their native results and do not wrap them in
        Time. Special Python operations such as indexing are not delegated.
        """
        native = object.__getattribute__(self, "__dict__").get("datetime")
        if native is None:
            raise AttributeError(f"{type(self).__name__!r} object has no attribute {name!r}")
        return getattr(native, name)

    @staticmethod
    def _convert(value, *, format=None, scale=None):
        """Select the pandas or scientific conversion path."""
        if scale is not None and scale not in AstroTime.SCALES:
            raise ValueError(f"unsupported time scale {scale!r}")
        if isinstance(value, Time):
            return value.datetime
        if isinstance(value, AstroTime):
            return Time._astropy_to_datetime(value)
        if value is None or value is pd.NaT:
            raise ValueError("time values must not contain None or NaT")
        native = (str, bytes, datetime.date, np.ndarray, np.generic,
                  pd.Series, pd.Index, int, float, complex, list, tuple, range, dict)
        if not isinstance(value, native):
            from spacepy.time import Ticktock

            if isinstance(value, Ticktock):
                return Time._astropy_to_datetime(value.APT)
            if isinstance(value, Iterable):
                value = list(value)
            else:
                raise TypeError(f"unsupported time input {type(value).__name__}")

        calendar = {None, "iso", "isot", "datetime", "datetime64"}
        unix = {"unix", "unix_ms", "unix_ns"}
        if format not in calendar | unix:
            return Time._astropy_to_datetime(AstroTime(value, format=format, scale=scale))
        result = Time._as_datetime(value, format=format)
        # Unix epochs denote instants independently of the display scale.
        numeric = np.asarray(value).dtype.kind in "iuf"
        if scale not in (None, "utc") and format not in unix and not (format is None and numeric):
            raw = result.to_datetime64() if isinstance(result, pd.Timestamp) else result.to_numpy()
            return Time._astropy_to_datetime(AstroTime(raw, format="datetime64", scale=scale))
        return result

    @staticmethod
    def _as_datetime(value, *, format=None):
        """Parse ordinary inputs and establish the pandas UTC/ns representation."""
        if isinstance(value, np.ma.MaskedArray):
            raise ValueError("masked time values are not supported")
        raw = np.asarray(value)
        if (isinstance(value, (list, tuple)) and raw.dtype.kind in "US"
                and any(not isinstance(item, (str, bytes)) for item in value)):
            raw = np.asarray(value, dtype=object)
        if raw.ndim > 1:
            raise ValueError("Time accepts only scalar or one-dimensional input")
        if raw.dtype.kind in "bc" or (raw.dtype.kind in "iuf" and not np.isfinite(raw).all()):
            raise ValueError("time values must be finite and cannot be boolean or complex")
        units = {"unix": "s", "unix_ms": "ms", "unix_ns": "ns"}
        numeric = raw.dtype.kind in "iuf"
        if numeric and format not in units and format is not None:
            raise ValueError(f"format {format!r} requires calendar datetime values")
        if format in units or numeric:
            unit = units.get(format, "s")
            if raw.size and raw.dtype.kind in "iu":
                factor = {"s": 1_000_000_000, "ms": 1_000_000, "ns": 1}[unit]
                if (int(raw.min()) * factor < pd.Timestamp.min.value
                        or int(raw.max()) * factor > pd.Timestamp.max.value):
                    raise ValueError("time is outside pandas' nanosecond range")
            result = pd.to_datetime(value, unit=unit, utc=True)
        elif raw.ndim == 1 and raw.size == 0:
            return pd.DatetimeIndex([], dtype="datetime64[ns]")
        elif (raw.ndim == 1 and raw.dtype.kind == "O"
              and not (isinstance(value, (pd.Series, pd.Index))
                       and pd.api.types.is_datetime64_any_dtype(value.dtype))):
            values = [Time._convert(item, format=format) for item in value]
            if any(not isinstance(item, pd.Timestamp) for item in values):
                raise ValueError("each time sequence element must be a scalar")
            result = pd.DatetimeIndex(values).tz_localize("UTC")
        else:
            result = pd.to_datetime(value, utc=True, format="mixed")
        result = pd.Timestamp(result) if raw.ndim == 0 else pd.DatetimeIndex(result).rename(None)
        if np.any(pd.isna(result)):
            raise ValueError("time values must not contain NaT")
        return result.tz_convert(None).as_unit("ns")

    @staticmethod
    def _astropy_to_datetime(value):
        """Convert a scientific instant to pandas through UTC datetime64."""
        if value.ndim > 1:
            raise ValueError("Time accepts only scalar or one-dimensional input")
        if np.any(value.mask):
            raise ValueError("masked time values are not supported")
        if value.size == 0:
            return pd.DatetimeIndex([], dtype="datetime64[ns]")
        utc = value.utc.replicate()
        utc.precision = 9
        try:
            raw = utc.to_value("datetime64")
        except ValueError as exc:
            raise ValueError("time cannot be represented by pandas datetime (including UTC leap seconds)") from exc
        return Time._as_datetime(raw)

    @property
    def isscalar(self):
        """Whether datetime is a single pandas Timestamp."""
        return isinstance(self.datetime, pd.Timestamp)

    def to_timestamp(self):
        """Return Unix seconds, retaining scalar or array form."""
        values = self.datetime.value if self.isscalar else self.datetime.asi8
        return values / 1e9

    def to_iso(self):
        """Return UTC ISO strings with nine fractional digits and a Z suffix."""
        def iso(value):
            return f"{value.strftime('%Y-%m-%dT%H:%M:%S')}.{value.microsecond:06d}{value.nanosecond:03d}Z"
        if self.isscalar:
            return iso(self.datetime)
        return np.asarray([iso(value) for value in self.datetime], dtype=str)

    def to_str(self, fmt="%Y-%m-%d %H:%M:%S"):
        """Format UTC times using pandas strftime (%f is microseconds)."""
        result = self.datetime.strftime(fmt)
        return result if self.isscalar else result.to_numpy(dtype=str)

    def to_numpy(self, unit="ns"):
        """Return datetime64 values at nanosecond or coarser precision."""
        if unit not in {"Y", "M", "W", "D", "h", "m", "s", "ms", "us", "ns"}:
            raise ValueError("unit must be a datetime64 unit with nanosecond or coarser precision")
        dtype = np.dtype(f"datetime64[{unit}]")
        if self.isscalar:
            return self.datetime.to_datetime64().astype(dtype)
        return self.datetime.to_numpy(dtype=dtype, copy=True)

    def to_astropy(self):
        """Create an Astropy Time on the UTC scale without losing nanoseconds."""
        return AstroTime(self.to_numpy(), format="datetime64", scale="utc", precision=9)

    def to_tick(self):
        """Create a SpacePy Ticktock using its Astropy adapter."""
        from spacepy.time import Ticktock

        if not self.isscalar and self.datetime.empty:
            result = Ticktock([datetime.datetime(1970, 1, 1)], "UTC")
            result.remove(0)
            return result
        return Ticktock(self.to_astropy(), "APT")

    def convert(self, format, *, scale="utc"):
        """Return values in an Astropy format and time scale."""
        if scale not in AstroTime.SCALES:
            raise ValueError(f"unsupported time scale {scale!r}")
        return getattr(self.to_astropy(), scale).to_value(format)

