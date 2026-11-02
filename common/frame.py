"""Shared pandas frame types, time-axis operations, and variable access.

This module provides :class:`TimeIndexedDataFrame` for shared row operations
and :class:`TimeSeriesDataFrame` for sampled physical quantities. GENA event
tables derive from the former; cross-mission sampled-data loaders use the
latter.

It establishes one normalized, timezone-naive UTC ``Datetime`` column,
keeps pandas operations type-aware, and supplies selection, interval splitting,
synchronization, concatenation, variable reading, and lazy plotting access.

The implementation defines table structures and operations. Loading,
mission-specific decoding, and scientific conversion remain outside this
module; callers provide decoded columns and this module handles their common
time-series behavior.

"""

from __future__ import annotations

from collections.abc import Sequence
from copy import deepcopy

import numpy as np

import pandas as pd

from common.time import Time

__all__ = [
    "TimeIndexedDataFrame",
    "TimeSeriesDataFrame",
    "concat",
    "read_variable",
    "split_time_intervals",
]


def split_time_intervals(times: pd.Series | pd.DatetimeIndex, *, max_gap="60min") -> list[tuple[pd.Timestamp, pd.Timestamp]]:
    """Group native pandas datetimes where adjacent gaps exceed ``max_gap``.

    Accept only a datetime Series or DatetimeIndex, with no format conversion.
    Sort the values and retain their timezone. Empty input returns []; NaT is
    rejected. A gap equal to the threshold remains in the same interval.
    """
    if not isinstance(times, (pd.Series, pd.DatetimeIndex)):
        raise TypeError("times must be a pandas datetime Series or DatetimeIndex")
    if not pd.api.types.is_datetime64_any_dtype(times.dtype):
        raise TypeError("times must have a pandas datetime dtype")
    dates = pd.DatetimeIndex(times)
    if dates.hasnans:
        raise ValueError("time values must not contain NaT")
    if len(dates) == 0:
        return []
    dates = dates.sort_values()
    gap_limit = pd.Timedelta(max_gap)
    intervals = []
    start = previous = dates[0]
    for current in dates[1:]:
        if current - previous > gap_limit:
            intervals.append((start, previous))
            start = current
        previous = current
    intervals.append((start, previous))
    return intervals


class TimeIndexedDataFrame(pd.DataFrame):
    """DataFrame with canonical UTC ``Datetime``/``Timestamp`` columns.
    The DataFrame contains standard UTC Datetime / Timestamp columns.

    The datetime uses pandas' timezone-naive `datetime64[ns]` type. Per project convention, this value represents UTC.

    Time-series operations read Datetime directly. This field must be supplied during construction and time normalization is performed via `common.time.Time`.

    Loaders parse raw data to generate Unix second-level Timestamp columns. The construction logic does **not** infer timestamps from other fields.
    """

    def __init__(self, *args, **kwargs):
        """Build a frame and normalize its required ``Datetime`` column.

        Accepts normal pandas DataFrame arguments. Exactly one ``Datetime`` column is required; its values are normalized as UTC.
        """
        super().__init__(*args, **kwargs)
        if list(self.columns).count("Datetime") != 1:
            raise ValueError("TimeIndexedDataFrame requires exactly one Datetime column")
        self["Datetime"] = Time(self["Datetime"]).datetime

    @property
    def _constructor(self):
        """Keep the concrete type only for results with a time column."""
        def construct(*args, **kwargs):
            frame = pd.DataFrame(*args, **kwargs)
            if list(frame.columns).count("Datetime") == 1:
                return type(self)(frame)
            return frame

        return construct

    """
    User defined methods for selection, selection by another frame, interval splitting, synchronization, and plotting.
    """

    def select(self, cond):
        """Select rows using a query string, callable, or boolean mask.

        The result has a reset index and retains this subclass when ``Datetime`` remains.
        """
        if isinstance(cond, str):
            result = self.query(cond)
        elif callable(cond):
            result = self[cond(self)]
        else:
            result = self[cond]
        return result.reset_index(drop=True)

    def select_by(
        self,
        other: TimeSeriesDataFrame,
        cond,
        *,
        method: str = "nearest",
        max_gap="auto",
    ):
        """Select original rows using a condition evaluated on another frame.

        ``other`` is synchronized to this frame's native timestamps before ``cond`` is evaluated.
        Make sure that two TSDs share same Datetime column and same indexing/numeber of rows.
        Otherwise, the selection may not be meaningful.

        The returned rows always come directly from this frame; synchronization never interpolates or replaces their data.

        ``cond`` accepts the same query string, callable, or boolean-mask forms as :meth:`select`.

        Nearest-neighbor synchronization is the default because this operation is intended for cross-instrument selection.  
        Pass ``method="linear"`` explicitly when the selector values should be interpolated.
        """
        if not isinstance(other, TimeSeriesDataFrame):
            raise TypeError("other must be a TimeSeriesDataFrame")

        selector = other.resample_to(self, method=method, max_gap=max_gap)
        if isinstance(cond, str):
            mask = selector.eval(cond, engine="python")
            if not isinstance(mask, pd.Series) or len(mask) != len(selector):
                raise ValueError("cond must evaluate to one boolean value per synchronized row")
            mask = pd.array(mask, dtype="boolean").fillna(False).to_numpy(dtype=bool)
            selected = selector.select(mask)
        else:
            selected = selector.select(cond)
        selected_dates = selected["Datetime"]
        return self.select(self["Datetime"].isin(selected_dates))

    def between(self, time_list):
        """Return reset-index frames for each consecutive open interval.

        ``time_list`` supplies timestamp-like boundaries and each result contains ``start < Datetime < end``.
        """
        boundaries = list(time_list)
        dates = self["Datetime"]
        frames = []
        for start, end in zip(boundaries[:-1], boundaries[1:]):
            t0, t1 = Time([start, end]).datetime
            result = self[(dates > t0) & (dates < t1)].reset_index(drop=True)
            frames.append(result)
        return frames

    def split(self, tgap=1800):
        """Split sorted samples at gaps larger than ``tgap`` seconds.

        Returns ``(count, intervals, frames)`` with inclusive interval endpoints and reset-index parts.
        """
        if len(self) == 0:
            return 0, [], []
        work = self.sort_values("Datetime").reset_index(drop=True)
        time_intervals = split_time_intervals(
            work["Datetime"],
            max_gap=pd.Timedelta(seconds=tgap),
        )
        frames = []
        intervals = []
        for start, end in time_intervals:
            selected = work[(work["Datetime"] >= start) & (work["Datetime"] <= end)]
            part = selected.reset_index(drop=True)
            frames.append(part)
            intervals.append([start, end])
        return len(frames), intervals, frames

    def tplot(self, key=None, trange=None, f=None, ax=None, **kwargs):
        """Plot one scalar or three-component time-series variable.

        The plotting import is intentionally lazy so the data contract does not
        depend on Matplotlib at import time.
        """
        from plots import tplot

        if f is not None:
            kwargs["f"] = f
        if ax is not None:
            kwargs["ax"] = ax
        return tplot(self, key=key, trange=trange, **kwargs)

class TimeSeriesDataFrame(TimeIndexedDataFrame):
    """Samples of time-dependent physical quantities, supporting interpolation.

    Particle event rows belong in TimeIndexedDataFrame subclasses instead.
    """

    def resample(self, tres, *, method: str = "linear", max_gap="auto"):
        """Interpolate physical quantities onto a regular UTC grid.

        tres is a positive pandas Timedelta-compatible duration, e.g. "4s".
        The grid is aligned inward to the source coverage. This returns a frame,
        not pandas' aggregation Resampler. For bin aggregation, use an ordinary
        pd.DataFrame view and its pandas resample method explicitly.
        method and max_gap follow resample_to's column and outage rules.
        """
        resolution = pd.Timedelta(tres)
        if pd.isna(resolution) or resolution <= pd.Timedelta(0):
            raise ValueError("tres must be positive")
        dates = self["Datetime"]
        if len(dates) == 0:
            target_dates = pd.DatetimeIndex([], dtype="datetime64[ns]")
        else:
            start = dates.min().ceil(resolution)
            end = dates.max().floor(resolution)
            target_dates = (pd.date_range(start, end, freq=resolution)
                            if start <= end else pd.DatetimeIndex([], dtype="datetime64[ns]"))
        return self.resample_at(target_dates, method=method, max_gap=max_gap)

    def resample_to(self, other: pd.DataFrame, *, method: str = "linear", max_gap="auto"):
        """Interpolate this table's physical quantities onto another table's times.

        other must be a DataFrame with a normalized UTC Datetime column; it may
        contain physical samples or event timestamps. Target times are sorted,
        deduplicated and restricted to the common coverage. Only this table's
        values are interpolated: the reference table is never modified.

        method is "linear" (default) or "nearest". Bool, integer, quality,
        string and array columns always use nearest samples. No extrapolation
        is performed. max_gap accepts "auto" (three median valid cadences per
        column), a positive Timedelta-compatible duration, or None. Linear
        checks the bracketing gap; nearest checks distance to the chosen sample.
        Returns the concrete frame subclass with shallow-copied attrs.
        """
        if not isinstance(other, pd.DataFrame):
            raise TypeError("other must be a pandas DataFrame")
        dates = self["Datetime"]
        other_dates = other["Datetime"]
        if len(dates) == 0 or len(other_dates) == 0:
            target_dates = pd.DatetimeIndex([], dtype="datetime64[ns]")
        else:
            start = max(dates.min(), other_dates.min())
            end = min(dates.max(), other_dates.max())
            target_dates = pd.DatetimeIndex(other_dates[(other_dates >= start) & (other_dates <= end)].unique()).sort_values()
        return self.resample_at(target_dates, method=method, max_gap=max_gap)

    def resample_at(self, times, *, method: str = "linear", max_gap="auto"):
        """Evaluate physical quantities at explicitly supplied UTC times.

        Parameters
        ----------
        times : one-dimensional array-like
            DatetimeIndex, Series, datetime/string sequence, or Unix-second
            numeric sequence. Values are normalized to timezone-naive UTC
            datetime64[ns]; invalid/NaT times are rejected. Input order and
            repeated times are preserved. Targets outside valid source coverage
            remain present with missing physical values; no extrapolation occurs.
        method : {"linear", "nearest"}, default "linear"
            Interpolate continuous numeric columns or use nearest samples.
            Discrete, quality and array-valued columns always use nearest.
        max_gap : str, Timedelta-compatible duration, or None, default "auto"
            Auto uses three median valid cadences per column. Linear limits
            the bracketing gap; nearest limits distance to the selected sample.
            None disables gap limits without permitting extrapolation.

        Returns
        -------
        TimeSeriesDataFrame
            Concrete frame subclass with one row per requested time, regenerated
            Unix-second Timestamp and shallow-copied attrs. Empty source data
            produces missing physical values on the requested axis.

        Use for custom or irregular evaluation times without constructing a
        reference table. resample and resample_to build their axes and delegate
        to this shared implementation. Independent particle rows cannot use it.
        """
        target_dates = Time(times).datetime
        if method not in {"linear", "nearest"}:
            raise ValueError("method must be 'linear' or 'nearest'")
        automatic_gap = isinstance(max_gap, str) and max_gap.lower() == "auto"
        gap_ns = None if max_gap is None or automatic_gap else pd.Timedelta(max_gap).value
        if gap_ns is not None and (pd.isna(pd.Timedelta(max_gap)) or gap_ns <= 0):
            raise ValueError("max_gap must be positive")

        source_dates = self["Datetime"]
        if len(target_dates) == 0:
            return self.iloc[0:0].copy()

        base = pd.DataFrame(self, copy=False).copy()
        base["Datetime"] = source_dates
        base = base.sort_values("Datetime", kind="stable").drop_duplicates("Datetime", keep="first")
        source_ns = pd.DatetimeIndex(base["Datetime"]).as_unit("ns").asi8
        target_ns = target_dates.as_unit("ns").asi8
        synchronized = {"Datetime": target_dates, "Timestamp": target_ns.astype(float) / 1e9}
        for column in base.columns:
            if column in {"Datetime", "Timestamp"}:
                continue
            series = base[column]
            column_name = str(column).lower()
            is_numeric = pd.api.types.is_numeric_dtype(series.dtype)
            is_discrete = (
                pd.api.types.is_bool_dtype(series.dtype)
                or pd.api.types.is_integer_dtype(series.dtype)
                or "quality" in column_name
            )
            column_method = (
                "linear"
                if method == "linear" and is_numeric and not is_discrete
                else "nearest"
            )
            column_source_ns, values = _clean_samples(
                source_ns,
                series.to_numpy(),
                dtype=float if column_method == "linear" else object,
            )
            column_gap_ns = (
                _automatic_gap_ns(column_source_ns) if automatic_gap else gap_ns
            )
            synchronized[column] = _sync_values(
                column_source_ns,
                values,
                target_ns,
                column_gap_ns,
                method=column_method,
            )

        result = self._constructor(synchronized)
        result.attrs = dict(self.attrs)
        return result


# ---------------------------  PRE SYNC ---------------------------
def _value_is_valid(value):
    """Return whether a scalar or array-valued cell is a usable sample; invalid scalar NaN/NA values are removed."""
    if isinstance(value, (np.ndarray, list, tuple)):
        return True
    if isinstance(value, (float, np.floating)) and not np.isfinite(value):
        return False
    missing = pd.isna(value)
    return not bool(missing) if np.ndim(missing) == 0 else True

def _clean_samples(source_ns, values, *, dtype):
    """Remove invalid cells while preserving aligned nanosecond source times and requested dtype."""
    source = np.asarray(source_ns, dtype=np.int64)
    data = np.asarray(values, dtype=object)
    if len(source) != len(data):
        raise ValueError("source times and values must have the same length")
    valid = np.fromiter(
        (_value_is_valid(value) for value in data), dtype=bool, count=len(data)
    )
    return source[valid], data[valid].astype(dtype, copy=False)


def _automatic_gap_ns(source_ns):
    """Compute the automatic three-median-cadence interpolation limit."""
    times = np.asarray(source_ns, dtype=np.int64)
    if len(times) < 2:
        return None
    positive = np.diff(times)
    positive = positive[positive > 0]
    return int(np.median(positive) * 3) if len(positive) else None


def _sync_values(source_ns, values, target_ns, max_gap_ns, *, method):
    """Map one cleaned column onto target times by linear or nearest sampling, respecting the gap limit."""
    if method == "linear":
        data = np.asarray(values, dtype=float)
        result = np.full(len(target_ns), np.nan, dtype=float)
    elif method == "nearest":
        data = np.asarray(values, dtype=object)
        result = np.full(len(target_ns), np.nan, dtype=object)
    else:
        raise ValueError("method must be 'linear' or 'nearest'")

    source = np.asarray(source_ns, dtype=np.int64)
    if len(source) == 0:
        return result

    target = np.asarray(target_ns, dtype=np.int64)
    if method == "linear":
        if len(source) == 1:
            result[target == source[0]] = data[0]
            return result

        inside = (target >= source[0]) & (target <= source[-1])
        origin = source[0]
        result[inside] = np.interp(
            (target[inside] - origin) / 1e9,
            (source - origin) / 1e9,
            data,
        )
        if max_gap_ns is not None:
            right = np.searchsorted(source, target, side="left")
            clipped = np.clip(right, 0, len(source) - 1)
            exact = (right < len(source)) & (source[clipped] == target)
            left = right - 1
            bracketed = (left >= 0) & (right < len(source))
            acceptable = exact.copy()
            acceptable[bracketed] |= (
                source[right[bracketed]] - source[left[bracketed]] <= max_gap_ns
            )
            result[~acceptable] = np.nan
        return result

    right = np.searchsorted(source, target, side="left")
    right_clipped = np.clip(right, 0, len(source) - 1)
    left = np.clip(right - 1, 0, len(source) - 1)
    choose_left = np.abs(target - source[left]) <= np.abs(source[right_clipped] - target)
    nearest = np.where(choose_left, left, right_clipped)
    available = (target >= source[0]) & (target <= source[-1])
    if max_gap_ns is not None:
        available &= np.abs(source[nearest] - target) <= max_gap_ns
    for index in np.where(available)[0]:
        value = data[nearest[index]]
        result[index] = value.copy() if isinstance(value, np.ndarray) else value
    return result


def concat(frames: Sequence[pd.DataFrame]):
    """Concatenate frames in input order, retaining all rows, and merge metadata.

    Rows are concatenated with pandas' column-union rules, without sorting or
    deduplicating timestamps, and the index is reset. The first
    frame's concrete constructor and attributes are preserved; all
    ``attrs["variables"]`` mappings are merged in input order. Empty input
    returns an empty :class:`TimeSeriesDataFrame`.
    """
    frames = list(frames)
    if not frames:
        return TimeSeriesDataFrame({
            "Datetime": pd.Series(dtype="datetime64[ns]"),
            "Timestamp": pd.Series(dtype=float),
        })
    tables = []
    for frame in frames:
        table = pd.DataFrame(frame, copy=False)
        table.attrs = {}
        tables.append(table)
    result = frames[0]._constructor(pd.concat(tables, ignore_index=True, sort=False))
    result.attrs = dict(frames[0].attrs)
    variables = {}
    for frame in frames:
        variables.update(frame.attrs.get("variables", {}))
    if variables:
        result.attrs["variables"] = variables
    return result


def read_variable(frame, name):
    """Return time, values, and named coordinates using explicit frame metadata."""
    info = frame.attrs.get("variables", {}).get(name)
    if info is None:
        values = frame[name].to_numpy()
        mask = ~pd.isna(values)
        return {
            "time": pd.DatetimeIndex(frame["Datetime"])[mask],
            "values": values[mask],
            "axes": {},
        }
    columns, kind = info["columns"], info["kind"]
    if kind == "scalar":
        values = frame[columns[0]].to_numpy()
        mask = ~pd.isna(values)
        values = values[mask]
    elif kind == "vector":
        values = frame[columns].to_numpy(dtype=float)
        mask = np.isfinite(values).any(axis=1)
        values = values[mask]
    elif kind == "array":
        cells = frame[columns[0]].to_numpy()
        mask = np.array(
            [isinstance(value, (np.ndarray, list, tuple)) for value in cells],
            dtype=bool,
        )
        try:
            values = (
                np.stack(cells[mask])
                if mask.any()
                else np.empty((0, *info.get("sample_shape", ())))
            )
        except ValueError as exc:
            raise ValueError(f"Variable {name!r} has irregular sample shapes") from exc
    else:
        raise ValueError(f"Variable {name!r} has unknown kind {kind!r}")
    axes = {}
    for key, spec in info.get("axes", {}).items():
        cells = frame[spec["column"]].to_numpy()[mask]
        try:
            samples = (
                np.stack(cells)
                if len(cells)
                else np.empty((0, *spec.get("sample_shape", ())))
            )
        except ValueError as exc:
            raise ValueError(
                f"Variable {name!r}, axis {key!r} has irregular sample shapes"
            ) from exc
        if spec.get("static") and not len(cells):
            axes[key] = np.asarray(spec["values"])
        elif spec.get("static") and all(
            np.array_equal(value, samples[0]) for value in samples
        ):
            axes[key] = samples[0].copy()
        else:
            axes[key] = samples
    return {
        "time": pd.DatetimeIndex(frame["Datetime"])[mask],
        "values": values,
        "axes": axes,
        "units": info.get("units"),
        "attrs": deepcopy(info.get("attrs", {})),
    }
