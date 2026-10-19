"""Shared helpers for loaders that accept one or more UTC intervals."""

from __future__ import annotations

import pandas as pd

from common.time import Time


def is_multi_interval(trange):
    """Return whether *trange* is a nested sequence of time intervals."""
    return (
        isinstance(trange, (list, tuple))
        and bool(trange)
        and isinstance(trange[0], (list, tuple))
    )


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
