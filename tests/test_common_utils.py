import copy
import pandas as pd
import pytest
from common import Time


def test_time_delegates_scalar_pandas_methods():
    time = Time("2025-10-10 13:21:30")
    assert time.floor("D") == pd.Timestamp("2025-10-10")
    assert time.ceil("D") == pd.Timestamp("2025-10-11")
    assert time.year == 2025
    assert time.strftime("%Y-%m-%d") == "2025-10-10"
    assert time.datetime.hour == 13


def test_time_delegates_index_methods_without_sorting_or_deduplicating():
    time = Time(["2025-10-11 13:00", "2025-10-10 12:00", "2025-10-11 01:00"])
    pd.testing.assert_index_equal(time.floor("D"), pd.DatetimeIndex(["2025-10-11", "2025-10-10", "2025-10-11"]).as_unit("ns"))
    assert time.shape == (3,)
    assert Time([]).floor("D").empty


def test_time_missing_attributes_and_copy():
    time = Time(0)
    for name in ("to_pandas", "normalize_days", "nonexistent"):
        with pytest.raises(AttributeError):
            getattr(time, name)
    assert copy.deepcopy(time).datetime == time.datetime
