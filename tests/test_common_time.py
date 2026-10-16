"""Time conversion boundaries, using in-memory data and no downloads."""

import builtins
import datetime

import numpy as np
import pandas as pd
import pytest
from astropy.time import Time as AstroTime

from common import Time
from common.frame import TimeSeriesDataFrame


@pytest.mark.parametrize("value", [
    "2000-01-01T12:00:00Z", "2000-01-01 12:00:00",
    "2000-01-01T20:00:00+08:00", datetime.datetime(2000, 1, 1, 12),
    datetime.datetime(2000, 1, 1, 20, tzinfo=datetime.timezone(datetime.timedelta(hours=8))),
    pd.Timestamp("2000-01-01T12:00:00Z"), np.datetime64("2000-01-01T12:00:00"),
    946728000, np.int64(946728000), np.array(946728000),
])
def test_scalar_inputs_share_utc_meaning(value):
    result = Time(value)
    assert result.isscalar and isinstance(result.datetime, pd.Timestamp)
    assert result.datetime == pd.Timestamp("2000-01-01T12:00:00")
    assert result.datetime.tzinfo is None
    assert result.datetime.unit == "ns"
    assert result.to_timestamp() == 946728000.0
    assert result.to_iso() == "2000-01-01T12:00:00.000000000Z"
    assert result.to_str("%Y/%m/%d %H:%M") == "2000/01/01 12:00"
    assert result.to_numpy() == np.datetime64("2000-01-01T12:00:00", "ns")


@pytest.mark.parametrize("value", [datetime.date(2000, 1, 1), "2000-01-01"])
def test_date_only_means_utc_midnight(value):
    assert Time(value).datetime == pd.Timestamp("2000-01-01")


@pytest.mark.parametrize("value", [[0, 1], (0, 1), np.array([0, 1]), pd.Series([0, 1]), pd.Index([0, 1])])
def test_array_operations_use_native_pandas(value):
    result = Time(value)
    assert not result.isscalar and result.datetime.shape == (2,)
    assert isinstance(result.datetime, pd.DatetimeIndex)
    assert str(result.datetime.dtype) == "datetime64[ns]"
    assert isinstance(result.datetime[0], pd.Timestamp)
    assert result.datetime[:1].shape == (1,)
    assert result.datetime[:0].empty
    np.testing.assert_array_equal(result.to_timestamp(), [0., 1.])
    np.testing.assert_array_equal(result.to_iso(), ["1970-01-01T00:00:00.000000000Z", "1970-01-01T00:00:01.000000000Z"])
    np.testing.assert_array_equal(result.to_str("%S"), ["00", "01"])
    assert result.to_numpy("s").dtype == np.dtype("datetime64[s]")
    with pytest.raises(ValueError):
        result.to_numpy("ps")


@pytest.mark.parametrize("value", [0, [0, 1]])
def test_time_is_a_converter_not_an_array_container(value):
    result = Time(value)
    for operation in (lambda: result[0], lambda: iter(result), lambda: len(result)):
        with pytest.raises(TypeError):
            operation()


@pytest.mark.parametrize("value", [[], (), np.array([], dtype="datetime64[ns]"), pd.DatetimeIndex([]), pd.Series([], dtype=float)])
def test_empty_arrays_keep_shape_and_conversion_units(value):
    result = Time(value)
    assert not result.isscalar and result.datetime.empty
    assert result.to_numpy().dtype == np.dtype("datetime64[ns]")
    assert result.to_timestamp().shape == result.to_iso().shape == result.to_str().shape == (0,)
    assert result.to_astropy().shape == (0,)
    assert result.convert("jd").shape == (0,)


@pytest.mark.parametrize("value", [
    pd.Timestamp("2026-10-12T08:00:00.123456789+08:00"),
    pd.DatetimeIndex(["2026-10-12T08:00:00.123456789+08:00"]),
    pd.Series(pd.DatetimeIndex(["2026-10-12T08:00:00.123456789+08:00"])),
    np.array(["2026-10-12T00:00:00.123456789"], dtype="datetime64[ns]"),
])
def test_nanoseconds_survive_native_and_astropy_roundtrips(value):
    result = Time(value)
    expected = np.datetime64("2026-10-12T00:00:00.123456789", "ns")
    np.testing.assert_array_equal(result.to_numpy(), expected if result.isscalar else [expected])
    restored = Time(result.to_astropy())
    difference = np.asarray(restored.to_numpy()).astype("int64") - np.asarray(result.to_numpy()).astype("int64")
    assert np.all(np.abs(difference) <= 1)
    scalar = result if result.isscalar else Time(result.datetime[0])
    assert scalar.to_iso() == "2026-10-12T00:00:00.123456789Z"
    np.testing.assert_array_equal(result.to_timestamp(), scalar.to_timestamp() if result.isscalar else [scalar.to_timestamp()])


@pytest.mark.parametrize("value, format, expected", [
    (2451545, None, "1970-01-29T08:59:05"),
    (2451545, "jd", "2000-01-01T12:00:00"),
    (51544.5, "mjd", "2000-01-01T12:00:00"),
    (946728000000, "unix_ms", "2000-01-01T12:00:00"),
    (946728000123456789, "unix_ns", "2000-01-01T12:00:00.123456789"),
])
def test_numeric_units_are_explicit(value, format, expected):
    assert Time(value, format=format).datetime == pd.Timestamp(expected)
    assert Time([value], format=format).datetime[0] == pd.Timestamp(expected)


def test_scientific_scales_are_converted_to_utc_not_relabelled():
    tai = AstroTime("2026-10-12T00:00:37", scale="tai")
    assert Time(tai).datetime == pd.Timestamp("2026-10-12")
    assert Time("2026-10-12T00:00:37", format="isot", scale="tai").datetime == pd.Timestamp("2026-10-12")
    time = Time(2451545, format="jd", scale="tt")
    assert time.datetime == pd.Timestamp("2000-01-01T11:58:55.816")
    assert time.convert("jd", scale="tt") == pytest.approx(2451545)
    assert time.to_astropy().scale == "utc"
    gps = time.convert("gps")
    assert abs(Time(gps, format="gps").to_timestamp() - time.to_timestamp()) < 1e-6
    assert time.convert("mjd", scale="tt") == pytest.approx(51544.5)


def test_scientific_input_converts_directly_without_reconstructing_time(monkeypatch):
    original_init = Time.__init__
    calls = []

    def tracked_init(self, value, **kwargs):
        calls.append(value)
        original_init(self, value, **kwargs)

    monkeypatch.setattr(Time, "__init__", tracked_init)
    result = Time(AstroTime("2026-10-13T00:00:37.123456789", scale="tai"))
    assert len(calls) == 1
    assert result.datetime == pd.Timestamp("2026-10-13T00:00:00.123456789")
    native = Time(pd.Timestamp("2026-10-13T00:00:37.123456789"), scale="tai")
    assert native.datetime == result.datetime


def test_mixed_inputs_and_existing_time():
    values = [0, "1970-01-01T00:00:01Z", np.datetime64("1970-01-01T00:00:02"),
              pd.Timestamp("1970-01-01T00:00:03"), AstroTime("1970-01-01T00:00:04"), Time(5)]
    np.testing.assert_array_equal(Time(values).to_timestamp(), np.arange(6))
    np.testing.assert_array_equal(Time([1, "1970-01-01T00:00:02Z"]).to_timestamp(), [1., 2.])
    assert Time(Time(0)).isscalar
    assert Time(Time([0])).datetime.shape == (1,)


@pytest.mark.parametrize("value", [
    None, pd.NaT, np.datetime64("NaT"), [None], [pd.NaT], [0, float("nan")],
    float("inf"), -float("inf"), True, 1j, [[0, 1]], np.ones((2, 2)),
    "not a date", "2500-01-01", np.array([0], dtype="datetime64[s]") + np.timedelta64(300000, "D"),
    np.ma.array([0], mask=[True]),
])
def test_invalid_or_unrepresentable_inputs_fail(value):
    with pytest.raises((TypeError, ValueError, OverflowError)):
        Time(value)


def test_bad_format_scale_mask_and_leap_seconds_fail():
    for construct in [lambda: Time(0, format="nonsense"), lambda: Time(0, scale="nonsense"),
                      lambda: Time(AstroTime(np.ma.array([0], mask=[True]), format="unix")),
                      lambda: Time(AstroTime("2016-12-31T23:59:60")),
                      lambda: Time("2016-12-31T23:59:60Z")]:
        with pytest.raises((TypeError, ValueError)):
            construct()


def test_pandas_path_does_not_import_scientific_packages(monkeypatch):
    original_import = builtins.__import__

    def guarded(name, *args, **kwargs):
        if name.startswith(("astropy", "spacepy")):
            raise AssertionError("ordinary pandas conversion imported a scientific package")
        return original_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", guarded)
    result = Time(pd.date_range("2026-10-12", periods=2))
    assert result.to_numpy().dtype == np.dtype("datetime64[ns]")
    assert Time([0, "1970-01-01T00:00:01Z"]).to_timestamp().tolist() == [0., 1.]


def test_dataframe_column_stays_native_and_supports_pandas_operations():
    times = Time(["2026-10-12T08:00:00+08:00", "2026-10-12T08:01:00+08:00"])
    frame = TimeSeriesDataFrame({"Datetime": times.datetime, "value": [1, 2]})
    assert str(frame["Datetime"].dtype) == "datetime64[ns]"
    assert frame["Datetime"].dt.hour.tolist() == [0, 0]
    assert frame[frame["Datetime"] > times.datetime[0]]["value"].tolist() == [2]
    np.testing.assert_array_equal(pd.DatetimeIndex(frame["Datetime"]).asi8, times.datetime.asi8)


def test_spacepy_conversion_handles_scalars_arrays_and_empty_input():
    from spacepy.time import Ticktock

    for value in [0, [0, 1], []]:
        result = Time(value)
        ticks = result.to_tick()
        assert isinstance(ticks, Ticktock) and len(ticks) == (1 if result.isscalar else len(result.datetime))
        restored = Time(ticks)
        np.testing.assert_allclose(restored.to_timestamp(), np.atleast_1d(result.to_timestamp()), atol=1e-6, rtol=0)
    leap = Ticktock(["2016-12-31T23:59:60"], "ISO")
    with pytest.raises(ValueError):
        Time(leap)


def test_day_tools_accept_time_objects():
    days = Time([0, 1, 86400]).floor("D").drop_duplicates().tz_localize("UTC")
    pd.testing.assert_index_equal(days, pd.date_range("1970-01-01", periods=2, tz="UTC").as_unit("ns"))


def test_integer_nanoseconds_do_not_overflow_or_wrap():
    assert Time(pd.Timestamp.max.value, format="unix_ns").datetime == pd.Timestamp.max
    assert Time(pd.Timestamp.min.value, format="unix_ns").datetime == pd.Timestamp.min
    for value in [np.uint64(2**64 - 1), np.uint64(2**63), [np.uint64(2**64 - 1)]]:
        with pytest.raises(ValueError, match="range"):
            Time(value, format="unix_ns")


def test_iterable_times_are_normalized_without_losing_numeric_meaning():
    np.testing.assert_array_equal(Time(range(3)).to_timestamp(), [0., 1., 2.])


@pytest.mark.parametrize("value, expected", [
    (Time(1), [1., 2., 3.]),
    (pd.Timestamp("1970-01-01T08:00:01+08:00"), [1., 2., 3.]),
    (pd.Series(pd.to_datetime(["1970-01-01T08:00:01+08:00", "1970-01-01T08:00:30+08:00"])),
     [[1., 30.], [2., 60.], [3., 90.]]),
])
def test_themis_position_uses_utc_seconds_for_every_input_shape(monkeypatch, value, expected):
    import importlib

    loader = importlib.import_module("load.load")

    calls = []

    def local_state(mission, trange, *args, **kwargs):
        calls.append(trange)
        return {"thb_pos_gsm": {"x": [0., 60.], "y": np.array([[0., 0., 0.], [60., 120., 180.]])}}

    monkeypatch.setattr(loader, "load_spedas", local_state)
    np.testing.assert_allclose(loader.get_themis_loc(value), expected)
    assert calls[0][0].startswith("1969-12-31T23:")
    assert calls[0][1].startswith("1970-01-01T01:")


def test_scene_time_normalization_keeps_native_pandas_and_numeric_seconds():
    from plots.v3d import SceneContext

    context = object.__new__(SceneContext)
    context.pose_timestamps = np.array([0., 60.])
    dates, seconds = context._normalize_times(Time([1, 30]))
    assert dates.tz is not None and str(dates.tz) == "UTC"
    np.testing.assert_array_equal(seconds, [1., 30.])
    with pytest.raises(ValueError, match="at least one"):
        context._normalize_times(Time([]))
