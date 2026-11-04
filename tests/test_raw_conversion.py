"""Scientific boundary checks for unmodified GENA event conversion."""

import numpy as np
import pandas as pd
import pytest

from gena.constants import init_param
from gena.raw import conversion as conv
from gena.raw import structure


def event_frame(probe="a", snid=1, tof=80.0):
    params = init_param(snid)
    row = {key: 0 for key in structure.df1a_type}
    row.update(UUID="event-1", Timestamp=1723180000.0, HV_C=-2500.0,
               HV_G=0.0, Tof_1=tof)
    for axis in ("X", "Y"):
        prefix = f"{probe.upper()}_{axis}"
        delay = params[f"{prefix}_SUM"]
        difference = -params[f"{probe.upper()}_B{axis}"] / params[f"{probe.upper()}_K{axis}"]
        row[f"T{axis}_1s"] = (delay + difference) / 2
        row[f"T{axis}_2s"] = (delay - difference) / 2
        for edge in (1, 2):
            row[f"T{axis}_{edge}e"] = row[f"T{axis}_{edge}s"] + tof
    return pd.DataFrame([row]), params


def test_ep_unsigned_codes_do_not_wrap():
    codes = pd.Series([0, 1, 2044, 2045, 4094, 4095], dtype="uint16")
    actual = conv.ep_by_index(codes)
    np.testing.assert_allclose(actual, [0, 2043 * 100 / 123, 0,
                                       -100 / 123, -2050 * 100 / 123, 0], rtol=1e-6)
    assert actual.dtype == np.float32


def test_invalid_ep_codes_remain_invalid():
    assert np.isnan(conv.ep_by_index(np.array([-1, 4096, np.inf, np.nan]))).all()


def test_main_tof_has_priority_and_is_not_overwritten():
    frame = pd.DataFrame({"Tof_M": [80.0], "Tof_X": [81.0], "Tof_Y": [82.0]})
    np.testing.assert_array_equal(conv.ENA_flytime(frame), [80.0])


def test_axis_tof_fallback_uses_absolute_or_relative_agreement():
    frame = pd.DataFrame({"Tof_M": [0.0, 1000.0],
                          "Tof_X": [20.0, 100.0], "Tof_Y": [23.0, 108.0]})
    np.testing.assert_array_equal(conv.ENA_flytime(frame), [20.0, 100.0])


@pytest.mark.parametrize("main,x,y", [
    (80, -10, -20), (80, 0, 80), (80, np.nan, 80),
    (80, 80, np.inf), (80, 80, 200), (0, 0, 0),
])
def test_unsupported_tof_is_not_fabricated(main, x, y):
    frame = pd.DataFrame({"Tof_M": [main], "Tof_X": [x], "Tof_Y": [y]})
    assert np.isnan(conv.ENA_flytime(frame)[0])


def test_invalid_tof_cannot_be_reclassified_by_ep():
    frame = pd.DataFrame({"Tof": [np.nan, np.inf, 0, -1, 20, 200, 80, 80, 80],
                          "EP": [0, 0, 0, 2000, np.nan, np.nan, 0, 2000, np.nan]})
    result = conv.ENA_spec(frame)
    np.testing.assert_array_equal(result, [3, 3, 3, 3, 1, 2, 1, 2, 0])
    assert result.dtype == np.uint8


def test_energy_requires_finite_positive_tof_geometry_and_known_species():
    frame = pd.DataFrame({"Tof": [100, 100, 0, -100, np.inf, 100, 100, 100],
                          "Spec": [1, 2, 1, 1, 1, 0, 3, 1],
                          "X1": [0, 0, 0, 0, 0, 0, 0, np.nan],
                          "X2": 0, "Y1": 0, "Y2": 0})
    result = conv.ENA_energy(frame)
    # A hydrogen atom travelling 50 mm in 100 ns has about 1.306 keV.
    assert result[0] == pytest.approx(1.306, rel=0.002)
    assert result[1] / result[0] == pytest.approx(15.8695, rel=1e-6)
    assert np.isnan(result[2:]).all()


def test_invalid_acceleration_has_no_finite_electron_correction():
    frame = pd.DataFrame({"HV_C": [0, np.nan, np.inf, -2500], "HV_G": 0})
    correction = conv.f_tof_ep(frame)
    assert np.isnan(correction[:3]).all()
    assert correction[3] == pytest.approx(1.794, rel=0.001)


@pytest.mark.parametrize("probe", ["a", "b"])
@pytest.mark.parametrize("snid", [1, 2])
def test_clean_events_preserve_inputs_and_calibration(probe, snid):
    frame, params = event_frame(probe, snid)
    original = frame.copy(deep=True)
    result = conv.l1a_2_l2(frame, params=params, prob=probe.upper())
    pd.testing.assert_frame_equal(frame, original)
    assert len(result) == 1
    assert result["UUID"].tolist() == ["event-1"]
    np.testing.assert_allclose(result[["X1", "X2", "Y1", "Y2"]], 0, atol=1e-5)
    np.testing.assert_allclose(result[["Tof_M", "Tof_X", "Tof_Y"]], 80, atol=1e-5)
    assert result["Tof"].iloc[0] == pytest.approx(81.794, abs=0.001)
    assert result["Spec"].iloc[0] == 1
    assert np.isfinite(result["Energy"]).all()


def test_explicit_ground_probe_b_calibration_is_used():
    frame, params = event_frame("b", 2)
    ground = conv.l1a_2_l2(frame, params=params, prob="b")
    orbit = conv.l1a_2_l2(frame, params=init_param(1), prob="b")
    assert ground["Y1"].iloc[0] == pytest.approx(0, abs=1e-5)
    assert abs(orbit["Y1"].iloc[0]) > 10


def test_invalid_l2_energy_does_not_remove_event_rows():
    frame, params = event_frame()
    frame["TY_1e"] += 200
    frame["TY_2e"] += 200
    result = conv.l1a_2_l2(frame, params=params)
    assert len(result) == len(frame)
    assert result["UUID"].tolist() == frame["UUID"].tolist()
    assert np.isnan(result["Tof"]).all()
    assert result["Spec"].tolist() == [3]
    assert np.isnan(result["Energy"]).all()


@pytest.mark.parametrize("probe,snid,mapping", [
    ("a", 1, [5, 7, 1, 3]), ("a", 2, [5, 7, 1, 3]),
    ("b", 1, [1, 3, 5, 7]), ("b", 2, [7, 5, 3, 1]),
])
def test_l0_channel_mapping_is_probe_and_satellite_specific(probe, snid, mapping):
    row = dict(UUID="raw-1", PID=1, BM=0, BX=0, Timestamp=1723180000.0,
               TimeInternal=0, EvtLevel=0, EP=np.uint16(2045),
               prob_flag=0 if probe == "a" else 1)
    for i in range(4):
        row[f"ch_tof_{i}"] = 100 + i
    for channel in (1, 3, 5, 7):
        for suffix in ("s", "e"):
            row[f"ch_{channel}{suffix}"] = 100 + channel
    result = conv.l0_2_l1(pd.DataFrame([row]), probe.upper(), snid)
    np.testing.assert_allclose(result[["TX_1s", "TX_2s", "TY_1s", "TY_2s"]],
                               [[(6 + channel) * 0.081 for channel in mapping]], rtol=1e-6)
    assert result["EP"].iloc[0] < 0


@pytest.mark.parametrize("probe,snid", [("c", 1), ("a", 3), ("b", 0)])
def test_invalid_probe_or_satellite_raises_library_error(probe, snid):
    with pytest.raises(ValueError):
        conv.l0_2_l1(pd.DataFrame(), probe, snid)


@pytest.mark.parametrize("method", [1, 2])
@pytest.mark.parametrize("empty", [False, True])
def test_legacy_repair_is_explicitly_unavailable(method, empty):
    frame, params = event_frame()
    if empty:
        frame = frame.iloc[:0]
    original = frame.copy(deep=True)
    with pytest.raises(NotImplementedError, match=f"method {method} is disabled"):
        conv.l1a_2_l1b(frame, params, "a", method)
    pd.testing.assert_frame_equal(frame, original)


def test_unknown_repair_method_raises_value_error():
    with pytest.raises(ValueError, match="method must"):
        conv.l1a_2_l1b(pd.DataFrame(), init_param(1), "a", 99)
