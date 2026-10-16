import numpy as np
import pandas as pd
import pytest

from common.frame import TimeSeriesDataFrame
from load.load import _time_columns, assemble, parse_spedas
from load.api.spedas import _variable_layouts
from libs.regions import classify_kiehas, classify_liuzzo, classify_runov


def _table(data):
    return TimeSeriesDataFrame({**data, **_time_columns(data["Datetime"])})


def _fixture_frame(data):
    layouts = _variable_layouts(data, vectors={name for name in data if name.endswith("_velocity_gsm")})
    return assemble(parse_spedas(data, layouts=layouts))


@pytest.mark.parametrize("probe", ["thb", "thc"])
def test_kiehas_accepts_timeseries_and_aligns_sparse_variables(probe):
    density_time = pd.date_range("2025-01-01", periods=4, freq="4s")
    temperature_time = density_time.delete(2)
    esa = _fixture_frame({
        f"{probe}_peir_density": {
            "x": density_time,
            "y": np.array([0.1, 0.1, 0.1, 0.1]),
        },
        f"{probe}_peir_avgtemp": {
            "x": temperature_time,
            "y": np.array([400.0, 500.0, 700.0]),
        },
    })

    result = classify_kiehas(esa)

    assert isinstance(result, TimeSeriesDataFrame)
    assert np.diff(result["Timestamp"]).min() == 4.0
    assert result["plasma_sheet"].all()


def test_liuzzo_classifies_lobe_from_aligned_field_and_position():
    time = pd.date_range("2025-01-01", periods=2701, freq="4s")
    esa = _table(
        {"Datetime": time, "thb_peir_density": np.full(len(time), 0.1)}
    )
    fgm = _table(
        {
            "Datetime": time,
            "thb_fgs_gsm_x": np.full(len(time), 10.0),
            "thb_fgs_gsm_y": np.zeros(len(time)),
            "thb_fgs_gsm_z": np.zeros(len(time)),
        }
    )

    def midnight_position(timestamps):
        return np.column_stack(
            (-60.0 * np.ones(len(timestamps)), np.zeros(len(timestamps)), np.zeros(len(timestamps)))
        )

    result = classify_liuzzo(esa, fgm, midnight_position)

    assert isinstance(result, TimeSeriesDataFrame)
    assert result["lobe"].all()
    assert (result["region_code"] == 4).all()


def test_liuzzo_uses_density_only_for_off_midnight_magnetosheath():
    time = pd.date_range("2025-01-01", periods=10, freq="4s")
    esa = _table(
        {"Datetime": time, "thb_peir_density": np.full(len(time), 2.0)}
    )
    fgm = _table(
        {
            "Datetime": time,
            "thb_fgs_gsm_x": np.ones(len(time)),
            "thb_fgs_gsm_y": np.zeros(len(time)),
            "thb_fgs_gsm_z": np.zeros(len(time)),
        }
    )

    def flank_position(timestamps):
        return np.column_stack(
            (-60.0 * np.ones(len(timestamps)), 60.0 * np.ones(len(timestamps)), np.zeros(len(timestamps)))
        )

    result = classify_liuzzo(esa, fgm, flank_position)

    assert result["magnetosheath"].all()
    assert (result["region_code"] == 2).all()


def test_runov_synchronizes_esa_and_sst_frames_before_classification():
    esa_time = pd.date_range("2025-01-01 00:00:00.1", periods=20, freq="4s")
    sst_time = pd.date_range("2025-01-01 00:00:00.2", periods=20, freq="4s")
    esa = _table(
        {
            "Datetime": esa_time,
            "thb_peir_density": 0.1,
            "thb_peir_avgtemp": 400.0,
            "thb_peer_density": 0.1,
            "thb_peer_avgtemp": 100.0,
        }
    )
    sst = _table(
        {
            "Datetime": sst_time,
            "thb_psif_density": 0.01,
            "thb_psif_avgtemp": 400.0,
        }
    )

    def position(timestamps):
        return np.column_stack(
            (-60.0 * np.ones(len(timestamps)), np.zeros(len(timestamps)), np.zeros(len(timestamps)))
        )

    result = classify_runov(esa, sst, position)

    assert isinstance(result, TimeSeriesDataFrame)
    assert result["plasma_sheet"].all()
    assert np.diff(result["Timestamp"]).min() == 4.0


def test_explicit_mms_columns_match_themis_classification():
    times = pd.date_range("2025-01-01", periods=4, freq="4s")
    themis = _table({
        "Datetime": times, "thb_peir_density": [0.1, 0.2, 0.7, 0.1],
        "thb_peir_avgtemp": [400., 200., 500., 600.],
    })
    mms = _table(themis.rename(columns={
        "thb_peir_density": "mms1_dis_numberdensity_fast",
        "thb_peir_avgtemp": "ion_temperature_ev",
    }))
    mms.attrs.update(source="mms", probe="1")
    result = classify_kiehas(mms, columns={
        "density": "mms1_dis_numberdensity_fast", "temperature": "ion_temperature_ev",
    })
    np.testing.assert_array_equal(result.plasma_sheet, classify_kiehas(themis).plasma_sheet)
    assert result.attrs["source"] == "mms"

    field = _table({"Datetime": times, "bx": 10., "by": 0., "bz": 0.})
    position = lambda t: np.tile([-60., 0., 0.], (len(t), 1))
    classified = classify_liuzzo(mms, field, position, columns={
        "density": "mms1_dis_numberdensity_fast", "magnetic": ["bx", "by", "bz"],
    })
    assert classified.lobe.all()
    assert classified.attrs["source"] == "mms"


def test_runov_accepts_explicit_moment_columns():
    time = pd.date_range("2025-01-01", periods=4, freq="4s")
    moments = _table({"Datetime": time, "ni": .1, "ti": 400., "ne": .1, "te": 100.})
    hot = _table({"Datetime": time, "nh": .01, "th": 400.})
    result = classify_runov(moments, hot, lambda t: np.tile([-60., 0., 0.], (len(t), 1)), columns={
        "ion_density": "ni", "ion_temperature": "ti", "electron_density": "ne",
        "electron_temperature": "te", "hot_density": "nh", "hot_temperature": "th",
    })
    assert result.plasma_sheet.all()
