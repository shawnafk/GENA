"""Offline boundary checks for explicit product layouts and plain construction."""
import importlib

import numpy as np
import pandas as pd
import pytest

from common.frame import TimeSeriesDataFrame, concat
from gena.frame import ParticleDataFrame
from load.api.spedas import _variable_layouts
from common.frame import read_variable
from load.api.public import _clip_variable_dict, _merge_variable_dicts, assemble
from load.api.spedas import parse_spedas


loader = importlib.import_module("load.load")


def spedas_dict_to_frame(data, *, layouts, **provenance):
    return assemble(parse_spedas(data, layouts=layouts), **provenance)


def _hro_to_frame(data, product="hro"):
    layouts = {
        name: {"kind": "scalar", "time_axis": 0, "columns": [name]}
        for name in data
    }
    return assemble(
        parse_spedas(data, layouts=layouts),
        source="omni",
        product=product,
        level="hro",
    )


@pytest.mark.parametrize("data", [{}, {"Timestamp": [1.]}, {"Date": ["2025-01-01"]}, {"Epoch": [1.]}])
def test_constructor_requires_datetime_without_alias_inference(data):
    with pytest.raises(ValueError, match="Datetime"):
        TimeSeriesDataFrame(data)


@pytest.mark.parametrize("times", [["2025-01-01T08:00:00+08:00"], [1735689600.], [np.datetime64("2025-01-01")]])
def test_constructor_uses_shared_normalization(times):
    frame = TimeSeriesDataFrame({"Datetime": times, "value": [1]}, index=[17])
    assert str(frame["Datetime"].dtype) == "datetime64[ns]"
    assert frame.loc[17, "Datetime"] == pd.Timestamp("2025-01-01")
    assert "Timestamp" not in frame  # The required input is Datetime alone.
    assert frame.to_dict(orient="list")["value"] == [1]
    assert not hasattr(frame, "to_variable_dict")


def test_constructor_rejects_invalid_or_duplicate_time_columns():
    with pytest.raises(ValueError):
        TimeSeriesDataFrame({"Datetime": ["not a time"]})
    with pytest.raises(ValueError, match="NaT"):
        TimeSeriesDataFrame({"Datetime": [pd.NaT]})
    with pytest.raises(ValueError, match="exactly one"):
        TimeSeriesDataFrame([[1, 2]], columns=["Datetime", "Datetime"])


@pytest.mark.parametrize("cls", [TimeSeriesDataFrame, ParticleDataFrame])
def test_pandas_results_with_and_without_time_keep_correct_types(cls):
    frame = cls({"Datetime": pd.date_range("2025-01-01", periods=3), "value": [1., 2., 3.]})
    frame.attrs = {"source": "test"}
    for result in (frame.iloc[:2], frame.copy(), frame.select(frame["value"] > 1), concat([frame, frame]), frame.iloc[:0]):
        assert type(result) is cls
        assert result.attrs == frame.attrs
    for result in (frame[["value"]], frame.drop(columns="Datetime")):
        assert type(result) is pd.DataFrame
    assert isinstance(frame["value"], pd.Series)
    assert "Datetime" in repr(frame)
    assert "Datetime" in frame._repr_html_()
    assert type(frame[["Datetime"]]) is cls



def test_explicit_three_channels_and_three_components_are_distinct():
    raw = {name: {"x": [1., 2.], "y": [[1, 2, 3], [4, 5, 6]]} for name in ("alpha", "beta")}
    raw["alpha"]["v"] = [10, 20, 30]
    raw["alpha"]["units"] = "counts"
    raw["alpha"]["attrs"] = {"DEPEND_0": "native_clock", "CUSTOM": "retained"}
    layouts = {
        "alpha": {"kind": "array", "columns": ["channels"], "time_axis": 0,
                  "axes": {"v": {"column": "frequency", "time_axis": None}}},
        "beta": {"kind": "vector", "columns": ["north", "east", "down"], "time_axis": 0},
    }
    frame = spedas_dict_to_frame(raw, layouts=layouts)
    assert "alpha_x" not in frame and "alpha_energy" not in frame
    np.testing.assert_array_equal(frame.loc[0, "channels"], [1, 2, 3])
    assert frame["north"].tolist() == [1, 4]
    restored = read_variable(frame, "alpha")
    assert restored["units"] == "counts"
    assert restored["attrs"] == raw["alpha"]["attrs"]
    np.testing.assert_array_equal(restored["axes"]["v"], [10, 20, 30])


@pytest.mark.parametrize("shape", [(2, 3), (2, 3, 4), (2, 3, 4, 5)])
def test_unknown_multidimensional_product_is_not_guessed_as_vector(shape):
    raw = {"unknown": {"x": [1., 2.], "y": np.arange(np.prod(shape)).reshape(shape)}}
    frame = spedas_dict_to_frame(raw, layouts=_variable_layouts(raw))
    assert frame.attrs["variables"]["unknown"]["kind"] == "array"
    assert frame.loc[0, "unknown"].shape == shape[1:]
    np.testing.assert_array_equal(read_variable(frame, "unknown")["values"], raw["unknown"]["y"])


def test_product_specific_vector_and_energy_rules():
    raw = {
        "tha_peir_velocity_gsm": {"x": [1., 2.], "y": np.ones((2, 3))},
        "tha_peir_en_eflux": {"x": [1., 2.], "y": np.ones((2, 3)), "v": [1, 2, 3]},
        "tha_peir_mode": {"x": [1., 2.], "y": np.zeros((2, 3))},
    }
    frame = spedas_dict_to_frame(raw, layouts=_variable_layouts(raw, vectors={"tha_peir_velocity_gsm"}, energy={"tha_peir_en_eflux"}))
    assert "tha_peir_velocity_gsm_x" in frame
    assert "tha_peir_en_eflux_energy" in frame
    assert "tha_peir_mode_x" not in frame
    mms = {"bulkv_gse_fast": raw["tha_peir_velocity_gsm"], "energyspectr_omni_fast": raw["tha_peir_en_eflux"]}
    frame = spedas_dict_to_frame(mms, layouts=_variable_layouts(mms, vectors={"bulkv_gse_fast"}, energy={"energyspectr_omni_fast"}))
    assert "bulkv_gse_fast_x" in frame and "energyspectr_omni_fast_energy" in frame


def test_time_axis_is_explicit_and_errors_include_variable_name():
    raw = {"alpha": {"x": [1., 2.], "y": np.arange(6).reshape(3, 2)}}
    spec = {"kind": "array", "columns": ["alpha"], "time_axis": 0}
    with pytest.raises(ValueError, match="alpha.*shape.*times"):
        spedas_dict_to_frame(raw, layouts={"alpha": spec})
    frame = spedas_dict_to_frame(raw, layouts={"alpha": {**spec, "time_axis": 1}})
    np.testing.assert_array_equal(frame.loc[0, "alpha"], [0, 2, 4])
    with pytest.raises(ValueError, match="alpha.*no product layout"):
        spedas_dict_to_frame(raw, layouts={})
    with pytest.raises(ValueError, match="alpha.*contain"):
        spedas_dict_to_frame({"alpha": {"x": [1]}}, layouts={})
    with pytest.raises(ValueError, match="alpha.*irregular"):
        spedas_dict_to_frame({"alpha": {"x": [1., 2.], "y": [[1], [2, 3]]}}, layouts={"alpha": spec})


def test_static_axis_length_equal_to_time_count_is_not_clipped():
    raw = {"s": {"x": [0., 1., 2.], "y": np.arange(9).reshape(3, 3), "v": np.array([10, 20, 30]), "v1": np.arange(6).reshape(3, 2)}}
    clipped = _clip_variable_dict(raw, [1., 1.])
    np.testing.assert_array_equal(clipped["s"]["v"], [10, 20, 30])
    np.testing.assert_array_equal(clipped["s"]["v1"], [[2, 3]])


def test_outer_join_duplicates_and_native_legacy_merge():
    a = {"density": {"x": [2., 1., 1.], "y": [20., 10., 99.]}}
    b = {"density": {"x": [2., 3.], "y": [200., np.nan]}}
    raw = _merge_variable_dicts([a, b])
    np.testing.assert_allclose(raw["density"]["y"], [10., 20., np.nan])
    assert raw["density"]["x"].size == 3
    data = {**a, "other": {"x": [3.], "y": [5.]}}
    frame = _hro_to_frame(data)
    assert frame["Timestamp"].tolist() == [1., 2., 3.]
    assert np.isnan(frame.loc[2, "density"])
    assert assemble({}).empty


def test_variable_reader_reads_only_requested_columns_and_concat_keeps_axes():
    def make(t, energy):
        raw = {"s": {"x": [t], "y": [[1., 2., 3.]], "v": energy, "attrs": {"array": np.arange(3)}}}
        return spedas_dict_to_frame(raw, layouts=_variable_layouts(raw))
    joined = concat([make(1., [1, 2, 3]), make(2., [4, 5, 6])])
    joined.attrs["variables"]["unrelated"] = {"kind": "scalar", "columns": ["absent"]}
    np.testing.assert_array_equal(read_variable(joined, "s")["axes"]["v"], [[1, 2, 3], [4, 5, 6]])


def test_loader_override_time_axis_is_used_during_clipping(monkeypatch):
    raw = {"custom": {"x": [0., 1., 2.], "y": np.arange(6).reshape(2, 3)}}
    layout = {"custom": {"kind": "array", "columns": ["custom"], "time_axis": 1}}
    monkeypatch.setattr(importlib.import_module("load.loaders.spedas"), "read_spedas", lambda *args, **kwargs: (raw, layout))
    frame = loader.load_spedas("themis", [0., 1.], "esd", variable_layouts=layout)
    assert frame["Timestamp"].tolist() == [0., 1.]
    np.testing.assert_array_equal(frame.loc[1, "custom"], [1, 4])
    legacy = loader.load_spedas("themis", [[0., 1.], [1., 2.]], "esd", variable_layouts=layout, return_type="dict")
    np.testing.assert_array_equal(legacy["custom"]["y"], raw["custom"]["y"])


def test_arbitrary_axis_names_and_sampled_coordinates():
    variable = {"time": [1., 2.], "values": np.ones((2, 3)), "kind": "array",
                "columns": ["power"], "time_axis": 0,
                "axes": {"frequency": {"column": "hz", "values": [[1, 2, 3], [4, 5, 6]], "time_axis": 0}}}
    frame = assemble({"signal": variable})
    np.testing.assert_array_equal(read_variable(frame, "signal")["axes"]["frequency"], [[1, 2, 3], [4, 5, 6]])
    variable["axes"]["frequency"]["values"] = [[1, 2, 3]]
    with pytest.raises(ValueError, match="signal.*frequency.*shape"):
        assemble({"signal": variable})


def test_legacy_merge_accepts_static_then_sampled_axis(monkeypatch):
    chunks = iter([
        {"custom": {"x": [0.], "y": [[1., 2.]], "v": [10., 20.]}},
        {"custom": {"x": [1.], "y": [[3., 4.]], "v": [[30., 40.]]}},
    ])
    def read_chunk(*args, **kwargs):
        raw = next(chunks)
        return raw, _variable_layouts(raw)

    monkeypatch.setattr(importlib.import_module("load.loaders.spedas"), "read_spedas", read_chunk)
    result = loader.load_spedas("themis", [[0., 0.], [1., 1.]], "esd", return_type="dict")
    np.testing.assert_array_equal(result["custom"]["v"], [[10., 20.], [30., 40.]])
