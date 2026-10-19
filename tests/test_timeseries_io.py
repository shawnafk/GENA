import importlib
import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pytest

from common.frame import TimeSeriesDataFrame, concat
from load.api.public import _time_columns, assemble
from load.api.spedas import _variable_layouts, parse_spedas
from common.frame import split_time_intervals
from common.time import Time
from gena.frame import ParticleDataFrame


def _table(data, **kwargs):
    time = data.get("Datetime", data.get("Timestamp", []))
    return TimeSeriesDataFrame({**data, **_time_columns(time)}, **kwargs)


def _ena_table(data):
    return ParticleDataFrame(_table(data))


def _fixture_frame(data, **kwargs):
    layouts = _variable_layouts(data)
    # Explicit descriptions for this synthetic test dataset.
    if "velocity_gsm" in data:
        layouts["velocity_gsm"].update(kind="vector", columns=["velocity_gsm_x", "velocity_gsm_y", "velocity_gsm_z"])
    if "peir_en_eflux" in data:
        layouts["peir_en_eflux"]["axes"]["v"].update(column="peir_en_eflux_energy", role="energy")
    return assemble(parse_spedas(data, layouts=layouts), **kwargs)


def test_shared_frame_outer_joins_and_round_trips_multidimensional_data():
    raw = {
        "density": {"x": [1.0, 2.0], "y": [10.0, 20.0]},
        "velocity_gsm": {"x": [1.0, 3.0], "y": [[1, 2, 3], [4, 5, 6]]},
        "peir_en_eflux": {
            "x": [2.0, 3.0],
            "y": [[100, 200], [300, 400]],
            "v": [[10, 20], [10, 20]],
        },
    }

    frame = _fixture_frame(raw, source="test")

    assert isinstance(frame, TimeSeriesDataFrame)
    assert str(frame["Datetime"].dtype) == "datetime64[ns]"
    assert frame["Timestamp"].tolist() == [1.0, 2.0, 3.0]
    assert ["velocity_gsm_x", "velocity_gsm_y", "velocity_gsm_z"] == [
        name for name in frame if name.startswith("velocity_gsm_")
    ]
    np.testing.assert_array_equal(frame.loc[1, "peir_en_eflux"], [100, 200])
    np.testing.assert_array_equal(frame.loc[1, "peir_en_eflux_energy"], [10, 20])

    from common.frame import read_variable
    restored_velocity = read_variable(frame, "velocity_gsm")
    restored_energy = read_variable(frame, "peir_en_eflux")
    np.testing.assert_allclose(restored_velocity["values"], [[1, 2, 3], [4, 5, 6]])
    np.testing.assert_allclose(restored_energy["axes"]["v"], [[10, 20], [10, 20]])


def test_time_inputs_are_normalized_to_naive_pandas_datetime():
    inputs = [
        ["2025-01-01T00:00:00Z"],
        [1735689600.0],
        [pd.Timestamp("2025-01-01", tz="UTC")],
        [np.datetime64("2025-01-01")],
    ]
    for values in inputs:
        frame = _fixture_frame({"value": {"x": values, "y": [1.0]}})
        assert str(frame["Datetime"].dtype) == "datetime64[ns]"
        assert frame.loc[0, "Datetime"] == pd.Timestamp("2025-01-01")


def test_naive_datetime_is_converted_to_unix_seconds_as_utc():
    assert Time(pd.Timestamp("2025-01-01")).to_timestamp() == 1735689600.0
    np.testing.assert_allclose(
        Time(pd.date_range("2025-01-01", periods=2, freq="min")).to_timestamp(),
        [1735689600.0, 1735689660.0],
    )


def test_eff_h_clips_energy_to_valid_range():
    from gena.constants import eff_H

    np.testing.assert_allclose(eff_H([-10.0, 0.0]), eff_H([1.0, 1.0]))
    np.testing.assert_allclose(eff_H([61.0, 100.0]), eff_H([60.0, 60.0]))


def test_eff_o_clips_energy_to_valid_range():
    from gena.constants import eff_O

    np.testing.assert_allclose(eff_O([-10.0, 4.0]), eff_O([5.0, 5.0]))
    np.testing.assert_allclose(eff_O([61.0, 100.0]), eff_O([60.0, 60.0]))


def test_interval_helpers_preserve_subclasses():
    ena = _ena_table({"Timestamp": [0.0, 10.0, 100.0], "value": [1, 2, 3]})
    selected = ena.select(ena["value"] > 1)
    between = ena.between([pd.Timestamp("1970-01-01"), pd.Timestamp("1970-01-01 00:01")])[0]
    count, _, split = ena.split(tgap=20)

    assert isinstance(selected, ParticleDataFrame)
    assert isinstance(between, ParticleDataFrame)
    assert count == 2
    assert all(isinstance(part, ParticleDataFrame) for part in split)


def test_concat_preserves_earlier_input_at_duplicate_times():
    first = _table(
        {
            "Timestamp": [2.0, 1.0, 0.0],
            "value": ["first-2", "first-1", "first-0"],
        }
    )
    repeated_times = [2.0, 1.0, 0.0] * 9
    later = _table(
        {
            "Timestamp": repeated_times,
            "value": [f"later-{index}" for index in range(len(repeated_times))],
        }
    )

    result = concat([first, later])

    assert result["Timestamp"].tolist() == [0.0, 1.0, 2.0]
    assert result["value"].tolist() == ["first-0", "first-1", "first-2"]


def test_sync_builds_regular_grid_and_aligns_sparse_numeric_columns():
    frame = _table(
        {
            "Datetime": pd.to_datetime([0.1, 0.2, 4.1, 4.2, 8.1, 8.2], unit="s"),
            "a": [0.1, np.nan, 4.1, np.nan, 8.1, np.nan],
            "b": [np.nan, 0.2, np.nan, 4.2, np.nan, 8.2],
        }
    )
    frame.attrs = {"source": "test"}

    result = frame.resample(tres="4s")

    assert isinstance(result, TimeSeriesDataFrame)
    assert result["Timestamp"].tolist() == [4.0, 8.0]
    np.testing.assert_allclose(result["a"], [4.0, 8.0])
    np.testing.assert_allclose(result["b"], [4.0, 8.0])
    assert result.attrs == {"source": "test"}


def test_sync_to_another_frame_uses_reference_cadence_and_common_coverage():
    source = _table(
        {"Timestamp": [0.0, 5.0, 10.0, 15.0], "value": [0.0, 5.0, 10.0, 15.0]}
    )
    reference = _table(
        {"Timestamp": [4.0, 8.0, 12.0, 16.0], "reference": [1, 2, 3, 4]}
    )

    result = source.resample_to(reference)

    assert result["Timestamp"].tolist() == [4.0, 8.0, 12.0]
    np.testing.assert_allclose(result["value"], [4.0, 8.0, 12.0])


def test_sync_method_parameter_selects_linear_or_nearest():
    source = _table(
        {"Timestamp": [0.0, 10.0], "value": [0.0, 10.0]}
    )
    reference = _table({"Timestamp": [4.0, 6.0]})

    linear = source.resample_to(reference, method="linear", max_gap=None)
    nearest = source.resample_to(reference, method="nearest", max_gap=None)

    np.testing.assert_allclose(linear["value"], [4.0, 6.0])
    assert nearest["value"].tolist() == [0.0, 10.0]


def test_sync_nearest_cleans_invalid_samples_before_matching():
    source = _table(
        {"Timestamp": [0.0, 5.0, 10.0], "label": ["early", None, "late"]}
    )
    reference = _table({"Timestamp": [6.0]})

    result = source.resample_to(reference, method="nearest", max_gap=None)

    assert result.loc[0, "label"] == "late"


def test_sync_requires_date_on_plain_dataframe_reference():
    source = _table({"Timestamp": [0.0, 10.0], "value": [0.0, 10.0]})
    reference = pd.DataFrame({"Timestamp": [2.0, 8.0]})

    with pytest.raises(KeyError, match="Datetime"):
        source.resample_to(reference)

    reference["Datetime"] = pd.to_datetime(reference["Timestamp"], unit="s")
    np.testing.assert_allclose(source.resample_to(reference)["value"], [2.0, 8.0])


def test_time_operations_preserve_alignment_with_nondefault_row_index():
    frame = _table(
        {"Timestamp": [0.0, 5.0, 10.0], "value": [0.0, 5.0, 10.0]},
        index=[30, 10, 20],
    )
    reference = pd.DataFrame({"Datetime": pd.to_datetime([2.0, 8.0], unit="s")}, index=[7, 4])

    np.testing.assert_allclose(frame.resample_to(reference)["value"], [2.0, 8.0])
    selected = frame.between(pd.to_datetime([0.0, 10.0], unit="s"))[0]
    assert selected["value"].tolist() == [5.0]


def test_select_by_aligns_condition_but_preserves_original_target_rows():
    esa = _table(
        {
            "Timestamp": [0.0, 2.0, 4.0, 6.0],
            "temperature": [500.0, 1500.0, 500.0, 1500.0],
        }
    )
    spectra = np.empty(3, dtype=object)
    spectra[:] = [np.array([10.0, 11.0]), np.array([20.0, 21.0]), np.array([30.0, 31.0])]
    sst = _table(
        {
            "Timestamp": [1.9, 4.1, 5.9],
            "spectrum": spectra,
            "sample_id": [101, 102, 103],
        }
    )
    sst.attrs = {"source": "sst"}

    result = sst.select_by(esa, "temperature > 1000")

    assert isinstance(result, TimeSeriesDataFrame)
    assert result["sample_id"].tolist() == [101, 103]
    assert result["Timestamp"].tolist() == sst.loc[[0, 2], "Timestamp"].tolist()
    np.testing.assert_array_equal(result.loc[0, "spectrum"], sst.loc[0, "spectrum"])
    np.testing.assert_array_equal(result.loc[1, "spectrum"], sst.loc[2, "spectrum"])
    assert result.attrs == {"source": "sst"}


def test_select_by_respects_selector_gap_limit():
    selector = _table(
        {"Timestamp": [0.0, 100.0], "hot": [True, True]}
    )
    target = _table(
        {"Timestamp": [0.0, 50.0, 100.0], "value": [1.0, 2.0, 3.0]}
    )

    result = target.select_by(selector, "hot", max_gap="10s")

    assert result["value"].tolist() == [1.0, 3.0]


def test_sync_does_not_bridge_a_gap_larger_than_max_gap():
    source = _table(
        {"Timestamp": [0.0, 100.0], "value": [0.0, 100.0]}
    )
    reference = _table({"Timestamp": [0.0, 50.0, 100.0]})

    result = source.resample_to(reference, max_gap="10s")

    np.testing.assert_allclose(result.loc[[0, 2], "value"], [0.0, 100.0])
    assert np.isnan(result.loc[1, "value"])


def test_sync_auto_gap_preserves_outage_longer_than_native_cadence():
    source = _table(
        {
            "Timestamp": [0.0, 4.0, 8.0, 100.0, 104.0, 108.0],
            "value": [0.0, 4.0, 8.0, 100.0, 104.0, 108.0],
        }
    )

    result = source.resample(tres="4s")

    assert np.isnan(result.loc[result["Timestamp"] == 52.0, "value"]).all()
    np.testing.assert_allclose(
        result.loc[result["Timestamp"].isin([0.0, 4.0, 104.0, 108.0]), "value"],
        [0.0, 4.0, 104.0, 108.0],
    )


def test_split_time_intervals_sorts_native_datetimes_and_uses_strict_gap_boundary():
    intervals = split_time_intervals(
        pd.DatetimeIndex([
            "2025-01-01 02:00:00",
            "2025-01-01 00:30:00",
            "2025-01-01 00:00:00",
            "2025-01-01 01:00:00",
        ]),
        max_gap="30min",
    )

    assert intervals == [
        (pd.Timestamp("2025-01-01 00:00:00"), pd.Timestamp("2025-01-01 01:00:00")),
        (pd.Timestamp("2025-01-01 02:00:00"), pd.Timestamp("2025-01-01 02:00:00")),
    ]
    assert split_time_intervals(pd.DatetimeIndex([]), max_gap="30min") == []


def test_wide_frame_repr_does_not_reconstruct_duplicate_time_columns():
    values = {f"value_{index}": np.arange(3) for index in range(24)}
    frame = _table({"Timestamp": [0.0, 1.0, 2.0], **values})

    rendered = repr(frame)
    html = frame._repr_html_()

    assert "value_23" in rendered
    assert "datetime64" in html or "Datetime" in html
    assert frame.columns.is_unique


def test_dataframe_has_builtin_tplot_for_scalar_and_vector_variables():
    frame = _fixture_frame(
        {
            "density": {"x": [0.0, 1.0], "y": [1.0, 2.0]},
            "velocity_gsm": {"x": [0.0, 1.0], "y": [[1, 2, 3], [4, 5, 6]]},
        }
    )

    scalar_figure, scalar_axis = frame.tplot("density")
    vector_figure, vector_axis = frame.tplot("velocity_gsm")

    assert len(scalar_axis.lines) == 1
    assert len(vector_axis.lines) == 3
    plt.close(scalar_figure)
    plt.close(vector_figure)


def test_tplot_has_one_public_implementation_outside_vis():
    from plots import tplot as top_level
    from plots.basic import tplot
    from plots.basic.timeseries import tplot as implementation
    from plots import vis

    assert tplot is implementation is top_level
    assert not hasattr(vis, "tplot")
    assert not hasattr(vis, "dataset_variable")
    assert not hasattr(vis, "resolve_dataset_key")


def _scalar_variables(names, start=0.0, end=60.0):
    return {name: {"x": np.array([start, end]), "y": np.ones((2, 3)) if name.endswith(("fgs_gsm", "velocity_gsm", "_flux", "bulkv_gse_fast")) else np.array([1.0, 2.0])} for name in names}


def test_themis_defaults_are_field_filtered_and_offline(monkeypatch):
    module = importlib.import_module("load.load")
    backend = importlib.import_module("pyspedas.projects.themis.load")
    calls = []

    def fake_load(**kwargs):
        calls.append(kwargs)
        raw = _scalar_variables(kwargs["varnames"])
        for name, variable in raw.items():
            if name.endswith(("velocity_gsm", "_flux", "bulkv_gse_fast")):
                variable["y"] = np.ones((2, 3))
        return raw

    monkeypatch.setattr(backend, "load", fake_load)
    monkeypatch.setattr(importlib.import_module("load.api.spedas"), "get_configured_paths", lambda: {"spedas_base": "/tmp/spedas"})

    result = module.load_spedas("themis", 
        ["1970-01-01", "1970-01-01 00:01"], "esa", probe="b"
    )

    assert isinstance(result, TimeSeriesDataFrame)
    assert len(calls) == 1
    assert calls[0]["no_update"] is True
    assert calls[0]["notplot"] is True
    assert len(calls[0]["varnames"]) == 12
    assert "thb_peir_en_eflux" in calls[0]["varnames"]
    assert "thb_peer_velocity_gsm" in calls[0]["varnames"]


def test_themis_esd_uses_native_names_and_preserves_three_dimensional_axes(monkeypatch):
    module = importlib.import_module("load.load")
    backend = importlib.import_module("pyspedas.projects.themis.load")
    calls = []

    def fake_load(**kwargs):
        calls.append(kwargs)
        time = np.array([1.0, 5.0, 9.0])
        cube = np.arange(24.0).reshape(3, 2, 4)
        return {
            "valid": {"x": time, "y": np.ones(3, dtype=np.uint8)},
            "data_quality": {"x": time, "y": np.zeros(3, dtype=np.int16)},
            "bins": {"x": time, "y": np.ones_like(cube, dtype=np.uint8), "v1": [0, 1], "v2": [0, 1, 2, 3]},
            "eff": {"x": time, "y": np.ones_like(cube), "v1": [0, 1], "v2": [0, 1, 2, 3]},
            "eflux": {"x": time, "y": cube, "v1": [0, 1], "v2": [0, 1, 2, 3]},
        }

    monkeypatch.setattr(backend, "load", fake_load)
    monkeypatch.setattr(importlib.import_module("load.api.spedas"), "get_configured_paths", lambda: {"spedas_base": "/tmp/spedas"})

    frame = module.load_spedas("themis", 
        ["1970-01-01", "1970-01-01 00:00:12"], "esd", probe="b", datatype="peif"
    )

    assert isinstance(frame, TimeSeriesDataFrame)
    assert calls[0]["instrument"] == "esd"
    assert calls[0]["datatype"] == "peif"
    assert calls[0]["varnames"] == ["valid", "data_quality", "bins", "eff", "eflux"]
    assert frame.attrs["datatype"] == "peif"
    np.testing.assert_array_equal(frame.loc[1, "eflux"], np.arange(8.0, 16.0).reshape(2, 4))
    np.testing.assert_array_equal(frame.loc[1, "eflux_v1"], [0, 1])
    np.testing.assert_array_equal(frame.loc[1, "eflux_v2"], [0, 1, 2, 3])
    from common.frame import read_variable
    restored = read_variable(frame, "eflux")
    assert restored["values"].shape == (3, 2, 4)
    assert restored["axes"]["v1"].shape == (2,)
    assert restored["axes"]["v2"].shape == (4,)


def test_themis_esd_rejects_non_l2_levels():
    module = importlib.import_module("load.load")

    with pytest.raises(ValueError, match="only at level='l2'"):
        module.load_spedas("themis", ["1970-01-01", "1970-01-01 00:01"], "esd", level="l1")


@pytest.mark.parametrize("probe", ["b", "c"])
def test_themis_esa_preserves_each_variables_native_samples(monkeypatch, probe):
    module = importlib.import_module("load.load")
    backend = importlib.import_module("pyspedas.projects.themis.load")
    root = f"th{probe}_"

    def fake_load(**kwargs):
        return {
            f"{root}peir_density": {
                "x": np.array([1.0, 5.0, 9.0]),
                "y": np.array([0.1, 0.2, 0.3]),
            },
            f"{root}peer_density": {
                "x": np.array([0.0, 4.0, 8.0, 12.0]),
                "y": np.array([0.0, 4.0, 8.0, 12.0]),
            },
        }

    monkeypatch.setattr(backend, "load", fake_load)
    monkeypatch.setattr(importlib.import_module("load.api.spedas"), "get_configured_paths", lambda: {"spedas_base": "/tmp/spedas"})

    result = module.load_spedas("themis", 
        ["1970-01-01", "1970-01-01 00:00:12"],
        "esa",
        probe=probe,
        variables=["peir_density", "peer_density"],
    )

    assert result["Timestamp"].tolist() == [0., 1., 4., 5., 8., 9., 12.]
    ion = result.loc[result[f"{root}peir_density"].notna()]
    electron = result.loc[result[f"{root}peer_density"].notna()]
    assert ion["Timestamp"].tolist() == [1., 5., 9.]
    assert electron["Timestamp"].tolist() == [0., 4., 8., 12.]
    np.testing.assert_array_equal(ion[f"{root}peir_density"], [0.1, 0.2, 0.3])
    np.testing.assert_array_equal(electron[f"{root}peer_density"], [0., 4., 8., 12.])
    assert "time_reference" not in result.attrs
    assert "time_alignment" not in result.attrs


def test_themis_local_first_miss_falls_back_once(monkeypatch):
    module = importlib.import_module("load.load")
    backend = importlib.import_module("pyspedas.projects.themis.load")
    calls = []

    def fake_load(**kwargs):
        calls.append(kwargs["no_update"])
        return {} if kwargs["no_update"] else _scalar_variables(kwargs["varnames"])

    monkeypatch.setattr(backend, "load", fake_load)
    monkeypatch.setattr(importlib.import_module("load.api.spedas"), "get_configured_paths", lambda: {"spedas_base": "/tmp/spedas"})
    module.load_spedas("themis", 
        ["1970-01-01", "1970-01-01 00:01"],
        "fgm",
        probe="b",
        download="local_first",
    )
    assert calls == [True, False]


def test_themis_default_local_miss_never_contacts_remote(monkeypatch):
    module = importlib.import_module("load.load")
    backend = importlib.import_module("pyspedas.projects.themis.load")
    calls = []

    def fake_load(**kwargs):
        calls.append(kwargs["no_update"])
        return {}

    monkeypatch.setattr(backend, "load", fake_load)
    monkeypatch.setattr(importlib.import_module("load.api.spedas"), "get_configured_paths", lambda: {"spedas_base": "/tmp/spedas"})
    module.load_spedas("themis", ["1970-01-01", "1970-01-01 00:01"], "fgm", probe="b")

    assert calls == [True]


@pytest.mark.parametrize("download_mode", [True, "download"])
def test_themis_download_mode_skips_offline_probe(monkeypatch, download_mode):
    module = importlib.import_module("load.load")
    backend = importlib.import_module("pyspedas.projects.themis.load")
    calls = []

    def fake_load(**kwargs):
        calls.append(kwargs)
        raw = _scalar_variables(kwargs["varnames"])
        for name, variable in raw.items():
            if name.endswith(("velocity_gsm", "_flux", "bulkv_gse_fast")):
                variable["y"] = np.ones((2, 3))
        return raw

    monkeypatch.setattr(backend, "load", fake_load)
    monkeypatch.setattr(importlib.import_module("load.api.spedas"), "get_configured_paths", lambda: {"spedas_base": "/tmp/spedas"})
    module.load_spedas("themis", 
        ["1970-01-01", "1970-01-01 00:01"],
        "fgm",
        probe="b",
        download=download_mode,
    )

    assert len(calls) == 1
    assert calls[0]["no_update"] is False
    assert calls[0]["force_download"] is False


def test_themis_variables_replaces_defaults(monkeypatch):
    module = importlib.import_module("load.load")
    backend = importlib.import_module("pyspedas.projects.themis.load")
    calls = []

    def fake_load(**kwargs):
        calls.append(kwargs["varnames"])
        raw = _scalar_variables(kwargs["varnames"])
        for name, variable in raw.items():
            if name.endswith(("velocity_gsm", "_flux", "bulkv_gse_fast")):
                variable["y"] = np.ones((2, 3))
        return raw

    monkeypatch.setattr(backend, "load", fake_load)
    monkeypatch.setattr(importlib.import_module("load.api.spedas"), "get_configured_paths", lambda: {"spedas_base": "/tmp/spedas"})

    module.load_spedas("themis", 
        ["1970-01-01", "1970-01-01 00:01"], "esa", probe="b",
        variables=["peir_density", "thb_peir_avgtemp"],
    )
    module.load_spedas("themis", 
        ["1970-01-01", "1970-01-01 00:01"], "esa", probe="b",
        variables=["peer_density"],
    )

    assert calls == [
        ["thb_peir_density", "thb_peir_avgtemp"],
        ["thb_peer_density"],
    ]


def test_mms_uses_local_backend_and_selected_defaults(monkeypatch):
    module = importlib.import_module("load.load")
    backend = importlib.import_module("pyspedas.projects.mms")
    calls = []

    def fake_fpi(**kwargs):
        calls.append(kwargs)
        raw = _scalar_variables(kwargs["varnames"])
        for name, variable in raw.items():
            if name.endswith(("velocity_gsm", "_flux", "bulkv_gse_fast")):
                variable["y"] = np.ones((2, 3))
        return raw

    monkeypatch.setattr(backend, "mms_load_fpi", fake_fpi)

    result = module.load_spedas("mms", 
        ["1970-01-01", "1970-01-01 00:01"], "fpi", probe=1, download=False
    )

    assert isinstance(result, TimeSeriesDataFrame)
    assert len(calls) == 1
    assert calls[0]["spdf"] is False
    assert calls[0]["no_update"] is True
    assert calls[0]["notplot"] is True
    assert calls[0]["varnames"] == [
        "mms1_dis_temppara_fast",
        "mms1_dis_tempperp_fast",
        "mms1_dis_numberdensity_fast",
        "mms1_dis_bulkv_gse_fast",
    ]
    assert "mms1_dis_numberdensity_fast" in result


def test_omni_and_indices_default_to_frames_with_legacy_escape_hatch(monkeypatch):
    module = importlib.import_module("load.load")
    captured = []

    def fake_omni(mission, trange, **kwargs):
        names = kwargs["variables"] or module.OMNI_VARIABLES
        captured.extend(names)
        raw = _scalar_variables(names)
        layouts = {
            name: {"kind": "scalar", "time_axis": 0, "columns": [name]}
            for name in raw
        }
        return raw, layouts

    monkeypatch.setattr(importlib.import_module("load.loaders.spedas"), "read_spedas", fake_omni)
    def fake_indices(trange, product=None, **kwargs):
        names = ["ae", "al", "au"] if product is None else [product]
        raw = {name: {"x": [0.0], "y": [-10.0]} for name in names}
        layouts = {name: {"kind": "scalar", "time_axis": 0, "columns": [name]} for name in names}
        return raw, layouts

    monkeypatch.setattr(importlib.import_module("load.loaders.indices"), "read_indices", fake_indices)

    omni_frame = module.load_spedas("omni", ["1970-01-01", "1970-01-01 00:01"])
    auroral_frame = module.load_indices(["1970-01-01", "1970-01-01 00:01"])
    auroral_dict = module.load_indices(
        ["1970-01-01", "1970-01-01 00:01"], return_type="dict"
    )
    dst_frame = module.load_indices(["1970-01-01", "1970-01-01 00:01"], "dst")
    dst_dict = module.load_indices(
        ["1970-01-01", "1970-01-01 00:01"], "dst", return_type="dict"
    )

    assert captured == ["BZ_GSM", "BY_GSM", "flow_speed", "proton_density", "Pressure"]
    assert isinstance(omni_frame, TimeSeriesDataFrame)
    assert isinstance(auroral_frame, TimeSeriesDataFrame)
    assert {"ae", "al", "au"}.issubset(auroral_frame.columns)
    assert set(auroral_dict) == {"ae", "al", "au"}
    assert isinstance(dst_frame, TimeSeriesDataFrame)
    assert "dst" in dst_frame
    assert set(dst_dict) == {"x", "y"}


def test_omni_download_uses_current_pyspedas_namespace(monkeypatch, tmp_path):
    module = importlib.import_module("load.load")
    backend = importlib.import_module("pyspedas.projects.omni.load")
    config = importlib.import_module("pyspedas.projects.omni.config")
    calls = []

    def fake_load(**kwargs):
        calls.append(kwargs)
        raw = _scalar_variables(kwargs["varnames"])
        for name, variable in raw.items():
            if name.endswith(("velocity_gsm", "_flux", "bulkv_gse_fast")):
                variable["y"] = np.ones((2, 3))
        return raw

    monkeypatch.setattr(backend, "load", fake_load)
    monkeypatch.setitem(config.CONFIG, "local_data_dir", config.CONFIG["local_data_dir"])

    result = module.load_spedas("omni", 
        ["1970-01-01", "1970-01-01 00:01"],
        prefix=str(tmp_path),
        download=True,
        variables=["BZ_GSM"],
        return_type="dict",
    )

    assert len(calls) == 1
    assert calls[0]["level"] == "hro"
    assert calls[0]["notplot"] is True
    assert calls[0]["varnames"] == ["BZ_GSM"]
    assert config.CONFIG["local_data_dir"] == str(tmp_path) + "/omni_data/"
    assert result["BZ_GSM"]["x"].tolist() == list(
        pd.to_datetime([0.0, 60.0], unit="s")
    )


def test_omni_download_bypasses_loopback_proxy(monkeypatch, tmp_path):
    module = importlib.import_module("load.load")
    backend = importlib.import_module("pyspedas.projects.omni.load")
    config = importlib.import_module("pyspedas.projects.omni.config")
    requests = importlib.import_module("requests")
    original_no_proxy = {key: os.environ.get(key) for key in ("NO_PROXY", "no_proxy")}
    calls = []

    def fake_load(**kwargs):
        calls.append({key: os.environ.get(key) for key in ("NO_PROXY", "no_proxy")})
        raw = _scalar_variables(kwargs["varnames"])
        for name, variable in raw.items():
            if name.endswith(("velocity_gsm", "_flux", "bulkv_gse_fast")):
                variable["y"] = np.ones((2, 3))
        return raw

    monkeypatch.setattr(backend, "load", fake_load)
    monkeypatch.setattr(
        requests.utils,
        "get_environ_proxies",
        lambda url: {"https": "http://127.0.0.1:8080"},
    )
    monkeypatch.setitem(config.CONFIG, "local_data_dir", config.CONFIG["local_data_dir"])

    result = module.load_spedas("omni", 
        ["1970-01-01", "1970-01-01 00:01"],
        prefix=str(tmp_path),
        download=True,
        variables=["BZ_GSM"],
        return_type="dict",
    )

    assert len(calls) == 1
    assert "spdf.gsfc.nasa.gov" in calls[0]["NO_PROXY"]
    assert "cdaweb.gsfc.nasa.gov" in calls[0]["no_proxy"]
    assert {key: os.environ.get(key) for key in original_no_proxy} == original_no_proxy
    assert result["BZ_GSM"]["y"].tolist() == [1.0, 2.0]


def test_event_rows_are_not_interpolated_but_can_be_selected_and_counted():
    events = _ena_table({"Timestamp": [0., 1., 1., 2.], "Energy": [1., 2., 3., 4.]})
    assert not isinstance(events, TimeSeriesDataFrame)
    assert not hasattr(events, "resample_to")
    assert not hasattr(events, "resample_at")
    # Inherited pandas resample aggregates existing events rather than making new rows.
    assert events.resample("1s", on="Datetime").size().sum() == 4
    selector = _table({"Timestamp": [0., 1., 2.], "hot": [False, True, False]})
    selected = events.select_by(selector, "hot", max_gap=None)
    assert isinstance(selected, ParticleDataFrame)
    assert selected["Energy"].tolist() == [2., 3.]


def test_resample_interfaces_have_distinct_targets():
    source = _table({"Timestamp": [0., 4.], "value": [0., 8.]})
    regular = source.resample("2s", max_gap=None)
    np.testing.assert_allclose(regular["value"], [0., 4., 8.])
    reference = _table({"Timestamp": [1., 3.]})
    np.testing.assert_allclose(source.resample_to(reference, max_gap=None)["value"], [2., 6.])
    with pytest.raises(TypeError):
        source.resample_to(tres="2s")


def test_resample_at_preserves_requested_axis_and_normalizes_utc():
    source = _table({"Timestamp": [0., 4.], "value": [0., 8.]})
    source.attrs["units"] = "test"
    result = source.resample_at([3., -1., 1., 3., 5.], max_gap=None)
    assert isinstance(result, TimeSeriesDataFrame)
    assert result["Timestamp"].tolist() == [3., -1., 1., 3., 5.]
    np.testing.assert_allclose(result["value"], [6., np.nan, 2., 6., np.nan])
    assert result.attrs == source.attrs
    zoned = source.resample_at(["1970-01-01T08:00:01+08:00"], max_gap=None)
    assert zoned["value"].tolist() == [2.]
    assert zoned["Datetime"].dt.tz is None
    with pytest.raises(ValueError):
        source.resample_at([pd.NaT])
    with pytest.raises(ValueError):
        source.resample_at([], max_gap="0s")
    with pytest.raises(ValueError):
        source.resample_at([], method="cubic")
    empty = source.iloc[:0].resample_at([1., 2.])
    assert empty["Timestamp"].tolist() == [1., 2.]
    assert empty["value"].isna().all()


def test_split_time_intervals_does_not_convert_through_time(monkeypatch):
    import common.frame as table_frames

    def unexpected_conversion(*args, **kwargs):
        raise AssertionError("pandas datetime grouping must not construct Time")

    monkeypatch.setattr(table_frames, "Time", unexpected_conversion)
    dates = pd.Series(pd.to_datetime([21., 0., 10., 10.], unit="s"))
    expected = [(pd.Timestamp(0, unit="s"), pd.Timestamp(10, unit="s")),
                (pd.Timestamp(21, unit="s"), pd.Timestamp(21, unit="s"))]
    assert table_frames.split_time_intervals(dates, max_gap="10s") == expected
    zoned = pd.DatetimeIndex(dates).tz_localize("UTC")
    intervals = table_frames.split_time_intervals(zoned, max_gap="10s")
    assert intervals[0][0].tzinfo is not None
    assert intervals[0][0] == expected[0][0].tz_localize("UTC")


@pytest.mark.parametrize("values", [[0., 1.], ["2025-01-01"], pd.Series([0., 1.]), pd.Series(["2025-01-01"])])
def test_split_time_intervals_requires_native_pandas_datetime(values):
    with pytest.raises(TypeError, match="pandas datetime"):
        split_time_intervals(values)
    with pytest.raises(ValueError, match="NaT"):
        split_time_intervals(pd.DatetimeIndex([pd.NaT]))


def test_moon_blocking_intervals_uses_native_pandas_grouping(monkeypatch):
    from gena.post import project

    mask = pd.DataFrame({"Datetime": pd.to_datetime([0., 10., 21.], unit="s"),
                         "partial_block": [True, True, True], "full_block": [True, False, True]})
    monkeypatch.setattr(project, "moon_blocking_mask", lambda *args, **kwargs: mask)
    partial = project.moon_blocking_intervals(pd.DataFrame(), max_gap="10s")
    assert partial == [(mask.loc[0, "Datetime"], mask.loc[1, "Datetime"]),
                       (mask.loc[2, "Datetime"], mask.loc[2, "Datetime"])]
    full = project.moon_blocking_intervals(pd.DataFrame(), max_gap="10s", full_disk=True)
    assert full == [(mask.loc[0, "Datetime"], mask.loc[0, "Datetime"]),
                    (mask.loc[2, "Datetime"], mask.loc[2, "Datetime"])]
