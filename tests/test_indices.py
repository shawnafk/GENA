"""Kyoto/GFZ API and index loader boundaries, with no network in tests."""
from copy import deepcopy
import importlib
import json

import numpy as np
import pandas as pd
import pytest
import requests

from common.frame import TimeSeriesDataFrame, read_variable
from load import load, load_indices, read_indices, parse_indices
from load.api.indices import INDICES_FIELDS

api = importlib.import_module("load.api.indices")
loader = importlib.import_module("load.loaders.indices")


def kyoto_payload(dataset="min_ae", times=None):
    fields = {"min_ae": ["aeValue", "alValue", "auValue", "aoValue"],
              "hour_dst": ["dstValue"],
              "min_asysym": ["asyDValue", "asyHValue", "symDValue", "symHValue"]}[dataset]
    if times is None:
        times = ["2025-01-01T00:29:30Z", "2025-01-01T01:29:30Z"] if dataset == "hour_dst" else ["2025-01-01T00:00:00Z", "2025-01-01T00:01:00Z"]
    params = [{"name": "Time", "units": "UTC"}] + [{"name": name, "units": "nT", "fill": "99999"} for name in fields] + [{"name": "versionCode"}]
    rows = [[time, *([-20.] if dataset == "hour_dst" else [100., -50., 50., 0.]), 20] for time in times]
    return {"status": {"code": 1200}, "parameters": params, "data": rows,
            "title": dataset, "cadence": "PT1H" if dataset == "hour_dst" else "PT1M",
            "timeStampLocation": "center", "resourceID": "doi:fixture"}


@pytest.fixture
def remote(monkeypatch):
    calls = []

    def get(url, *, params, timeout):
        calls.append((url, deepcopy(params), timeout))
        if "kyoto" in url:
            payload = kyoto_payload(params["id"])
        else:
            payload = {"datetime": ["2025-01-01T00:00:00Z", "2025-01-01T03:00:00Z"],
                       params["index"]: [100., -1.], "meta": {"source": "GFZ Potsdam"}, "status": ["def", "def"]}
        class Response:
            def raise_for_status(self):
                pass

            def json(self):
                return deepcopy(payload)
        return Response()

    monkeypatch.setattr(requests, "get", get)
    return calls


def write_cache(root, source, dataset, payload, month="202501"):
    path = root / "indices" / source / dataset / f"{month}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload))
    return path


def test_ae_default_group_uses_one_request_and_preserves_layout(remote, tmp_path):
    data, layouts = read_indices(["2025-01-01", "2025-01-01T00:01"], prefix=tmp_path)
    assert list(data) == ["ae", "al", "au"]
    assert len(remote) == 1
    assert remote[0][1]["id"] == "min_ae"
    assert remote[0][1]["include"] == "header"
    assert data["al"]["y"].tolist() == [-50., -50.]
    assert layouts["ae"] == {"kind": "scalar", "time_axis": 0, "columns": ["ae"]}
    assert parse_indices(data, layouts)["ae"]["attrs"]["processing_versions_in_source"] == [20]
    # A different AE component reuses the same full dataset cache.
    read_indices(["2025-01-01", "2025-01-02"], "ao", prefix=tmp_path, download=False)
    assert len(remote) == 1


def test_dst_native_centers_and_missing_fill(tmp_path):
    payload = kyoto_payload("hour_dst")
    payload["data"][1][1] = 99999
    write_cache(tmp_path, "kyoto", "hour_dst", payload)
    frame = load_indices(["2025-01-01", "2025-01-01T02:00"], "dst", prefix=tmp_path, download=False)
    assert isinstance(frame, TimeSeriesDataFrame)
    assert frame["Datetime"].tolist() == list(pd.to_datetime(["2025-01-01T00:29:30", "2025-01-01T01:29:30"]))
    assert frame.loc[0, "dst"] == -20
    assert np.isnan(frame.loc[1, "dst"])
    variable = read_variable(frame, "dst")
    assert variable["units"] == "nT"
    assert variable["attrs"]["timeStampLocation"] == "center"


def test_mixed_sources_and_cadences_outer_join(remote, tmp_path):
    frame = load_indices(["2025-01-01", "2025-01-01T04:00"], ["ae", "dst", "kp", "f107"], prefix=tmp_path)
    assert len(remote) == 4
    assert set(frame.columns) == {"Datetime", "Timestamp", "ae", "dst", "kp", "f107"}
    assert frame["Datetime"].is_monotonic_increasing
    at_dst = frame[frame["Datetime"] == pd.Timestamp("2025-01-01T00:29:30")].iloc[0]
    assert at_dst["dst"] == -20 and np.isnan(at_dst["ae"])
    assert frame.attrs["variables"]["f107"]["units"] == "sfu"
    assert frame.attrs["variables"]["kp"]["attrs"]["source"] == "gfz"
    assert frame["kp"].notna().sum() == 1


@pytest.mark.parametrize("product,field", [("kp", "Kp"), ("ap", "ap"), ("ap_daily", "Ap"), ("sn", "SN"), ("f107", "Fobs"), ("f107adj", "Fadj")])
def test_gfz_field_selection_and_native_time(remote, tmp_path, product, field):
    data, _ = read_indices(["2025-01-01", "2025-01-01T03:00"], product, prefix=tmp_path)
    assert remote[0][1]["index"] == field
    assert data[product]["x"][0] == pd.Timestamp("2025-01-01")
    assert np.isnan(data[product]["y"][1])
    if product in {"ap", "ap_daily"}:
        assert data[product]["units"] == "2 nT"


def test_sym_asy_share_kyoto_request(remote, tmp_path):
    products = ["sym_h", "sym_d", "asy_h", "asy_d"]
    data, _ = read_indices(["2025-01-01", "2025-01-02"], products, prefix=tmp_path)
    assert len(remote) == 1 and remote[0][1]["id"] == "min_asysym"
    assert data["sym_h"]["y"].tolist() == [0., 0.]
    assert set(data) == set(products)


@pytest.mark.parametrize("download", [False, "offline", "local"])
def test_offline_miss_never_requests_or_creates_directory(remote, tmp_path, download):
    with pytest.raises(FileNotFoundError, match="cache"):
        read_indices(["2025-01-01", "2025-01-02"], "dst", prefix=tmp_path, download=download)
    assert remote == []
    assert not (tmp_path / "indices").exists()


def test_cache_reuse_and_refresh(remote, tmp_path):
    args = (["2025-01-01", "2025-01-02"], "dst")
    read_indices(*args, prefix=tmp_path)
    read_indices(*args, prefix=tmp_path, download=True)
    read_indices(*args, prefix=tmp_path, download=False)
    assert len(remote) == 1
    read_indices(*args, prefix=tmp_path, download="refresh")
    assert len(remote) == 2


def test_bad_refresh_preserves_existing_cache(monkeypatch, tmp_path):
    path = write_cache(tmp_path, "kyoto", "hour_dst", kyoto_payload("hour_dst"))
    original = path.read_bytes()
    class Response:
        def raise_for_status(self):
            pass
        def json(self):
            return {"status": {"code": 1400, "message": "bad request"}}
    monkeypatch.setattr(requests, "get", lambda *args, **kwargs: Response())
    with pytest.raises(ValueError, match="HAPI error"):
        read_indices(["2025-01-01", "2025-01-02"], "dst", prefix=tmp_path, download="refresh")
    assert path.read_bytes() == original
    assert not list(path.parent.glob("*.tmp"))


def test_http_failure_is_not_cached(monkeypatch, tmp_path):
    def fail(*args, **kwargs):
        raise requests.Timeout("fixture timeout")
    monkeypatch.setattr(requests, "get", fail)
    with pytest.raises(requests.Timeout):
        read_indices(["2025-01-01", "2025-01-02"], "dst", prefix=tmp_path)
    assert not (tmp_path / "indices").exists()


def test_multiple_months_and_disjoint_intervals_keep_native_samples(tmp_path):
    write_cache(tmp_path, "kyoto", "hour_dst", kyoto_payload("hour_dst", ["2025-01-31T23:29:30Z"]), "202501")
    write_cache(tmp_path, "kyoto", "hour_dst", kyoto_payload("hour_dst", ["2025-02-01T00:29:30Z", "2025-02-01T01:29:30Z"]), "202502")
    trange = [["2025-01-31T23:00", "2025-02-01T00:30"], ["2025-02-01T00:00", "2025-02-01T02:00"]]
    frame = load_indices(trange, "dst", prefix=tmp_path, download=False)
    assert len(frame) == 3 and frame["Datetime"].is_unique
    legacy = load_indices(trange, "dst", prefix=tmp_path, download=False, return_type="dict")
    assert set(legacy) >= {"x", "y"} and len(legacy["x"]) == 3


def test_legacy_single_multi_and_top_level_dispatch(remote, tmp_path):
    trange = ["2025-01-01", "2025-01-02"]
    one = load_indices(trange, "dst", prefix=tmp_path, return_type="dict")
    assert one["y"].tolist() == [-20., -20.]
    multiple = load_indices(trange, prefix=tmp_path, return_type="dict")
    assert set(multiple) == {"ae", "al", "au"}
    frame = load("indices", trange, "dst", prefix=tmp_path, download=False)
    assert frame["dst"].tolist() == one["y"].tolist()


def test_aliases_and_duplicate_products_use_canonical_name(remote, tmp_path):
    data, _ = read_indices(["2025-01-01", "2025-01-02"], ["F107OBS", "f107"], prefix=tmp_path)
    assert list(data) == ["f107"] and len(remote) == 1


@pytest.mark.parametrize("product", [[], "unknown", ["ae", "unknown"]])
def test_invalid_product_rejected_before_io(remote, tmp_path, product):
    with pytest.raises(ValueError):
        load_indices(["2025-01-01", "2025-01-02"], product, prefix=tmp_path)
    assert remote == []


@pytest.mark.parametrize("trange", ["2025-01-01", [], [0], [0, 1, 2], [1, 0]])
def test_invalid_range_rejected_before_io(remote, tmp_path, trange):
    with pytest.raises(ValueError):
        read_indices(trange, "dst", prefix=tmp_path)
    assert remote == []


def test_empty_kyoto_response_retains_typed_columns(tmp_path):
    payload = kyoto_payload("hour_dst", [])
    payload["status"]["code"] = 1201
    write_cache(tmp_path, "kyoto", "hour_dst", payload)
    frame = load_indices(["2025-01-01", "2025-01-02"], "dst", prefix=tmp_path, download=False)
    assert frame.empty and "dst" in frame.columns
    assert str(frame["Datetime"].dtype) == "datetime64[ns]"


def test_parsers_reject_shape_errors_and_copy_metadata():
    item = {"x": [0.], "y": [3.], "units": "nT", "attrs": {"source": "fixture"}}
    layout = {"dst": {"kind": "scalar", "time_axis": 0, "columns": ["dst"]}}
    parsed = parse_indices({"dst": item}, layout)
    parsed["dst"]["attrs"]["source"] = "changed"
    assert item["attrs"]["source"] == "fixture"
    with pytest.raises(ValueError, match="timestamps"):
        parse_indices({"dst": {**item, "y": [[3.]]}}, layout)
    with pytest.raises(ValueError, match="rows"):
        api._kyoto_records({**kyoto_payload(), "data": [["2025-01-01"]]}, ["ae"])
    with pytest.raises(ValueError, match="timestamps"):
        api._gfz_records({"datetime": [], "Kp": [1.]}, ["kp"])
