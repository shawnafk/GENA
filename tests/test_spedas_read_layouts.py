"""Read boundaries return data and matching layouts without network access."""
import importlib
import sys
from types import ModuleType

import pytest

spedas = importlib.import_module("load.api.spedas")


def mock_backend(monkeypatch, mission, raw):
    """Replace pySPEDAS modules, leaving the mission read functions exercised."""
    names = ["pyspedas", "pyspedas.projects", f"pyspedas.projects.{mission}"]
    config_name = f"pyspedas.projects.{mission}." + ("mms_config" if mission == "mms" else "config")
    names.append(config_name)
    if mission != "mms":
        names.append(f"pyspedas.projects.{mission}.load")
    modules = {name: ModuleType(name) for name in names}
    for name, module in modules.items():
        monkeypatch.setitem(sys.modules, name, module)
        parent, _, child = name.rpartition(".")
        if parent in modules:
            setattr(modules[parent], child, module)
    modules[config_name].CONFIG = {"remote_data_dir": "https://example.invalid/"}
    calls = []

    def read(**kwargs):
        calls.append(kwargs)
        return raw

    backend = modules[f"pyspedas.projects.{mission}"]
    if mission == "mms":
        backend.mms_load_fpi = read
        backend.mms_load_hpca = read
    else:
        modules[f"pyspedas.projects.{mission}.load"].load = read
    return calls


@pytest.mark.parametrize("mission,product,probe,name,backend", [
    ("themis", "esa", "b", "thb_peir_velocity_gsm", "_read_themis"),
    ("mms", "fpi", 2, "mms2_dis_bulkv_gse_fast", "_read_mms"),
    ("omni", "hro", None, "BZ_GSM", "_read_omni"),
])
def test_read_returns_matching_layouts(monkeypatch, mission, product, probe, name, backend):
    vector = mission != "omni"
    raw = {name: {"x": [0., 1.], "y": [[1, 2, 3], [4, 5, 6]] if vector else [1, 2]}}
    mock_backend(monkeypatch, mission, raw)
    data, layouts = spedas.read_spedas(mission, [0, 1], product=product, probe=probe, variables=[name], prefix="/tmp/spedas", download=mission == "omni")
    assert data == raw
    assert layouts[name]["kind"] == ("vector" if vector else "scalar")
    assert layouts[name]["columns"] == ([name + "_" + c for c in "xyz"] if vector else [name])
    assert spedas.parse_spedas(data, layouts=layouts)[name]["values"] is raw[name]["y"]


def test_read_applies_override_before_parsing(monkeypatch):
    raw = {"custom": {"x": [0., 1.], "y": [[1, 2], [3, 4], [5, 6]]}}
    mock_backend(monkeypatch, "themis", raw)
    override = {"custom": {"kind": "array", "time_axis": 1, "columns": ["custom"]}}
    _, layouts = spedas.read_spedas("themis", [0, 1], product="esd", prefix="/tmp/spedas", variable_layouts=override)
    assert layouts["custom"] == override["custom"]


def test_mms_loader_preserves_backend_names(monkeypatch):
    loader = importlib.import_module("load.loaders.spedas")
    name = "mms2_dis_bulkv_gse_fast"
    raw = {name: {"x": [0.], "y": [[1, 2, 3]]}}
    mock_backend(monkeypatch, "mms", raw)
    frame = loader.load_spedas("mms", [0, 1], "fpi", variables=["bulkv_gse_fast"], probe=2, prefix="/tmp/spedas")
    assert frame[name + "_x"].tolist() == [1]
    assert frame.attrs["variables"][name]["columns"] == [name + "_" + c for c in "xyz"]


@pytest.mark.parametrize("mission,product,probe,name", [
    ("themis", "esa", "c", "thc_peir_en_eflux"),
    ("themis", "sst", "a", "tha_psif_en_eflux"),
    ("mms", "fpi", 3, "mms3_dis_energyspectr_omni_fast"),
])
def test_energy_and_unknown_axes(monkeypatch, mission, product, probe, name):
    raw = {
        name: {"x": [0, 1], "y": [[1, 2], [3, 4]], "v": [10, 20]},
        "unknown": {"x": [0, 1], "y": [[1, 2, 3], [4, 5, 6]], "v1": [[1, 2, 3], [3, 4, 5]]},
    }
    calls = mock_backend(monkeypatch, mission, raw)
    data, layouts = spedas.read_spedas(
        mission, [0, 1], product=product, probe=probe,
        variables=list(raw), prefix="/tmp/spedas",
    )
    assert layouts[name]["axes"]["v"]["role"] == "energy"
    assert layouts[name]["axes"]["v"]["time_axis"] is None
    if mission == "themis":
        assert layouts["unknown"]["kind"] == "array"
        assert layouts["unknown"]["axes"]["v1"]["time_axis"] == 0
    assert calls[0]["no_update"] is (mission == "themis")


@pytest.mark.parametrize("product,name", [("fgm", "thd_fgs_gsm"), ("mom", "thd_peim_velocity_gsm")])
def test_themis_other_product_vectors(monkeypatch, product, name):
    raw = {name: {"x": [0], "y": [[1, 2, 3]]}}
    mock_backend(monkeypatch, "themis", raw)
    _, layouts = spedas.read_spedas("themis", [0, 1], product=product, probe="d", prefix="/tmp/spedas")
    assert layouts[name]["kind"] == "vector"


def test_mms_hpca_vector_and_empty_result(monkeypatch):
    name = "mms4_hpca_hplus_ion_bulk_velocity_GSM"
    raw = {name: {"x": [0], "y": [[1, 2, 3]]}}
    mock_backend(monkeypatch, "mms", raw)
    _, layouts = spedas.read_spedas("mms", [0, 1], product="hpca", probe=4, variables=[name], prefix="/tmp/spedas")
    assert layouts[name]["kind"] == "vector"
    mock_backend(monkeypatch, "mms", {})
    assert spedas.read_spedas("mms", [0, 1], product="hpca", prefix="/tmp/spedas") == ({}, {})
