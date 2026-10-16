"""The public facade follows the current API/loader split without doing I/O."""
import importlib

import pytest

package = importlib.import_module("load")
facade = importlib.import_module("load.load")


def test_public_exports_exist_and_point_to_current_owners():
    for module in (package, facade):
        assert all(hasattr(module, name) for name in module.__all__)
    assert facade.load_spedas.__module__ == "load.loaders.spedas"
    assert facade.load_indices.__module__ == "load.loaders.indices"
    assert facade.load_gena.__module__ == "load.loaders.gena"
    assert facade.read_spedas.__module__ == "load.api.spedas"
    assert facade.read_indices.__module__ == "load.api.indices"
    assert facade.get_themis_loc.__module__ == "load.load"
    assert not {"esa_mode2_angles", "themis_decompress", "read_themis_cdf"}.intersection(facade.__all__)


@pytest.mark.parametrize("mission", ["themis", "mms", "omni"])
def test_spedas_dispatch_forwards_options_once(monkeypatch, mission):
    calls = []
    sentinel = object()
    def read(*args, **kwargs):
        calls.append((args, kwargs))
        return sentinel
    monkeypatch.setattr(facade, "load_spedas", read)
    assert facade.load(mission.upper(), [0, 1], "product", probe=2, variables=["value"], download=False) is sentinel
    assert calls == [((mission, [0, 1], "product"), {"probe": 2, "variables": ["value"], "download": False})]


@pytest.mark.parametrize("mission", ["indices", "index", "kyoto"])
def test_indices_dispatch_uses_current_loader(monkeypatch, mission):
    sentinel = object()
    def read(trange, **kwargs):
        assert trange == [0, 1]
        assert kwargs == {"product": "dst", "download": False}
        return sentinel
    monkeypatch.setattr(facade, "load_indices", read)
    assert facade.load(mission, [0, 1], key="dst", download=False) is sentinel


def test_gena_dispatch(monkeypatch):
    def read(**kwargs):
        assert kwargs == {"trange": [0, 1], "product": "l2_a", "columns": ["Energy"]}
        return "gena"
    monkeypatch.setattr(facade, "load_gena", read)
    assert facade.load("gena", [0, 1], "l2_a", columns=["Energy"]) == "gena"


def test_unknown_source_is_rejected():
    with pytest.raises(ValueError, match="mission"):
        facade.load("unknown", [0, 1])
