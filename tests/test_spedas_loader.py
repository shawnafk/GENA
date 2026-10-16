"""Offline checks for the consolidated mission-loading pipeline."""
import importlib
import sys
from types import ModuleType

import numpy as np
import pandas as pd
import pytest

from load.api.spedas import OMNI_FIELDS, default_variables, _variable_layouts

loader = importlib.import_module("load.loaders.spedas")


@pytest.mark.parametrize("mission,product,probe", [("themis", "fgm", "b"), ("mms", "fpi", 2), ("omni", "hro", None)])
@pytest.mark.parametrize("return_type", ["dataframe", "dict"])
def test_shared_loader_clips_and_merges_intervals(monkeypatch, mission, product, probe, return_type):
    raw = {"value": {"x": [-1., 0., 1., 2., 3., 4.], "y": [9., 10., 11., 12., 13., 14.]}}
    calls = []

    def read(*args, **kwargs):
        calls.append(kwargs)
        return raw, _variable_layouts(raw)

    monkeypatch.setattr(loader, "read_spedas", read)
    result = loader.load_spedas(mission, [[0, 1], [2, 3]], product, probe=probe, return_type=return_type)
    assert len(calls) == 2
    if return_type == "dataframe":
        assert result["Timestamp"].tolist() == [0., 1., 2., 3.]
        assert result["value"].tolist() == [10., 11., 12., 13.]
        assert result.attrs["source"] == mission
    else:
        np.testing.assert_array_equal(result["value"]["y"], [10., 11., 12., 13.])


def test_omni_defaults_are_declared():
    assert default_variables("omni", "hro", None) == loader.OMNI_VARIABLES
    assert OMNI_FIELDS["hro"]["kind"] == "scalar"


def test_local_omni_path_uses_shared_assembly(monkeypatch, tmp_path):
    folder = tmp_path / "omni_data" / "hro_1min" / "1970"
    folder.mkdir(parents=True)
    (folder / "omni_hro_1min_19700101_v01.cdf").touch()

    class Field:
        attrs = {"UNITS": "nT", "FILLVAL": -999.}

        def __array__(self, dtype=None, copy=None):
            return np.array([2., -999., 4.], dtype=dtype)

    class CDF(dict):
        def __init__(self, path):
            super().__init__(Epoch=pd.date_range("1970-01-01", periods=3, freq="s"), BZ_GSM=Field())

        def __enter__(self):
            return self

        def __exit__(self, *args):
            return False

    spacepy = ModuleType("spacepy")
    spacepy.pycdf = ModuleType("spacepy.pycdf")
    spacepy.pycdf.CDF = CDF
    monkeypatch.setitem(sys.modules, "spacepy", spacepy)
    result = loader.load_spedas("omni", [0, 1], variables=["BZ_GSM"], prefix=tmp_path, download=False)
    assert result["Timestamp"].tolist() == [0., 1.]
    assert result.loc[0, "BZ_GSM"] == 2.
    assert np.isnan(result.loc[1, "BZ_GSM"])
