"""Range consumers validate endpoints before invoking I/O or plotting."""

import importlib

import pandas as pd
import pytest


@pytest.mark.parametrize("loader, product, backend", [
    ("themis", "fgm", "read_spedas"),
    ("mms", "fpi", "read_spedas"),
    ("omni", "hro", "read_spedas"),
    ("load_indices", "ae", "read_indices"),
])
@pytest.mark.parametrize("trange", ["2026-04-07", [], ["2026-04-07"], [0, 1, 2]])
def test_loaders_reject_wrong_endpoint_count_before_io(monkeypatch, loader, product, backend, trange):
    module = importlib.import_module("load.load")

    def unexpected_io(*args, **kwargs):
        raise AssertionError("invalid time range reached the data backend")

    if backend == "read_spedas":
        monkeypatch.setattr(importlib.import_module("load.loaders.spedas"), backend, unexpected_io)
        call = lambda: module.load_spedas(loader, trange, product)
    else:
        monkeypatch.setattr(importlib.import_module("load.loaders.indices"), backend, unexpected_io)
        call = lambda: getattr(module, loader)(trange, product=product)
    with pytest.raises(ValueError, match="exactly two"):
        call()


@pytest.mark.parametrize("trange", ["2026-04-07", [], [0], [0, 1, 2]])
def test_scene_ranges_require_two_endpoints(trange):
    from plots.v3d import SceneContext
    context = object.__new__(SceneContext)
    context.pose_timestamps = pd.Index([0., 1., 2.]).to_numpy()
    with pytest.raises(ValueError, match="exactly two"):
        context.make_timerange(trange)


@pytest.mark.parametrize("trange", ["2026-04-07", [], [0], [0, 1, 2]])
def test_overview_range_validation_precedes_figure_creation(trange):
    import matplotlib.pyplot as plt
    from plots import Panel, overview
    from plots.basic import tplot

    before = plt.get_fignums()
    with pytest.raises(ValueError, match="exactly two"):
        overview([Panel(tplot)], trange=trange)
    assert plt.get_fignums() == before
