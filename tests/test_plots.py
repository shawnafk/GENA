"""Shared plotting primitives and offline cross-mission composition."""

import importlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from plots import vis
from plots.basic import tplot
from common.frame import TimeSeriesDataFrame
from gena.frame import ParticleDataFrame
from plots import Panel, env, overview
from plots.basic import gena, themis, mms, omni


def test_region_strip_uses_shared_classification_codes():
    from plots import plot_region_code
    result = TimeSeriesDataFrame({"Datetime": pd.date_range("2025-01-01", periods=3), "region_code": [4, 5, 6]})
    figure, axes = overview([Panel(plot_region_code, kwargs={"result": result, "legend": True})])
    try:
        assert axes[0].collections[0].get_array().size == 3
        figure.canvas.draw()
    finally:
        plt.close(figure)


@pytest.mark.parametrize("mission", ["gena", "themis", "mms", "omni"])
def test_basic_modules_do_not_reexport_tplot(mission):
    module = importlib.import_module(f"plots.basic.{mission}")
    assert not hasattr(module, "tplot")


def test_joint_figure_uses_four_missions_and_preserves_native_samples():
    t = pd.date_range("2026-04-07", periods=8, freq="min")
    events = ParticleDataFrame({"Datetime": t, "Timestamp": t.as_unit("ns").asi8 / 1e9})
    field = {"fgs_btotal": {"x": t[:3], "y": [1., 2., 3.]},
             "fgs_gsm": {"x": t[:3], "y": np.ones((3, 3))}}
    plasma = {
        name: {"x": t[::2], "y": np.ones((4, 3)) if name == "bulkv_gse_fast" else np.arange(1., 5.)}
        for name in ("numberdensity_fast", "tempperp_fast", "temppara_fast", "bulkv_gse_fast")
    }
    solar_wind = {"BZ_GSM": {"x": t[::3], "y": [-2., 0., 2.]}}
    figure, axes = overview([
        Panel(gena.counts, ([events],), {"window": "2min"}),
        Panel(themis.fgm_bfield, (field,)),  # This function has no f parameter.
        Panel(mms.fpi, (plasma,), nrows=3),
        Panel(omni.plot_panel, (solar_wind, "BZ_GSM")),
    ], trange=[t[0], t[-1]])
    try:
        assert len(axes) == 6
        assert len(figure.axes) == 6
        assert all(axis.get_shared_x_axes().joined(axes[0], axis) for axis in axes[1:])
        assert [len(axes[i].lines[0].get_ydata()) for i in (1, 2, 5)] == [3, 4, 3]
        np.testing.assert_array_equal(axes[2].lines[0].get_xdata(), t[::2].to_numpy())
        np.testing.assert_allclose(axes[-1].get_xlim(), mdates.date2num([t[0], t[-1]]))
        assert axes[-1].get_xlabel() == "UTC"
        figure.canvas.draw()
    finally:
        plt.close(figure)


def test_overview_reuses_existing_axes_and_rejects_mixed_figures():
    first, axes = plt.subplots(2, 1)
    second, other = plt.subplots()
    data = {"x": pd.date_range("2026-04-07", periods=2), "y": [1., 2.]}
    panels = [Panel(tplot, (data,)), Panel(tplot, (data,))]
    try:
        result, actual = overview(panels, ax=axes)
        assert result is first
        assert actual == list(axes)
        with pytest.raises(ValueError, match="same|belong"):
            overview(panels, ax=[axes[0], other])
        with pytest.raises(ValueError, match="2 axes"):
            overview(panels, ax=axes[0])
    finally:
        plt.close(first)
        plt.close(second)


def test_overview_rejects_invalid_calls_before_creating_figure():
    before = plt.get_fignums()
    with pytest.raises(ValueError, match="nonempty"):
        overview([])
    with pytest.raises(ValueError, match="start before end"):
        overview([Panel(tplot)], trange=["2026-04-08", "2026-04-07"])
    with pytest.raises(TypeError):
        overview([Panel(tplot)])
    assert plt.get_fignums() == before
