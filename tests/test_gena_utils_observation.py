import matplotlib
import numpy as np
import pandas as pd
import pytest

matplotlib.use("Agg")

from common.constants import RE
from common.frame import TimeSeriesDataFrame
from gena import utils
from plots.advance import env
from gena.frame import ParameterDataFrame
from libs import ephemeris


def test_env_uses_cols_to_create_and_group_dynamic_axes():
    times = pd.date_range("2026-04-07", periods=4, freq="min")
    omni = TimeSeriesDataFrame(
        {
            "Datetime": times,
            "flow_speed": [400.0, 410.0, 420.0, 430.0],
            "BZ_GSM": [-2.0, -1.0, 1.0, 2.0],
        }
    )
    indices = TimeSeriesDataFrame(
        {
            "Datetime": times,
            "dst": [-10.0, -11.0, -12.0, -13.0],
            "ae": [100.0, 110.0, 120.0, 130.0],
            "al": [-80.0, -90.0, -100.0, -110.0],
            "au": [20.0, 20.0, 20.0, 20.0],
        }
    )

    figure, axes = env(
        omni,
        indices,
        cols=["dst", ("ae", "al", "au"), "flow_speed", "BZ_GSM"],
    )

    assert len(axes) == 4
    assert [len(axis.lines) for axis in axes] == [2, 4, 1, 2]
    assert axes[-1].get_xlabel() == "UTC"
    figure.canvas.draw()


def test_env_rejects_wrong_number_of_supplied_axes():
    import matplotlib.pyplot as plt

    figure, axes = plt.subplots(2, 1)
    with pytest.raises(ValueError, match="3 axes"):
        env(cols=["dst", "al", "ae"], f=figure, ax=axes)


def test_env_accepts_one_column_as_a_string():
    times = pd.date_range("2026-04-07", periods=2, freq="min")
    indices = TimeSeriesDataFrame({"Datetime": times, "dst": [-10.0, -11.0]})

    figure, axis = env(indices, cols="dst")

    assert not isinstance(axis, list)
    assert len(axis.lines) == 2
    figure.canvas.draw()


def test_pointing_geometry_returns_forward_plane_intersections(monkeypatch):
    dates = pd.date_range("2026-04-07", periods=2, freq="min")
    pdf = ParameterDataFrame(
        {"Datetime": pd.to_datetime([value.timestamp() for value in dates], unit="s"), 
            "Timestamp": [value.timestamp() for value in dates],
            "MJ2000_Star_X": np.full(2, 10.0 * RE),
            "MJ2000_Star_Y": np.full(2, 20.0 * RE),
            "MJ2000_Star_Z": np.full(2, 30.0 * RE),
            "MJ2000_q1": np.full(2, 0.1),
            "MJ2000_q2": np.full(2, 0.1),
            "MJ2000_q3": np.full(2, 0.1),
        }
    )

    class FakeCoord:
        def __init__(self, dates, values, coord_type, rotation):
            self.values = np.asarray(values, dtype=float)
            self.coord_type = coord_type
            self.result = None

        def convert(self, coord):
            if self.coord_type == "ECI2000":
                self.result = self.values
            else:
                self.result = np.tile([-1.0, -1.0, -1.0], (len(self.values), 1))

        def __getitem__(self, coord):
            return self.result

    monkeypatch.setattr(utils, "Coord", FakeCoord)
    monkeypatch.setattr(ephemeris, "moon", lambda dates: np.zeros((len(dates), 3)))

    geometry = utils.pointing_geometry(pdf, probe="A", coord="GSM")

    np.testing.assert_allclose(geometry[["xy_x", "xy_y", "xy_z"]], [[-20, -10, 0]] * 2)
    np.testing.assert_allclose(geometry[["xz_x", "xz_y", "xz_z"]], [[-10, 0, 10]] * 2)
    np.testing.assert_allclose(geometry[["yz_x", "yz_y", "yz_z"]], [[0, 10, 20]] * 2)


def test_orbit_plots_precomputed_geometry_without_recomputing_pose():
    times = pd.date_range("2026-04-07", periods=3, freq="min")
    geometry = pd.DataFrame(
        {
            "Datetime": times,
            "spacecraft_x": [-60.0, -59.9, -59.8],
            "spacecraft_y": [0.0, 0.1, 0.2],
            "spacecraft_z": [1.0, 1.1, 1.2],
            "xy_x": [-40.0, -40.1, -40.2],
            "xy_y": [2.0, 2.1, 2.2],
            "xz_x": [-60.0, -60.0, -60.0],
            "xz_z": [0.0, 0.1, 0.2],
            "yz_y": [1.0, 1.1, 1.2],
            "yz_z": [0.0, 0.1, 0.2],
        }
    )

    figure, axes = utils.orbit(geometry=geometry)

    assert len(axes) == 4
    assert [len(axis.lines) for axis in axes] == [3, 2, 2, 2]
    figure.canvas.draw()
