"""Solar-wind loading combines native interval samples without evaluation."""

import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

from common.frame import TimeSeriesDataFrame, read_variable
from gena.frame import ParameterDataFrame
from load import load_sw
from load.api.variable import assemble


def scalar_frame(dates, values, **provenance):
    return assemble({
        name: {"time": dates, "values": np.asarray(samples, dtype=float),
               "time_axis": 0, "kind": "scalar", "columns": [name], "units": "test-unit"}
        for name, samples in values.items()
    }, **provenance)


class SolarWindLoadingTests(unittest.TestCase):
    def setUp(self):
        self.trange = ["2025-03-10T00:00:00", "2025-03-10T01:00:00"]
        self.omni_dates = pd.date_range(self.trange[0], periods=4, freq="20min")
        self.dst_dates = pd.to_datetime(["2025-03-10T00:29:30", "2025-03-10T00:59:30"])
        self.omni = scalar_frame(self.omni_dates, {
            "Pressure": [1, 4, 9, 16], "BY_GSM": [0, 2, -3, 4], "BZ_GSM": [-2, -4, -8, -16],
            "Beta": [1, np.nan, 3, 4], "Mgs_mach_num": [5, 6, 7, 8], "SYM_H": [-10, -12, -13, -15],
        }, source="omni", product="hro")
        self.dst = scalar_frame(self.dst_dates, {"dst": [-25, -40]}, source="indices", product="dst")

    def test_combines_native_samples_using_load_facade(self):
        with patch("load.facade.load_data", side_effect=[self.omni, self.dst]) as loader:
            frame = load_sw(self.trange, prefix="/local/data", download=False)
        self.assertIs(type(frame), TimeSeriesDataFrame)
        self.assertEqual(len(frame), 6)
        np.testing.assert_array_equal(frame["Datetime"], self.omni_dates.union(self.dst_dates))
        for name in ("Pressure", "BZ_GSM", "SYM_H"):
            variable = read_variable(frame, name)
            np.testing.assert_array_equal(variable["time"], self.omni_dates)
            np.testing.assert_array_equal(variable["values"], self.omni[name])
        np.testing.assert_array_equal(read_variable(frame, "dst")["time"], self.dst_dates)
        np.testing.assert_array_equal(read_variable(frame, "dst")["values"], [-25, -40])
        self.assertTrue(frame.loc[frame["Datetime"].isin(self.dst_dates), "Pressure"].isna().all())
        self.assertTrue(frame.loc[frame["Datetime"].isin(self.omni_dates), "dst"].isna().all())
        self.assertTrue(frame.loc[frame["Datetime"] == self.omni_dates[1], "Beta"].isna().all())
        self.assertEqual(frame.attrs["sources"]["indices"]["product"], "dst")
        self.assertEqual(frame.attrs["variables"]["Pressure"]["units"], "test-unit")
        self.assertEqual([call.args[0] for call in loader.call_args_list], ["omni", "indices"])
        for call in loader.call_args_list:
            self.assertEqual(call.args[1], self.trange)
            self.assertEqual(call.kwargs["prefix"], "/local/data")
            self.assertIs(call.kwargs["download"], False)

    def test_missing_fields_stay_nan_without_fallback_constants(self):
        empty = assemble({})
        with patch("load.facade.load_data", side_effect=[self.omni, empty]):
            frame = load_sw(self.trange, download=False)
        self.assertEqual(len(frame), 4)
        self.assertTrue(frame["dst"].isna().all())
        with patch("load.facade.load_data", side_effect=[empty, self.dst]):
            frame = load_sw(self.trange, download=False)
        self.assertTrue(frame["Pressure"].isna().all())
        np.testing.assert_array_equal(frame["dst"], [-25, -40])

    def test_empty_sources_return_empty_timeseries(self):
        with patch("load.facade.load_data", side_effect=[assemble({}), assemble({})]):
            frame = load_sw(self.trange, download=False)
        self.assertIs(type(frame), TimeSeriesDataFrame)
        self.assertTrue(frame.empty)
        self.assertIn("Pressure", frame)
        self.assertIn("dst", frame)

    def test_invalid_ranges_are_rejected_before_loading(self):
        for interval in (self.trange[0], [], self.trange[:1], self.trange * 2, self.trange[::-1]):
            with self.subTest(interval=interval), patch("load.facade.load_data") as loader:
                with self.assertRaisesRegex(ValueError, "trange"):
                    load_sw(interval)
                loader.assert_not_called()

    def test_loading_errors_propagate(self):
        for results in ([OSError("missing OMNI")], [self.omni, OSError("missing Dst")]):
            with self.subTest(results=results), patch("load.facade.load_data", side_effect=results):
                with self.assertRaisesRegex(OSError, "missing"):
                    load_sw(self.trange, download=False)

    def test_model_averages_native_physical_columns(self):
        from libs.space import get_obj_ave

        with patch("load.facade.load_data", side_effect=[self.omni, self.dst]):
            weather = load_sw(self.trange, download=False)
        with patch("libs.space.load_sw", return_value=weather) as loader, patch(
            "libs.space.Bowshock",
        ) as bowshock, patch("libs.space.Magnetopause"), patch("libs.space.T96_fieldlines") as model:
            get_obj_ave(self.trange)
        loader.assert_called_once_with(self.trange)
        bowshock.assert_called_once_with(-7.5, 8 / 3, 6.5, 7.5)
        midpoint = pd.Timestamp("2025-03-10T00:30:00").timestamp()
        model.assert_called_once_with(midpoint, 7.5, -12.5, 0, -7.5)

    def test_scene_explicitly_resamples_weather_to_cache_times(self):
        from plots.v3d import SceneContext

        with patch("load.facade.load_data", side_effect=[self.omni, self.dst]):
            weather = load_sw(self.trange, download=False)
        pose = ParameterDataFrame({
            "Datetime": self.omni_dates, "Timestamp": self.omni_dates.asi8 / 1e9,
            "MJ2000_Star_X": 4e8, "MJ2000_Star_Y": 1e6, "MJ2000_Star_Z": 1e6,
            "MJ2000_q1": 0.1, "MJ2000_q2": 0.1, "MJ2000_q3": 0.1,
        })
        with patch("plots.v3d.load_data", return_value=pose), patch(
            "plots.v3d.load_sw", return_value=weather,
        ) as loader:
            scene = SceneContext.load(
                self.trange, satellites=(), cache_cadence="20min",
                include_field_lines=False, include_neutral_sheet=False,
            )
        self.assertEqual(len(loader.call_args.args[0]), 2)
        np.testing.assert_array_equal(scene.space_weather["Pressure"], [1, 4, 9, 16])
        self.assertEqual(len(scene.space_weather["dst"]), 4)
        self.assertTrue(np.isfinite(scene.space_weather["dst"][2]))


if __name__ == "__main__":
    unittest.main()
