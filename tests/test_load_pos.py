"""Position loading retains native interval samples; interpolation is explicit."""

import importlib
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd
from spacepy import pycdf

from common.frame import TimeSeriesDataFrame, read_variable
from gena.frame import ParameterDataFrame, ParticleDataFrame
from load import load_data, load_pos


class PositionLoadingTests(unittest.TestCase):
    def setUp(self):
        self.dates = pd.date_range("2025-03-10T00:00:00", periods=7, freq="10s")
        # Deliberately nonlinear: resampling the endpoints cannot reproduce this.
        self.positions = np.arange(21, dtype=float).reshape(7, 3) ** 2 + 1
        self.trange = [self.dates[1], self.dates[5]]

    def test_themis_reads_all_native_cdf_samples(self):
        config = importlib.import_module("pyspedas.projects.themis.config").CONFIG
        with tempfile.TemporaryDirectory() as root:
            path = Path(root) / "themis_data/thb/l1/state/2025/thb_l1_state_20250310.cdf"
            path.parent.mkdir(parents=True)
            with pycdf.CDF(str(path), "") as cdf:
                cdf["Epoch"] = self.dates.to_pydatetime()
                cdf["Epoch"].attrs["VAR_TYPE"] = "support_data"
                cdf["thb_pos_gsm"] = self.positions
                cdf["thb_pos_gsm"].attrs.update(
                    DEPEND_0="Epoch", VAR_TYPE="data", UNITS="km",
                )
            with patch.dict(config), patch(
                "requests.sessions.Session.request",
                side_effect=AssertionError("offline loading must not use the network"),
            ):
                frame = load_pos("themis", self.trange, "b", prefix=root, download=False)
                direct = load_data(
                    "themis", self.trange, "state", probe="b",
                    variables=["pos_gsm"], prefix=root, download=False,
                )
            self.assertIsInstance(frame, TimeSeriesDataFrame)
            self.assertEqual(len(frame), 5)
            np.testing.assert_array_equal(frame["Datetime"], self.dates[1:6])
            np.testing.assert_array_equal(
                read_variable(frame, "thb_pos_gsm")["values"], self.positions[1:6],
            )
            pd.testing.assert_frame_equal(frame, direct)
            self.assertEqual(frame.attrs["level"], "l1")
            self.assertEqual(frame.attrs["probe"], "b")
            sampled = frame.resample_at([self.dates[2] + pd.Timedelta("5s")])
            np.testing.assert_allclose(
                read_variable(sampled, "thb_pos_gsm")["values"],
                [(self.positions[2] + self.positions[3]) / 2],
            )

    def test_gena_reads_position_and_particle_rows_without_interpolation(self):
        columns = ["MJ2000_Star_X", "MJ2000_Star_Y", "MJ2000_Star_Z"]
        source = pd.DataFrame(self.positions, columns=columns)
        source["Timestamp"] = self.dates.as_unit("ns").asi8 / 1e9
        with tempfile.TemporaryDirectory() as root:
            source.to_parquet(Path(root) / "20250310_pdf.pqt")
            particles = source.assign(
                Energy=np.arange(7, dtype=float) ** 2, PID=np.arange(7), X1=1.0,
            )
            particles.to_parquet(Path(root) / "20250310_l2_a.pqt")
            with patch("load.loaders.gena.get_configured_paths", return_value={"pqt_base": root}):
                for mission in ("gena", "queqiao", "qq"):
                    frame = load_pos(mission, self.trange)
                    self.assertIs(type(frame), TimeSeriesDataFrame)
                    self.assertEqual(len(frame), 5)
                    np.testing.assert_array_equal(frame["Datetime"], self.dates[1:6])
                    np.testing.assert_array_equal(frame[columns], self.positions[1:6])
                    self.assertEqual(frame.attrs["units"], "m")
                    self.assertEqual(frame.attrs["source"], "gena")
                    sampled = frame.resample_at([self.dates[2] + pd.Timedelta("5s")])
                    self.assertIs(type(sampled), TimeSeriesDataFrame)
                    np.testing.assert_allclose(
                        sampled[columns], [(self.positions[2] + self.positions[3]) / 2],
                    )
                self.assertIsInstance(load_data("gena", self.trange, "pdf"), ParameterDataFrame)
                events = load_data("gena", self.trange, "l2_a")
            self.assertIsInstance(events, ParticleDataFrame)
            np.testing.assert_array_equal(events["Energy"], np.arange(1, 6) ** 2)
            np.testing.assert_array_equal(events[columns], self.positions[1:6])

    def test_mms_and_goes_clip_native_vectors(self):
        cases = [
            ("mms", "2", "mms2_mec_r_gsm", "pyspedas.projects.mms.mms_load_mec"),
            ("goes", "16", "g16_orbit_XYZ_GSM", "pyspedas.projects.goes.load_orbit"),
        ]
        for mission, probe, variable, target in cases:
            with self.subTest(mission=mission), patch(target, return_value={
                variable: {"x": self.dates.as_unit("ns").asi8 / 1e9, "y": self.positions},
            }) as reader:
                frame = load_pos(mission, self.trange, probe, download=False, prefix="/unused")
                self.assertIsInstance(frame, TimeSeriesDataFrame)
                self.assertEqual(len(frame), 5)
                np.testing.assert_array_equal(read_variable(frame, variable)["values"], self.positions[1:6])
                self.assertTrue(reader.call_args.kwargs["no_update"])
                np.testing.assert_array_equal(
                    pd.to_datetime(reader.call_args.kwargs["trange"]), pd.DatetimeIndex(self.trange),
                )

    def test_invalid_ranges_fail_before_io(self):
        for trange in (self.dates[0], [], self.dates[:1], self.dates[:3], self.trange[::-1]):
            with self.subTest(trange=trange), self.assertRaisesRegex(ValueError, "trange"):
                load_pos("themis", trange)

    def test_empty_position_products_still_return_frames(self):
        with patch("load.api.spedas._read_themis", return_value=({}, {})):
            frame = load_pos("themis", self.trange, "b", download=False, prefix="/unused")
        self.assertIsInstance(frame, TimeSeriesDataFrame)
        self.assertTrue(frame.empty)
        with tempfile.TemporaryDirectory() as root, patch(
            "load.loaders.gena.get_configured_paths", return_value={"pqt_base": root},
        ):
            frame = load_pos("gena", self.trange)
        self.assertIs(type(frame), TimeSeriesDataFrame)
        self.assertTrue(frame.empty)


if __name__ == "__main__":
    unittest.main()
