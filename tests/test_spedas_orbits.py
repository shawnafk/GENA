"""Registered orbit products use the shared SPEDAS source and frame pipeline."""

import importlib
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd
from spacepy import pycdf

from common.frame import TimeSeriesDataFrame, read_variable
from load import load_data, load_pos
from load.loaders.spedas import load_spedas


class RegisteredOrbitTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.dates = pd.date_range("2025-03-10", periods=7, freq="10s")
        self.trange = [self.dates[1].isoformat(), self.dates[5].isoformat()]
        self.positions = np.arange(21, dtype=float).reshape(7, 3) ** 2
        for module in ("pyspedas.projects.mms.mms_config", "pyspedas.projects.goes.config"):
            restore = patch.dict(importlib.import_module(module).CONFIG)
            restore.start()
            self.addCleanup(restore.stop)
        restore = patch.dict(os.environ)
        restore.start()
        self.addCleanup(restore.stop)
        offline = patch("requests.sessions.Session.request", side_effect=AssertionError("unexpected network I/O"))
        offline.start()
        self.addCleanup(offline.stop)

    def record(self):
        return {"x": self.dates.as_unit("ns").asi8 / 1e9, "y": self.positions, "units": "km"}

    def test_mms_mec_registered_defaults_and_position_entry_are_identical(self):
        variable = "mms2_mec_r_gsm"
        with patch("pyspedas.projects.mms.mms_load_mec", return_value={
            variable: self.record(), "unused": self.record(),
        }) as backend:
            direct = load_data("mms", self.trange, "mec", probe=2, prefix=self.root, download=False)
            positions = load_pos("mms", self.trange, "2", prefix=self.root, download=False)
        pd.testing.assert_frame_equal(direct, positions)
        self.assertIs(type(direct), TimeSeriesDataFrame)
        np.testing.assert_array_equal(read_variable(direct, variable)["values"], self.positions[1:6])
        self.assertNotIn("unused", direct)
        self.assertEqual(direct.attrs["product"], "mec")
        self.assertEqual(direct.attrs["probe"], "2")
        self.assertEqual(direct.attrs["level"], "l2")
        call = backend.call_args.kwargs
        self.assertEqual(call["varnames"], [variable])
        self.assertEqual(call["data_rate"], "srvy")
        self.assertEqual(call["datatype"], "epht89q")
        self.assertTrue(call["no_update"])
        self.assertFalse(call["spdf"])
        config = importlib.import_module("pyspedas.projects.mms.mms_config").CONFIG
        self.assertEqual(config["local_data_dir"], str(self.root / "mms_data"))
        self.assertEqual(os.environ["MMS_DATA_DIR"], str(self.root / "mms_data"))

    def test_mec_selection_options_and_dictionary_output(self):
        variable = "mms4_mec_r_gse"
        with patch("pyspedas.projects.mms.mms_load_mec", return_value={variable: self.record()}) as backend:
            for selection in ("r_gse", variable):
                raw = load_data(
                    "mms", self.trange, "mec", probe="4", key=selection,
                    data_rate="srvy", datatype="custom", spdf=True,
                    download="local", prefix=self.root, return_type="dict",
                )
                np.testing.assert_array_equal(raw[variable]["y"], self.positions[1:6])
                self.assertEqual(backend.call_args.kwargs["varnames"], [variable])
                self.assertEqual(backend.call_args.kwargs["datatype"], "custom")
                self.assertTrue(backend.call_args.kwargs["spdf"])

    def test_goes_reads_local_cdf_with_unprefixed_selection(self):
        path = self.root / "goes_data/goes16/orbit/2025/goes16_ephemeris_ssc_20250101_v01.cdf"
        path.parent.mkdir(parents=True)
        with pycdf.CDF(str(path), "") as cdf:
            cdf["Epoch"] = self.dates.to_pydatetime()
            cdf["Epoch"].attrs["VAR_TYPE"] = "support_data"
            for variable in ("XYZ_GSM", "XYZ_GSE"):
                cdf[variable] = self.positions
                cdf[variable].attrs.update(DEPEND_0="Epoch", VAR_TYPE="data", UNITS="km")
        direct = load_data("goes", self.trange, probe=16, prefix=self.root, download=False)
        positions = load_pos("goes", self.trange, 16, prefix=self.root, download=False)
        pd.testing.assert_frame_equal(direct, positions)
        np.testing.assert_array_equal(read_variable(direct, "g16_orbit_XYZ_GSM")["values"], self.positions[1:6])
        self.assertEqual(direct.attrs["product"], "orbit")
        self.assertEqual(direct.attrs["probe"], "16")
        self.assertIsNone(direct.attrs["level"])
        for selection in ("XYZ_GSE", "g16_orbit_XYZ_GSE"):
            frame = load_data("goes-r", self.trange, "orbit", probe="016", key=selection,
                              prefix=self.root, download="offline")
            self.assertEqual(set(frame.attrs["variables"]), {"g16_orbit_XYZ_GSE"})
            np.testing.assert_array_equal(frame["Datetime"], self.dates[1:6])

    def test_goes_download_modes_and_configured_cache(self):
        for mode, offline, refresh in ((False, True, False), ("local", True, False),
                                       ("local_first", False, False), ("refresh", False, True)):
            with self.subTest(mode=mode), patch(
                "pyspedas.projects.goes.load_orbit", return_value={"g15_orbit_XYZ_GSM": self.record()},
            ) as backend:
                frame = load_spedas("goes", self.trange, "orbit", probe="15", download=mode, prefix=self.root)
                call = backend.call_args.kwargs
                self.assertEqual(call["varnames"], ["XYZ_GSM"])
                self.assertIs(call["no_update"], offline)
                self.assertIs(call["force_download"], refresh)
                self.assertTrue(call["notplot"])
                self.assertEqual(len(frame), 5)
                config = importlib.import_module("pyspedas.projects.goes.config").CONFIG
                self.assertEqual(config["local_data_dir"], str(self.root / "goes_data"))

    def test_existing_mms_products_keep_their_backend_dispatch(self):
        cases = [("fpi", "numberdensity_fast", "mms1_dis_numberdensity_fast", "mms_load_fpi"),
                 ("hpca", "hplus_number_density", "mms1_hpca_hplus_number_density", "mms_load_hpca")]
        for product, selection, variable, target in cases:
            with self.subTest(product=product), patch(
                f"pyspedas.projects.mms.{target}",
                return_value={variable: {"x": self.dates.as_unit("ns").asi8 / 1e9, "y": np.arange(7, dtype=float)}},
            ) as backend:
                frame = load_data("mms", self.trange, product, key=selection, download=False, prefix=self.root)
                self.assertEqual(backend.call_args.kwargs["varnames"], [variable])
                self.assertTrue(backend.call_args.kwargs["no_update"])
                np.testing.assert_array_equal(frame[variable], np.arange(1, 6))
                if product == "fpi":
                    self.assertEqual(backend.call_args.kwargs["datatype"], ["dis-moms"])

    def test_empty_or_invalid_orbit_requests(self):
        for mission, product, probe, target in (
            ("mms", "mec", "1", "mms.mms_load_mec"), ("goes", "orbit", "16", "goes.load_orbit"),
        ):
            with self.subTest(mission=mission), patch(f"pyspedas.projects.{target}", return_value={}):
                frame = load_data(mission, self.trange, product, probe=probe, prefix=self.root, download=False)
                self.assertIs(type(frame), TimeSeriesDataFrame)
                self.assertTrue(frame.empty)
                self.assertIn("empty_reason", frame.attrs)
        for probe in (None, 0, "x", -1, 16.5):
            with self.subTest(probe=probe), self.assertRaisesRegex(ValueError, "GOES"):
                load_data("goes", self.trange, "orbit", probe=probe, prefix=self.root, download=False)
        with self.assertRaisesRegex(ValueError, "no level"):
            load_data("goes", self.trange, "orbit", probe="16", level="l2", prefix=self.root, download=False)


if __name__ == "__main__":
    unittest.main()
