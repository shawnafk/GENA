"""Event geometry preserves invalid rows when one attitude array is reused."""

import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

from common.coord import transform_vectors, fov_to_vector
from common.constants import RE
from gena.frame import ParticleDataFrame, ParameterDataFrame
from gena.utils import velocity
from gena.post.project import (
    event_axes, events_rays, position_as_coord, projected_fov_axis_angle,
    plane_intersections, ray_sphere_intersections, obj_in_view,
    project_events, moon_blocking_mask, moon_blocking_intervals,
)


class EventGeometryTests(unittest.TestCase):
    def setUp(self):
        dates = pd.date_range("2025-03-10", periods=4, freq="1s")
        self.events = ParticleDataFrame({
            "Datetime": dates, "Timestamp": dates.as_unit("ns").asi8 / 1e9,
            "MJ2000_q1": [0., np.nan, 2., 0.],
            "MJ2000_q2": [0., 0., 0., 0.],
            "MJ2000_q3": [0., 0., 0., np.sin(np.pi / 4)],
            "MJ2000_Star_X": 0., "MJ2000_Star_Y": 0., "MJ2000_Star_Z": 0.,
            "X1": [0., 1., 2., 3.], "Y1": [1., 2., 3., 4.],
            "X2": [10., 11., 12., 13.], "Y2": [-4., -3., -2., -1.],
            "Tof": [100., 100., 100., 100.],
        }, index=[7, 3, 5, 1])
        self.events.attrs = {"product": "l2_a"}
        self.valid = [0, 3]
        self.invalid = [1, 2]

    def test_pose_accessors_preserve_rows_units_and_invalid_attitude(self):
        positions = np.array([[RE, 0., 0.], [0., RE, np.nan], [0., 0., 0.], [2 * RE, RE, 0.]])
        for frame_type in (ParticleDataFrame, ParameterDataFrame):
            with self.subTest(frame_type=frame_type.__name__):
                frame = frame_type(self.events)
                frame["Datetime"] = frame["Datetime"].iloc[[3, 0, 0, 1]].to_numpy()
                frame[["MJ2000_Star_X", "MJ2000_Star_Y", "MJ2000_Star_Z"]] = positions
                original = frame.copy(deep=True)
                pos, rot = frame.pos(), frame.rot()
                np.testing.assert_allclose(pos, positions)
                np.testing.assert_allclose(rot[self.valid], [
                    [0., 0., 0., 1.], [0., 0., np.sin(np.pi / 4), np.cos(np.pi / 4)],
                ])
                self.assertTrue(np.isnan(rot[self.invalid]).all())
                pos[:] = 0.
                rot[:] = 0.
                pd.testing.assert_frame_equal(frame, original)
                self.assertEqual(frame.iloc[:0].pos().shape, (0, 3))
                self.assertEqual(frame.iloc[:0].rot().shape, (0, 4))

    def test_ray_direction_and_velocity_keep_invalid_rows_and_input(self):
        original = self.events.copy(deep=True)
        for probe in ("A", "B"):
            with self.subTest(probe=probe):
                ray = events_rays(self.events, probe=probe, coord="ECI2000")["ECI2000"]
                vectors = velocity(self.events, probe=probe, coord="ECI2000")
                direction = vectors[["DX", "DY", "DZ"]].to_numpy()
                speeds = vectors[["VX", "VY", "VZ"]].to_numpy()
                self.assertTrue(np.isnan(ray[self.invalid]).all())
                self.assertTrue(np.isnan(direction[self.invalid]).all())
                self.assertTrue(np.isnan(speeds[self.invalid]).all())
                expected = transform_vectors(
                    np.tile([10., -5., 50.], (2, 1)), source=probe, target="ECI2000",
                    attitude=[[0., 0., 0., 1.], [0., 0., np.sin(np.pi / 4), np.cos(np.pi / 4)]],
                )
                np.testing.assert_allclose(ray[self.valid], expected)
                np.testing.assert_allclose(direction[self.valid], expected / np.linalg.norm(expected, axis=1)[:, None])
                np.testing.assert_allclose(speeds[self.valid], expected * 10.)
                self.assertEqual(list(vectors.index), list(self.events.index))
                self.assertIsInstance(vectors, ParticleDataFrame)
                self.assertEqual(vectors.attrs["product"], "l2_a")
                self.assertEqual(vectors.attrs["vector_coordinate"], "ECI2000")
        pd.testing.assert_frame_equal(self.events, original)
        self.assertEqual(self.events.attrs, original.attrs)

    def test_velocity_keeps_directions_when_tof_is_invalid(self):
        events = self.events.copy()
        events[["MJ2000_q1", "MJ2000_q2", "MJ2000_q3"]] = 0.
        events["Tof"] = [100., 0., -1., np.nan]
        result = velocity(events, probe="QQ", coord="ECI2000", vscale=2.)
        expected_direction = np.array([10., -5., 50.]) / np.linalg.norm([10., -5., 50.])
        np.testing.assert_allclose(result[["DX", "DY", "DZ"]], np.tile(expected_direction, (4, 1)))
        np.testing.assert_allclose(result[["VX", "VY", "VZ"]].iloc[0], [200., -100., 1000.])
        self.assertTrue(result[["VX", "VY", "VZ"]].iloc[1:].isna().all().all())
        empty = velocity(events.iloc[:0], coord="ECI2000")
        self.assertTrue(empty.empty)
        self.assertTrue(set(("DX", "DY", "DZ", "VX", "VY", "VZ")).issubset(empty.columns))

    def test_velocity_inplace_supports_custom_columns_and_per_row_scaling(self):
        events = self.events.iloc[self.valid].copy()
        events["Tof"] = [100., 200.]
        expected_ray = transform_vectors(
            np.tile([10., -5., 50.], (2, 1)), source="A", target="ECI2000",
            attitude=[[0., 0., 0., 1.], [0., 0., np.sin(np.pi / 4), np.cos(np.pi / 4)]],
        )
        result = velocity(events, probe="A", coord="ECI2000", inplace=True,
                          vscale=lambda frame: np.array([2., 3.]),
                          direction_cols=("d1", "d2", "d3"), velocity_cols=("v1", "v2", "v3"))
        self.assertIs(result, events)
        np.testing.assert_allclose(result[["d1", "d2", "d3"]],
                                   expected_ray / np.linalg.norm(expected_ray, axis=1)[:, None])
        np.testing.assert_allclose(result[["v1", "v2", "v3"]], expected_ray * [[20.], [15.]])
        self.assertNotIn("DX", result)
        self.assertNotIn("VX", result)

    def test_all_three_axes_retain_invalid_attitude_masks(self):
        for probe in ("A", "B"):
            with self.subTest(probe=probe):
                axes = event_axes(self.events, probe=probe, coord="ECI2000")
                components = np.stack([axis["ECI2000"] for axis in axes], axis=1)
                self.assertTrue(np.isnan(components[self.invalid]).all())
                products = components[self.valid] @ components[self.valid].transpose(0, 2, 1)
                np.testing.assert_allclose(products, np.broadcast_to(np.eye(3), (2, 3, 3)), atol=1e-14)

    def test_fov_projection_preserves_nan_rows_on_both_rays(self):
        with patch("gena.post.project.position_as_coord", return_value={"A": np.zeros((4, 3))}):
            angles, average = projected_fov_axis_angle(
                self.events, pos_offset=np.zeros((4, 3)), probe="A", coord="A", plane_normal=(0., 0., 1.),
                plane_point=(0., 0., 100.), plot_axes=(0, 1),
            )
        self.assertTrue(np.isnan(angles[self.invalid]).all())
        np.testing.assert_allclose(angles[self.valid], 0., atol=1e-14)
        self.assertAlmostEqual(average, 0.)

    def test_positions_use_explicit_moon_coordinates_without_loading(self):
        events = self.events.copy()
        columns = ["MJ2000_Star_X", "MJ2000_Star_Y", "MJ2000_Star_Z"]
        events[columns] = np.tile([RE, 0., 0.], (4, 1))
        moon = np.tile([60., 0., 0.], (4, 1))
        with patch("libs.ephemeris.moon", side_effect=AssertionError("must use prepared Moon positions")):
            position = position_as_coord(events["Datetime"], events[columns].to_numpy() / RE,
                                         pos_offset=moon, res_coord="ECI2000")
            projected = project_events(events, pos_offset=moon, probe="A", coord="ECI2000",
                                       plane_normal=(0, 0, 1), plane_point=(0, 0, 10))
        np.testing.assert_allclose(position["ECI2000"], np.tile([61., 0., 0.], (4, 1)))
        self.assertEqual(position.unit, "Re")
        np.testing.assert_allclose(projected[self.valid, 2], 10.)
        self.assertTrue(np.isnan(projected[self.invalid]).all())

    def test_positions_apply_optional_offset_in_source_frame_without_loading(self):
        locations = np.tile([1., 2., 3.], (4, 1))
        original = locations.copy()
        with patch("libs.ephemeris.moon", side_effect=AssertionError("unexpected ephemeris loading")):
            for offset in (None, np.array([10., 20., 30.])):
                with self.subTest(offset=offset):
                    position = position_as_coord(self.events["Datetime"], locations,
                                                 pos_offset=offset, ref_coord="GSM", res_coord="AGSM")
                    expected = transform_vectors(
                        locations if offset is None else locations + offset,
                        source="GSM", target="AGSM", times=self.events["Datetime"],
                    )
                    np.testing.assert_allclose(position["AGSM"], expected)
                    np.testing.assert_allclose(position["GSM"], locations if offset is None else locations + offset)
                    self.assertEqual(position.unit, "Re")
        np.testing.assert_array_equal(locations, original)

    def test_plane_projection_has_explicit_forward_semantics(self):
        origins = [[0., 0., 1.], [0., 0., 1.], [0., 0., 1.]]
        directions = [[0., 0., 1.], [0., 0., -1.], [1., 0., 0.]]
        lines = plane_intersections(origins, directions, plane_normal=(0, 0, 1))
        rays = plane_intersections(origins, directions, plane_normal=(0, 0, 1), forward=True)
        np.testing.assert_allclose(lines[:2], 0.)
        self.assertTrue(np.isnan(lines[2]).all())
        self.assertTrue(np.isnan(rays[[0, 2]]).all())
        np.testing.assert_allclose(rays[1], 0.)

    def test_sphere_intersections_choose_the_nearest_forward_hit(self):
        hits = ray_sphere_intersections(
            [[-3., 0., 0.], [0., 0., 0.], [3., 0., 0.], [-3., 0., 0.]],
            [[1., 0., 0.], [1., 0., 0.], [1., 0., 0.], [0., 0., 0.]], radius=1.,
        )
        np.testing.assert_allclose(hits[:2], [[-1., 0., 0.], [1., 0., 0.]])
        self.assertTrue(np.isnan(hits[2:]).all())

    def test_instrument_snapshot_uses_prepared_pose_and_preserves_angular_geometry(self):
        with patch("load.load_pos", side_effect=AssertionError("unexpected data loading")), patch(
            "load.load_rot", side_effect=AssertionError("unexpected data loading"),
        ), patch("load.load_sw", side_effect=AssertionError("unexpected data loading")):
            image = obj_in_view("2025-03-10", [0., 0., 10.], [0., 0., 0., 1.],
                                {"earth": {"kind": "sphere", "center": [0., 0., 0.], "radius": 1.}},
                                probe="QQ", samples=33)["earth"]
        np.testing.assert_allclose(image["center_px_py"], [0., 0.])
        self.assertEqual(image["px_py"].shape, (33, 2))
        look = fov_to_vector(*image["px_py"].T, direction="look")
        np.testing.assert_allclose(-look[:, 2], np.sqrt(1 - 1 / 100))
        self.assertTrue(image["all_in_view"])
        behind = obj_in_view("2025-03-10", [0., 0., -10.], [0., 0., 0., 1.],
                             {"earth": {"kind": "sphere", "center": [0., 0., 0.], "radius": 1.}},
                             probe="QQ")["earth"]
        self.assertTrue(np.isnan(behind["px_py"]).all())
        self.assertFalse(behind["any_in_view"])

    def test_moon_blocking_selects_positions_and_attitudes_from_the_same_rows(self):
        dates = pd.date_range("2025-03-10", periods=5, freq="1min")
        pose = ParameterDataFrame({
            "Datetime": dates, "MJ2000_Star_X": [0., np.nan, 0., 0., 0.],
            "MJ2000_Star_Y": 0., "MJ2000_Star_Z": 10 * RE,
            "MJ2000_q1": [0., 0., np.nan, 0., 0.],
            "MJ2000_q2": [0., 0., 0., 0., 1.], "MJ2000_q3": 0.,
            "PosePresent": [True, True, True, False, True],
        })
        for probe in ("A", "B", "QQ"):
            with self.subTest(probe=probe):
                mask = moon_blocking_mask(pose, probe=probe, gap=1)
                np.testing.assert_array_equal(mask["Datetime"], dates[[0, 4]])
                self.assertEqual(mask["partial_block"].tolist(), [True, False])
                self.assertEqual(mask["full_block"].tolist(), [True, False])
        intervals = moon_blocking_intervals(pose, probe="QQ", gap=1)
        self.assertEqual(intervals, [(dates[0], dates[0])])

    def test_export_uses_position_arrays_with_prepared_offset(self):
        from gena.post.export import dump

        events = self.events.copy()
        events[["MJ2000_Star_X", "MJ2000_Star_Y", "MJ2000_Star_Z"]] = np.tile([RE, 0., 0.], (4, 1))
        with patch("libs.ephemeris.moon", side_effect=AssertionError("unexpected ephemeris loading")):
            result = dump(events, coord="ECI2000", pos_offset=[60., 0., 0.])
        np.testing.assert_allclose(result["QQ_X_ECI2000"], 61.)
        self.assertEqual(len(result), len(events))


if __name__ == "__main__":
    unittest.main()
