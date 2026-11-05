"""Object geometry and viewing share one projection path without data loading."""

import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd
from scipy.spatial.transform import Rotation

from common.constants import RE, RM
from common.coord import fov_to_vector, transform_vectors, vector_to_fov
from gena.frame import ParameterDataFrame, ParticleDataFrame
from gena.post.project import dipole_fieldlines, moon_blocking_mask, obj_in_view, sphere_limb


class ObjectViewTests(unittest.TestCase):
    def test_sphere_limb_is_on_surface_and_tangent_to_observer(self):
        center, observer = np.array([3., 4., 5.]), np.array([10., 7., 20.])
        limb = sphere_limb(center, 2., observer, samples=65)
        np.testing.assert_allclose(np.linalg.norm(limb - center, axis=1), 2.)
        np.testing.assert_allclose(np.einsum("ij,ij->i", limb - center, observer - limb), 0., atol=1e-13)
        np.testing.assert_allclose(limb[0], limb[-1])

    def test_mixed_objects_reuse_one_pose_and_preserve_definitions(self):
        position = np.array([2., 3., 10.])
        targets = {
            "earth": {"kind": "sphere", "center": position + [0., 0., -10.], "radius": 1.},
            "moon": {"kind": "sphere", "center": position + [0., 0., -5.], "radius": RM / RE},
            "satellites": {"kind": "points", "points": position + [[0.1, 0.05, -2.], [0., 0., 2.], [np.nan, 0., 0.]]},
            "tail": {"kind": "curve", "points": position + [[-3., 0., -10.], [3., 0., -10.]], "label": "prepared tail"},
        }
        for probe in ("A", "B"):
            with self.subTest(probe=probe), patch("libs.ephemeris.moon", side_effect=AssertionError("unexpected ephemeris")):
                attitude = Rotation.from_euler("z", 90., degrees=True)
                view = obj_in_view("2025-03-10", position, attitude, targets, probe=probe)
                sight = transform_vectors([[0.1, 0.05, -2.]], source="ECI2000", target=probe,
                                          attitude=attitude)
                np.testing.assert_allclose(view["satellites"]["px_py"][0], np.array(vector_to_fov(sight[0], direction="look")))
                np.testing.assert_array_equal(view["satellites"]["in_view"], [True, False, False])
                self.assertTrue(view["earth"]["any_in_view"])
                self.assertFalse(view["earth"]["all_in_view"])
                self.assertTrue(view["moon"]["all_in_view"])
                self.assertEqual(view["tail"]["label"], "prepared tail")
        np.testing.assert_allclose(targets["tail"]["points"], position + [[-3., 0., -10.], [3., 0., -10.]])
        self.assertNotIn("px_py", targets["earth"])

    def test_curve_crossing_is_clipped_even_with_both_samples_outside(self):
        targets = {"tail": {"kind": "curve", "points": [[-10., 0., -10.], [10., 0., -10.]]}}
        line = obj_in_view("2025-03-10", [0., 0., 0.], Rotation.identity(), targets,
                           probe="QQ", ax_deg=10., ay_deg=10.)["tail"]
        self.assertFalse(line["in_view"].any())
        self.assertTrue(line["any_in_view"])
        self.assertFalse(line["all_in_view"])
        np.testing.assert_allclose(line["segments_px_py"], [[[-10., 0.], [10., 0.]]], atol=1e-12)

    def test_large_disk_covers_view_without_any_limb_sample_inside(self):
        targets = {"body": {"kind": "sphere", "center": [0., 0., -10.], "radius": 7.}}
        for samples in (4, 181):
            body = obj_in_view("2025-03-10", [0., 0., 0.], Rotation.identity(), targets,
                               probe="QQ", ax_deg=10., ay_deg=10., samples=samples)["body"]
            self.assertFalse(body["in_view"].any())
            self.assertTrue(body["any_in_view"])
            self.assertFalse(body["all_in_view"])
            self.assertTrue(body["covers_view"])

    def test_sphere_overlap_handles_edges_and_corners_independently_of_sampling(self):
        targets = {
            "edge": {"kind": "sphere", "center": 10 * fov_to_vector(15., 0., direction="look"),
                     "radius": 10 * np.sin(np.deg2rad(6.))},
            "corner_miss": {"kind": "sphere", "center": 10 * fov_to_vector(14., 14., direction="look"),
                            "radius": 10 * np.sin(np.deg2rad(4.5))},
            "corner_hit": {"kind": "sphere", "center": 10 * fov_to_vector(14., 14., direction="look"),
                           "radius": 10 * np.sin(np.deg2rad(6.))},
        }
        for samples in (4, 181):
            view = obj_in_view("2025-03-10", [0., 0., 0.], Rotation.identity(), targets,
                               probe="QQ", ax_deg=10., ay_deg=10., samples=samples)
            self.assertTrue(view["edge"]["any_in_view"])
            self.assertFalse(view["edge"]["all_in_view"])
            self.assertFalse(view["corner_miss"]["any_in_view"])
            self.assertTrue(view["corner_hit"]["any_in_view"])

    def test_dipole_targets_keep_sm_geometry_and_convert_at_view_time(self):
        targets = dipole_fieldlines(l_shells=(4,), mlt_hours=(0,), samples=33)
        obj = targets["dipole_L4_MLT0"]
        original = obj["points"].copy()
        np.testing.assert_allclose(original[16], [-4., 0., 0.], atol=1e-14)
        xyz = transform_vectors(original, source="SM", target="ECI2000", times="2025-03-10")
        expected = np.column_stack(vector_to_fov(xyz - [0., 0., 10.], direction="look"))
        view = obj_in_view("2025-03-10", [0., 0., 10.], Rotation.identity(), targets, probe="QQ")
        np.testing.assert_allclose(view["dipole_L4_MLT0"]["px_py"], expected)
        np.testing.assert_array_equal(obj["points"], original)

    def test_moon_mask_reports_disk_containment_and_full_view_coverage_separately(self):
        pose = ParameterDataFrame({
            "Datetime": ["2025-03-10"], "MJ2000_Star_X": [0.], "MJ2000_Star_Y": [0.],
            "MJ2000_Star_Z": [1.1 * RM], "MJ2000_q1": [0.], "MJ2000_q2": [0.],
            "MJ2000_q3": [0.], "PosePresent": [True],
        })
        mask = moon_blocking_mask(pose, probe="QQ", gap=1, samples=4)
        self.assertTrue(mask["partial_block"].iloc[0])
        self.assertFalse(mask["full_block"].iloc[0])
        self.assertTrue(mask["covers_view"].iloc[0])

    def test_fov_plot_accepts_generic_curves_and_point_targets(self):
        from matplotlib import pyplot as plt
        from plots.gena import fov

        targets = {
            "tail": {"kind": "curve", "points": [[-10., 0., -10.], [10., 0., -10.]],
                     "style": {"color": "cyan"}},
            "satellite": {"kind": "points", "points": [0., 0., -2.]},
        }
        view = obj_in_view("2025-03-10", [0., 0., 0.], Rotation.identity(), targets, probe="QQ")
        events = ParticleDataFrame({"Datetime": pd.date_range("2025-03-10", periods=2, freq="1s"),
                                    "PX": [0., 1.], "PY": [0., 1.]})
        fig, ax = fov(events, projection=view)
        try:
            fig.canvas.draw()
            self.assertEqual(len(ax.collections), 2)
            self.assertEqual(len(ax.lines), 1)
        finally:
            plt.close(fig)


if __name__ == "__main__":
    unittest.main()
