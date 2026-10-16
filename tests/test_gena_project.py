import numpy as np
import pandas as pd
import pytest

from gena.frame import ParticleDataFrame
from gena.post import project


@pytest.mark.parametrize("probe", ["A", "B"])
def test_event_vectors_are_non_mutating_for_both_probes(probe):
    frame = ParticleDataFrame(
        {
            "Timestamp": [1735689600.0],
            "Datetime": pd.to_datetime([1735689600.0], unit="s"),
            "Energy": [10.0],
            "Spec": [1],
            "X1": [0.0],
            "Y1": [0.0],
            "X2": [1.0],
            "Y2": [1.0],
            "Tof": [100.0],
            "MJ2000_q1": [0.1],
            "MJ2000_q2": [0.1],
            "MJ2000_q3": [0.1],
        }
    )

    directed = project.event_vectors(frame, coord="GSM", probe=probe)

    assert not {"DX", "VX"}.intersection(frame.columns)
    assert {"DX", "DY", "DZ", "VX", "VY", "VZ"}.issubset(directed.columns)
    assert np.linalg.norm(directed.loc[0, ["DX", "DY", "DZ"]].to_numpy(float)) == pytest.approx(1)
    assert directed.attrs["probe"] == probe


def test_projection_is_invariant_to_source_versus_flight_direction():
    start = np.array([[12.0, 8.0, -3.0], [15.0, -4.0, 2.0]])
    source_look = np.array([[2.0, 1.0, -50.0], [-3.0, 2.0, -50.0]])
    particle_flight = -source_look

    projected_from_look = project.project(
        start, source_look, plane_normal=[1, 0, 0], plane_point=[-65, 0, 0]
    )
    projected_from_flight = project.project(
        start, particle_flight, plane_normal=[1, 0, 0], plane_point=[-65, 0, 0]
    )

    assert projected_from_flight == pytest.approx(projected_from_look)


def test_axial_mean_angle_treats_180_degree_separated_axes_as_equal():
    mean = project.axial_mean_angle([89.0, -89.0])

    assert abs(mean) == pytest.approx(90.0)


def test_axial_mean_angle_supports_event_weights():
    mean = project.axial_mean_angle([-10.0, 20.0], weights=[3.0, 1.0])
    expected = 0.5 * np.rad2deg(
        np.arctan2(
            3.0 * np.sin(np.deg2rad(-20.0)) + np.sin(np.deg2rad(40.0)),
            3.0 * np.cos(np.deg2rad(-20.0)) + np.cos(np.deg2rad(40.0)),
        )
    )

    assert mean == pytest.approx(expected)
