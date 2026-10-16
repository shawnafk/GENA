"""Stateless FOV representations and angle conventions."""

import numpy as np
import pytest

from common.coord import angle_between, angle_to_plane, fov_to_vector, vector_to_fov


@pytest.mark.parametrize('direction', ['look', 'propagation'])
def test_fov_conventions_round_trip_and_sample_basis(direction):
    angles = np.array([-10., 0., 12., 20.])
    vectors = fov_to_vector(angles, angles / 2, direction=direction)
    px, py = vector_to_fov(vectors, direction=direction)
    assert px == pytest.approx(angles)
    assert py == pytest.approx(angles / 2)
    basis = np.tile(np.diag([1., 1., -1. if direction == 'look' else 1.]), (4, 1, 1))
    assert vector_to_fov(vectors, basis)[0] == pytest.approx(angles)
    assert fov_to_vector(0, 0, direction=direction).shape == (3,)
    with pytest.raises(ValueError):
        fov_to_vector(90, 0)


def test_angle_units_signs_and_invalid_vectors():
    assert angle_between([1, 0, 0], [0, 1, 0]) == pytest.approx(90.)
    assert angle_between([1, 0, 0], [0, 1, 0], degrees=False) == pytest.approx(np.pi / 2)
    assert angle_to_plane([0, 0, -1], [0, 0, 1]) == pytest.approx(-90.)
    assert angle_to_plane([0, 0, -1], [0, 0, 1], signed=False) == pytest.approx(90.)
    with pytest.raises(ValueError, match='nonzero'):
        angle_between([0, 0, 0], [1, 0, 0])




@pytest.mark.parametrize('single', [True, False])
def test_rotation_matrix_and_quaternion_apply_same_inverse(single):
    from scipy.spatial.transform import Rotation
    from common.coord import rotate_vectors, quaternion_array

    rotation = Rotation.from_euler('z', 35., degrees=True)
    vectors = np.array([1., 2., 3.]) if single else np.array([[1., 2., 3.], [4., 5., 6.]])
    expected = rotation.apply(vectors)
    actual = rotate_vectors(vectors, rotation.as_matrix())
    assert actual.shape == vectors.shape
    assert actual == pytest.approx(expected)
    assert rotate_vectors(vectors, rotation) == pytest.approx(expected)
    assert rotate_vectors(vectors, quaternion_array(rotation)) == pytest.approx(expected)
    assert rotate_vectors(actual, rotation.as_matrix(), inverse=True) == pytest.approx(vectors)
    assert rotate_vectors(actual, rotation.as_quat(), inverse=True) == pytest.approx(vectors)


def test_quaternion_interpolation_and_validation():
    from scipy.spatial.transform import Rotation
    from common.coord import quaternion_array, interpolate_quaternions

    rotation = Rotation.from_euler('z', [0., 90.], degrees=True)
    result = interpolate_quaternions([0., 10.], rotation, [10., 5., 0., 5.])
    assert Rotation.from_quat(result).apply([1., 0., 0.]) == pytest.approx(
        Rotation.from_euler('z', [90., 45., 0., 45.], degrees=True).apply([1., 0., 0.]))
    assert interpolate_quaternions([0., 10.], [0., 0., 0., 2.], [5.]) == pytest.approx(np.array([[0., 0., 0., 1.]]))
    with pytest.raises(ValueError, match='nonzero'):
        quaternion_array([0., 0., 0., 0.])
    with pytest.raises(ValueError, match='coverage'):
        interpolate_quaternions([0., 10.], rotation, [11.])
    with pytest.raises(ValueError, match='sample count'):
        quaternion_array(rotation, count=3)
