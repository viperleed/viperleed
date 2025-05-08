import numpy as np
import pytest

from viperleed_jax.files.displacements.direction import (
    Direction,
    UnsupportedDirectionError,
)


def test_basis_vector_single():
    d = Direction('x')
    assert d.dof == 1
    np.testing.assert_allclose(d.vectors[0], [1, 0, 0])


def test_basis_vector_plane():
    d = Direction('xy')
    assert d.dof == 2
    np.testing.assert_allclose(d.vectors[0], [1, 0, 0])
    np.testing.assert_allclose(d.vectors[1], [0, 1, 0])


def test_full_space():
    d = Direction('xyz')
    assert d.dof == 3
    expected = np.eye(3)
    np.testing.assert_allclose(d.vectors, expected)


def test_vector_direction_xy():
    d = Direction('xy[1 1]')
    assert d.dof == 1
    expected = np.array([1, 1, 0]) / np.sqrt(2)
    np.testing.assert_allclose(d.vectors[0], expected)


def test_vector_direction_xyz():
    d = Direction('xyz[1 2 3]')
    norm = np.linalg.norm([1, 2, 3])
    expected = np.array([1, 2, 3]) / norm
    np.testing.assert_allclose(d.vectors[0], expected)


def test_invalid_mismatch_components():
    with pytest.raises(ValueError):
        Direction('xy[1 2 3]')


def test_invalid_labels():
    with pytest.raises(ValueError):
        Direction('xzq')


def test_invalid_vector_format():
    with pytest.raises(ValueError):
        Direction('xy[1,2]')  # invalid separator


def test_zero_vector():
    with pytest.raises(ValueError):
        Direction('xyz[0 0 0]')

@pytest.mark.parametrize(
    'direction_str',
    ['a', 'b', 'c', 'ab', 'abc', 'ab[1 0]', 'abc[1 0 0]', 'ab[1 2]']
)
def test_invalid_fractional(direction_str):
    """Test that fractional directions are not accepted."""
    with pytest.raises(UnsupportedDirectionError):
        Direction(direction_str)

@pytest.mark.parametrize(
    'direction_str',
    [
        'azi(xy[1 1])',
        'azi(ab[1 0])',
        'r(xy[1 1])',
        'r(ab[1 1])',
    ],
)
def test_invalid_radial_azimuthal(direction_str):
    """Test that the old radial and azimuthal directions are not accepted."""
    with pytest.raises(UnsupportedDirectionError):
        Direction(direction_str)
