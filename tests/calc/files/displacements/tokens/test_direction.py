import numpy as np
import pytest

from viperleed.calc.files.new_displacements.tokens.direction import (
    DirectionToken,
    DirectionTokenParserError,
)


def test_basis_vector_single():
    d = DirectionToken('x')
    assert d.dof == 1
    assert d.vectors[0] == pytest.approx([0, 1, 0])


def test_basis_vector_plane():
    d = DirectionToken('xy')
    assert d.dof == 2
    assert d.vectors[0] == pytest.approx([0, 1, 0])
    assert d.vectors[1] == pytest.approx([0, 0, 1])


def test_full_space():
    d = DirectionToken('xyz')
    assert d.dof == 3
    expected = np.eye(3)[:, [2, 0, 1]] # zxy order!
    assert d.vectors == pytest.approx(expected)


def test_vector_direction_xy():
    d = DirectionToken('xy[1 1]')
    assert d.dof == 1
    expected = np.array([0, 1, 1]) / np.sqrt(2)
    assert d.vectors[0] == pytest.approx(expected)


def test_vector_direction_xyz():
    d = DirectionToken('xyz[1 2 3]')
    norm = np.linalg.norm([3, 1, 2])
    expected = np.array([3, 1, 2]) / norm
    assert d.vectors[0] == pytest.approx(expected)


def test_invalid_mismatch_components():
    with pytest.raises(ValueError):
        DirectionToken('xy[1 2 3]')


def test_invalid_labels():
    with pytest.raises(ValueError):
        DirectionToken('xzq')


def test_invalid_vector_format():
    with pytest.raises(ValueError):
        DirectionToken('xy[1,2]')  # invalid separator


def test_zero_vector():
    with pytest.raises(ValueError):
        DirectionToken('xyz[0 0 0]')

@pytest.mark.parametrize(
    'direction_str',
    ['a', 'b', 'c', 'ab', 'abc', 'ab[1 0]', 'abc[1 0 0]', 'ab[1 2]']
)
def test_invalid_fractional(direction_str):
    """Test that fractional directions are not accepted."""
    with pytest.raises(DirectionTokenParserError):
        DirectionToken(direction_str)

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
    with pytest.raises(DirectionTokenParserError):
        DirectionToken(direction_str)
