"""Tests for module viperleed.calc.files.new_displacements.tokens.direction."""

__authors__ = ('Alexandra Mia Imre (@alexmiame)',)
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2024-10-14'
__license__ = 'GPLv3+'

import numpy as np
import pytest

from viperleed.calc.files.new_displacements.tokens.direction import (
    CartesianDirectionToken,
    DirectionTokenParserError,
    FractionalDirectionToken,
    make_direction_token,
)


def test_basis_vector_single():
    d = CartesianDirectionToken('x')
    assert d.dof == 1
    assert d.vectors_zxy[0] == pytest.approx([0, 1, 0])
    assert d.vectors_xyz[0] == pytest.approx([1, 0, 0])


def test_basis_vector_plane():
    d = CartesianDirectionToken('xy')
    assert d.dof == 2
    assert d.vectors_zxy[0] == pytest.approx([0, 1, 0])
    assert d.vectors_zxy[1] == pytest.approx([0, 0, 1])


def test_full_space():
    d = CartesianDirectionToken('xyz')
    assert d.dof == 3
    expected = np.eye(3)[:, [2, 0, 1]]  # zxy order!
    assert d.vectors_zxy == pytest.approx(expected)


def test_vector_direction_xy():
    d = CartesianDirectionToken('xy[1 1]')
    assert d.dof == 1
    expected = np.array([0, 1, 1]) / np.sqrt(2)
    assert d.vectors_zxy[0] == pytest.approx(expected)


def test_vector_direction_xyz():
    d = CartesianDirectionToken('xyz[1 2 3]')
    norm = np.linalg.norm([3, 1, 2])
    expected = np.array([3, 1, 2]) / norm
    assert d.vectors_zxy[0] == pytest.approx(expected)


def test_invalid_mismatch_components():
    with pytest.raises(ValueError):
        CartesianDirectionToken('xy[1 2 3]')


def test_invalid_labels():
    with pytest.raises(ValueError):
        CartesianDirectionToken('xzq')


def test_invalid_vector_format():
    with pytest.raises(ValueError):
        CartesianDirectionToken('xy[1,2]')  # invalid separator


def test_zero_vector():
    with pytest.raises(ValueError):
        CartesianDirectionToken('xyz[0 0 0]')


@pytest.mark.parametrize(
    ('direction_str', 'match'),
    [
        ('c', 'c unit-cell vector'),
        ('abc', 'c unit-cell vector'),
        ('bc', 'c unit-cell vector'),
        ('ac', 'c unit-cell vector'),
        ('abc[1 0 0]', 'c unit-cell vector'),
        ('az', 'mix Cartesian'),
        ('xb', 'mix Cartesian'),
        ('xa[1 0]', 'mix Cartesian'),
        ('azi(ab[1 0])', 'not yet supported'),
        ('r([1 0])', 'not yet supported'),
    ],
)
def test_invalid_fractional(direction_str, match):
    """Test that unsupported fractional directions raise clearly."""
    with pytest.raises(DirectionTokenParserError, match=match):
        make_direction_token(direction_str)


@pytest.mark.parametrize(
    'direction_str',
    ['ab[1 2 3]', 'ab[1]', 'ab[0 0]', 'ab[1,2]', 'ba[1 0]', 'x'],
)
def test_invalid_fractional_token(direction_str):
    """Test FractionalDirectionToken rejects invalid or other input."""
    with pytest.raises(DirectionTokenParserError):
        FractionalDirectionToken(direction_str)


@pytest.mark.parametrize('direction_str', ['a', 'b', 'ab', 'ab[1 0]', '[1 0]'])
def test_cartesian_rejects_fractional(direction_str):
    """Test that CartesianDirectionToken rejects fractional input."""
    with pytest.raises(DirectionTokenParserError, match='fractional'):
        CartesianDirectionToken(direction_str)


class TestFractional:
    """Tests for FractionalDirectionToken and its dispatching."""

    @pytest.mark.parametrize(
        ('direction_str', 'expect_type'),
        [
            ('x', CartesianDirectionToken),
            ('xy[1 0]', CartesianDirectionToken),
            ('a', FractionalDirectionToken),
            ('ab', FractionalDirectionToken),
            ('ab[1 2]', FractionalDirectionToken),
            ('[1 2]', FractionalDirectionToken),
            (' [0.5 1] ', FractionalDirectionToken),
        ],
    )
    def test_dispatch(self, direction_str, expect_type):
        """Test that make_direction_token picks the right class."""
        assert isinstance(make_direction_token(direction_str), expect_type)

    @pytest.mark.parametrize(
        ('first', 'second'),
        [
            ('a', 'ab[1 0]'),
            ('a', '[1 0]'),
            ('ab[1 0]', '[1 0]'),
            ('b', 'ab[0 1]'),
            ('b', '[0 1]'),
            ('ab[2 4]', '[1 2]'),
            ('ab[0.5 1]', '[1 2]'),
        ],
    )
    def test_equal(self, first, second):
        """Test equality of equivalent fractional directions."""
        assert make_direction_token(first) == make_direction_token(second)

    @pytest.mark.parametrize(
        ('first', 'second'),
        [
            ('a', 'b'),
            ('a', 'ab'),
            ('a', '[-1 0]'),
            ('a', 'x'),
            ('ab', 'xy'),
        ],
    )
    def test_not_equal(self, first, second):
        """Test inequality of different directions."""
        assert make_direction_token(first) != make_direction_token(second)

    @pytest.mark.parametrize(
        ('direction_str', 'dof'),
        [('a', 1), ('b', 1), ('ab', 2), ('ab[0.5 1]', 1), ('[1 -1]', 1)],
    )
    def test_dof(self, direction_str, dof):
        """Test the number of degrees of freedom."""
        assert FractionalDirectionToken(direction_str).dof == dof

    def test_float_components(self):
        """Test that floating-point components are parsed."""
        token = FractionalDirectionToken('ab[0.5 1]')
        expected = np.array([0.5, 1]) / np.linalg.norm([0.5, 1])
        assert token.components[0] == pytest.approx(expected)

    _hex_ab_cell = np.array([[2.77, -1.385], [0, 2.399]])  # a, b columns

    @pytest.mark.parametrize('comps', [(1, 2), (1, 0), (0, 1), (2, -1)])
    def test_hexagonal_matches_old_reader(self, comps):
        """Check the Cartesian vector against the old reader."""
        token = FractionalDirectionToken(f'ab[{comps[0]} {comps[1]}]')
        # As in calc.files.displacements.readDISPLACEMENTS_block
        abst = self._hex_ab_cell.T
        old = comps[0]*abst[0] + comps[1]*abst[1]
        old = old / np.linalg.norm(old)
        vectors = token.get_vectors_xyz(self._hex_ab_cell)
        assert vectors.shape == (1, 3)
        assert vectors[0] == pytest.approx([*old, 0])
        assert token.get_vectors_zxy(self._hex_ab_cell)[0] == pytest.approx(
            [0, *old]
        )

    def test_plane_vectors(self):
        """Check that 'ab' gives an orthonormal basis of the plane."""
        vectors = FractionalDirectionToken('ab').get_vectors_xyz(
            self._hex_ab_cell
        )
        assert vectors.shape == (2, 3)
        assert vectors @ vectors.T == pytest.approx(np.eye(2))
        assert vectors[:, 2] == pytest.approx([0, 0])


@pytest.mark.parametrize(
    ('first', 'second'),
    [('x', 'xy[1 0]'), ('y', 'xy[0 1]'), ('z', 'xyz[0 0 1]')],
)
def test_cartesian_shorthands_equal(first, second):
    """Test #413: 'x' == 'xy[1 0]' and 'y' == 'xy[0 1]'."""
    assert CartesianDirectionToken(first) == CartesianDirectionToken(second)


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
        CartesianDirectionToken(direction_str)
