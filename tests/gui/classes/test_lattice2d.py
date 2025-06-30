"""Tests for module lattice2d of viperleed.gui.classes."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2020-01-11'
__license__ = 'GPLv3+'

import numpy as np
import pytest
from pytest_cases import parametrize

from viperleed.calc.lib.math_utils import angle as angle_between
from viperleed.gui.classes import planegroup
from viperleed.gui.classes.lattice2d import Lattice2D

_HEX_ACUTE = (1, 0), (0.5, 3**0.5/2)
_HEX_OBTUSE = (1, 0), (-0.5, 3**0.5/2)
_OBLIQUE_ACUTE = (1, 0), (1.2, 3.5)
_OBLIQUE_OBTUSE = (1, 0), (-1.2, 3.5)
_RECT = (1, 0), (0, 3.5)
_RHOMBIC_ACUTE = (2, 1), (2, -1)
_RHOMBIC_OBTUSE = (-2, 1), (2, 1)
_SQUARE = (1, 0), (0, 1)


def test_str_repr():
    """Check expected outcome of __str__ and __repr__ methods."""
    basis = _SQUARE
    lattice = Lattice2D(basis)
    _repr = "Lattice2D([[1,0], [0,1]], space='real', group='p1', limit=1)"
    _str = "Square, real-space Lattice2D([[1,0], [0,1]], group=p1)"
    assert repr(lattice) == _repr
    assert str(lattice) == _str


class TestProperties:
    """Collection of tests for accessing properties."""

    _basis = {
        'hex': _HEX_OBTUSE,
        'oblique': _OBLIQUE_OBTUSE,
        'rect': _RECT,
        'rhombic': _RHOMBIC_OBTUSE,
        'square': _SQUARE,
        }

    @parametrize(basis=_basis.values(), ids=_basis)
    def test_basis(self, basis):
        """Check correct setting of lattice basis."""
        lattice = Lattice2D(basis)
        assert lattice.space == 'real'
        assert np.all(lattice.basis == basis)

    _cell = {
        _HEX_ACUTE: 'Hexagonal', _HEX_OBTUSE: 'Hexagonal',
        _OBLIQUE_ACUTE: 'Oblique', _OBLIQUE_OBTUSE: 'Oblique',
        _RECT: 'Rectangular',
        _RHOMBIC_ACUTE: 'Rhombic', _RHOMBIC_OBTUSE: 'Rhombic',
        _SQUARE: 'Square',
        }

    @parametrize('basis,shape', _cell.items(), ids=(str(c) for c in _cell))
    def test_cell_shape(self, basis, shape):
        """Check correct identification of the cell shape."""
        assert Lattice2D(basis).cell_shape == shape

    def test_group_upon_basis_change(self, subtests):
        """Check that the group becomes p1 when a new basis shape is given."""
        with subtests.test('upon construction'):
            lattice = Lattice2D(_RECT, group='pmm')
            assert lattice.group == 'pmm'
        with subtests.test('rect->square'):
            lattice.basis = _SQUARE
            assert lattice.group == 'pmm'
        with subtests.test('square->rhombic'):
            lattice.basis = _RHOMBIC_OBTUSE
            assert lattice.group == 'p1'

    _group = {
        None: 'p1',  # Default
        'pgg': 'pgg'
        }

    @parametrize('group,expect', _group.items(), ids=_group)
    def test_group_valid(self, group, expect):
        """Check assignment of a valid plane group."""
        lattice = Lattice2D(_RECT)
        lattice.group = group
        assert lattice.group == expect

    _params = {
        (_HEX_ACUTE, 'real'): (1, 1, 120),
        (_HEX_OBTUSE, 'real'): (1, 1, 120),
        (_HEX_ACUTE, 'reciprocal'): (1, 1, 60),
        (_HEX_OBTUSE, 'reciprocal'): (1, 1, 60),
        (_OBLIQUE_ACUTE, 'real'): (1, 3.7, 71.075),
        (_OBLIQUE_ACUTE, 'reciprocal'): (1, 3.7, 71.075),
        (_OBLIQUE_OBTUSE, 'real'): (1, 3.7, 108.925),
        (_OBLIQUE_OBTUSE, 'reciprocal'): (1, 3.7, 108.925),
        (_RECT, 'real'): (1, 3.5, 90),
        (_RECT, 'reciprocal'): (1, 3.5, 90),
        (_RHOMBIC_ACUTE, 'real'): (5**0.5, 5**0.5, 126.87),
        (_RHOMBIC_OBTUSE, 'real'): (5**0.5, 5**0.5, 126.87),
        (_RHOMBIC_ACUTE, 'reciprocal'): (5**0.5, 5**0.5, 53.13),
        (_RHOMBIC_OBTUSE, 'reciprocal'): (5**0.5, 5**0.5, 53.13),
        (_SQUARE, 'real'): (1, 1, 90),
        (_SQUARE, 'reciprocal'): (1, 1, 90),
        }

    @parametrize('args,expect', _params.items(), ids=(str(p) for p in _params))
    def test_lattice_parameters(self, args, expect):
        """Check correct values of vector lengths and angle."""
        basis, space = args
        lat = Lattice2D(basis, space=space)
        assert lat.lattice_parameters == pytest.approx(expect, abs=1e-3)
        assert np.linalg.det(basis) == np.linalg.det(lat.basis)

    _pts = {'real': ('n_points', 19),
            'reciprocal': ('n_beams', 7)}

    @parametrize('space,expect', _pts.items(), ids=_pts)
    def test_lattice_pts(self, space, expect):
        """Check correct number of lattice points/beams."""
        lattice = Lattice2D(_RECT, space=space, limit=3)
        attr, value = expect
        assert getattr(lattice, attr) == value

    _real_rec = {  # (self.basis, real_basis, reciprocal_basis)
        (_HEX_ACUTE, 'real'): (
            ((-0.5, -3**0.5/2), (1, 0)),  # OBTUSE, swap a'=-b, b'=a
            ((-0.5, -3**0.5/2), (1, 0)),  # OBTUSE, swap a'=-b, b'=a
            2*np.pi * np.array(((0, -2/3**0.5), (1, -1/3**0.5))),
            ),
        (_HEX_OBTUSE, 'real'): (
            _HEX_OBTUSE,
            _HEX_OBTUSE,
            2*np.pi * np.array(((1, 1/3**0.5), (0, 2/3**0.5))),
            ),
        (_HEX_ACUTE, 'reciprocal'): (
            _HEX_ACUTE,
            2*np.pi * np.array(((1, -1/3**0.5), (0, 2/3**0.5))),
            2*np.pi * np.array(((1, -1/3**0.5), (0, 2/3**0.5))),
            ),
        (_HEX_OBTUSE, 'reciprocal'): (
            ((0.5, -3**0.5/2), (1, 0)),  # OBTUSE, swap a'=-b, b'=a
            2*np.pi * np.array(((0, -2/3**0.5), (1, 1/3**0.5))),
            2*np.pi * np.array(((0, -2/3**0.5), (1, 1/3**0.5))),
            ),
        (_OBLIQUE_ACUTE, 'real'): (
            _OBLIQUE_ACUTE,
            _OBLIQUE_ACUTE,
            2*np.pi/3.5 * np.array(((3.5, -1.2), (0, 1))),
            ),
        (_OBLIQUE_ACUTE, 'reciprocal'): (
            _OBLIQUE_ACUTE,
            2*np.pi/3.5 * np.array(((3.5, -1.2), (0, 1))),
            2*np.pi/3.5 * np.array(((3.5, -1.2), (0, 1))),
            ),
        (_OBLIQUE_OBTUSE, 'real'): (
            _OBLIQUE_OBTUSE,
            _OBLIQUE_OBTUSE,
            2*np.pi/3.5 * np.array(((3.5, 1.2), (0, 1))),
            ),
        (_OBLIQUE_OBTUSE, 'reciprocal'): (
            _OBLIQUE_OBTUSE,
            2*np.pi/3.5 * np.array(((3.5, 1.2), (0, 1))),
            2*np.pi/3.5 * np.array(((3.5, 1.2), (0, 1))),
            ),
        (_RECT, 'real'): (
            _RECT,
            _RECT,
            2*np.pi * np.array(((1, 0), (0, 1/3.5))),
            ),
        (_RECT, 'reciprocal'): (
            _RECT,
            2*np.pi * np.array(((1, 0), (0, 1/3.5))),
            2*np.pi * np.array(((1, 0), (0, 1/3.5))),
            ),
        (_RHOMBIC_ACUTE, 'real'): (
            _RHOMBIC_OBTUSE,
            _RHOMBIC_OBTUSE,
            2*np.pi/4 *np.array(((-1, 2), (1, 2))),
            ),
        (_RHOMBIC_OBTUSE, 'real'): (
            _RHOMBIC_OBTUSE,
            _RHOMBIC_OBTUSE,
            2*np.pi/4 *np.array(((-1, 2), (1, 2))),
            ),
        (_RHOMBIC_ACUTE, 'reciprocal'): (
            _RHOMBIC_ACUTE,
            2*np.pi/4 *np.array(((1, 2), (1, -2))),
            2*np.pi/4 *np.array(((1, 2), (1, -2))),
            ),
        (_RHOMBIC_OBTUSE, 'reciprocal'): (
            ((-2, -1), (-2, 1)),
            2*np.pi/4 *np.array(((-1, -2), (-1, 2))),
            2*np.pi/4 *np.array(((-1, -2), (-1, 2))),
            ),
        (_SQUARE, 'real'): (
            _SQUARE,
            _SQUARE,
            2*np.pi * np.array(_SQUARE),
            ),
        (_SQUARE, 'reciprocal'): (
            _SQUARE,
            2*np.pi * np.array(_SQUARE),
            2*np.pi * np.array(_SQUARE),
            ),
        }

    @parametrize('args,expect', _real_rec.items(),
                 ids=(str(p) for p in _real_rec))
    def test_real_reciprocal_bases(self, args, expect, subtests):
        basis, space = args
        expect_self, expect_real, expect_rec = expect
        lattice = Lattice2D(basis, space=space)
        with subtests.test('self basis'):
            assert np.allclose(lattice.basis, expect_self)
        with subtests.test('real basis'):
            assert np.allclose(lattice.real_basis, expect_real)
        with subtests.test('reciprocal basis'):
            assert np.allclose(lattice.reciprocal_basis, expect_rec)

    _acute = {
        (_HEX_ACUTE, 'real') : True,
        (_HEX_OBTUSE, 'real'): False,
        (_HEX_ACUTE, 'reciprocal'): False,
        (_HEX_OBTUSE, 'reciprocal'): True,
        (_OBLIQUE_ACUTE, 'real'): False,  # Only hex & rhombic!
        (_OBLIQUE_OBTUSE, 'real'): False,
        (_OBLIQUE_ACUTE, 'reciprocal'): False,
        (_OBLIQUE_OBTUSE, 'reciprocal'): False,
        (_RECT, 'real'): False,
        (_RECT, 'reciprocal'): False,
        (_RHOMBIC_ACUTE, 'real'): True,
        (_RHOMBIC_OBTUSE, 'real'): False,
        (_RHOMBIC_ACUTE, 'reciprocal'): False,
        (_RHOMBIC_OBTUSE, 'reciprocal'): True,
        (_SQUARE, 'real'): False,
        (_SQUARE, 'reciprocal'): False,
        }

    @parametrize('args,expect', _acute.items(), ids=(str(a) for a in _acute))
    def test_was_acute(self, args, expect):
        """Check correct value of property was_acute."""
        basis, space = args
        lat = Lattice2D(basis, space=space)
        assert lat.was_acute == expect


class TestRaises:
    """Collection of tests for complaints by Lattice2D objects."""

    _init = {
        'space': (_SQUARE, {'space': 'invalid_space'}, ValueError),
        'limit': (_SQUARE, {'limit': 'not_a_number'}, TypeError),
        'basis': ('not_a_list', {}, TypeError),
        'basis shape': (np.eye(3), {}, ValueError),
        'basis singular': (np.zeros((2, 2)), {}, ValueError),
        'group for shape': (_HEX_OBTUSE, {'group': 'p4'}, ValueError),
        'not a group': (_HEX_OBTUSE, {'group': 'invalid_group'}, ValueError),
        'group type': (_HEX_OBTUSE, {'group': 123}, ValueError),
        }

    @parametrize('basis,kwargs,exc', _init.values(), ids=_init)
    def test_init_invalid(self, basis, kwargs, exc):
        """Check complaints for making a Lattice2D with invalid arguments."""
        with pytest.raises(exc):
            Lattice2D(basis, **kwargs)

    _attr = {
        'n_beams': Lattice2D(_SQUARE, space='real'),
        'n_points': Lattice2D(_SQUARE, space='reciprocal')
        }

    @parametrize('attr,lattice', _attr.items(), ids=_attr)
    def test_attribute_error(self, lattice, attr):
        """Check complaints when accessing a property for the wrong space."""
        with pytest.raises(AttributeError):
            _ = getattr(lattice, attr)


class TestSpecialDirections:
    """Tests for determination of mirror directions."""

    @staticmethod
    def _any_mismatch(results, expected):
        """Check that results == expected, but items may be None or arrays."""
        if len(results) != len(expected):
            return f'length mismatch: {len(results)} vs. {len(expected)}'
        for i, (res_, exp_) in enumerate(zip(results, expected)):
            if res_ is None and exp_ is None:
                continue
            if any(v is None for v in (res_, exp_)):
                # Only one is None
                return f'mismatch at ind={i}: {res_}, {exp_}'
            if not np.allclose(res_, exp_):
                return f'mismatch at ind={i}: {res_}, {exp_}'
        return False

    _angles = {
        (_HEX_ACUTE, 'p6m'): [0, 120, 90, 30, 60, 150],
        (_HEX_OBTUSE, 'p6m'): [0, 120, 90, 30, 60, 150],
        (_OBLIQUE_ACUTE, 'p2'): [],
        (_OBLIQUE_OBTUSE, 'p2'): [],
        (_RECT, 'pmm'): [0, 90],
        (_RHOMBIC_ACUTE, 'cmm'): [116.565, 26.565],
        (_RHOMBIC_OBTUSE, 'cmm'): [116.565, 26.565],
        (_SQUARE, 'p4m'): [0, 90, 45, 135],
        }

    @parametrize('args,expect', _angles.items(), ids=(str(k) for k in _angles))
    def test_angles(self, args, expect):
        """Check correct orientation of the special directions."""
        basis, group = args
        lattice = Lattice2D(basis, group=group)
        angles = [np.degrees(angle_between(lattice.basis[0], d)) % 180
                  for d in lattice.special_directions
                  if d is not None]
        assert angles == pytest.approx(expect, abs=1e-3)

    _directions = {  # Enough to check the highest-symmetry group
        (_HEX_ACUTE, 'p6m'): [
            None, (-0.5, -3**0.5/2), (1, 0), (3**0.5/2, -0.5), (0, -1),
            (0.5, -3**0.5/2), (3**0.5/2, 0.5), None, None, None,
            None, None
            ],
        (_HEX_OBTUSE, 'p6m'): [
            None, (1, 0), (-0.5, 3**0.5/2), (0, 1), (3**0.5/2, 0.5),
            (0.5, 3**0.5/2), (-3**0.5/2, 0.5), None, None, None,
            None, None
            ],
        (_OBLIQUE_ACUTE, 'p2'): [None, None],
        (_OBLIQUE_OBTUSE, 'p2'): [None, None],
        (_RECT, 'pmm'): [None, (1, 0), (0, 1), None],
        (_RHOMBIC_ACUTE, 'cmm'): [None, None, (0, 1), (1, 0)],
        (_RHOMBIC_OBTUSE, 'cmm'): [None, None, (0, 1), (1, 0)],
        (_SQUARE, 'p4m'): [None, (1, 0), (0, 1), (2**-0.5, 2**-0.5),
                           (-2**-0.5, 2**-0.5), None, None, None],
        }

    @parametrize('args,expect', _directions.items(),
                 ids=(str(k) for k in _directions))
    def test_directions(self, args, expect):
        """Check exact match of all the directions."""
        basis, group = args
        lattice = Lattice2D(basis, group=group)
        directions = lattice.special_directions
        print(directions)
        assert not self._any_mismatch(directions, expect)


# Prepare some linear combinations for testing high_symm_transform.
LINEAR_COMBINATIONS = np.random.randint(-20, 20, size=20)
LINEAR_COMBINATIONS[LINEAR_COMBINATIONS==0] += 1


class TestTransforms:
    """Collection of tests for transformations of a Lattice2D."""

    _rot_basis = {
        'square 90': (_SQUARE, 90, ((0, 1), (-1, 0))),
        'rhomb -90': (_RHOMBIC_OBTUSE, -90, ((1, 2), (1, -2))),
        'hex 60': (_HEX_OBTUSE, 60, ((0.5, 3**0.5/2), (-1, 0))),
        }

    @parametrize('basis,angle,expect', _rot_basis.values(), ids=_rot_basis)
    def test_rotated_basis(self, basis, angle, expect):
        """Check correct result of get_rotated_basis(angle)."""
        lattice = Lattice2D(basis)
        rotated_basis = lattice.get_rotated_basis(angle)
        assert np.allclose(rotated_basis, expect)

    def test_rotated_lattice(self):
        """Check correct result of get_rotated_lattice_points(angle)."""
        lattice = Lattice2D(_RECT, space='real', limit=3)
        pts = lattice.lattice
        expect_pts = [[-4,   0. ], [-3,   0. ], [-2,  -3.5], [-2,   0. ],
                      [-2,   3.5], [-1,  -3.5], [-1,   0. ], [-1,   3.5],
                      [ 0,  -3.5], [ 0,   0. ], [ 0,   3.5], [ 1,  -3.5],
                      [ 1,   0. ], [ 1,   3.5], [ 2,  -3.5], [ 2,   0. ],
                      [ 2,   3.5], [ 3,   0. ], [ 4,   0. ]]
        expect_rot_pts = [[ 0.0, -4], [ 0.0, -3], [ 3.5, -2], [ 0.0, -2],
                          [-3.5, -2], [ 3.5, -1], [ 0.0, -1], [-3.5, -1],
                          [ 3.5,  0], [ 0.0,  0], [-3.5,  0], [ 3.5,  1],
                          [ 0.0,  1], [-3.5,  1], [ 3.5,  2], [ 0.0,  2],
                          [-3.5,  2], [ 0.0,  3], [ 0.0,  4]]
        assert np.allclose(pts, expect_pts)
        rot_pts = lattice.get_rotated_lattice_points(90)
        assert np.allclose(rot_pts, expect_rot_pts)

    def test_transform(self):
        """Check correct transformation of a lattice basis."""
        lattice = Lattice2D(_RECT)
        # Notice that here we should only use transformation
        # matrices that make sense for the current shape!
        lattice.transform(planegroup.Mx)
        assert np.array_equal(lattice.basis, ((1, 0), (0, -3.5)))
        assert lattice.cell_shape == 'Rectangular'

    _transformed = {
        planegroup.Cm6: ((0.5, -3**0.5/2), (0.5, 3**0.5/2)),
        planegroup.M10: ((1, 0), (-0.5, -3**0.5/2)),
        }

    @parametrize('transform,expect_basis', _transformed.items())
    def test_transformed(self, transform, expect_basis):
        """Check generation of a transformed lattice."""
        lattice = Lattice2D(_HEX_OBTUSE)
        # Notice that here we should only use transformation
        # matrices that make sense for the current shape!
        new_lattice = lattice.transformed(transform)
        assert new_lattice is not lattice
        assert np.array_equal(lattice.basis, _HEX_OBTUSE)
        assert new_lattice.cell_shape == 'Hexagonal'
        assert np.array_equal(new_lattice.basis, expect_basis)

    _bases = (_HEX_OBTUSE, _RECT, _RHOMBIC_OBTUSE, _SQUARE)

    @parametrize(basis=_bases, lin_comb=LINEAR_COMBINATIONS)
    def test_high_sym(self, basis, lin_comb, subtests):
        """Check correct value of Lattice2D.is_high_symmetry()."""
        lattice = Lattice2D(basis)
        # Use only transformations that preserve the area,
        # as the algorithm does not try to reduce it
        transform = np.eye(2)
        transform[0, 1] = lin_comb
        self._check_high_symmetry(lattice, transform, subtests, 'modify 1st')
        transform = np.eye(2)
        transform[1, 0] = lin_comb
        self._check_high_symmetry(lattice, transform, subtests, 'modify 2nd')

    @staticmethod
    def _check_high_symmetry(lattice, transform, subtests, txt):
        """Check high-symmetry for a lattice after transforming it."""
        shape = lattice.cell_shape
        new_lattice = lattice.transformed(transform)
        with subtests.test(f'{txt} - is_high_symmetry'):
            expect = new_lattice.cell_shape != 'Oblique'
            assert new_lattice.is_high_symmetry() == expect
        with subtests.test(f'{txt} - make_high_symmetry'):
            new_lattice.make_high_symmetry()
            assert new_lattice.cell_shape == shape
