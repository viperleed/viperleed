"""Tests for module coordinates of viperleed.calc.lib."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-08-25'
__license__ = 'GPLv3+'

from collections import namedtuple

import numpy as np
import pytest
from pytest_cases import parametrize

from viperleed.calc.lib.coordinates import COLLAPSE_EPS as EPS
from viperleed.calc.lib.coordinates import add_edges_and_corners
from viperleed.calc.lib.coordinates import collapse
from viperleed.calc.lib.coordinates import collapse_fractional


class TestAddEdgesAndCorners:
    """Tests for the add_edges_and_corners function."""

    _add_args = namedtuple('_add_args',
                           ('cart', 'frac', 'releps', 'ucell', 'props'))
    _add_expect = namedtuple('_expect', ('cart', 'props'))

    _valid = {
        '2D': (
            _add_args(cart=np.array([[0.1, 0.1], [0.9, 0.9]]),
                      frac=np.array([[0.0, 0.0], [1.0, 1.0]]),
                      releps=[0.1, 0.1],
                      ucell=np.eye(2),
                      props=['A', 'B']),
            _add_expect(cart=np.array([[0.1, 0.1], [0.9, 0.9],
                                      [1.1, 0.1], [0.1, 1.1], [1.1, 1.1],
                                      [-0.1, 0.9], [0.9, -0.1], [-0.1, -0.1]]),
                        props=['A', 'B', 'A', 'A', 'A', 'B', 'B', 'B'])
            ),
        '3D': (
            _add_args(cart=np.array([[0.1, 0.1, 0.1], [0.9, 0.9, 0.9]]),
                      frac=np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]]),
                      releps=[0.1, 0.1, 0.1],
                      ucell=np.eye(3),
                      props=['A', 'B']),
            _add_expect(
                cart=np.array([
                    [0.1, 0.1, 0.1], [0.9, 0.9, 0.9],
                    [1.1, 0.1, 0.1], [0.1, 1.1, 0.1], [0.1, 0.1, 1.1],
                    [1.1, 1.1, 0.1], [1.1, 0.1, 1.1], [0.1, 1.1, 1.1],
                    [1.1, 1.1, 1.1],
                    [-0.1, 0.9, 0.9],  [0.9, -0.1, 0.9],  [0.9, 0.9, -0.1],
                    [-0.1, -0.1, 0.9], [-0.1, 0.9, -0.1], [0.9, -0.1, -0.1],
                    [-0.1, -0.1, -0.1],
                    ]),
                props=['A', 'B', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'B', 'B',
                       'B', 'B', 'B', 'B', 'B'],
                )
            ),
        '2d replicas of 3d, via fractional': (
            _add_args(cart=np.array([[0.1, 0.1, 0.1], [0.9, 0.9, 0.9]]),
                      frac=np.array([[0.0, 0.0], [1.0, 1.0]]),
                      releps=[0.1, 0.1, 0.1],
                      ucell=np.eye(3),
                      props=None),
            _add_expect(
                cart=np.array([
                    [0.1, 0.1, 0.1], [0.9, 0.9, 0.9],
                    [1.1, 0.1, 0.1], [0.1, 1.1, 0.1], [1.1, 1.1, 0.1],
                    [-0.1, 0.9, 0.9],  [0.9, -0.1, 0.9], [-0.1, -0.1, 0.9],
                    ]),
                props=None,
                )
            ),
        '2d replicas of 3d, via releps': (
            _add_args(cart=np.array([[0.1, 0.1, 0.1], [0.9, 0.9, 0.9]]),
                      frac=np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]]),
                      releps=[0.1, 0.1],
                      ucell=np.eye(3),
                      props=None),
            _add_expect(
                cart=np.array([
                    [0.1, 0.1, 0.1], [0.9, 0.9, 0.9],
                    [1.1, 0.1, 0.1], [0.1, 1.1, 0.1], [1.1, 1.1, 0.1],
                    [-0.1, 0.9, 0.9],  [0.9, -0.1, 0.9], [-0.1, -0.1, 0.9],
                    ]),
                props=None,
                )
            ),
        '2d replicas of 3d, via ucell': (
            _add_args(cart=np.array([[0.1, 0.1, 0.1], [0.9, 0.9, 0.9]]),
                      frac=np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]]),
                      releps=[0.1, 0.1, 0.1],
                      ucell=np.eye(3)[:2],
                      props=None),
            _add_expect(
                cart=np.array([
                    [0.1, 0.1, 0.1], [0.9, 0.9, 0.9],
                    [1.1, 0.1, 0.1], [0.1, 1.1, 0.1], [1.1, 1.1, 0.1],
                    [-0.1, 0.9, 0.9],  [0.9, -0.1, 0.9], [-0.1, -0.1, 0.9],
                    ]),
                props=None,
                )
            ),
        'nothing to replicate': (
            _add_args(cart=np.array([[0.1, 0.1, 0.1]]),
                      frac=np.array([[0.5, 0.4, 0.7]]),
                      releps=[0.1, 0.1, 0.1],
                      ucell=np.eye(3),
                      props=None),
            _add_expect(cart=np.array([[0.1, 0.1, 0.1]]), props=None)
            ),
        }

    @parametrize('args,expect', _valid.values(), ids=_valid)
    def test_valid(self, args, expect):
        """Check correct outcome of adding edges/corners."""
        extended_cartesian, extended_props = add_edges_and_corners(*args)
        np.testing.assert_allclose(extended_cartesian, expect.cart)
        assert extended_props == expect.props

    _raises = {
        # NB: shape mismatches on fractional and releps are silently
        # ignored currently to prevent excess checks in the function.
        # No checking is done on the length of props before execution
        # for the same reason.
        'ucell shape mismatch': (
            _add_args(cart=np.array([[0.0, 0.0], [1.0, 1.0]]),
                      frac=np.array([[0.0, 0.0], [1.0, 1.0]]),
                      releps=np.array([1e-5, 1e-5]),
                      ucell=np.array([[1.0, 0.0, 0.0]]),
                      props=None),
            ValueError,
            'operands could not be broadcast',
            ),
        'props length mismatch': (
            _add_args(cart=np.array([[0.0, 0.0], [1.0, 1.0]]),
                      frac=np.array([[0.0, 0.0], [1.0, 1.0]]),
                      releps=np.array([1e-5, 1e-5]),
                      ucell=np.eye(2),
                      props=[]),
            IndexError,
            None,
            ),
        }

    @parametrize('args,exc,exc_txt', _raises.values(), ids=_raises)
    def test_raises(self, args, exc, exc_txt):
        """Check complaints for invalid arguments."""
        with pytest.raises(exc, match=exc_txt):
            add_edges_and_corners(*args)


class TestCollapseCartesian:
    """Tests for the collapse function."""

    _collapse_args = namedtuple('_collapse_args', ('cart', 'ucell'))
    _collapse_expect = namedtuple('_collapse_expect', ('cart', 'frac'))

    _floor = {
        '2d': (
            _collapse_args(
                cart=np.array([[1.5, 1.5], [0.5, 0.5], [0.1, 0.1]]),
                ucell=np.eye(2),
                ),
            _collapse_expect(
                cart=np.array([[0.5, 0.5], [0.5, 0.5], [0.1, 0.1]]),
                frac=np.array([[0.5, 0.5], [0.5, 0.5], [0.1, 0.1]]),
                ),
            ),
        }

    @parametrize('args,expect', _floor.values(), ids=_floor)
    def test_floor(self, args, expect):
        """Check correct collapsing of Cartesian coordinates by flooring."""
        collapsed_cartesians, fractional_coordinates = collapse(*args,
                                                                method='floor')
        np.testing.assert_allclose(collapsed_cartesians, expect.cart)
        np.testing.assert_allclose(fractional_coordinates, expect.frac)

    _round = {
        '2d': (
            _collapse_args(
                cart=np.array([[3.9, 0.9], [0.499, 0.5], [0.501, -0.5]]),
                ucell=np.eye(2),
                ),
            _collapse_expect(
                cart=np.array([[-0.1, -0.1], [0.499, 0.5], [-0.499, -0.5]]),
                frac=np.array([[-0.1, -0.1], [0.499, 0.5], [-0.499, -0.5]]),
                ),
            ),
        }

    @parametrize('args,expect', _round.values(), ids=_round)
    def test_round(self, args, expect):
        """Check correct collapsing of Cartesian coordinates by rounding."""
        collapsed_cartesians, fractional_coordinates = collapse(*args,
                                                                method='round')
        np.testing.assert_allclose(collapsed_cartesians, expect.cart)
        np.testing.assert_allclose(fractional_coordinates, expect.frac)


class TestCollapseFractional:
    """Tests for the collapse_fractional function."""

    _floor = {
        'midway': (np.array([[1.5, 1.5], [0.5, 0.5], [0.1, 0.1]]),
                   np.array([[0.5, 0.5], [0.5, 0.5], [0.1, 0.1]])),
        'just below integer': (np.array([1 - EPS, -EPS, 2 - EPS, -1 - EPS]),
                               -EPS),
        'just above integer': (
            np.array([1 + EPS, -1 + EPS, 2 + EPS, -2 + EPS, EPS]),
            EPS,
            ),
        'exact integer': (np.array([0.0, -1.0, 2.0, -3.0]), 0),
        'close to zero': (np.array([EPS, -EPS, -0.1*EPS, -1e-10*EPS]),
                          np.array([EPS, -EPS, -0.1*EPS, -1e-10*EPS])),
        'empty': (np.array([]).reshape(0, 2), np.array([]).reshape(0, 2)),
        }

    @parametrize('fractional,expect', _floor.values(), ids=_floor)
    def test_floor(self, fractional, expect):
        """Check correct collapsing of fractional coordinates by flooring."""
        collapsed = collapse_fractional(fractional, method='floor')
        np.testing.assert_allclose(collapsed, expect)
        assert collapsed is not fractional

    @parametrize('fractional,expect', _floor.values(), ids=_floor)
    def test_in_place(self, fractional, expect):
        """Check correct in-place collapsing of fractional coordinates."""
        fractional = fractional.copy()  # Not to mess with other tests
        collapsed = collapse_fractional(fractional,
                                        method='floor',
                                        in_place=True)
        np.testing.assert_allclose(fractional, expect)
        assert collapsed is fractional

    _raises = {
        'invalid method': (np.array([[1.5, 1.5], [0.5, 0.5], [0.1, 0.1]]),
                           {'method': 'invalid'},
                           ValueError),
        }

    @parametrize('fractional,kwargs,exc', _raises.values(), ids=_raises)
    def test_raises(self, fractional, kwargs, exc):
        """Check complaints for invalid arguments."""
        with pytest.raises(exc):
            collapse_fractional(fractional, **kwargs)

    _round = {
        'midway': (np.array([[1.8, -2.3], [0.5, 0.5], [0.1, 0.1]]),
                   np.array([[-0.2, -0.3], [0.5, 0.5], [0.1, 0.1]])),
        'just below integer': (np.array([1 - EPS, -EPS, 2 - EPS, -1 - EPS]),
                               -EPS),
        'just above integer': (
            np.array([1 + EPS, -1 + EPS, 2 + EPS, -2 + EPS, EPS]),
            EPS,
            ),
        'exact integer': (np.array([0.0, -1.0, 2.0, -3.0]), 0),
        'close to zero': (np.array([EPS, -EPS, -0.1*EPS, -1e-10*EPS]),
                          np.array([EPS, -EPS, -0.1*EPS, -1e-10*EPS])),
        }

    @parametrize('fractional,expect', _round.values(), ids=_round)
    def test_round(self, fractional, expect):
        """Check correct collapsing of fractional coordinates by flooring."""
        collapsed = collapse_fractional(fractional, method='round')
        np.testing.assert_allclose(collapsed, expect)

    @parametrize(method=('floor', 'round'))
    def test_large_values(self, method):
        """Check collapsing of large close-to-integer values."""
        fractional = np.array([1e8 + 1 + EPS, -1e8 - EPS])
        expect = np.array([EPS, -EPS])
        collapsed = collapse_fractional(fractional, method=method)
        np.testing.assert_allclose(collapsed, expect, atol=0.5*EPS)
