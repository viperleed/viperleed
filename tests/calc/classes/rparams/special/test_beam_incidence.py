"""Tests for module beam_incidence of calc.classes.rparams.special."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-24'
__license__ = 'GPLv3+'

import copy
from numbers import Real
import pickle

import pytest
from pytest_cases import parametrize

from viperleed.calc.classes.rparams.rparams import Rparams
from viperleed.calc.classes.rparams.special.base import SpecialParameter
from viperleed.calc.classes.rparams.special.beam_incidence import (
    BeamIncidence,
    )


class TestBeamIncidenceValid:
    """Collection of tests for valid inputs to BeamIncidence objects."""

    valid = {
        'default': ((), (0.0, 0.0)),
        'ints': ((30, 50), (30.0, 50.0)),
        'floats': ((19.3, 138.0), (19.3, 138.0)),
        'theta 90': ((90, 10), (90.0, 10.0)),
        'phi 360 kept': ((10, 360), (10.0, 360.0)),
        'phi above 360': ((10, 370), (10.0, 10.0)),
        'phi negative': ((10, -30), (10.0, 330.0)),
        'negative theta': ((-30, 50), (30.0, 230.0)),
        'negative theta, phi wraps': ((-30, 270), (30.0, 90.0)),
        'negative theta, phi 180': ((-30, 180), (30.0, 0.0)),
        'negative theta, phi above 360': ((-30, 400), (30.0, 220.0)),
        'theta -90': ((-90, 0), (90.0, 180.0)),
        }

    @parametrize('args,expect', valid.values(), ids=valid)
    def test_angles(self, args, expect):
        """Check normalised values of .theta and .phi."""
        incidence = BeamIncidence(*args)
        assert (incidence.theta, incidence.phi) == expect
        assert all(isinstance(v, Real)
                   for v in (incidence.theta, incidence.phi))

    _types = {
        'ints': ((30, 50), (int, int)),
        'floats': ((19.3, 138.0), (float, float)),
        'mixed': ((30, 138.0), (int, float)),
        }

    @parametrize('args,expect', _types.values(), ids=_types)
    def test_types_kept(self, args, expect):
        """Check that real numbers are stored with their type."""
        incidence = BeamIncidence(*args)
        assert (type(incidence.theta), type(incidence.phi)) == expect

    _normal = {
        'zero': (0, 0, True),
        'zero, phi': (0, 123, True),
        'at eps': (1e-4, 0, True),
        'above eps': (2e-4, 0, False),
        'negative small': (-1e-5, 0, True),
        'large': (20, 0, False),
        }

    @parametrize('theta,phi,expect', _normal.values(), ids=_normal)
    def test_is_normal(self, theta, phi, expect):
        """Check the .is_normal property."""
        assert BeamIncidence(theta, phi).is_normal is expect

    def test_equality(self):
        """Check equality and hashing of BeamIncidence objects."""
        incidence = BeamIncidence(30, 50)
        assert incidence == BeamIncidence(30.0, 50.0)
        assert incidence == BeamIncidence(-30, 230)
        assert incidence != BeamIncidence(30, 51)
        assert incidence != (30, 50)
        assert hash(incidence) == hash(BeamIncidence(-30, 230))

    def test_repr(self):
        """Check the representation string."""
        expect = 'BeamIncidence(theta=30.0, phi=230.0)'
        assert repr(BeamIncidence(-30.0, 50)) == expect

    @parametrize(value=((12, 34), [12, 34], BeamIncidence(12, 34)))
    def test_from_value(self, value):
        """Check creation from a (theta, phi) pair or an instance."""
        incidence = BeamIncidence.from_value(value)
        assert incidence == BeamIncidence(12, 34)

    def test_rparams_default(self):
        """Check that Rparams has normal incidence by default."""
        rpars = Rparams()
        assert rpars.BEAM_INCIDENCE == BeamIncidence(0, 0)
        assert isinstance(rpars.BEAM_INCIDENCE, BeamIncidence)

    def test_registered(self):
        """Check that BeamIncidence is registered for its parameter."""
        cls = SpecialParameter.get_subclass('BEAM_INCIDENCE')
        assert cls is BeamIncidence

    @parametrize(attr=('theta', 'phi', '_theta', 'other'))
    def test_immutable(self, attr):
        """Check that attributes cannot be set or deleted."""
        incidence = BeamIncidence(10, 20)
        with pytest.raises(AttributeError):
            setattr(incidence, attr, 3)
        with pytest.raises(AttributeError):
            delattr(incidence, attr)
        assert incidence == BeamIncidence(10, 20)

    @parametrize(copier=(copy.copy,
                         copy.deepcopy,
                         lambda x: pickle.loads(pickle.dumps(x))))
    def test_copy(self, copier):
        """Check that copies and pickles are equal to the original."""
        incidence = BeamIncidence(-12.5, 340.25)
        copied = copier(incidence)
        assert copied == incidence
        assert isinstance(copied, BeamIncidence)


class TestBeamIncidenceInvalid:
    """Tests for invalid inputs to BeamIncidence objects."""

    _invalid = {
        'theta too large': ((90.5, 0), ValueError),
        'theta too negative': ((-91, 0), ValueError),
        'theta str': (('10', 0), TypeError),
        'phi str': ((10, 'abc'), TypeError),
        'theta None': ((None, 0), TypeError),
        }

    @parametrize('args,exc', _invalid.values(), ids=_invalid)
    def test_raises(self, args, exc):
        """Check complaints for invalid input."""
        with pytest.raises(exc):
            BeamIncidence(*args)
