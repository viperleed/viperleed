"""Tests for module checker of viperleed.calc.files.parameters."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-10-25'
__license__ = 'GPLv3+'

import pytest
from pytest_cases import fixture, parametrize

from viperleed.calc.classes.rparams import Rparams
from viperleed.calc.files.parameters.checker import ParametersChecker
from viperleed.calc.files.parameters.errors import ParameterConflictError

from ....helpers import not_raises


@fixture(name='rpars_with_attrs', scope='session')
def factory_rpars_with_attrs():
    """Return an Rparams with specific attributes set."""
    def _make(**attribute_defs):
        rpars = Rparams()
        for name, value in attribute_defs.items():
            setattr(rpars, name, value)
        return rpars
    return _make


@fixture(name='checker')
def fixture_checker():
    """Return a fresh ParametersChecker object."""
    return ParametersChecker()


class TestElementConflicts:
    """Tests for element-name conflict checks."""

    def test_no_collision(self, checker, rpars_with_attrs):
        """Check that no exceptions are raised if there's no collision."""
        rpars = rpars_with_attrs(ELEMENT_RENAME={'A': 'B'},
                                 ELEMENT_MIX={'C': 'D'})
        with not_raises(Exception):
            checker.check_parameter_conflicts(rpars)

    def test_collision(self, checker, rpars_with_attrs):
        """Ensure exceptions are raised if element name collisions exist."""
        rpars = rpars_with_attrs(ELEMENT_RENAME={'A': 'B', 'C': 'D'},
                                 ELEMENT_MIX={'C': 'E'})
        with pytest.raises(ParameterConflictError):
            checker.check_parameter_conflicts(rpars)

    missing = {
        'no mix': {'ELEMENT_RENAME': {'A': 'B'}},
        'no rename': {'ELEMENT_MIX': {'C': 'D'}},
        'neither': {},
        }

    @parametrize(rpars_kwargs=missing.values(), ids=missing)
    def test_missing_parameters(self, rpars_kwargs, checker, rpars_with_attrs):
        """Check no complaints when some values are at their default."""
        rpars = rpars_with_attrs(**rpars_kwargs)
        with not_raises(Exception):
            checker.check_parameter_conflicts(rpars)
