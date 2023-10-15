"""Tests for module _known_parameters of viperleed.tleedmlib.files.parameters.

Created on 2023-10-15

@author: Michele Riva (@michele-riva)
"""

import pytest
from pytest_cases import parametrize

from viperleed.tleedmlib.files.parameters._known_parameters import (
    from_alias
    )
from viperleed.tleedmlib.files.parameters import errors


class TestFromAlias:
    """Collection of tests for the from_alias function."""
    _known = {        # (arg_to_from_alias, result)
        'correct_name': ('LOG_LEVEL', 'LOG_LEVEL'),
        'no_underscore': ('symmetryeps', 'SYMMETRY_EPS'),
        'alias': ('ivplot', 'PLOT_IV'),
        'lowercase': ('tl_version', 'TL_VERSION'),
        'mix_case': ('SymmetryFix', 'SYMMETRY_FIX'),
        }

    @parametrize('alias,expected', _known.values(), ids=_known)
    def test_known(self, alias, expected):
        """Check correct return value for known parameters."""
        param = from_alias(alias)
        assert param == expected

    _invalid = {
        'unknown': ('unknown_parameter', errors.ParameterNotRecognizedError),
        'empty': ('', errors.ParameterNotRecognizedError),
        'spaces': ('   ', errors.ParameterNotRecognizedError),
        'non-string': (123, TypeError),
        }

    @parametrize('invalid_arg,exc', _invalid.values(), ids=_invalid)
    def test_unknown(self, invalid_arg, exc):
        """Check complaints for unknown parameters."""
        with pytest.raises(exc):
            from_alias(invalid_arg)

    def test_two_aliases(self):
        """Check that two aliases of the same parameter match."""
        param1 = from_alias('fdoptimize')
        param2 = from_alias('fdoptimization')
        assert param1 == param2 == 'OPTIMIZE'
