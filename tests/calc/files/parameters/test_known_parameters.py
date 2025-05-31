"""Tests for module known_parameters of viperleed.calc.files.parameters."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-10-15'
__license__ = 'GPLv3+'

import logging

import pytest
from pytest_cases import parametrize

from viperleed.calc.files.parameters import errors
from viperleed.calc.files.parameters.known_parameters import KNOWN_PARAMS
from viperleed.calc.files.parameters.known_parameters import _PARAM_ALIAS
from viperleed.calc.files.parameters.known_parameters import did_you_mean
from viperleed.calc.files.parameters.known_parameters import from_alias
from viperleed.calc.files.parameters.known_parameters import warn_if_deprecated


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


class TestDidYouMean:
    """Collection of tests for function did_you_mean."""

    def test_known_params(self, subtests):
        """Check that the close match of known parameters is itself."""
        for known_param in KNOWN_PARAMS:
            with subtests.test(known_param):
                assert did_you_mean(known_param) == known_param

    def test_aliases(self, subtests):
        """Check that the close match of an alias is the known parameter."""
        for alias, known_param in _PARAM_ALIAS.items():
            with subtests.test(alias):
                assert did_you_mean(alias) == known_param

    valid = {
        'sotp': 'STOP',
        'bulkliek': 'BULK_LIKE_BELOW',
        'bluk_bellow': 'BULK_LIKE_BELOW',
        'cmprssn': 'ZIP_COMPRESSION_LEVEL',
        'plt': 'PLOT_IV',
        'T_deb': 'T_DEBYE',
        'Tensot': 'TENSOR_OUTPUT',
        'level': 'LOG_LEVEL',
        'vers': 'TL_VERSION',
        }
    invalid = ('zip', 'zi', 'suprel', 'super', 'asdf', 'unknown_param')

    @parametrize('typo,known', valid.items(), ids=valid)
    def test_close_typos(self, typo, known):
        """Check correct identification of a parameter from a typo."""
        assert did_you_mean(typo) == known

    @parametrize(typo=invalid)
    def test_far_typos(self, typo):
        """Ensure far typos raise exceptions."""
        with pytest.raises(errors.ParameterNotRecognizedError):
            did_you_mean(typo)

    def test_type_error(self):
        """Check complaints when passing a non-string."""
        with pytest.raises(TypeError):
            did_you_mean(5)


class TestDeprecated:
    """Collection of tests for the reporting of parameter deprecation."""

    not_deprecated = {  # parameter, current version
        'param': ('FORTRAN_COMP', '0.10.0'),
        'will be deprecated': ('ParabolaFit', '0.10.0'),
        'unknown': ('invalid_param', '0.11.0'),
        }

    @parametrize('param,version', not_deprecated.values(), ids=not_deprecated)
    def test_not_deprecated(self, param, version, caplog):
        """Check a non-deprecated parameter."""
        with caplog.at_level(logging.WARNING):
            is_deprecated = warn_if_deprecated(param, version)
        assert not is_deprecated
        assert 'deprecated' not in caplog.text

    deprecated = {  # parameter, current version
        'current': ('PARABOLA_FIT', '0.11.0'),
        'was deprecated': ('PARABOLAfit', '1.0.0'),
        }

    @parametrize('param,version', deprecated.values(), ids=deprecated)
    def test_deprecated(self, param, version, caplog):
        """Check a non-deprecated parameter."""
        with caplog.at_level(logging.WARNING):
            is_deprecated = warn_if_deprecated(param, version)
        assert is_deprecated
        assert all(s in caplog.text for s in ('deprecated', param, version))

    @pytest.mark.skip(reason='No parameter has been revived yet')
    def test_revived(self, caplog):
        """Check that a revived PARAMETER is not deprecarted."""
        self.test_not_deprecated('PARABOLA_FIT', '2.3.0', caplog)
