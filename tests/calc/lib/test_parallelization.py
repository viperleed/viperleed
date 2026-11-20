"""Tests for module parallelization of viperleed.calc.lib."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-20'
__license__ = 'GPLv3+'

import logging
import os

from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.calc.lib.parallelization import BLAS_THREAD_VARIABLES
from viperleed.calc.lib.parallelization import blas_thread_env
from viperleed.calc.lib.parallelization import log_blas_thread_env


_DEFAULT = '1'  # BLAS thread limit if not set by the user
_UNRELATED = {'PATH': '/some/path', 'VIPERLEED_SOME_VARIABLE': 'value'}
_USER_VARIABLES = (
    'OPENBLAS_NUM_THREADS',
    'GOTO_NUM_THREADS',
    'MKL_NUM_THREADS',
    'OMP_NUM_THREADS',
    )
_BLANK = ('', '  ')  # Values that count as not set
_LIMITS = {  # Set if no variable is set by the user
    'OPENBLAS_NUM_THREADS': _DEFAULT,
    'MKL_NUM_THREADS': _DEFAULT,
    'OMP_NUM_THREADS': _DEFAULT,
    }
# User environment (without OMP_NUM_THREADS) and variables
# expected to be set automatically
_PER_LIBRARY = {
    'MKL': ({'MKL_NUM_THREADS': '2'},
            {'OPENBLAS_NUM_THREADS': _DEFAULT, 'OMP_NUM_THREADS': _DEFAULT}),
    'OPENBLAS': ({'OPENBLAS_NUM_THREADS': '4'},
                 {'MKL_NUM_THREADS': _DEFAULT, 'OMP_NUM_THREADS': _DEFAULT}),
    'GOTO': ({'GOTO_NUM_THREADS': '4'},
             {'MKL_NUM_THREADS': _DEFAULT, 'OMP_NUM_THREADS': _DEFAULT}),
    'OPENBLAS+GOTO': ({'OPENBLAS_NUM_THREADS': '4', 'GOTO_NUM_THREADS': '2'},
                      {'MKL_NUM_THREADS': _DEFAULT,
                       'OMP_NUM_THREADS': _DEFAULT}),
    'MKL+OPENBLAS': ({'MKL_NUM_THREADS': '2', 'OPENBLAS_NUM_THREADS': '0'},
                     {'OMP_NUM_THREADS': _DEFAULT}),
    'MKL+GOTO': ({'MKL_NUM_THREADS': 'garbage', 'GOTO_NUM_THREADS': '4'},
                 {'OMP_NUM_THREADS': _DEFAULT}),
    }


@fixture(name='patch_environ')
def fixture_patch_environ(mocker):
    """Return a callable that replaces os.environ for one test.

    The callable takes keyword arguments for the variables to set
    in addition to some unrelated ones. No other variable is set.
    """
    def _patch(**extra):
        mocker.patch.dict(os.environ, {**_UNRELATED, **extra}, clear=True)
    return _patch


class TestBlasThreadEnv:
    """Tests for the blas_thread_env function."""

    @parametrize(blank=(None, *_BLANK))
    def test_sets_limit(self, blank, patch_environ):
        """Check that all variables are set to 1 if none is set."""
        if blank is None:
            patch_environ()
        else:  # Blank counts as not set
            patch_environ(**dict.fromkeys(_USER_VARIABLES, blank))
        env = blas_thread_env()
        assert all(env[var] == _DEFAULT for var in BLAS_THREAD_VARIABLES)

    @parametrize(value=('4', '0', 'garbage'))
    @parametrize(other=({}, {'MKL_NUM_THREADS': '2'},
                        {'OPENBLAS_NUM_THREADS': '3'}))
    def test_manual_mode(self, value, other, patch_environ):
        """Check that nothing is set if the user set OMP_NUM_THREADS."""
        user_env = {'OMP_NUM_THREADS': value, **other}
        patch_environ(**user_env)
        env = blas_thread_env()
        assert env == {**_UNRELATED, **user_env}

    @parametrize(blank=_BLANK)
    def test_manual_mode_with_blank(self, blank, patch_environ):
        """Check that blank variables are passed on in manual mode."""
        user_env = {'OMP_NUM_THREADS': '8', 'OPENBLAS_NUM_THREADS': blank}
        patch_environ(**user_env)
        env = blas_thread_env()
        assert env == {**_UNRELATED, **user_env}

    @parametrize('user_env,expect', tuple(_PER_LIBRARY.values()),
                 ids=tuple(_PER_LIBRARY))
    def test_per_library(self, user_env, expect, patch_environ):
        """Check that only libraries not set by the user get a limit."""
        patch_environ(**user_env)
        env = blas_thread_env()
        assert env == {**_UNRELATED, **user_env, **expect}

    @parametrize(user_set=_USER_VARIABLES)
    @parametrize(blank=_BLANK)
    def test_blank_is_unset(self, user_set, blank, patch_environ):
        """Check that a blank variable is treated as not set."""
        patch_environ(**{user_set: blank})
        env = blas_thread_env()
        expect = {**_UNRELATED, user_set: blank, **_LIMITS}
        assert env == expect

    @parametrize(blank=_BLANK)
    def test_blank_with_other_library(self, blank, patch_environ):
        """Check that blank variables get a limit next to a set one."""
        patch_environ(MKL_NUM_THREADS='2', OPENBLAS_NUM_THREADS=blank,
                      OMP_NUM_THREADS=blank)
        env = blas_thread_env()
        assert env == {**_UNRELATED, 'MKL_NUM_THREADS': '2',
                       'OPENBLAS_NUM_THREADS': _DEFAULT,
                       'OMP_NUM_THREADS': _DEFAULT}

    @parametrize(user_env=({}, {'OMP_NUM_THREADS': ''},
                           {'MKL_NUM_THREADS': '2'}))
    def test_environ_unchanged(self, user_env, patch_environ):
        """Check that os.environ is not modified."""
        patch_environ(**user_env)
        before = os.environ.copy()
        env = blas_thread_env()
        assert env is not os.environ
        assert os.environ == before

    def test_passes_unrelated(self, patch_environ):
        """Check that other environment variables are passed through."""
        patch_environ()
        env = blas_thread_env()
        assert all(env[var] == value for var, value in _UNRELATED.items())


class TestLogBlasThreadEnv:
    """Tests for the log_blas_thread_env function."""

    @parametrize(user_env=({}, {'OMP_NUM_THREADS': ' '}))
    def test_limit(self, user_env, patch_environ, caplog):
        """Check the debug message with the limit set automatically."""
        patch_environ(**user_env)
        before = os.environ.copy()
        caplog.set_level(logging.DEBUG)
        log_blas_thread_env('some calculations')
        expect_log = ('BLAS thread limit for some calculations: '
                      'OPENBLAS_NUM_THREADS=1, OMP_NUM_THREADS=1, '
                      'MKL_NUM_THREADS=1')
        assert expect_log in caplog.text
        assert all(r.levelno == logging.DEBUG for r in caplog.records)
        assert os.environ == before

    @parametrize('user_env,user_set,automatic', (
        ({'MKL_NUM_THREADS': '2'}, 'MKL_NUM_THREADS=2',
         'OPENBLAS_NUM_THREADS=1, OMP_NUM_THREADS=1'),
        ({'OPENBLAS_NUM_THREADS': '4'}, 'OPENBLAS_NUM_THREADS=4',
         'OMP_NUM_THREADS=1, MKL_NUM_THREADS=1'),
        ({'GOTO_NUM_THREADS': '4'}, 'GOTO_NUM_THREADS=4',
         'OMP_NUM_THREADS=1, MKL_NUM_THREADS=1'),
        ), ids=('MKL', 'OPENBLAS', 'GOTO'))
    def test_per_library(self, user_env, user_set, automatic,
                         patch_environ, caplog):
        """Check the info message if the user set one library only."""
        patch_environ(**user_env, OMP_NUM_THREADS='')
        before = os.environ.copy()
        caplog.set_level(logging.INFO)
        log_blas_thread_env('some calculations')
        expect_log = ('Number of BLAS threads for some calculations set '
                      f'by the user: {user_set}. Setting the others '
                      f'automatically: {automatic}.')
        assert expect_log in caplog.text
        assert all(r.levelno == logging.INFO for r in caplog.records)
        assert os.environ == before

    def test_manual_mode(self, patch_environ, caplog):
        """Check the info message naming the variables set by a user."""
        patch_environ(OMP_NUM_THREADS='8', GOTO_NUM_THREADS='2',
                      MKL_NUM_THREADS='')
        before = os.environ.copy()
        caplog.set_level(logging.INFO)
        log_blas_thread_env('some calculations')
        expect_log = ('Number of BLAS threads for some calculations set '
                      'by the user: GOTO_NUM_THREADS=2, OMP_NUM_THREADS=8. '
                      'Not setting a limit automatically.')
        assert expect_log in caplog.text
        assert all(r.levelno == logging.INFO for r in caplog.records)
        assert os.environ == before
