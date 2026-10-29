"""Test configuration for viperleed.tests.

Defines fixtures and fixture factories used in multiple tests, and
loads the plugins with the options and hooks of the test suite:
check_slow.py (--check-slow), quick_option.py (--quick), rng_seed.py
(--rng-seed) and slow_marks.py (slow and tensorleed marks, stale
SLOW_CASES guard).

Random-input tests (those requesting `draws`) check N_DRAWS draws
seeded from the test id, so that outcomes are reproducible and known
failures can be listed (xfail). The opt-in option --rng-seed=<int> or
--rng-seed=random derives the seeds from a different base, to explore
other draws locally; it is not meant for CI. With the option, the
xfail marks of random-input tests are removed (rng_seed.py), and
`draws.fresh` is True: tests that expect a known failure (e.g., via
may_fail or inline pytest.xfail) must then just assert, so that the
real outcome is reported. With `random`, the base is picked once by
the controller (and passed on to pytest-xdist workers) and is printed
in the report header and in the terminal summary, to reproduce a run.

Fixtures
--------
check_log_records (factory)
    Raise unless caplog records are exactly as expected.
data_path
    Path to the top-level folder containing test data.
draws
    N_DRAWS seeded random generators for one random-input test.
    Its `fresh` attribute tells whether --rng-seed was given.
mock_import_module (factory)
    Replace importlib.import_module to return selected fake modules.
    Retain normal behavior for all the others.
first_case
    The first of the current pytest-cases.
re_match (factory)
    Return a match object from a pattern and a string.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    'Alexandra Mia Imre (@alexmiame)',
)
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2023-02-28'
__license__ = 'GPLv3+'

import collections
from importlib import import_module
import re
import zlib

import numpy as np
from pytest_cases import fixture

from . import optional_dependencies
from .helpers import TEST_DATA

pytest_plugins = (
    'pytester',  # For testing the options themselves
    'tests.check_slow',
    'tests.quick_option',
    'tests.rng_seed',
    'tests.slow_marks',
    )

# Fixture of optional_dependencies that needs no wrapping
fixture_missing_dependency = optional_dependencies.fixture_missing_dependency

# Number of seeded draws a random-input test checks, all within one
# test (one id). This is the only place where the number is defined.
N_DRAWS = 5

Draw = collections.namedtuple('Draw', ('index', 'seed', 'rng'))


class Draws(tuple):
    """A tuple of Draw objects, knowing whether the base seed is fresh."""

    fresh = False


def pytest_addoption(parser):
    """Add the options of the test suite."""
    optional_dependencies.pytest_addoption(parser)


@fixture
def check_log_records(caplog):
    """Raise unless log records are exactly as expected."""
    def _check(expected_records):
        logged = tuple(r.getMessage() for r in caplog.records)
        assert len(logged) == len(expected_records)
        for log, expect in zip(logged, expected_records):
            if isinstance(expect, str):
                assert log == expect
            else:
                assert expect.fullmatch(log)
    return _check


@fixture(scope='session')
def re_match():  # This is actually a fixture factory
    """Return a re.match object from a pattern and a string."""
    def _match(pattern, string):
        return re.match(pattern, string)
    return _match


@fixture(scope='session')
def data_path():
    """Return the Path to the top-level folder containing test data."""
    return TEST_DATA


@fixture(name='first_case')
def first_case(current_cases):
    """Return the first of the current cases."""
    def _find_case(cases_dict):
        for value in cases_dict.values():
            if isinstance(value, dict):
                try:
                    return _find_case(value)
                except ValueError:
                    pass
            try:
                value.id
            except AttributeError:
                pass
            else:
                return value
        raise ValueError('No case found')
    return _find_case(current_cases)


@fixture
def draws(request):
    """Return N_DRAWS Draw tuples, seeded from the test id."""
    # Seed from the test id (without the file path), so that the
    # random values of a test do not depend on which other tests
    # run, in which order, or in which process: outcomes are then
    # the same in every run. Not hash(): it depends on PYTHONHASHSEED.
    _, test_id = request.node.nodeid.split('::', 1)
    base = request.config.option.rng_seed
    prefix = '' if base is None else f'{base}:'
    seeds = (
        zlib.crc32(f'{prefix}{test_id}#{index}'.encode())
        for index in range(N_DRAWS)
        )
    draws = Draws(
        Draw(index, seed, np.random.default_rng(seed))
        for index, seed in enumerate(seeds)
        )
    draws.fresh = base is not None
    return draws


@fixture
def mock_import_module(mocker):
    """Make importlib.import_module return fake modules."""
    def _patch(**fake_modules):
        def _fake_import(module_name, *args, **kwargs):
            try:
                return fake_modules[module_name]
            except KeyError:
                return import_module(module_name, *args, **kwargs)
        mocker.patch('importlib.import_module', side_effect=_fake_import)
        mocker.patch.dict('sys.modules', fake_modules)
    return _patch
