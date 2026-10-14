"""Test configuration for viperleed.tests.

Defines fixtures and fixture factories used in multiple tests.

Fixtures
--------
check_log_records (factory)
    Raise unless caplog records are exactly as expected.
data_path
    Path to the top-level folder containing test data.
draws
    N_DRAWS seeded random generators for one random-input test.
mock_import_module (factory)
    Replace importlib.import_module to return selected fake modules.
    Retain normal behavior for all the others.
first_case
    The first of the current pytest-cases.
re_match (factory)
    Return a match object from a pattern and a string.
rng
    A random-number generator seeded from the test id.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    'Alexandra Mia Imre (@alexmiame)',
)
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-02-28'
__license__ = 'GPLv3+'

import collections
from importlib import import_module
import re
import zlib

import numpy as np
from pytest_cases import fixture

from .helpers import TEST_DATA

# Number of seeded draws a random-input test checks, all within one
# test (one id). This is the only place where the number is defined.
N_DRAWS = 5

Draw = collections.namedtuple('Draw', ('index', 'seed', 'rng'))


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


@fixture(name='rng')
def fixture_rng(request):
    """Return a numpy random generator seeded from the test id."""
    # Seed from the test id (without the file path), so that the
    # random values of a test do not depend on which other tests
    # run, in which order, or in which process: outcomes are then
    # the same in every run. Not hash(): it depends on PYTHONHASHSEED.
    _, test_id = request.node.nodeid.split('::', 1)
    return np.random.default_rng(zlib.crc32(test_id.encode()))


@fixture
def draws(request):
    """Return N_DRAWS Draw tuples, seeded from the test id."""
    # Seed from the test id (without the file path), so that the
    # random values of a test do not depend on which other tests
    # run, in which order, or in which process: outcomes are then
    # the same in every run. Not hash(): it depends on PYTHONHASHSEED.
    _, test_id = request.node.nodeid.split('::', 1)
    seeds = (
        zlib.crc32(f'{test_id}#{index}'.encode())
        for index in range(N_DRAWS)
        )
    return tuple(
        Draw(index, seed, np.random.default_rng(seed))
        for index, seed in enumerate(seeds)
        )


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
