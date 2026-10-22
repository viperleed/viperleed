"""Test configuration for viperleed.tests.

Defines fixtures and fixture factories used in multiple tests.

Random-input tests (those requesting `draws`) check N_DRAWS draws
seeded from the test id, so that outcomes are reproducible and known
failures can be listed (xfail). The opt-in option --rng-seed=<int> or
--rng-seed=random derives the seeds from a different base, to explore
other draws locally; it is not meant for CI. With the option, the
xfail marks of random-input tests are removed, and `draws.fresh` is
True: tests that expect a known failure (e.g., via may_fail or inline
pytest.xfail) must then just assert, so that the real outcome is
reported. With `random`, the base is picked once by the controller
(and passed on to pytest-xdist workers) and is printed in the report
header and in the terminal summary, to reproduce a run.

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
import random
import re
import zlib

import numpy as np
import pytest
from pytest_cases import fixture

from . import optional_dependencies
from .helpers import TEST_DATA

pytest_plugins = ('pytester',)  # For testing the options themselves

# Fixture of optional_dependencies that needs no wrapping
fixture_missing_dependency = optional_dependencies.fixture_missing_dependency

# Number of seeded draws a random-input test checks, all within one
# test (one id). This is the only place where the number is defined.
N_DRAWS = 5

Draw = collections.namedtuple('Draw', ('index', 'seed', 'rng'))

RNG_SEED_OPTION = '--rng-seed'
WORKER_SEED_KEY = 'rng_seed'  # In pytest-xdist's workerinput

class Draws(tuple):
    """A tuple of Draw objects, knowing whether the base seed is fresh."""

    fresh = False


def pytest_addoption(parser):
    """Add the options of the test suite."""
    optional_dependencies.pytest_addoption(parser)
    parser.addoption(
        RNG_SEED_OPTION,
        dest='rng_seed',
        metavar='INT|random',
        default=None,
        help=('Derive the seeds of random-input tests from this base '
              'integer, or from a random one if "random" (printed at '
              'the end of the run, to reproduce it). Known failures '
              '(xfail) of random-input tests then report their real '
              'outcome. For local exploration, not for CI.'),
        )


def pytest_configure(config):
    """Resolve --rng-seed to an integer, once, in the controller.

    A pytest-xdist worker does not resolve it: it takes the base that
    the controller resolved and handed over (pytest_configure_node).
    """
    workerinput = getattr(config, 'workerinput', None)
    if workerinput is not None:
        config.option.rng_seed = workerinput[WORKER_SEED_KEY]
        return
    value = config.option.rng_seed
    if value is None:
        return  # No option
    if value == 'random':
        value = random.SystemRandom().randrange(2**32)
    try:
        config.option.rng_seed = int(value)
    except ValueError:
        raise pytest.UsageError(
            f'{RNG_SEED_OPTION} must be an integer or "random", '
            f'not {value!r}'
            ) from None


@pytest.hookimpl(optionalhook=True)  # Only exists with pytest-xdist
def pytest_configure_node(node):
    """Hand the resolved base seed over to a pytest-xdist worker."""
    node.workerinput[WORKER_SEED_KEY] = node.config.option.rng_seed


def _rng_seed_message(config):
    """Return a line with the base seed, or None if not used."""
    base = config.option.rng_seed
    if base is None:
        return None
    return f'rng-seed base: {base} (reproduce with --rng-seed={base})'


def pytest_report_header(config):
    """Report the base seed, if used."""
    return _rng_seed_message(config)


def pytest_terminal_summary(terminalreporter, config):
    """Report the base seed, if used, also with -q."""
    message = _rng_seed_message(config)
    if message:
        terminalreporter.write_line(message)


def pytest_collection_modifyitems(items):
    """Do not expect failures of random-input tests with fresh seeds."""
    for item in items:
        if (item.config.option.rng_seed is not None
                and 'draws' in item.fixturenames):
            item.own_markers[:] = [
                mark for mark in item.own_markers if mark.name != 'xfail'
                ]


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
