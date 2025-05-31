"""Test configuration for viperleed.tests.

Defines fixtures and fixture factories used in multiple tests.

Fixtures
--------
check_log_records (factory)
    Raise unless caplog records are exactly as expected.
data_path
    Path to the top-level folder containing test data.
first_case
    The first of the current pytest-cases.
re_match (factory)
    Return a match object from a pattern and a string.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    'Alexander M. Imre (@amimre)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-02-28'
__license__ = 'GPLv3+'

import re

from pytest_cases import fixture

from .helpers import TEST_DATA


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
