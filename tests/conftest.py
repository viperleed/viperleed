"""Test configuration for viperleed.tests.

Defines fixtures and fixture factories used in multiple tests.

Fixtures
--------
data_path
    Path to the top-level folder containing test data.
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
