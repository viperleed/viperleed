"""Cases for tests/calc/bookkeeper/history."""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-07-19'
__license__ = 'GPLv3+'


import pytest
from pytest_cases import case

from viperleed.calc.bookkeeper.history.constants import HISTORY_INFO_SEPARATOR

from ..tag import BookkeeperTag as Tag
from .entry.cases_entry import MOCK_TIME_GERMAN
from .entry.cases_entry import MOCK_TIME_ISO
from .entry.cases_entry import NOTES_TEST_CONTENT


class CasesHistoryInfo:
    """Collection of test cases for bookkeeper.history tests."""

    @case(tags=(Tag.BOOKKEEPER, Tag.EMPTY, Tag.CANT_FIX))
    def case_no_history_file(self):
        """Return None as a marker for 'no file.'"""
        return None

    @case(tags=(Tag.BOOKKEEPER, Tag.HISTORY, Tag.EMPTY, Tag.CANT_FIX))
    def case_empty_history_info(self):
        """Return the contents of an empty file."""
        return ''

    @case(tags=(Tag.HISTORY, Tag.MULTI_ENTRY, Tag.CANT_FIX))
    @pytest.mark.skip(reason='needs implementation of file fix')
    def case_entries_mixed_time_formats(self):                                  # TODO: 
        """Return the contents of two consecutive valid entries."""
        return f'''\
# TENSORS   1, 2, 3
# JOB ID    1, 2, 3
# TIME      {MOCK_TIME_GERMAN}
# FOLDER    t003.r001_010203-040506
Notes: {NOTES_TEST_CONTENT}
{HISTORY_INFO_SEPARATOR}
# TENSORS   3
# JOB ID    4
# TIME      {MOCK_TIME_ISO}
# FOLDER    t003.r004_010203-040506
Notes: {NOTES_TEST_CONTENT}'''
