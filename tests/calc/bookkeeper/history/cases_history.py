"""Cases for tests/calc/bookkeeper/history."""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-07-19'
__license__ = 'GPLv3+'


from pytest_cases import case
from pytest_cases import parametrize

from viperleed.calc.bookkeeper.history.constants import HISTORY_INFO_SEPARATOR

from ..tag import BookkeeperTag as Tag
from .entry import cases_entry
from .entry.cases_entry import MOCK_TIME_GERMAN
from .entry.cases_entry import MOCK_TIME_ISO
from .entry.cases_entry import NOTES_TEST_CONTENT

_NO_NOTES = cases_entry.CasesInfoEntryCorrect().case_no_notes()

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

    @case(tags=(Tag.HISTORY, Tag.MULTI_ENTRY, Tag.AUTO_FIX))
    def case_mixed_time_formats(self):
        """Return two valid entries with mixed time-stamp formats."""
        return f'''\
# TENSORS   1, 2, 3
# JOB ID    1, 2, 3
# TIME      {MOCK_TIME_ISO}
# FOLDER    t003.r001_010203-040506
Notes:
{HISTORY_INFO_SEPARATOR}
# TENSORS   3
# JOB ID    4
# TIME      {MOCK_TIME_GERMAN}
# FOLDER    t003.r004_010203-040506
Notes:'''

    @case(tags=(Tag.HISTORY, Tag.MULTI_ENTRY, Tag.AUTO_FIX))
    def case_mixed_time_formats_with_notes(self):
        """Return two valid entries with mixed time-stamp formats."""
        return f'''\
# TENSORS   1, 2, 3
# JOB ID    1, 2, 3
# TIME      {MOCK_TIME_ISO}
# FOLDER    t003.r001_010203-040506
Notes:{NOTES_TEST_CONTENT}
{HISTORY_INFO_SEPARATOR}
# TENSORS   3
# JOB ID    4
# TIME      {MOCK_TIME_GERMAN}
# FOLDER    t003.r004_010203-040506
Notes:{NOTES_TEST_CONTENT}'''


    @case(tags=Tag.NEEDS_NO_FIX)
    def case_multi_entry(self):
        """Return the contents of a few entries."""
        cases = cases_entry.CasesInfoEntryCorrect()
        entries = (case()
                   for attr, case in cases.__dict__.items()
                   if attr.startswith('case'))
        return HISTORY_INFO_SEPARATOR.join(entries)

    @case(tags=(Tag.HISTORY, Tag.MULTI_ENTRY, Tag.AUTO_FIX))
    def case_old_time_formats(self):
        """Return two valid entries with old-style time-stamp formats."""
        return f'''\
# TENSORS   1, 2, 3
# JOB ID    1, 2, 3
# TIME      {MOCK_TIME_GERMAN}
# FOLDER    t003.r001_010203-040506
Notes:
{HISTORY_INFO_SEPARATOR}
# TENSORS   3
# JOB ID    4
# TIME      {MOCK_TIME_GERMAN}
# FOLDER    t003.r004_010203-040506
Notes:'''

    @case(tags=(Tag.HISTORY, Tag.EMPTY))
    @parametrize(char=' \n')
    def case_whitespace_only(self, char):
        """Return the contents of a history.info file with only white space."""
        return char*5

    @case(tags=Tag.HISTORY)
    def case_separator_at_end(self):
        """Return the contents of a history.info file with only white space."""
        return _NO_NOTES + HISTORY_INFO_SEPARATOR

    @case(tags=Tag.HISTORY)
    def case_separator_at_end_edited(self):
        """Return the contents of an history file with an edited separator."""
        return _NO_NOTES + HISTORY_INFO_SEPARATOR.rstrip()

    @case(tags=(Tag.HISTORY, Tag.EMPTY))
    @parametrize(whitespace=('', ' '*5, '\n'*5))
    def case_separator_only(self, whitespace):
        """Return the contents of a history.info file with only white space."""
        return (whitespace + HISTORY_INFO_SEPARATOR) * 3 + whitespace
