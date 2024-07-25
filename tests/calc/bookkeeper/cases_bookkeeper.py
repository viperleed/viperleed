"""Cases for tests/calc/bookkeeper.

This module results from a refactoring of bookkeeper/conftest.py.
Collects tagged cases to test the functionality of the bookkeeper
utility of viperleed.calc.
"""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-07-19'
__license__ = 'GPLv3+'

from dataclasses import fields
from enum import IntEnum, auto
import logging

from pytest_cases import case
from pytest_cases import parametrize

from viperleed.calc.bookkeeper.history import _DISCARDED
from viperleed.calc.bookkeeper.history import EntrySyntaxError
from viperleed.calc.bookkeeper.history import HISTORY_INFO_SEPARATOR
from viperleed.calc.bookkeeper.history import HistoryInfoEntry
from viperleed.calc.lib.dataclass_utils import is_optional_field


class BookkeeperTag(IntEnum):
    """Tags for cases for the bookkeeper utility."""
    AUTO_FIX = auto()        # Entry can be auto-fixed
    CANT_FIX = auto()        # Entry cannot be auto-fixed
    EMPTY = auto()           # history.info file is empty
    ENTRY = auto()           # Use only for HistoryInfoEntry
    BOOKKEEPER = auto()      # Use only for bookkeeper.bookkeeper tests
    HISTORY = auto()         # Use only for history.info file tests
    MISS_MANDATORY = auto()  # Some important entry fields are missing
    MULTI_ENTRY = auto()     # Contains more than one entry
    NEEDS_NO_FIX = auto()    # Entry fields do not require any fixing
    RAISES = auto()          # Entry causes exceptions when parsed

    @classmethod
    def all_targets(cls):
        """Return all bookkeeper tags."""
        return cls.BOOKKEEPER, cls.HISTORY


Tag = BookkeeperTag


NOTES_TEST_CONTENT = 'This is a test note.\n   Over multiple lines.'
MOCK_TIME_ISO = '2003-02-01 04:03:06'
MOCK_TIME_GERMAN = '03.02.01 04:03:06'


class CasesHistoryInfo:
    """Collection of test cases for bookkeeper.history tests."""

    @case(tags=(Tag.BOOKKEEPER, Tag.EMPTY, Tag.CANT_FIX))
    def case_no_history_file(self):
        """Return None as a marker for 'no file.'"""
        return None

    @case(tags=(*Tag.all_targets(), Tag.EMPTY, Tag.CANT_FIX))
    def case_empty_history_info(self):
        """Return the contents of an empty file."""
        return ''


class CasesHistoryInfoFunny:
    """Test cases for contents of history.info likely edited by users."""

    @case(tags=(Tag.HISTORY, Tag.AUTO_FIX))
    def case_mixed_tensors_separators(self):
        """Return entry contents with mixed separators for TENSORS & JOB ID."""
        return f'''\
# TENSORS   1, 2  29
# JOB ID    33 24, 12
# RUN       1 2 3
# TIME      {MOCK_TIME_ISO}
# FOLDER    t003.r001_010203-040506
Notes:'''

    @case(tags=(Tag.HISTORY, Tag.CANT_FIX))
    def case_r_factor_out_of_range(self):
        """Return entry contents with R REF outside the expected range."""
        return f'''\
# TENSORS   1, 2, 3
# JOB ID    4, 5, 6
# TIME      {MOCK_TIME_ISO}
# R REF     3.1234
# FOLDER    t003.r001_010203-040506
Notes:'''

    @case(tags=(Tag.HISTORY, Tag.MULTI_ENTRY, Tag.CANT_FIX))
    def _todo_case_entries_mixed_time_formats(self):                            # TODO
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

    @case(tags=(Tag.RAISES, Tag.HISTORY, Tag.CANT_FIX))
    @parametrize(with_notes=(True, False))
    def case_multi_entries_no_separator(self, with_notes):
        """Return contents of more than one entry, with Notes removed."""
        one_entry = self.case_r_factor_out_of_range()
        if not with_notes:
            one_entry = one_entry.replace('Notes:', '')
        return '\n'.join(one_entry for _ in range(3)), EntrySyntaxError

    @case(tags=(Tag.RAISES, Tag.ENTRY, Tag.MULTI_ENTRY))
    @parametrize(with_notes=(True, False))
    def case_multi_entries_separator(self, with_notes):
        """Return contents of more than one entry, with Notes removed."""
        one_entry = f'''\
# TENSORS   1, 2, 3
# JOB ID    1, 2, 3
# TIME      {MOCK_TIME_ISO}
# FOLDER    t003.r001_010203-040506
'''
        if with_notes:
            one_entry += f'Notes: {NOTES_TEST_CONTENT}\n'
        multi_entry = HISTORY_INFO_SEPARATOR.join(one_entry for _ in range(4))
        return multi_entry, ValueError


class CasesInfoEntryCommented:
    """Collection of history.info contents with some empty fields."""

    @case(tags=(Tag.HISTORY, Tag.AUTO_FIX))
    def case_comment_folder(self):
        """Return entry contents where FOLDER has been edited."""
        return f'''\
# TENSORS   1, 2, 29
# JOB ID    24, 37, 99
# RUN       1 2 3
# TIME      {MOCK_TIME_ISO}
# FOLDER    t003.r001_010203-040506  # This is a comment
Notes:'''

    @case(tags=(Tag.HISTORY, Tag.CANT_FIX))
    def case_comment_r_factor(self):
        """Return entry contents where R REF was has user comments."""
        return f'''\
# TENSORS   1, 2, 3
# JOB ID    4, 5, 6
# TIME      {MOCK_TIME_ISO}
# R REF     0.1234  # This is quite good!
# FOLDER    t003.r001_010203-040506
Notes:'''

    @case(tags=(Tag.HISTORY, Tag.CANT_FIX))
    def case_comment_tensors(self):
        """Return entry contents with comments next to TENSORS."""
        return f'''\
# TENSORS   1, 2, 29  # Actually, 29 was 92
# JOB ID    2, 4, 12
# RUN       1 2 3
# TIME      {MOCK_TIME_ISO}
# FOLDER    t003.r001_010203-040506
Notes:'''

    @case(tags=(Tag.HISTORY, Tag.CANT_FIX))
    def case_comment_timestamp(self):
        """Return entry contents with comments next to TENSORS."""
        return f'''\
# TENSORS   1, 2, 29
# JOB ID    2, 4, 12
# RUN       1 2 3
# TIME      {MOCK_TIME_ISO}  # Comment
# FOLDER    t003.r001_010203-040506
Notes:'''


class CasesInfoEntryCorrect:
    """Collection of cases for non-user-edited history.info entries."""

    @case(tags=(Tag.HISTORY, Tag.NEEDS_NO_FIX))
    def case_no_notes(self):
        """Return one full history.info entry without notes."""
        return f'''\
# TENSORS   1, 2, 29
# JOB ID    24, 37, 99
# RUN       1 2 3
# TIME      {MOCK_TIME_ISO}
# FOLDER    t003.r001_010203-040506
Notes:'''

    @case(tags=(Tag.HISTORY, Tag.AUTO_FIX))
    def case_german_datetime(self):
        """Return one full history.info entry without notes."""
        return f'''\
# TENSORS   1, 2, 29
# JOB ID    24, 37, 99
# RUN       1 2 3
# TIME      {MOCK_TIME_GERMAN}
# FOLDER    t003.r001_010203-040506
Notes:'''

    @case(tags=(Tag.HISTORY, Tag.NEEDS_NO_FIX))
    def case_jobname_and_notes(self):
        """Return one full entry with notes and a specific job name."""
        return f'''\
# TENSORS   1, 2, 3
# JOB ID    4, 5, 6
# JOB NAME  test_jobname
# TIME      {MOCK_TIME_ISO}
# FOLDER    t003.r001_010203-040506_test_jobname
Notes: {NOTES_TEST_CONTENT}
'''

    @case(tags=(Tag.HISTORY, Tag.NEEDS_NO_FIX))
    def case_discarded(self):
        """Return the contents of an entry marked as DISCARDED."""
        return f'''\
# TENSORS   1, 2, 3
# JOB ID    4, 5, 6
# JOB NAME  test_jobname
# TIME      {MOCK_TIME_ISO}
# FOLDER    t003.r001_010203-040506_test_jobname
Notes: {NOTES_TEST_CONTENT}
{_DISCARDED}
'''

    @case(tags=(Tag.HISTORY, Tag.NEEDS_NO_FIX))
    def case_notes(self):
        """Return a full entry with user notes."""
        return f'''\
# TENSORS   1, 2, 3
# JOB ID    4, 5, 6
# TIME      {MOCK_TIME_ISO}
# FOLDER    t003.r001_010203-040506
Notes: {NOTES_TEST_CONTENT}
'''

    @case(tags=(Tag.HISTORY, Tag.NEEDS_NO_FIX))
    def case_newline_notes(self):
        """Return a full entry with a newline-only note."""
        return f'''\
# TENSORS   1, 2, 3
# JOB ID    4, 5, 6
# TIME      {MOCK_TIME_ISO}
# FOLDER    t003.r001_010203-040506
Notes:
'''

    @case(tags=(Tag.HISTORY, Tag.NEEDS_NO_FIX))
    def case_r_ref(self):
        """Return a full entry with REFCALC R-factor info."""
        return f'''\
# TENSORS   1, 2, 3
# JOB ID    4, 5, 6
# TIME      {MOCK_TIME_ISO}
# R REF     0.1234
# FOLDER    t003.r001_010203-040506
Notes:'''

    @case(tags=(Tag.HISTORY, Tag.NEEDS_NO_FIX))
    def case_r_super(self):
        """Return a full entry with SUPERPOS R-factor info."""
        return f'''\
# TENSORS   1, 2, 3
# JOB ID    4, 5, 6
# TIME      {MOCK_TIME_ISO}
# R SUPER   0.1234
# FOLDER    t003.r001_010203-040506
Notes:'''

    @case(tags=(Tag.HISTORY, Tag.NEEDS_NO_FIX))
    def case_run_info(self):
        """Return a full entry with RUN information."""
        return f'''\
# TENSORS   1, 2, 3
# JOB ID    4, 5, 6
# RUN       1 2 3
# TIME      {MOCK_TIME_ISO}
# FOLDER    t003.r001_010203-040506
Notes:'''

    @case(tags=(*Tag.all_targets(), Tag.NEEDS_NO_FIX))
    def case_tensors_none(self):
        """Return an entry where no TENSORS were used."""
        return f'''\
# TENSORS   None
# JOB ID    1
# RUN       1 2 3
# TIME      {MOCK_TIME_ISO}
# FOLDER    t000.r001_010203-040506
Notes:'''

    # This one cannot be used as a single entry, and would be
    # labelled as "needs_fixing" if used as a single string,
    # but cannot be auto-fixed
    @case(tags=(Tag.HISTORY, Tag.MULTI_ENTRY))
    def case_two_entries(self):
        """Return the contents of two consecutive valid entries."""
        return f'''\
# TENSORS   1, 2, 3
# JOB ID    1, 2, 3
# TIME      {MOCK_TIME_ISO}
# FOLDER    t003.r001_010203-040506
Notes: {NOTES_TEST_CONTENT}
{HISTORY_INFO_SEPARATOR}
# TENSORS   3
# JOB ID    4
# TIME      {MOCK_TIME_ISO}
# FOLDER    t003.r004_010203-040506
Notes: {NOTES_TEST_CONTENT}'''


class CasesInfoEntryEmpty:
    """Collection of history.info contents with some empty fields."""

    @case(tags=(Tag.HISTORY, Tag.CANT_FIX))
    def case_empty_folder(self):
        """Return entry contents where FOLDER has no contents."""
        return f'''\
# TENSORS   1, 2, 29
# JOB ID    24, 37, 99
# RUN       1 2 3
# TIME      {MOCK_TIME_ISO}
# FOLDER
Notes:'''

    @case(tags=(Tag.HISTORY, Tag.CANT_FIX))
    def case_empty_folder_trailing_spaces(self):
        """Return entry contents where FOLDER has no contents."""
        return f'''\
# TENSORS   1, 2, 29
# JOB ID    24, 37, 99
# RUN       1 2 3
# TIME      {MOCK_TIME_ISO}
# FOLDER    {' '*3}
Notes:'''

    @case(tags=(Tag.HISTORY, Tag.CANT_FIX))
    def case_empty_job_nums(self):
        """Return entry contents with an empty JOB ID field."""
        return f'''\
# TENSORS   1, 2, 29
# JOB ID
# RUN       1 2 3
# TIME      {MOCK_TIME_ISO}
# FOLDER    t003.r001_010203-040506
Notes:'''

    @case(tags=(Tag.HISTORY, Tag.CANT_FIX))
    def case_empty_r_factor(self):
        """Return a full entry with REFCALC R-factor info."""
        return f'''\
# TENSORS   1, 2, 3
# JOB ID    4, 5, 6
# TIME      {MOCK_TIME_ISO}
# R REF
# FOLDER    t003.r001_010203-040506
Notes:'''

    @case(tags=(Tag.HISTORY, Tag.CANT_FIX))
    def case_empty_run(self):
        """Return a full entry with RUN information."""
        return f'''\
# TENSORS   1, 2, 3
# JOB ID    4, 5, 6
# RUN
# TIME      {MOCK_TIME_ISO}
# FOLDER    t003.r001_010203-040506
Notes:'''

    @case(tags=(Tag.HISTORY, Tag.CANT_FIX))
    def case_empty_timestamp(self):
        """Return a full entry with RUN information."""
        return '''\
# TENSORS   1, 2, 3
# JOB ID    4, 5, 6
# TIME
# FOLDER    t003.r001_010203-040506
Notes:'''


_missing = tuple(f.name
                 for f in fields(HistoryInfoEntry)
                 if f.init and not is_optional_field(f))
_missing += (
    'notes',  # Optional with empty default, but field is mandatory
    )

@case(tags=(Tag.HISTORY, Tag.MISS_MANDATORY, Tag.CANT_FIX))
@parametrize(field=_missing)
def case_missing_field(field, caplog):
    """Return the contents of an entry with a missing mandatory field."""
    with caplog.at_level(logging.CRITICAL):
        dummy = HistoryInfoEntry.from_string(
            CasesInfoEntryCorrect().case_notes()
            )
    kwargs = {f.name: getattr(dummy, f.name)
              for f in fields(dummy)
              if f.init and f.name != field}
    with caplog.at_level(logging.CRITICAL):
        entry_str = str(HistoryInfoEntry(**kwargs))
    # Reason for disable: We would need an Enum of known field tags
    # pylint: disable-next=magic-value-comparison
    if field == 'notes':
        entry_str = '\n'.join(line for line in entry_str.splitlines()
                              if field not in line.lower())
    return entry_str


class CasesInfoEntryPureComment:
    """Entries containing also pure-comment values."""

    @case(tags=(Tag.HISTORY, Tag.ENTRY))
    def case_pure_comment_entry(self):
        """Return the contents of a single, comment-only entry."""
        return '''
    Some text
        with indentation too

    and empty lines, that is a pure-comment block.
      It cannot be discarded.
    '''

    @case(tags=(Tag.HISTORY, Tag.MULTI_ENTRY))
    def case_multi_pure_comments(self):
        """An entry with more comments."""
        one_entry = self.case_pure_comment_entry()
        return HISTORY_INFO_SEPARATOR.join(one_entry for _ in range(2))

    comment = object()
    entry = object()

    @parametrize(which_ones=(
        (comment, entry),
        (entry, comment),
        (entry, comment, entry),
        (comment, entry, comment)
        ))
    @case(tags=(Tag.HISTORY, Tag.MULTI_ENTRY))
    def case_correct_and_comment(self, which_ones):
        """Return a mix of correct and comment-only entries."""
        comment = self.case_pure_comment_entry()
        entry = CasesInfoEntryCorrect().case_no_notes()
        return HISTORY_INFO_SEPARATOR.join(
            comment if which is self.comment else '\n' + entry
            for which in which_ones
            )


class CasesInfoEntryReplaced:
    """Contents of HistoryInfoEntry with some user-replaced fields."""

    @case(tags=(Tag.HISTORY, Tag.CANT_FIX))
    def case_replaced_folder(self):
        """Return entry contents where FOLDER has been edited."""
        return f'''\
# TENSORS   1, 2, 29
# JOB ID    24, 37, 99
# RUN       1 2 3
# TIME      {MOCK_TIME_ISO}
# FOLDER    this is not really a folder
Notes:'''

    @case(tags=(Tag.HISTORY, Tag.CANT_FIX))
    def case_replaced_run(self):
        """Return entry contents where RUN was user edited."""
        return  f'''\
# TENSORS   1, 2, 3
# JOB ID    4, 5, 6
# RUN       someone has edit this
# TIME      {MOCK_TIME_ISO}
# FOLDER    t003.r001_010203-040506
Notes:'''
