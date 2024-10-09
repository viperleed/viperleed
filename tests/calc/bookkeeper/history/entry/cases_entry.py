"""Cases for tests/calc/bookkeeper/history/entry."""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-07-19'
__license__ = 'GPLv3+'

from dataclasses import fields as data_fields
import re

from pytest_cases import case
from pytest_cases import parametrize

from viperleed.calc.bookkeeper.history.constants import HISTORY_INFO_SEPARATOR
from viperleed.calc.bookkeeper.history.entry.entry import HistoryInfoEntry
from viperleed.calc.bookkeeper.history.entry.enums import FieldTag
from viperleed.calc.bookkeeper.history.entry.field import FieldBase
from viperleed.calc.bookkeeper.history.entry.notes_field import _DISCARDED
from viperleed.calc.lib.log_utils import logging_silent

from .....helpers import with_case_tags
from ...tag import BookkeeperTag as Tag


NOTES_TEST_CONTENT = 'This is a test note.\n   Over multiple lines.'
MOCK_TIME_ISO = '2003-02-01 04:03:06'
MOCK_TIME_GERMAN = '03.02.01 04:03:06'


class CasesHistoryInfoFunny:
    """Test cases for contents of history.info likely edited by users."""

    @case(tags=Tag.AUTO_FIX_ENTRY)
    def case_extra_lines(self):
        """Return an entry with intermixed comment lines."""
        return f'''\
# TENSORS   None
This is a comment line without a hash
# This one another with a leading hash
# JOB ID    1
# TIME      {MOCK_TIME_ISO}
-- Here another comment line
# FOLDER    t003.r001_010203-040506
Notes:'''

    @case(tags=(Tag.HISTORY, Tag.AUTO_FIX, Tag.AUTO_FIX_ENTRY))
    def case_field_and_entry_fix(self):
        """Return an entry with issues in both a field and the whole entry."""
        return f'''\
# JOB ID    4 5, 6
# TENSORS   1, 2, 3
# R REF     0.1234
# FOLDER    t003.r001_010203-040506
Notes:
# TIME      {MOCK_TIME_ISO}'''

    @case(tags=Tag.CANT_FIX)
    def case_missing_and_fewer_extra_than_n_fields(self):
        """Return an entry with intermixed comment lines and missing fields."""
        return '''\
This entry has 7 _raw_fields, all 5 missing mandatory, and all 4
optional ones. In total n_fields == 5+4 == 9 > len(_raw_fields)
# JOB NAME  name
# RUN       0
# R REF     0.123
# R SUPER   1.234
-- Here another comment line.'''

    @case(tags=Tag.CANT_FIX)
    def case_missing_and_more_extra_than_n_fields(self):
        """Return an entry with intermixed comment lines and missing fields."""
        return '''\
This entry has 12 _raw_fields, 3 filled mandatory ones, 2 missing mandatory,
and all 4 optional ones. In total n_fields == 3+2+4 == 9 < len(_raw_fields)
This is a comment line without a hash. TENSORS is missing
# This one another with a leading hash
# JOB ID    1
# JOB NAME  name
# RUN       0
# R REF     0.123
# R SUPER   1.234
-- Here another comment line. TIME is also missing
# FOLDER    t003.r001_010203-040506_name
Notes:'''

    @case(tags=Tag.CANT_FIX)
    def case_missing_and_same_extra_as_n_fields(self):
        """Return an entry with intermixed comment lines and missing fields."""
        return '''\
This entry has 9 _raw_fields, all 5 missing mandatory,
and all 4 optional ones. In total n_fields == 5+4 == 9 == len(_raw_fields)
This is a comment line without a hash. TENSORS is missing
# This one another with a leading hash
# JOB NAME  name
# RUN       0
# R REF     0.123
# R SUPER   1.234
-- Here another comment line. TIME is also missing'''

    _multi = {
        'non-empty notes': 'case_notes',
        'empty notes': 'case_no_notes',
        }

    @case(tags=(Tag.HISTORY, Tag.AUTO_FIX_ENTRY))
    @parametrize(with_notes=(True, False))
    @parametrize(one_entry_name=_multi.values(), ids=_multi)
    def case_multi_entries_no_separator(self, one_entry_name, with_notes):
        """Return contents of more than one entry, without a separator."""
        one_entry = getattr(CasesInfoEntryCorrect(), one_entry_name)()
        if not with_notes:
            one_entry = re.sub(r'Notes:\s*(.*\n)*', '', one_entry)
        return '\n'.join(one_entry for _ in range(3))

    @case(tags=(Tag.RAISES, Tag.ENTRY, Tag.MULTI_ENTRY))
    @parametrize(with_notes=(True, False))
    def case_multi_entries_separator(self, with_notes):
        """Return contents of more than one entry, with a separator."""
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

    @case(tags=(Tag.HISTORY, Tag.AUTO_FIX_ENTRY))
    def case_multiple_notes_fields(self):
        """Return entry contents with duplicate identical fields."""
        return f'''\
# TENSORS   1, 2, 3
# JOB ID    4, 5, 6
# TIME      {MOCK_TIME_ISO}
# R REF     0.1234
# FOLDER    t003.r001_010203-040506
Notes: {NOTES_TEST_CONTENT}
Notes: Another note with {NOTES_TEST_CONTENT}'''

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

    @case(tags=(Tag.HISTORY, Tag.AUTO_FIX_ENTRY))
    def case_repeated_fields_identical(self):
        """Return entry contents with duplicate identical fields."""
        return f'''\
# TENSORS   1, 2, 3
# JOB ID    4, 5, 6
# JOB ID    4, 5, 6
# TIME      {MOCK_TIME_ISO}
# R REF     0.1234
# FOLDER    t003.r001_010203-040506
Notes:'''

    @case(tags=(Tag.HISTORY, Tag.CANT_FIX))
    def case_repeated_fields_different(self):
        """Return entry contents with duplicate fields and different values."""
        return f'''\
# TENSORS   1, 2, 3
# JOB ID    4, 5, 6
# JOB ID    14, 15, 16
# TIME      {MOCK_TIME_ISO}
# R REF     0.1234
# FOLDER    t003.r001_010203-040506
Notes:'''

    @case(tags=(Tag.HISTORY, Tag.AUTO_FIX_ENTRY))
    def case_sorting_messed_up(self):
        """Return entry contents with duplicate fields and different values."""
        return f'''\
# JOB ID    4, 5, 6
# TENSORS   1, 2, 3
# R REF     0.1234
# FOLDER    t003.r001_010203-040506
Notes:
# TIME      {MOCK_TIME_ISO}'''

    @case(tags=(Tag.MULTI_ENTRY, Tag.RAISES, Tag.CANT_FIX))
    def case_two_entries(self):
        """Return the contents of two consecutive valid entries."""
        contents = f'''\
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
        return contents, ValueError


@with_case_tags(Tag.AUTO_FIX)
class CasesInfoEntryAutoFixFields:
    """Collection of entries that can be automatically fixed."""

    @case(tags=Tag.OLD)
    def case_german_datetime(self):
        """Return one full history.info entry without notes."""
        return f'''\
# TENSORS   1, 2, 29
# JOB ID    24, 37, 99
# RUN       1 2 3
# TIME      {MOCK_TIME_GERMAN}
# FOLDER    t003.r001_010203-040506
Notes:'''

    @case(tags=Tag.HISTORY)
    def case_mixed_tensors_separators(self):
        """Return entry contents with mixed separators for TENSORS & JOB ID."""
        return f'''\
# TENSORS   1, 2  29
# JOB ID    33 24, 12
# RUN       1 2 3
# TIME      {MOCK_TIME_ISO}
# FOLDER    t003.r001_010203-040506
Notes:'''

    @case(tags=Tag.CANT_FIX)
    def case_missing_and_autofix(self):
        """Return entry contents with a fixable field and missing TIME."""
        return '''\
# TENSORS   1, 2  29
# JOB ID    33 24, 12
# RUN       1 2 3
# FOLDER    t003.r001_010203-040506
Notes:'''


@with_case_tags(Tag.NO_ISSUES, Tag.HISTORY)
class CasesInfoEntryCommented:
    """Collection of history.info contents with some empty fields."""

    def case_comment_folder(self):
        """Return entry contents where FOLDER has been edited."""
        return f'''\
# TENSORS   1, 2, 29
# JOB ID    24, 37, 99
# RUN       1 2 3
# TIME      {MOCK_TIME_ISO}
# FOLDER    t003.r001_010203-040506  # This is a comment
Notes:'''

    def case_comment_r_factor(self):
        """Return entry contents where R REF was has user comments."""
        return f'''\
# TENSORS   1, 2, 3
# JOB ID    4, 5, 6
# TIME      {MOCK_TIME_ISO}
# R REF     0.1234  # This is quite good!
# FOLDER    t003.r001_010203-040506
Notes:'''

    def case_comment_tensors(self):
        """Return entry contents with comments next to TENSORS."""
        return f'''\
# TENSORS   1, 2, 29  # Actually, 29 was 92
# JOB ID    2, 4, 12
# RUN       1 2 3
# TIME      {MOCK_TIME_ISO}
# FOLDER    t003.r001_010203-040506
Notes:'''

    def case_comment_timestamp(self):
        """Return entry contents with comments next to TENSORS."""
        return f'''\
# TENSORS   1, 2, 29
# JOB ID    2, 4, 12
# RUN       1 2 3
# TIME      {MOCK_TIME_ISO}  # Comment
# FOLDER    t003.r001_010203-040506
Notes:'''


@with_case_tags(Tag.HISTORY, Tag.NEEDS_NO_FIX)
class CasesInfoEntryCorrect:
    """Collection of cases for non-user-edited history.info entries."""

    def case_no_notes(self):
        """Return one full history.info entry without notes."""
        return f'''\
# TENSORS   1, 2, 29
# JOB ID    24, 37, 99
# RUN       1 2 3
# TIME      {MOCK_TIME_ISO}
# FOLDER    t003.r001_010203-040506
Notes:'''

    def case_jobname_and_notes(self):
        """Return one full entry with notes and a specific job name."""
        return f'''\
# TENSORS   1, 2, 3
# JOB ID    4, 5, 6
# JOB NAME  test_jobname
# TIME      {MOCK_TIME_ISO}
# FOLDER    t003.r001_010203-040506_test_jobname
Notes: {NOTES_TEST_CONTENT}'''

    def case_discarded(self):
        """Return the contents of an entry marked as DISCARDED."""
        return f'''\
# TENSORS   1, 2, 3
# JOB ID    4, 5, 6
# JOB NAME  test_jobname
# TIME      {MOCK_TIME_ISO}
# FOLDER    t003.r001_010203-040506_test_jobname
Notes: {NOTES_TEST_CONTENT}
{_DISCARDED}'''

    def case_notes(self):
        """Return a full entry with user notes."""
        return f'''\
# TENSORS   1, 2, 3
# JOB ID    4, 5, 6
# TIME      {MOCK_TIME_ISO}
# FOLDER    t003.r001_010203-040506
Notes: {NOTES_TEST_CONTENT}'''

    def case_newline_notes(self):
        """Return a full entry with a newline-only note."""
        return f'''\
# TENSORS   1, 2, 3
# JOB ID    4, 5, 6
# TIME      {MOCK_TIME_ISO}
# FOLDER    t003.r001_010203-040506
Notes:
'''

    def case_r_ref(self):
        """Return a full entry with REFCALC R-factor info."""
        return f'''\
# TENSORS   1, 2, 3
# JOB ID    4, 5, 6
# TIME      {MOCK_TIME_ISO}
# R REF     0.1234
# FOLDER    t003.r001_010203-040506
Notes:'''

    def case_r_super(self):
        """Return a full entry with SUPERPOS R-factor info."""
        return f'''\
# TENSORS   1, 2, 3
# JOB ID    4, 5, 6
# TIME      {MOCK_TIME_ISO}
# R SUPER   0.1234
# FOLDER    t003.r001_010203-040506
Notes:'''

    def case_run_info(self):
        """Return a full entry with RUN information."""
        return f'''\
# TENSORS   1, 2, 3
# JOB ID    4, 5, 6
# RUN       1 2 3
# TIME      {MOCK_TIME_ISO}
# FOLDER    t003.r001_010203-040506
Notes:'''

    @case(tags=(Tag.BOOKKEEPER, Tag.ENTRY))
    def case_tensors_none(self):
        """Return an entry where no TENSORS were used."""
        return f'''\
# TENSORS   None
# JOB ID    1
# RUN       1 2 3
# TIME      {MOCK_TIME_ISO}
# FOLDER    t000.r001_010203-040506
Notes:'''


@with_case_tags(Tag.HISTORY, Tag.CANT_FIX)
class CasesInfoEntryEmpty:
    """Collection of history.info contents with some empty fields."""

    def case_empty_folder(self):
        """Return entry contents where FOLDER has no contents."""
        return f'''\
# TENSORS   1, 2, 29
# JOB ID    24, 37, 99
# RUN       1 2 3
# TIME      {MOCK_TIME_ISO}
# FOLDER
Notes:'''

    @case(tags=Tag.NOT_PRESERVED)  # We can't handle trailing spaces
    def case_empty_folder_trailing_spaces(self):
        """Return entry contents where FOLDER has no contents."""
        return f'''\
# TENSORS   1, 2, 29
# JOB ID    24, 37, 99
# RUN       1 2 3
# TIME      {MOCK_TIME_ISO}
# FOLDER    {' '*3}
Notes:'''

    def case_empty_job_nums(self):
        """Return entry contents with an empty JOB ID field."""
        return f'''\
# TENSORS   1, 2, 29
# JOB ID
# RUN       1 2 3
# TIME      {MOCK_TIME_ISO}
# FOLDER    t003.r001_010203-040506
Notes:'''

    def case_empty_r_factor(self):
        """Return a full entry with REFCALC R-factor info."""
        return f'''\
# TENSORS   1, 2, 3
# JOB ID    4, 5, 6
# TIME      {MOCK_TIME_ISO}
# R REF
# FOLDER    t003.r001_010203-040506
Notes:'''

    def case_empty_run(self):
        """Return a full entry with RUN information."""
        return f'''\
# TENSORS   1, 2, 3
# JOB ID    4, 5, 6
# RUN
# TIME      {MOCK_TIME_ISO}
# FOLDER    t003.r001_010203-040506
Notes:'''

    def case_empty_timestamp(self):
        """Return a full entry with RUN information."""
        return '''\
# TENSORS   1, 2, 3
# JOB ID    4, 5, 6
# TIME
# FOLDER    t003.r001_010203-040506
Notes:'''


_missing = tuple(t
                 for t in FieldTag
                 if FieldBase.for_tag(t).is_mandatory)

@case(tags=(Tag.HISTORY, Tag.MISS_MANDATORY, Tag.CANT_FIX))
@parametrize(field_tag=_missing)
def case_missing_field(field_tag):
    """Return the contents of an entry with a missing mandatory field."""
    with logging_silent():
        dummy = HistoryInfoEntry.from_string(
            CasesInfoEntryCorrect().case_notes()
            )
        kwargs = {f.name: getattr(dummy, f.name)
                  for f in data_fields(dummy)
                  if f.init and f.type.tag is not field_tag}
        entry_str = str(HistoryInfoEntry(**kwargs))
    if field_tag is FieldTag.NOTES:
        entry_str = '\n'.join(line for line in entry_str.splitlines()
                              if field_tag.value not in line)
    return entry_str


@with_case_tags(Tag.HISTORY)
class CasesInfoEntryPureComment:
    """Entries containing also pure-comment values."""

    @case(tags=Tag.ENTRY)
    def case_pure_comment_entry(self):
        """Return the contents of a single, comment-only entry."""
        return '''
    Some text
        with indentation too

    and empty lines, that is a pure-comment block.
      It cannot be discarded.
    '''

    @case(tags=Tag.MULTI_ENTRY)
    def case_multi_pure_comments(self):
        """An entry with more comments."""
        one_entry = self.case_pure_comment_entry()
        return HISTORY_INFO_SEPARATOR.join(one_entry for _ in range(2))

    comment = '_comment_'
    entry = '_entry_'

    @parametrize(which_ones=(
        (comment, entry),
        (entry, comment),
        (entry, comment, entry),
        (comment, entry, comment)
        ))
    @case(tags=Tag.MULTI_ENTRY)
    def case_correct_and_comment(self, which_ones):
        """Return a mix of correct and comment-only entries."""
        comment = self.case_pure_comment_entry()
        entry = CasesInfoEntryCorrect().case_no_notes()
        return HISTORY_INFO_SEPARATOR.join(
            comment if which is self.comment else '\n' + entry
            for which in which_ones
            )


@with_case_tags(Tag.HISTORY, Tag.CANT_FIX)
class CasesInfoEntryReplaced:
    """Contents of HistoryInfoEntry with some user-replaced fields."""

    def case_replaced_folder(self):
        """Return entry contents where FOLDER has been edited."""
        return f'''\
# TENSORS   1, 2, 29
# JOB ID    24, 37, 99
# RUN       1 2 3
# TIME      {MOCK_TIME_ISO}
# FOLDER    this is not really a folder
Notes:'''

    def case_replaced_run(self):
        """Return entry contents where RUN was user edited."""
        return  f'''\
# TENSORS   1, 2, 3
# JOB ID    4, 5, 6
# RUN       someone has edit this
# TIME      {MOCK_TIME_ISO}
# FOLDER    t003.r001_010203-040506
Notes:'''
