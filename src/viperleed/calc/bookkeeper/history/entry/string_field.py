"""Module string_field of viperleed.calc.bookkeeper.history.entry.

Defines the StringField base class for handling fields whose value
is string-only, and its subclasses FolderField, and JobNameField.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2024-07-29'
__license__ = 'GPLv3+'

import re

from viperleed.calc.lib.dataclass_utils import frozen
from viperleed.calc.lib.dataclass_utils import non_init_field
from viperleed.calc.lib.dataclass_utils import set_frozen_attr

from ..errors import EntrySyntaxError
from ..errors import FixableSyntaxError
from .enums import FieldTag
from .field import NoneIsEmptyField
from .field import MissingField


@frozen
class StringField(NoneIsEmptyField):
    """A string-only field."""

    value: str = MissingField

    @property
    def _value_str(self):
        """Return the string version of this field."""
        # It's always identical to value, unless this is missing
        return None if self.is_missing else self.value

    def _check_str_value(self):
        """Check a string value."""
        super()._check_str_value()
        # Clean up string to avoid extra spaces at the left
        # pylint: disable-next=no-member    # Can't infer it's a string
        set_frozen_attr(self, 'value', self.value.lstrip())


_JOB_NAME_PATTERN = r'[ \w]+'
_FOLDER_NAME_RE = re.compile(
    r'(?P<head>'
    r't\d+'          # Tensor number
    r'\.r\d+'        # Run number, given by bookkeeper
    r'(?:\.\d+)?'    # Optional search number, from workhistory number
    '_(?:moved-)?'   # Optional, if bookkeeper runs twice
    r'\d{6}-\d{6})'  # Date and time, with log-name format
    rf'(?:_(?P<job_name>{_JOB_NAME_PATTERN}))?'
    r'(?P<tail>_moved-\d{6}-\d{6})?'  # Ran twice within one second
    )


@frozen
class FolderField(StringField, tag=FieldTag.FOLDER, mandatory=True):
    """A field containing the name of the subfolder(s) of history."""

    # The job_name is automatically extracted from self.value
    job_name: str = non_init_field(repr=True)

    def check_has_job_name(self, job_name):
        """Raise if self.job_name is inconsistent with `job_name`."""
        if not isinstance(job_name, JobNameField):
            job_name = JobNameField(job_name)
        try:
            job_name.check_value()
        except EntrySyntaxError as exc:
            raise ValueError(f'Invalid job_name {job_name!r}') from exc
        if job_name.is_missing:
            raise ValueError('Cannot check_has_job_name '
                             'with a missing job_name')
        with self._register_errors():
            self._check_job_name_consistent(job_name.value)

    def _check_job_name_consistent(self, job_name):
        """Raise if self.job_name does not match `job_name`."""
        if self.is_missing or self.is_empty or not self.was_understood:
            return
        if self.job_name == job_name:
            return
        match_ = _FOLDER_NAME_RE.fullmatch(self.value)
        head, tail = (match_[g] for g in ('head', 'tail'))
        tail = tail or ''
        suffix = f'_{job_name}{tail}'
        reason = (f'is inconsistent with {FieldTag.JOB_NAME.stripped} '
                  f'value {job_name!r}. Should end with {suffix!r}.')
        raise FixableSyntaxError(reason=reason,
                                 fixed_value=f'{head}{suffix}')

    def _check_str_value(self):
        """Check a string value."""
        super()._check_str_value()
        # About the disable: pylint can't infer that self.value is
        # a string by now, as we don't call this if it is not.
        # pylint: disable-next=no-member
        match_ = _FOLDER_NAME_RE.fullmatch(self.value.strip())
        if not match_:
            raise EntrySyntaxError(
                'Does not have format '
                'tTTT.rRRR_[moved-]yymmdd-HHMMSS[_job_name]'
                )
        set_frozen_attr(self, 'value', match_.group())  # Stripped
        set_frozen_attr(self, 'job_name', match_['job_name'])


@frozen
class JobNameField(StringField, tag=FieldTag.JOB_NAME):
    """A field containing a job name."""

    def _check_str_value(self):
        """Complain if a string value contains non-word characters."""
        super()._check_str_value()
        if not re.fullmatch(_JOB_NAME_PATTERN, self.value):
            raise EntrySyntaxError('Can only contain alphanumeric characters, '
                                   'spaces, and underscores')
