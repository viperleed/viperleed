"""Module notes_field of viperleed.calc.bookkeeper.history.entry.

Defines the NotesField class that handles user notes as well as
discarding.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-07-29'
__license__ = 'GPLv3+'

import copy

from viperleed.calc.lib.dataclass_utils import frozen
from viperleed.calc.lib.dataclass_utils import non_init_field
from viperleed.calc.lib.dataclass_utils import set_frozen_attr

from .enums import FieldTag
from .field import MultiLineField
from .field import CommentLessField
from .field import UnknownField


_DISCARDED = 'DISCARDED'    # For entries marked via --discard


# !! IMPORTANT !!
# Notice that we never set _value_str for NotesField. This is on
# purpose, so we don't have to fiddle with updating it depending on
# the status of .is_discarded. This is also the reason why NotesField
# is not a StringField subclass, as those always return self.value for
# _value_str.


@frozen
class NotesField(MultiLineField, CommentLessField,
                 tag=FieldTag.NOTES, mandatory=True):
    """A multi-line field containing user notes."""

    # is_discarded is set via __post_init__
    is_discarded: bool = non_init_field(default=False, repr=True)

    def __post_init__(self):
        """See if these notes were discarded and mark them accordingly."""
        self._update_discarded()

    def __add__(self, notes):
        """Add more notes to this one."""
        lines = (
            notes.value if isinstance(notes, (NotesField, UnknownField))
            else notes
            )
        try:
            lines += '' if lines.endswith('\n') else '\n'
        except AttributeError:
            return NotImplemented
        value = self.value if not self.is_missing else ''
        cls = type(self)
        instance = cls(f'{value}\n{lines}')
        set_frozen_attr(instance, 'is_discarded',
                        instance.is_discarded or self.is_discarded)
        return instance

    __iadd__ = __add__
    __radd__ = __add__

    def __bool__(self):
        """Return whether there are any user notes."""
        return not self.is_missing and bool(self.value)

    def __str__(self):
        """Return a string version of these notes."""
        if not self.is_discarded:
            return super().__str__()
        # When the notes are discarded, automatically
        # 'fix' them, even if they used to be missing
        notes = self._get_string_value().strip()
        discarded = f'{_DISCARDED}' if self.is_discarded else ''
        if discarded and notes:
            notes += '\n'
        return self._format_string_value(notes + discarded)

    def as_discarded(self):
        """Return a version of these notes marked as DISCARDED."""
        if self.is_discarded:
            return self
        discarded = copy.deepcopy(self)
        set_frozen_attr(discarded, 'is_discarded', True)
        return discarded

    def _check_not_empty(self):
        """Don't complain, as notes can be empty."""
        # The parent class would raise an EntrySyntaxError otherwise.

    def _check_str_value(self):
        """Mark these notes as discarded as needed."""
        super()._check_str_value()
        self._update_discarded()

    def _format_string_value(self, value_str):
        """Return a formatted and tagged version of `value_str`."""
        # Notes don't have the extra white space after the colon
        return f'{self.tag.value} {value_str}'

    def _update_discarded(self):
        """Update the .is_discarded attribute depending on the value."""
        if self.is_missing:
            return
        # About the pylint disables below: pylint cannot infer that
        # the previous check guarantees that self.value is a string
        # pylint: disable-next=no-member
        notes = self.value.strip()
        is_discarded = self.is_discarded or notes.endswith(_DISCARDED)
        clean_notes = notes.rstrip(_DISCARDED).rstrip()
        set_frozen_attr(self, 'is_discarded', is_discarded)
        set_frozen_attr(self, 'value', clean_notes)
