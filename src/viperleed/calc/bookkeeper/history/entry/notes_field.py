"""Module notes_field of viperleed.calc.bookkeeper.history.entry.

Defines the NotesField class that handles user notes as well as
discarding.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2024-07-29'
__license__ = 'GPLv3+'

from collections import namedtuple
import copy
import re
from typing import ClassVar
from typing import Tuple

from viperleed.calc.bookkeeper.history.entry.enums import FieldTag
from viperleed.calc.bookkeeper.history.entry.field import MissingField
from viperleed.calc.bookkeeper.history.entry.field import MultiLineField
from viperleed.calc.bookkeeper.history.entry.field import UnknownField
from viperleed.calc.bookkeeper.history.errors import EntrySyntaxError
from viperleed.calc.bookkeeper.history.errors import FixableSyntaxError
from viperleed.calc.lib.dataclass_utils import frozen
from viperleed.calc.lib.dataclass_utils import non_init_field
from viperleed.calc.lib.dataclass_utils import set_frozen_attr


_DISCARDED = 'DISCARDED'    # For entries marked via --discard
_DiscardedInfo = namedtuple('_DiscardedInfo', ('pos', 'line'))
_DISCARDED_RE = re.compile(rf'(\s*{_DISCARDED}\s*)')


# A note that contains _DISCARDED in its string value is considered
# indeed discarded only if _DISCARDED appears on a line by itself
# (excluding white spaces). It does not matter where this line is
# among the various lines that the note field may be composed of.
# If multiple lines are _DISCARDED, the last one is the one that
# "counts". The "normal" way is to have _DISCARDED on the last
# non-empty line, but having it somewhere else does not cause the
# field to be flagged as 'wrongly formatted'. This is to prevent
# situations in which a user may first --discard, then add more
# comments after the _DISCARDED line. If lines containing both
# _DISCARDED and other non-white-space text are found, the field
# is marked as "to be fixed" (but not as DISCARDED).


@frozen
class NotesField(MultiLineField, tag=FieldTag.NOTES, mandatory=True):
    """A multi-line field containing user notes."""

    value: [Tuple[str], str] = MissingField

    # is_discarded is set via __post_init__
    is_discarded: bool = non_init_field(default=False, repr=True)

    # The index in self.value where the _DISCARDED line is placed,
    # and the full raw value of that line (excluding trailing spaces)
    # to ensure user-edited lines are preserved.
    _discarded_info: _DiscardedInfo = non_init_field(
        default=_DiscardedInfo(None, None)
        )
    # Make sure to skip is_discarded when fixing an instance, as it
    # is updated via __post_init__ and we don't want to override it
    _skip_when_fixing: ClassVar[set] = {
        *MultiLineField._skip_when_fixing,
        'is_discarded',
        }

    def __post_init__(self):
        """See if these notes were discarded and mark them accordingly."""
        self._update_discarded(noisy=False)

    def __add__(self, notes):
        """Add more `notes` to this one as new lines."""
        return self._concatenate(self, notes)

    def __radd__(self, notes):
        """Return the reverse addition of `notes` and this one."""
        return self._concatenate(notes, self)

    __iadd__ = __add__

    def __bool__(self):
        """Return whether there are any user notes."""
        return not self._is_empty_like and any(self.value)

    def __str__(self):
        """Return a string version of these notes."""
        if not self.is_discarded:
            return super().__str__()
        # When the notes are discarded, automatically 'fix' them, even
        # if they used to be missing. This is why we don't use super()
        # here, which would return an empty string for missing notes
        # even if they are discarded.
        return self._format_string_value(self._get_string_value())

    def as_discarded(self):
        """Return a version of these notes marked as DISCARDED."""
        if self.is_discarded:
            return self
        discarded = copy.deepcopy(self)
        set_frozen_attr(discarded, 'is_discarded', True)
        return discarded

    @property
    def _is_empty_like(self):
        """Return whether these notes are missing or empty."""
        return self.is_missing or self.is_empty

    def _check_not_empty(self):
        """Don't complain, as notes can be empty."""
        try:
            super()._check_not_empty()
        except EntrySyntaxError:  # Notes can be empty
            pass

    def _check_str_value(self):
        """Mark these notes as discarded as needed."""
        super()._check_str_value()
        # NB: the following line may seem at risk of infinite recursion
        # (since check_value is called in _update_discarded) however,
        # super()._check_str_value converts the string to a tuple. This
        # means that the call to check_value in _update_discarded will
        # call _check_tuple_value rather than _check_str_value again.
        self._update_discarded()

    @classmethod
    def _collect_lines_from_notes(cls, notes):
        """Return a tuple of lines from acceptable `notes`."""
        is_field = isinstance(notes, (NotesField, UnknownField))
        if is_field and notes.is_missing:
            return ()
        if isinstance(notes, NotesField):
            discarded_line = (
                () if not notes.is_discarded
                # pylint: disable-next=protected-access     # Same type
                else (notes._discarded_info.line or _DISCARDED,)
                )
            if notes.is_discarded and notes.is_empty:
                # It's a single line
                return discarded_line
            lines = ('',) if notes.is_empty else notes.value
            lines += discarded_line
            return lines
        if is_field and notes.is_empty:
            return ('',)
        lines_str = notes.value if is_field else notes
        if not lines_str:
            return ('',)
        return cls._convert_str_to_tuple(lines_str)

    @classmethod
    def _concatenate(cls, *operands):
        """Return a new NotesField with concatenated lines from `operands`."""
        collected_lines = (cls._collect_lines_from_notes(item)
                           for item in operands)
        try:
            return cls(sum(collected_lines, tuple()))
        except AttributeError:  # Operands are unacceptable
            return NotImplemented

    def _find_discarded_line(self):
        """Return the index of the last DISCARDED line in self.value."""
        try:  # pylint: disable=too-many-try-statements
            # Here we make sure to realize whether there are issues
            # with the value. This may mean that .check_value was not
            # called appropriately before. We cannot place the "with"
            # statement outside the try...except block, as we would
            # miss the EntrySyntaxError that we want to register.
            with self._register_errors():
                self._check_tuple_value()
        except EntrySyntaxError:
            raise TypeError('Cannot _find_discarded_line when '
                            'value is not a tuple of strings.') from None
        if self._is_empty_like:
            # May be turned to EmptyField by _check_tuple_value
            return None, None

        indices_and_lines = list(enumerate(self.value))
        indices_and_lines.reverse()
        try:
            return next((ind, line.rstrip())
                        for (ind, line) in indices_and_lines
                        if _DISCARDED_RE.fullmatch(line))
        except StopIteration:
            pass
        # See if by chance there is a line that contains _DISCARDED
        # but contains other stuff too. This is considered fixable.
        try:
            discarded_info = next(
                (ind, line.rstrip())
                for (ind, line) in indices_and_lines
                if _DISCARDED_RE.search(line)
                )
        except StopIteration:  # None of them
            return None, None
        # Found one. Fixing corresponds to splitting up the value
        # at the first _DISCARDED, like so:
        # before:  (..., 'DISCARDED line with DISCARDED stuff', ...)
        # after:   (..., 'DISCARDED', 'line with DISCARDED stuff', ...)
        ind, line = discarded_info
        new_lines = (s.rstrip() for s in _DISCARDED_RE.split(line, maxsplit=1))
        new_lines = (s for s in new_lines if s)
        fixed_value = list(self.value)
        fixed_value[ind:ind+1] = new_lines
        reason = (f'Found {_DISCARDED!r} in a line containing '
                  f'additional text (line {ind+1})')
        raise FixableSyntaxError(reason=reason, fixed_value=tuple(fixed_value))

    def _prepare_lines_for_str(self):
        """Return a tuple of lines ready to '\n'-join. None otherwise."""
        # The behavior differs from the one of a MultiLineField only in
        # case the notes are discarded. In that case, we need to place
        # the _DISCARDED line in the right spot and "automatically fix"
        # missing notes. It is important to call the base-class method
        # as it may potentially update self.value to empty if needed.
        # This should happen whether is_discarded or not.
        super_lines = super()._prepare_lines_for_str()
        if not self.is_discarded:
            return super_lines

        lines = [] if self._is_empty_like else list(self.value)

        # Now place _DISCARDED in the right spot
        # pylint: disable-next=unpacking-non-sequence    # Pylint #7437
        ind, line = self._discarded_info
        if ind is None:
            ind, line = len(lines), _DISCARDED
        lines.insert(ind, line)
        return tuple(lines)

    def _update_discarded(self, noisy=True):
        """Update the .is_discarded attribute depending on the value."""
        try:
            self.check_value()
        except EntrySyntaxError:
            return
        except FixableSyntaxError:
            pass

        if self._is_empty_like:
            return

        # At this point self.value should be a tuple of strings.
        # Figure out which line (if any) is _DISCARDED.
        try:
            ind, line = self._find_discarded_line()
        except FixableSyntaxError:
            # _DISCARDED is in a line with other stuff
            if noisy:
                raise
            return
        if ind is None:  # No line in self.value is _DISCARDED
            return

        # There is one _DISCARDED line. Remember this, and clean up
        # self.value to exclude that line. This way, a note that only
        # contains a discarded line is considered to be free of extra
        # user notes.
        set_frozen_attr(self, 'is_discarded', True)
        set_frozen_attr(self, '_discarded_info', _DiscardedInfo(ind, line))

        clean_value = list(self.value)
        clean_value.pop(ind)
        set_frozen_attr(self, 'value', tuple(clean_value))
