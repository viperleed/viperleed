"""Module enums of viperleed.calc.bookkeeper.history.entry.

Defines enumeration classes used in multiple spots.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-08-01'
__license__ = 'GPLv3+'

from enum import Enum
from enum import auto


class FaultyLabel(Enum):
    """An enumeration of labels used to mark faulty entries/fields."""

    # For fields
    EDITED  = 'edited?'
    FIXABLE = 'fixable'
    MISSING = 'missing'
    OK = ''

    # For fix actions of entries
    DUPLICATE = 'duplicate'
    EXTRA = 'extra'
    SORTING = 'scrambled'

    def __str__(self):
        """Return a formatted label for this FaultyLabel."""
        label = self.value
        if label:
            label += _FAULTY_MARK_SEP
        return f'{label:<{_FAULTY_MARK_SPACING}}'

    @classmethod
    def for_action(cls, action, field):
        """Return a label adequate for performing a fix `action` on `field`."""
        _action_to_label = {
            FixAction.FIX_FIELDS: cls.for_field(field),
            FixAction.MERGE_NOTES: cls.DUPLICATE,
            FixAction.MOVE_EXTRAS_TO_NOTES: cls.EXTRA,
            FixAction.NO_ACTION: cls.OK,
            FixAction.REMOVE_DUPLICATES: cls.DUPLICATE,
            FixAction.SORTING: cls.SORTING,
            FixAction.UNFIXABLE: cls.EDITED,
            }
        try:
            return _action_to_label[action]
        except KeyError:
            raise ValueError(f'Unknown action {action}') from None

    @classmethod
    def for_field(cls, field):
        """Return a label adequate for `field`."""
        _label_for_attr = {
            # !! IMPORTANT: sorting here is relevant, as
            # we go through these field attributes in order!
            'is_missing': cls.MISSING,
            'needs_fixing': cls.FIXABLE,
            'was_understood': cls.OK,
            }
        for attr, label in _label_for_attr.items():
            if getattr(field, attr):
                return label
        return cls.EDITED


_FAULTY_MARK_SEP = ' -> '
_FAULTY_MARK_SPACING = len(_FAULTY_MARK_SEP) + max(
    len(f.value) for f in FaultyLabel
    )
_HISTORY_INFO_SPACING = 12  # For FieldTag, i.e., the leftmost bit


class DuplicateType(Enum):
    """Types of possible duplicates. Values are message prefixes."""

    DIFFERENT = 'Problematic lines'
    IDENTICAL = 'Problematic lines with identical contents'
    NONE = ''
    NOTES = 'Duplicate notes'


class FieldTag(Enum):
    """The left side of a history.info field."""

    # !! VERY IMPORTANT !! The values below are sorted as expected
    # in an 'untouched' entry. UNKNOWN is the only exception, as it
    # may happen anywhere.
    TENSOR_NUMS = '# TENSORS'
    JOB_NUMS = '# JOB ID'
    JOB_NAME = '# JOB NAME'
    RUN_INFO = '# RUN'
    TIMESTAMP = '# TIME'
    R_REF = '# R REF'
    R_SUPER = '# R SUPER'
    FOLDER = '# FOLDER'
    NOTES = 'Notes:'

    UNKNOWN = ''  # Unrecognized line

    def __str__(self):
        """Return a string version of this tag."""
        if self is FieldTag.UNKNOWN:
            # No added space at all here
            return self.value
        if self is FieldTag.NOTES:
            # One white space to (potential) text
            return self.value + ' '
        return f'{self.value:<{_HISTORY_INFO_SPACING}}'

    @property
    def stripped(self):
        """Return the bare tag value, without leading '#'."""
        return self.value.replace('# ', '').replace(':', '').strip()


class FixAction(Enum):
    """An action required to fix the contents of an entry."""

    # !! IMPORTANT !!
    # Sorting matters, as this is the order in which fix actions
    # are performed!

    FIX_FIELDS = auto()            # Taken care of by FieldBase
    MERGE_NOTES = auto()           # Multiple notes in one entry
    MOVE_EXTRAS_TO_NOTES = auto()  # Move these into the notes
    NO_ACTION = auto()
    REMOVE_DUPLICATES = auto()     # Identical lines found
    SORTING = auto()               # Fields are in the wrong spots
    UNFIXABLE = auto()             # Something can't be fixed

    @property
    def needs_action(self):
        """Whether this action needs some activity to be fixed."""
        return self not in {FixAction.UNFIXABLE, FixAction.NO_ACTION}

    @classmethod
    def for_duplicates(cls, duplicate_):
        """Return a FixAction for a DuplicateType."""
        _map = {
            DuplicateType.DIFFERENT: cls.UNFIXABLE,
            DuplicateType.IDENTICAL: cls.REMOVE_DUPLICATES,
            DuplicateType.NONE: cls.NO_ACTION,
            DuplicateType.NOTES: cls.MERGE_NOTES,
            }
        try:
            return _map[duplicate_]
        except KeyError:
            raise ValueError(f'{duplicate_} is not a DuplicateType') from None
