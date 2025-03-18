"""Module entry of viperleed.calc.bookkeeper.history.entry.

Defines classes that handle the contents of a single 'block'
of the history.info file.
"""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2020-01-30'
__license__ = 'GPLv3+'

from collections import defaultdict
from contextlib import AbstractContextManager
from contextlib import nullcontext
from dataclasses import InitVar
from dataclasses import fields as data_fields
import logging
from typing import ClassVar
from typing import Dict

from viperleed.calc.lib.collections_utils import IdentitySet
from viperleed.calc.lib.dataclass_utils import frozen
from viperleed.calc.lib.dataclass_utils import non_init_field
from viperleed.calc.lib.dataclass_utils import replace_values
from viperleed.calc.lib.dataclass_utils import set_frozen_attr
from viperleed.calc.lib.log_utils import logging_silent

from ...mode import BookkeeperMode as Mode
from ..constants import HISTORY_INFO_NAME
from ..constants import HISTORY_INFO_SEPARATOR
from ..errors import EntrySyntaxError
from ..errors import FieldsScrambledError
from ..errors import FixableSyntaxError
from ..errors import FixFailedError
from .enums import DuplicateType
from .enums import FaultyLabel
from .enums import FixAction
from .field import FieldBase
from .field import UnknownField
from .field_collection import FieldList
from .list_of_int_field import JobIdsField
from .list_of_int_field import RunInfoField
from .list_of_int_field import TensorNumsField
from .notes_field import NotesField
from .rfactor_field import RRefField
from .rfactor_field import RSuperField
from .string_field import FolderField
from .string_field import JobNameField
from .time_field import TimestampField
from .time_field import TimestampFormat


LOGGER = logging.getLogger(__name__)


# Format string for FixableSyntaxError raised because issues
# are found in the whole entry when it is read from_string
_FIX_RAW_MSG = 'Found {faulty} lines in entry'


# !! IMPORTANT !!
# The arguments of a HistoryInfoEntry must be sorted the same way as
# in .enums.FieldTag, i.e., the same order as they should appear in
# an entry in the history.info file.

@frozen
class HistoryInfoEntry:
    """A container for information in a single "block" of history.info.

    All attributes below, when given as a string, accept an optional
    comment on the right side. Having comments marks this entry as
    non-removable. The same is true for any non-fixed or non-understood
    fields. Not all fields are always present. The optional ones will
    exist as attributes, but will be missing. You can test this via
    the .is_missing attribute of each field.

    Attributes
    ----------
    tensor_nums : TensorNumsField
        The progressive identifiers of Tensors used during this run.
        String, or list/tuple of non-negative integers are also
        tolerated. A string containing a comma- or space-separated list
        of non-negative integers is also accepted at object creation.
        A space-separated list will cause the entry to be marked as
        'needs fixing'. A string that cannot be converted to a list of
        integers is accepted and retained, but causes the entry to be
        flagged as 'not understood' (and thus not fixable). None,
        'None', [None], and (None,) are accepted at object creation,
        and interpreted as if no Tensors were used or produced. The
        same result is obtained with CalcSection.INITIALIZATION.value
        (== 0) instead of None.
    job_name : JobNameField
        A name assigned to this run by the user. This field is optional
        and may be missing.
    job_nums : JobIdsField
        The progressive identifiers of the executions that used certain
        Tensors during this run. String, or list/tuple of non-negative
        integers are also tolerated. A string containing a comma- or
        space-separated list of non-negative integers is also accepted
        at object creation. A space-separated list will cause the entry
        to be marked as 'needs fixing'. A string that cannot be
        converted to a list of non-negative integers is accepted, but
        causes the entry to be flagged as 'not understood' (and thus
        not fixable).
    run_info : RunInfoField
        The sequence of segments that were executed. Expected format:
        space-separated list of integers. Any non-conforming string
        marks the field as 'not understood' (and thus not fixable).
        This field is optional and may be missing.
    timestamp : TimestampField
        The time when this run was started. A string representing
        a date-time is also accepted. 'moved-' is discarded from
        the string before parsing. If a valid string is given, but
        does not conform with TimestampFormat.DEFAULT, the value
        is accepted but marked as 'needs fixing'. A string that
        cannot be turned into a date-time is also accepted, but
        causes the entry to be flagged as 'not understood' (and
        thus not fixable).
    r_ref : RRefField
        The R factor resulting from a reference-calculation run.
        Numbers outside the interval [0, 2] are flagged as 'not
        understood'. A string is also accepted at object creation.
        In this case, it is considered acceptable if it can be
        converted to a valid R factor. An invalid string is also
        retained, but will cause the field to be marked as 'not
        understood'. This field is optional and may be missing.
    r_super : RSuperField
        The R factor resulting from a superpos-calculation run.
        Numbers outside the interval [0, 2] are flagged as 'not
        understood'. A string is also accepted at object creation.
        In this case, it is considered acceptable if it can be
        converted to a valid R factor. An invalid string is also
        retained, but will cause the field to be marked as 'not
        understood'. This field is optional and may be missing.
    folder_name : FolderField
        The name of the **main** history folder that corresponds
        to this run (i.e., not one of the workhistory folders).
    notes : NotesField
        The notes that users have added for this run. It excludes the
        DISCARDED tag. May span multiple lines. The default is an empty
        string.
    is_discarded : bool
        Whether the user decided to mark this run as DISCARDED. Default
        is False. This value can be set at creation of an entry via the
        discarded_ keyword.
    """

    tensor_nums: TensorNumsField = TensorNumsField()
    job_nums: JobIdsField = JobIdsField()
    job_name: JobNameField = JobNameField()
    run_info: RunInfoField = RunInfoField()
    timestamp: TimestampField = TimestampField()
    r_ref: RRefField = RRefField()
    r_super: RSuperField = RSuperField()
    folder_name: FolderField = FolderField()
    notes: NotesField = NotesField()
    discarded_: InitVar[bool] = False

    # The following init-only keyword should only be used internally.
    # It will trigger the emission of a logging.WARNING at the end of
    # __post_init__ if the arguments given at initialization cause the
    # entry to be marked as 'to be fixed' or 'not understood'. Notice
    # that logging messages for the individual fields are ALWAYS logged
    # irrespective of the value given for this one.
    _should_log_init_warning: InitVar[bool] = True

    # _field_to_attr is a {type(field): attr} map at class level
    _field_to_attr: ClassVar[Dict[FieldBase, str]] = {}

    # _raw_fields is only set when this instance is read from_string.
    # It contains lines in the string, as identified when reading.
    _raw_fields: FieldList = non_init_field(default_factory=FieldList)

    # _fix_todos is a set of actions to be taken for fixing the
    # overall entry (e.g., duplicate fields, sorting, ...). The
    # FieldBase subclasses take care of their own fixing. We use
    # an IdentitySet for the factory rather than a normal set as
    # sets would use hash AND __eq__ to test for item membership.
    # However, we really need to make sure we collect all FieldBase
    # objects, irrespective of their values.
    _fix_todos: Dict[FixAction, IdentitySet] = non_init_field(
        default_factory=lambda: defaultdict(IdentitySet),
        )

    def __post_init__(self, discarded_, _should_log_init_warning):
        """Determine whether arguments make sense and if they need fixing."""
        for attr, dataclass_field in self._iter_fields():
            self._ensure_field_subclass(attr, dataclass_field.type)
            self._check_field_value(getattr(self, attr))
        if discarded_:
            set_frozen_attr(self, 'notes', self.notes.as_discarded())
        self._check_consistency()
        self._log_if_entry_problematic(_should_log_init_warning)

    def __str__(self):
        """Return a string version of this entry, ready for writing to file."""
        # pylint: disable-next=not-an-iterable       # See pylint #7437
        txt = '\n'.join(str(f) for f in self._raw_fields)
        if not txt:  # Not read from_string
            formatted = (str(field) for field in self)
            txt = '\n'.join(f for f in formatted if f.rstrip())
        return '\n' + txt

    def __iter__(self):
        """Yield field values one at a time."""
        yield from self._iter_values()

    @property
    def can_be_removed(self):
        """Return whether this entry can be removed from history.info."""
        return (self.was_understood
                and not self.needs_fixing
                and not self.has_notes
                and not self.has_comments)

    @property
    def is_discarded(self):
        """Return whether this entry was marked as discarded."""
        return self.notes.is_discarded

    @property
    def is_only_outdated(self):
        """Return whether the only fix reason is a legacy-format issue."""
        if not self.needs_fixing:
            return False
        if self._needs_fix_for_entry:
            return False
        other_reasons = (
            not self.was_understood
            or self.has_notes
            or self.has_comments
            )
        if other_reasons:
            return False
        # Only fields need fixing
        fields_to_fix = {f for f in self if f.needs_fixing}
        assert fields_to_fix

        # The only out-of-date field that we currently have is
        # timestamp, with a non-default time format
        out_of_date_fields = {self.timestamp}
        also_other_fields = any(f in fields_to_fix
                                for f in self
                                if f not in out_of_date_fields)
        if also_other_fields:
            return False
        return True

    @property
    def has_comments(self):
        """Return whether any field was edited by adding comments."""
        return any(f.has_comments for f in self)

    @property
    def has_notes(self):
        """Return whether this entry has user notes."""
        return bool(self.notes)

    @property
    def misses_mandatory_fields(self):
        """Return whether this entry has missing fields."""
        return any(f.is_missing for f in self if f.is_mandatory)

    @property
    def needs_fixing(self):
        """Return whether any of the entry fields are wrongly formatted."""
        return self._needs_fix_for_fields or self._needs_fix_for_entry

    @property
    def time_format(self):
        """Return the format used for self.timestamp."""
        return self.timestamp.time_format or TimestampFormat.DEFAULT

    @property
    def was_understood(self):
        """Return whether all fields in this entry were recognized."""
        fields_understood = (f.was_understood for f in self)
        # pylint: disable-next=unsupported-membership-test        #7437
        entry_understood = FixAction.UNFIXABLE not in self._fix_todos
        return all(fields_understood) and entry_understood

    def as_discarded(self):
        """Return a discarded version of this entry."""
        if self.is_discarded:
            return self
        # replace_values goes via __init__ + __post_init__, and would
        # log the same messages as for this instance. Mute logger.
        new_notes = self.notes.as_discarded()
        return self._with_replaced_fields(silent=True, notes=new_notes)

    def as_fixed(self):
        """Return a version of this entry with fixable fields replaced."""
        if not self.needs_fixing:
            return self
        kwargs = {
            attr: f.as_fixed()
            for (attr, _), f in zip(self._iter_fields(), self)
            if f.was_understood  # Fix only those that were understood
            }
        fixed = self._with_replaced_fields(silent=True, **kwargs)
        fixed.fix_entry_issues()
        return fixed

    def fix_entry_issues(self):
        """Fix problems at the level of the whole entry."""
        if self._needs_fix_for_fields:  # Safeguard for the future
            raise FixFailedError('Failed to fix fields for entry:\n'
                                 + self.format_problematic_fields())
        # pylint: disable-next=no-member             # See pylint #7437
        self._fix_todos.pop(FixAction.FIX_FIELDS, None)
        if not self._needs_fix_for_entry:
            return
        for todo in FixAction:
            if not todo.needs_action:
                continue
            try:
                fixer = getattr(self, f'_do_fix_action_{todo.name.lower()}')
            except AttributeError:  # No fixer
                continue
            fixer()
        # pylint: disable-next=no-member             # See pylint #7437
        self._fix_todos.pop(FixAction.NO_ACTION, None)
        # By now, we should be left only with UNFIXABLE stuff
        if self._needs_fix_for_entry:
            raise FixFailedError('Failed to fix entry:\n'
                                 + self.format_problematic_fields())

    def format_problematic_fields(self):
        """Return a string with problematic fields highlighted."""
        if not self.needs_fixing and self.was_understood:
            return ''
        lines = (
            field.format_faulty(with_label=label)
            for field, label in self._get_labels_for_problematic_fields()
            )
        return '\n'.join(line for line in lines if line)

    @classmethod
    def from_string(cls, entry_str):
        """Return an entry object from its string version.

        Parameters
        ----------
        entry_str : str
            The string version of this entry.

        Returns
        -------
        entry: HistoryInfoEntry or PureCommentEntry
            The entry. A PureCommentEntry is returned when entry_str
            does not contain any of the expected fields but only text.

        Raises
        ------
        TypeError
            If `entry_str` is not a string.
        ValueError
            If `entry_str` contains multiple separated entries.
        """
        if not isinstance(entry_str, str):
            raise TypeError('Invalid entry_str. Expected \'str\', '
                            f'Found {type(entry_str).__name__!r}.')
        if HISTORY_INFO_SEPARATOR.rstrip() in entry_str:
            raise ValueError('Found multiple entries. Split them at '
                             'HISTORY_INFO_SEPARATOR beforehand')
        cleaned_str = entry_str.lstrip()
        lines = cleaned_str.splitlines()
        if cleaned_str.endswith('\n'):
            # Make sure not to loose a potential last newline
            lines.append('')
        fields = (FieldBase.from_string(line) for line in lines)
        raw_fields = FieldList(*cls._join_adjacent_note_lines(fields))
        if not raw_fields or raw_fields.has_only_pure_comments:
            return PureCommentEntry(entry_str)

        # Collect the first occurrence of each field as value
        cls._collect_field_to_attr_map()  # Does stuff once per class
        kwargs = {attr: raw_fields.first_by_type(field_cls)
                  for field_cls, attr in cls._field_to_attr.items()}
        kwargs = {k: v for k, v in kwargs.items() if v is not None}
        # Log only at _check_raw_fields that the "entry has issues"
        kwargs['_should_log_init_warning'] = False
        instance = cls(**kwargs)
        set_frozen_attr(instance, '_raw_fields', raw_fields)
        instance._check_raw_fields()
        return instance

    def with_time_format(self, fmt):
        """Return a version of this entry with a specific timestamp format.

        Parameters
        ----------
        fmt : str or TimestampFormat
            The new date-time format to be applied.

        Returns
        -------
        new_entry : HistoryInfoEntry
            This entry, with the new date-time format applied.

        Raises
        ------
        EntrySyntaxError
            If there is no timestamp for which a format can
            be applied.
        """
        new_timestamp = self.timestamp.with_format(fmt)
        return self._with_replaced_fields(silent=True, timestamp=new_timestamp)

    @property
    def _needs_fix_for_entry(self):
        """Return whether the entry as a whole requires a fix."""
        return any(todo.needs_action
                   # pylint: disable-next=not-an-iterable  # Issue 7437
                   for todo in self._fix_todos
                   if todo is not FixAction.FIX_FIELDS)

    @property
    def _needs_fix_for_fields(self):
        """Return whether any of the fields requires a fix."""
        return any(f.needs_fixing for f in self)

    def _check_consistency(self):
        """Check that fields are consistent with one another."""
        checkers = {  # field: (checker_method_name, args, kwargs, skip)
            self.folder_name: (
                'check_has_job_name',
                (self.job_name,),
                {},
                self.job_name.is_missing,  # Can't check a missing one
                ),
            }
        for field, (checker_name, args, kwargs, skip) in checkers.items():
            if skip:
                continue
            with self._log_syntax_errors(for_field=field):
                checker = getattr(field, checker_name)
                # About the next disable: pylint can't infer that
                # kwargs is indeed a mapping. Probably a bug.
                # pylint: disable-next=not-a-mapping
                checker(*args, **kwargs)

    def _check_field_value(self, field):
        """Check that the value of a `field` is OK. Log otherwise."""
        with self._log_syntax_errors(for_field=field):
            field.check_value()

    def _check_raw_fields(self):
        """Log issues found in the raw lines read from string."""
        # First of all, check the values of the raw fields,
        # skipping the ones already stored in self, as we
        # have checked them already in __post_init__.
        checked_already = set(self)
        for field in self._raw_fields:  # pylint: disable=E1133   #7437
            if field in checked_already:
                continue
            self._check_field_value(field)

        # Then check the overall entry
        checkers = (
            self._check_raw_fields_extra_lines,
            self._check_raw_fields_duplicates,
            self._check_raw_fields_sorting,
            )
        for checker in checkers:
            with self._log_syntax_errors():
                checker()
        self._log_if_entry_problematic()

    def _check_raw_fields_duplicates(self):
        """Raise if there is any duplicate raw field."""
        duplicates = defaultdict(list)
        for _, field in self._iter_fields():
            # pylint: disable-next=no-member         # See pylint #7437
            which, duplicate_fields = self._raw_fields.duplicates(field.type)
            if which is DuplicateType.NONE:
                continue
            duplicates[which].extend(duplicate_fields)
        if not duplicates:
            return
        exc = FixableSyntaxError
        for type_ in DuplicateType:
            faulty = duplicates[type_]
            if not faulty:
                continue
            action = FixAction.for_duplicates(type_)
            if action is FixAction.UNFIXABLE:
                exc = EntrySyntaxError
            # pylint: disable-next=E1136             # See pylint #7437
            self._fix_todos[action].update(faulty)
        raise exc(_FIX_RAW_MSG.format(faulty='duplicate'))

    def _check_raw_fields_extra_lines(self):
        """Raise if there are any extra raw lines."""
        # pylint: disable-next=no-member             # See pylint #7437
        extras = self._raw_fields.select_type(UnknownField)
        if not extras:
            return
        # pylint: disable-next=E1136                 # See pylint #7437
        self._fix_todos[FixAction.MOVE_EXTRAS_TO_NOTES].update(extras)
        raise FixableSyntaxError(_FIX_RAW_MSG.format(faulty='unexpected'))

    def _check_raw_fields_sorting(self):
        """Raise if the sorting of the raw lines is not as expected."""
        try:
            # pylint: disable-next=no-member         # See pylint #7437
            self._raw_fields.check_sorted()
        except FieldsScrambledError as exc:
            reason = _FIX_RAW_MSG.format(faulty='incorrectly sorted')
            reason += f' {exc}'
            action = FixAction.SORTING
            raise FixableSyntaxError(reason, action=action) from None

    @classmethod
    def _collect_field_to_attr_map(cls):
        """Produce the internal class attribute mapping types and names."""
        if cls._field_to_attr:
            return
        field_to_attr = {f.type: f_name
                         for f_name, f in cls._iter_fields()}
        setattr(cls, '_field_to_attr', field_to_attr)

    def _do_fix_action_merge_notes(self):
        """Merge duplicate notes fields into one."""
        # pylint: disable-next=no-member             # See pylint #7437
        to_merge = self._fix_todos.pop(FixAction.MERGE_NOTES, None)
        if not to_merge:
            return
        # Filter out the one we're storing already as a field
        to_merge = [n for n in to_merge if n is not self.notes]
        set_frozen_attr(self, 'notes', sum(to_merge, self.notes))
        # pylint: disable-next=no-member             # See pylint #7437
        self._raw_fields.remove_fields(*to_merge)

    def _do_fix_action_move_extras_to_notes(self):
        """Move lines found in intermediate positions to the notes."""
        # pylint: disable-next=no-member             # See pylint #7437
        to_merge = self._fix_todos.pop(FixAction.MOVE_EXTRAS_TO_NOTES, None)
        if not to_merge:
            return
        set_frozen_attr(self, 'notes', sum(to_merge, self.notes))
        # pylint: disable-next=no-member             # See pylint #7437
        self._raw_fields.remove_fields(*to_merge)

    def _do_fix_action_remove_duplicates(self):
        """Remove identical duplicate lines."""
        # pylint: disable-next=no-member             # See pylint #7437
        duplicates = self._fix_todos.pop(FixAction.REMOVE_DUPLICATES, None)
        if not duplicates:
            return
        # Keep all the duplicates that we're storing as field values
        self_values = IdentitySet(self)
        duplicates = (d for d in duplicates if d not in self_values)
        # pylint: disable-next=no-member             # See pylint #7437
        self._raw_fields.remove_fields(*duplicates)

    def _do_fix_action_sorting(self):
        """Fix sorting issues in the raw fields."""
        # pylint: disable-next=no-member             # See pylint #7437
        to_sort = self._fix_todos.pop(FixAction.SORTING, None)
        if not to_sort:
            return
        raw = self._raw_fields
        # pylint: disable-next=no-member             # See pylint #7437
        raw.sort()
        for field_cls, attr in self._field_to_attr.items():
            # pylint: disable-next=no-member         # See pylint #7437
            fixed_field = raw.first_by_type(field_cls)
            if fixed_field is None:
                continue
            set_frozen_attr(self, attr, fixed_field)

    def _ensure_field_subclass(self, attr, field_cls):
        """Make sure self.<attr> is the right FieldBase subclass."""
        value = getattr(self, attr)
        if isinstance(value, field_cls):
            return
        if isinstance(value, FieldBase):
            raise TypeError(f'Invalid type {type(value).__name__!r} '
                            f'for {attr!r}')
        set_frozen_attr(self, attr, field_cls(value))

    def _log_if_entry_problematic(self, should_warn=True):
        """Emit a warning if `should_warn` and this entry has issues."""
        should_warn &= self.needs_fixing or not self.was_understood
        if should_warn:
            LOGGER.warning(f'{HISTORY_INFO_NAME}: Faulty entry is\n\n%s',
                           self.format_problematic_fields())

    def _get_labels_for_problematic_fields(self):
        """Yield fields and labels to display their problematic values."""
        # pylint: disable-next=no-member     # _raw_fields is FieldList
        fields = self._raw_fields.copy() or FieldList(*self)
        if self._raw_fields:
            # Make sure to include the missing ones. All the
            # others should already be there from _raw_fields
            missing = (f for f in self if f.is_missing)
            for field in missing:
                fields.insert_sorted(field)
        # pylint: disable-next=unsupported-membership-test       # 7437
        is_unsorted = FixAction.SORTING in self._fix_todos
        for field in fields:
            try:
                fix_action = next(
                    action
                    # pylint: disable-next=E1101    #7437
                    for action, action_fields in self._fix_todos.items()
                    if field in action_fields
                    )
            except StopIteration:  # Nothing bad with this one
                label = FaultyLabel.for_field(field)
            else:
                # Notice that, since we use an IdentitySet, we do not
                # have to bother with specially checking for unhashable
                # fields (as we would with a set). IdentitySet handles
                # them gracefully.
                label = FaultyLabel.for_action(fix_action, field)
            if label is FaultyLabel.OK and is_unsorted:
                label = FaultyLabel.SORTING
            yield field, label

    @classmethod
    def _iter_fields(cls):
        """Yield fields that are a subclass of FieldBase."""
        fields = ((f.name, f) for f in data_fields(cls)
                  if f.init and issubclass(f.type, FieldBase))
        yield from fields

    def _iter_values(self):
        """Yield attribute values of self that are a subclass of FieldBase."""
        yield from (getattr(self, f) for f, _ in self._iter_fields())

    @staticmethod
    def _join_adjacent_note_lines(fields):
        """Yield fields with adjacent note lines joined."""
        fields = iter(fields)  # Just to be sure
        for notes_field in fields:
            if not isinstance(notes_field, NotesField):
                yield notes_field
                continue
            for next_field in fields:
                if isinstance(next_field, NotesField):
                    # Yield the previous one, then go on with this one
                    yield notes_field
                    notes_field = next_field
                    continue
                if not isinstance(next_field, UnknownField):
                    yield notes_field  # The notes from before
                    yield next_field   # The one we already pulled
                    break
                notes_field += next_field
            else:
                # All lines after a note belonged to the note
                yield notes_field

    def _log_syntax_errors(self, for_field=None):
        """Return a context for catching and logging entry syntax errors."""
        return SyntaxErrorLogger(self._fix_todos, for_field=for_field)

    def _replace_raw_fields_from(self, other, fields):
        """Replace this instance's _raw_fields with those from `other`.

        Parameters
        ----------
        other : HistoryInfoEntry
            The entry from which _raw_fields should be pulled.
        fields : dict
            Keys are names of attributes of HistoryInfoEntry,
            values replace other.<key> in the new _raw_fields
            of this instance. If other.<key> is missing, the
            value is inserted in _raw_fields in the correct
            order.

        Returns
        -------
        None.
        """
        # Concerning the disable: other must be an instance of the same
        # type, from which self has just been created. We don't really
        # care about accessing the private attribute here. It is also
        # somewhat of a safeguard for other to be an HistoryInfoEntry.
        # pylint: disable-next=protected-access
        new_raw = other._raw_fields.copy()
        if not new_raw:  # Nothing to do
            return
        missing = []  # Will insert these in the right order later
        for field_name, new_value in fields.items():
            old_value = getattr(other, field_name)
            if old_value is new_value:
                continue
            if old_value.is_missing and old_value not in new_raw:
                missing.append(new_value)
                continue
            new_raw.replace(old_value, new_value)
        # Now fill in those that were not there
        for missing_value in missing:
            new_raw.insert_sorted(missing_value)

        # Set the new raw fields, and check them to update _fix_todos
        set_frozen_attr(self, '_raw_fields', new_raw)
        self._check_raw_fields()

    def _with_replaced_fields(self, silent=True, **fields):
        """Return a new instance with different field values."""
        context = logging_silent if silent else nullcontext
        with context():
            kwargs = {
                'skip': {
                    '_raw_fields',  # Needs special care
                    '_fix_todos',   # __post_init__ + _check_raw_fields
                    },
                # Py3.7 complains even if InitVar has a default,
                # so discarded_ must be specified explicitly
                'discarded_': self.is_discarded,
                # The next one may be overridden by passing a
                # specific kwarg in **fields
                '_should_log_init_warning': True,
                **fields,
                }
            new_instance = replace_values(self, **kwargs)
            # About the disable: don't care to access a protected one,
            # as we've just created new_instance as a 'copy' of self.
            # pylint: disable-next=protected-access
            new_instance._replace_raw_fields_from(self, fields)
        return new_instance


# Not great to do it this way, but it guarantees that the ClassVar is
# set up correctly as soon as this module is imported. Making the
# method public does not really make much sense, considering that
# it does stuff only once per class.
# pylint: disable-next=protected-access
HistoryInfoEntry._collect_field_to_attr_map()


@frozen
class PureCommentEntry:
    """An entry with only comments."""

    raw_comment: str = ''

    @property
    def can_be_removed(self):
        """Return that a pure-comment entry can never be removed."""
        return False

    def __str__(self):
        """Return the whole raw content of this comment-only entry."""
        return self.raw_comment


class SyntaxErrorLogger(AbstractContextManager):
    """A context manager that logs information upon entry-syntax issues.

    Whenever a syntax error occurs (and is silenced) it also updates
    a dictionary of fixing TODOs with the appropriate fixing action.
    """

    def __init__(self, todos, for_field=None):
        """Initialize logger.

        Parameters
        ----------
        todos : MutableMapping[FixAction, Set]
            Fixing TODOs that will be updated upon syntax errors.
            Keys are FixAction enumerations. When a FixableSyntaxError
            occurs, values are updated with `for_field`, if hashable
            or if Set can handle non-hashable items, `repr(for_field)`
            otherwise. When an EntrySyntaxError occurs, None is added
            to the value for the FixAction.UNFIXABLE key. Values may
            be any data structure with an .add method.
        for_field : FieldBase, optional
            The field for which syntax errors are caught and logged.
            If not given (or None), errors are assumed to be related
            to the entry as whole. Default is None.

        Returns
        -------
        None.
        """
        self.todos = todos
        self.field = for_field
        self.tag = None if not for_field else self.field.tag.stripped

    def __enter__(self):
        """Begin context."""
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        """Catch entry-syntax exceptions."""
        if exc_type is FixableSyntaxError:
            self._handle_fixable(exc_value)
            return True  # Suppressed
        if exc_type is EntrySyntaxError:
            self._handle_unfixable(exc_value)
            return True
        return super().__exit__(exc_type, exc_value, traceback)

    def _add_fixable_todo(self, action):
        """Add the correct fixing TODO."""
        action = (FixAction.FIX_FIELDS if self.field
                  else action or FixAction.NO_ACTION)
        try:
            self.todos[action].add(self.field)
        except TypeError:  # Something unhashable
            self.todos[action].add(repr(self.field))

    def _get_fixable_log_msg(self, problem):
        """Return a message for the logger for a fixable problem."""
        if not self.field:
            return problem
        msg = f'Found entry with {problem}.'
        if self.tag not in msg:
            msg += f' Field: {self.tag}.'
        return msg + ' '

    def _get_unfixable_log_msg(self):
        """Return a message for the logger for an unfixable problem."""
        if not self.field:
            return 'entry'
        value_msg = ('' if self.field.is_missing or self.field.is_empty
                     else f' with value {self.field.value!r}')
        return f'{self.tag} field{value_msg}'

    def _handle_fixable(self, exc):
        """Update TODOs, then log a message for a fixable problem."""
        self._add_fixable_todo(exc.action)
        LOGGER.warning(
            f'{HISTORY_INFO_NAME}: {self._get_fixable_log_msg(exc.reason)}'
            f'Consider running \'bookkeeper {Mode.FIX.long_flag}\'.'
            )

    def _handle_unfixable(self, exc):
        """Update TODOs, then log a message for an unfixable problem."""
        self.todos[FixAction.UNFIXABLE].add(None)
        LOGGER.error(
            f'{HISTORY_INFO_NAME}: Could not understand '
            f'{self._get_unfixable_log_msg()}. Reason: {exc}'
            )
