"""Tests for module enums of viperleed.calc.bookkeeper.history.entry."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-08-31'
__license__ = 'GPLv3+'


import pytest
from pytest_cases import parametrize

from viperleed.calc.bookkeeper.history.entry.enums import _HISTORY_INFO_SPACING
from viperleed.calc.bookkeeper.history.entry.enums import DuplicateType
from viperleed.calc.bookkeeper.history.entry.enums import FaultyLabel
from viperleed.calc.bookkeeper.history.entry.enums import FieldTag
from viperleed.calc.bookkeeper.history.entry.enums import FixAction



class MockFieldOK:  # pylint: disable=too-few-public-methods
    """A fake field that was understood."""

    is_missing = False
    needs_fixing = False
    was_understood = True


class MockFieldMissing:  # pylint: disable=too-few-public-methods
    """A fake missing field."""

    is_missing = True
    needs_fixing = False
    was_understood = False


class MockFieldToFix:  # pylint: disable=too-few-public-methods
    """A fake field that needs fixing."""

    is_missing = False
    needs_fixing = True
    was_understood = True


class MockFieldNotUnderstood:  # pylint: disable=too-few-public-methods
    """A fake field that was not understood."""

    is_missing = False
    needs_fixing = False
    was_understood = False


class _TestEnumBase:
    """Base class for testing enumerations."""

    enum = None
    values = {}

    def get_member(self, member):
        """Return a member of the enumeration under test."""
        return getattr(self.enum, member)

    def test_values(self):
        """Check values of a enumeration members."""
        if self.values is None:  # We don't want to check this enum
            return
        for member, expect in self.values.items():
            assert self.get_member(member).value == expect

    def test_skipped_no_member(self):
        """Check that we haven't missed any member."""
        if self.values is None:  # We don't want to check this enum
            return
        members = {e.name for e in self.enum}
        assert members == set(self.values)


class TestDuplicateType(_TestEnumBase):
    """Tests for the DuplicateType enumeration."""

    enum = DuplicateType
    values = {
        'DIFFERENT': 'Problematic lines',
        'IDENTICAL': 'Problematic lines with identical contents',
        'NONE': '',
        'NOTES': 'Duplicate notes',
        }


class TestFaultyLabel(_TestEnumBase):
    """Tests for the FaultyLabel enumeration."""

    enum = FaultyLabel
    values = {
        'EDITED': 'edited?',
        'FIXABLE': 'fixable',
        'MISSING': 'missing',
        'OK': '',
        'DUPLICATE': 'duplicate',
        'EXTRA': 'extra',
        'SORTING': 'scrambled',
        }

    _str = {
        'EDITED': '  edited? -> ',
        'OK': ' ' * 13,
        }

    @parametrize('member,expect', _str.items(), ids=_str)
    def test_str(self, member, expect):
        """Test the string version of FaultyLabel."""
        assert str(self.get_member(member)) == expect

    _actions = {
        'no action': (FixAction.NO_ACTION, MockFieldOK, 'OK'),
        'missing': (FixAction.FIX_FIELDS, MockFieldMissing, 'MISSING'),
        }

    @parametrize('action,field_cls,member', _actions.values(), ids=_actions)
    def test_for_action(self, action, field_cls, member):
        """Test the for_action method of FaultyLabel."""
        expect = self.get_member(member)
        assert FaultyLabel.for_action(action, field_cls()) is expect

    def test_for_action_invalid(self):
        """Check complaints when giving an invalid action."""
        with pytest.raises(ValueError):
            FaultyLabel.for_action('INVALID_ACTION', MockFieldOK())

    _fields = {
        MockFieldMissing: 'MISSING',
        MockFieldNotUnderstood: 'EDITED',
        MockFieldOK: 'OK',
        MockFieldToFix: 'FIXABLE',
        }

    @parametrize('field_cls,member', _fields.items(), ids=_fields.values())
    def test_faulty_label_for_field(self, field_cls, member):
        """Test the for_field method of FaultyLabel."""
        expect = self.get_member(member)
        assert FaultyLabel.for_field(field_cls()) is expect


class TestFieldTag(_TestEnumBase):
    """Tests for the FieldTag enumeration."""

    enum = FieldTag
    values = {
        'TENSOR_NUMS': '# TENSORS',
        'JOB_NUMS': '# JOB ID',
        'JOB_NAME': '# JOB NAME',
        'RUN_INFO': '# RUN',
        'TIMESTAMP': '# TIME',
        'R_REF': '# R REF',
        'R_SUPER': '# R SUPER',
        'FOLDER': '# FOLDER',
        'NOTES': 'Notes:',
        'UNKNOWN': '',
        }

    _str = {k: f'{v:<{_HISTORY_INFO_SPACING}}' for k, v in values.items()}
    _str['UNKNOWN'] = ''
    _str['NOTES'] = 'Notes: '
    _stripped = {
        'TENSOR_NUMS': 'TENSORS',
        'UNKNOWN': '',
        'NOTES': 'Notes',
        }

    @parametrize('member,expect', _str.items(), ids=_str)
    def test_str(self, member, expect):
        """Test the string version of FaultyLabel."""
        assert str(self.get_member(member)) == expect

    @parametrize('member,expect', _stripped.items(), ids=_stripped)
    def test_stripped(self, member, expect):
        """Check correct stripping of extra characters."""
        assert self.get_member(member).stripped == expect

    def test_sorting(self):
        """Check the expected sorting of tags."""
        assert tuple(e.name for e in self.enum) == tuple(self.values)


class TestFixAction(_TestEnumBase):
    """Tests for the FixAction enumeration."""

    enum = FixAction
    values = None  # Skip value checks, as we use auto()

    _needs = {
        'NO_ACTION': False,
        'UNFIXABLE': False,
        }
    for action in FixAction:
        _needs.setdefault(action.name, True)

    @parametrize('member,expect', _needs.items(), ids=_needs)
    def test_needs_action(self, member, expect):
        """Check the needs_action property."""
        assert self.get_member(member).needs_action == expect

    _dupes = {
        DuplicateType.DIFFERENT: 'UNFIXABLE',
        DuplicateType.IDENTICAL: 'REMOVE_DUPLICATES',
        DuplicateType.NONE: 'NO_ACTION',
        DuplicateType.NOTES: 'MERGE_NOTES',
        }

    @parametrize('duplicate,expect', _dupes.items(), ids=_dupes.values())
    def test_for_duplicates(self, duplicate, expect):
        """Test the for_duplicates method of FixAction."""
        member = self.get_member(expect)
        assert FixAction.for_duplicates(duplicate) is member

    def test_for_duplicates_invalid(self):
        """Check complaints for an invalid duplicate."""
        with pytest.raises(ValueError):
            FixAction.for_duplicates('INVALID_DUPLICATE_TYPE')

    def test_sorting(self):
        """Check the expected sorting of tags."""
        expect = (
            'FIX_FIELDS',
            'MERGE_NOTES',
            'MOVE_EXTRAS_TO_NOTES',
            'NO_ACTION',
            'REMOVE_DUPLICATES',
            'SORTING',
            'UNFIXABLE',
            )
        assert tuple(e.name for e in self.enum) == expect
