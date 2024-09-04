"""Tests for module field of viperleed.calc.bookkeeper.history.entry."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-08-30'
__license__ = 'GPLv3+'

from dataclasses import FrozenInstanceError
import functools
import re
from typing import Union

import pytest
from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.calc.bookkeeper.history.errors import EntrySyntaxError
from viperleed.calc.bookkeeper.history.errors import FixableSyntaxError
from viperleed.calc.bookkeeper.history.errors import HistoryInfoError
from viperleed.calc.bookkeeper.history.entry.enums import FaultyLabel
from viperleed.calc.bookkeeper.history.entry.enums import FieldTag
from viperleed.calc.bookkeeper.history.entry.field import CommonRegex
from viperleed.calc.bookkeeper.history.entry.field import CommentLessField
from viperleed.calc.bookkeeper.history.entry.field import DefaultMessage
from viperleed.calc.bookkeeper.history.entry.field import EmptyField
from viperleed.calc.bookkeeper.history.entry.field import FieldBase
from viperleed.calc.bookkeeper.history.entry.field import FixedFieldValue
from viperleed.calc.bookkeeper.history.entry.field import MissingField
from viperleed.calc.bookkeeper.history.entry.field import MultiLineField
from viperleed.calc.bookkeeper.history.entry.field import UnknownField
from viperleed.calc.lib.dataclass_utils import frozen
from viperleed.calc.lib.dataclass_utils import set_frozen_attr

from .....helpers import not_raises
from .conftest import MockFieldTag


class TestCommentlessField:
    """Tests for (concrete subclasses of) CommentLessField."""

    @fixture(name='commentless')
    def fixture_commentless(self, make_concrete_field, make_and_check_field):
        """Return a concrete subclass of CommentLessField."""
        field_cls = make_concrete_field(CommentLessField,
                                        tag=MockFieldTag.TAG_1)
        return functools.partial(make_and_check_field, field_cls)

    _init = {
        '': {'is_empty': True, 'was_understood': False, 'value': EmptyField},
        'value # comment': {'is_empty': False, 'was_understood': True,
                            'value': 'value # comment'},
        '! only comment': {'is_empty': False, 'was_understood': True,
                           'value': '! only comment'},
        '   leading': {'is_empty': False, 'was_understood': True,
                       'value': '   leading'},
        'trailing   ': {'is_empty': False,  'was_understood': True,
                        'value': 'trailing'},
        '#': {'is_empty': False, 'was_understood': True, 'value': '#'},
        }

    @parametrize('value,expect', _init.items(), ids=_init)
    def test_init(self, value, expect, commentless):
        """Check attributes of a comment-less field with a value."""
        field = commentless(value)
        assert not field.has_comments
        assert not field.needs_fixing
        for attr, expected in expect.items():
            assert getattr(field, attr) == expected

    def test_strip_comments(self, commentless):
        """Check correct (non-)stripping of comments."""
        value_with_comment = 'value with # comment'
        field = commentless(value_with_comment)
        assert field.value == value_with_comment
        # pylint: disable-next=protected-access           # OK in tests
        assert field._get_string_value() == value_with_comment
        assert not field.has_comments


class TestCommonRegex:
    """Tests for the regular-expression patterns in CommonRegex."""

    _match = {
        'comma space': ('123,9, 17, 456', 'COMMA_OR_SPACE_SEPARATED_INTS'),
        'any': ('$@123sncbASD_89FCZVJBPA4#€ \t', 'ANY'),
        'one line': ('$@123sncbASD_89FCZVJBPA4#€ \t', 'ONE_LINE'),
        'multiple lines': ('line\n'*10 + '  ', 'MULTILINE'),
        }

    @parametrize('string,member', _match.values(), ids=_match)
    @parametrize(matcher=(re.match, re.fullmatch))
    def test_match(self, member, string, matcher):
        """Check regex patterns."""
        pattern = getattr(CommonRegex, member).value
        assert matcher(pattern, string)

    def test_match_fails(self):
        """Check failure or regex match."""
        string = '123,, , 456'  # No ints between some commas
        pattern = CommonRegex.COMMA_OR_SPACE_SEPARATED_INTS.value
        assert re.fullmatch(pattern, string) is None


class TestFieldBase:
    """Tests for the FieldBase class."""

    _init = {  # field_args, attrs_to_check
        'missing': ((UnknownField,),
                    {'is_missing': True, 'is_empty': False,
                     'was_understood': True, 'needs_fixing': False,
                     'has_comments': False, '_value_str': None}),
        'empty string': ((UnknownField, ''),
                         {'is_missing': False, 'is_empty': True,
                          'was_understood': False, 'needs_fixing': False,
                          'has_comments': False, '_value_str': None}),
        'empty list': ((UnknownField, []),
                       {'is_missing': False, 'is_empty': False,
                        'was_understood': False, 'needs_fixing': False,
                        'has_comments': False, '_value_str': None}),
        'non-empty': ((UnknownField, 'test_value'),
                      {'is_missing': False, 'is_empty': False,
                       'was_understood': True, 'needs_fixing': False,
                       'has_comments': False, '_value_str': None,
                       'value': 'test_value'}),
        'minimal': ((UnknownField, 'x'),
                    {'is_missing': False, 'is_empty': False,
                     'was_understood': True, 'needs_fixing': False,
                     'has_comments': False, '_value_str': None,
                     'value': 'x'}),
        'field value': ((UnknownField, UnknownField('test_value')),
                        {'is_missing': False, 'is_empty': False,
                         'was_understood': True, 'needs_fixing': False,
                         'has_comments': False, '_value_str': None,
                         'value': 'test_value'}),
        'white space only': ((UnknownField, '    '),
                             {'is_missing': False, 'is_empty': True,
                              'was_understood': False, 'needs_fixing': False,
                              'has_comments': False, '_value_str': None}),
        'rstrip': ((UnknownField, 'test_value   \t     '),
                   {'is_missing': False, 'is_empty': False,
                   'was_understood': True, 'needs_fixing': False,
                   'has_comments': False, '_value_str': None,
                   'value': 'test_value'}),
        'special chars': ((UnknownField, 'special chars: \n\t©∆'),
                          {'value': 'special chars: \n\t©∆'})
        }

    @parametrize('args,expect', _init.values(), ids=_init)
    def test_init(self, args, expect, make_and_check_field):
        """Check initialization of a field."""
        field = make_and_check_field(*args)
        for attr, value in expect.items():
            assert getattr(field, attr) == value

    _invalid_value = {
        'None': ((UnknownField, None), DefaultMessage.NOT_STRING),
        'invalid type hash': ((UnknownField, {}), DefaultMessage.NOT_STRING),
        'invalid type int': ((UnknownField, 123), DefaultMessage.NOT_STRING),
        'invalid type list': ((UnknownField, [1, 2, 3]),
                              DefaultMessage.NOT_STRING),
        'multiple lines': ((UnknownField, 'String\n Over  \n Many\n Lines'),
                           'Unexpected format'),
        }
    for key in ('empty string', 'white space only'):
        _invalid_value[key] = (_init[key][0], DefaultMessage.EMPTY)
    for key in ('empty list',):
        _invalid_value[key] = (_init[key][0], DefaultMessage.NOT_STRING)

    @parametrize('args,reason', _invalid_value.values(), ids=_invalid_value)
    def test_invalid_value(self, args, reason, make_field):
        """Check complaints when an invalid value is given."""
        field = make_field(*args)
        try:
            reason = reason.value
        except AttributeError:  # Likely a string, not DefaultMessage
            pass
        with pytest.raises(EntrySyntaxError, match=reason):
            field.check_value()

    def test_valid_value(self, make_field):
        """Check there are no complaints when checking a valid value."""
        args, *_ = self._init['non-empty']
        field = make_field(*args)
        with not_raises(Exception):
            field.check_value()

    _big = {
        'large numeric string': str(10**100),
        'long string': 'abcdefghijklmnopqrstuvwxzy' * 1000,
        'large number': 10**100,
        }

    @parametrize(value=_big.values(), ids=_big)
    def test_big_value(self, value, make_and_check_field):
        """Check handling of large values."""
        field = make_and_check_field(UnknownField, value)
        assert field.value == value

    _str = {
        'missing': '',
        'empty string': '',
        'non-empty': 'test_value',
        'white space only': ''
        }

    @parametrize('key,expect', _str.items(), ids=_str)
    def test_str(self, key, expect, make_and_check_field):
        """Check string version of a field."""
        args, *_ = self._init[key]
        field = make_and_check_field(*args)
        assert str(field) == expect

    @parametrize(key=_init)
    def test_as_fixed_no_fix_needed(self, key, make_and_check_field):
        """Check "fixing" a field that does not need fixing."""
        args, *_ = self._init[key]
        field = make_and_check_field(*args)
        assert field.as_fixed() is field

    def test_as_fixed_fix_needed(self):
        """Check fixing of a field that needs it."""
        field = UnknownField(value='valid_string')
        fixed_value = FixedFieldValue('fix_reason', 'fixed_value')
        set_frozen_attr(field, '_needs_fix', fixed_value)
        fixed_field = field.as_fixed()
        assert fixed_field is not field
        assert fixed_field != field
        assert fixed_field.value == fixed_value.value
        assert not fixed_field.needs_fixing

    def test_as_fixed_no_fix_value(self):
        """Check complaints when trying to fix something with no fix value."""
        field = UnknownField(value='valid_string')
        fixed_value = FixedFieldValue('fix_reason', None)
        set_frozen_attr(field, '_needs_fix', fixed_value)
        with pytest.raises(HistoryInfoError):
            field.as_fixed()

    def test_for_tag(self):
        """Check returning a field subclass from a valid tag."""
        assert FieldBase.for_tag(FieldTag.UNKNOWN) is UnknownField

    _invalid_tag = (
        'INVALID_TAG',
        MockFieldTag.TAG_1,
        )

    @parametrize(tag=_invalid_tag)
    def test_for_tag_invalid(self, tag):
        """Check complaints when an invalid tag is requested."""
        with pytest.raises(ValueError, match='No subclass to handle'):
            FieldBase.for_tag(tag)

    _faulty = {  # value, expected
        'valid': ('   valid_value', ' '*13 + '   valid_value'),
        'missing': (MissingField, ''),
        'empty': ('', str(FaultyLabel.EDITED)),
        }

    @parametrize('value,expect', _faulty.values(), ids=_faulty)
    def test_format_faulty_non_mandatory(self, value, expect,
                                         make_and_check_field):
        """Check formatting of faulty field."""
        field = make_and_check_field(UnknownField, value=value)
        assert field.format_faulty() == expect

    _faulty_mandatory = {
        'missing': ((UnknownField,), str(FaultyLabel.MISSING) + '???'),
        'empty': ((UnknownField, ''), str(FaultyLabel.EDITED)),
        'wrong': ((UnknownField, []), str(FaultyLabel.EDITED)),
        }

    @parametrize('args,expect',
                 _faulty_mandatory.values(),
                 ids=_faulty_mandatory)
    def test_format_faulty_mandatory(self, args, expect, monkeypatch,
                                     make_and_check_field):
        """Check formatting of a faulty, mandatory field."""
        field_cls, *_ = args
        monkeypatch.setattr(field_cls, 'is_mandatory', True)
        field = make_and_check_field(*args)
        assert field.format_faulty() == expect

    def test_immutable(self, make_field):
        """Check that instances can't be modified."""
        field = make_field(UnknownField)
        with pytest.raises(FrozenInstanceError):
            field.value = 3

    def test_invalid_type_hint(self, make_concrete_field,
                               make_and_check_field):
        """Check complaints when an invalid type is given."""
        @frozen
        class _Dummy(FieldBase):
            value: str = MissingField
        field_cls = make_concrete_field(_Dummy, MockFieldTag.TAG_1)
        field = make_and_check_field(field_cls, 123)
        assert not field.was_understood

    def test_register_errors_fixable(self, make_and_check_field):
        """Check correct registering of fixable errors."""
        field = make_and_check_field(UnknownField)
        reason, fixed_value = 'Fixable', 'fixed'
        # pylint: disable-next=protected-access           # OK in tests
        with pytest.raises(FixableSyntaxError), field._register_errors():
            raise FixableSyntaxError(reason, fixed_value=fixed_value)
        assert field.needs_fixing
        assert field.was_understood
        # pylint: disable-next=protected-access           # OK in tests
        assert field._needs_fix.reason == reason
        # pylint: disable-next=protected-access           # OK in tests
        assert field._needs_fix.value == fixed_value
        fixed_field = field.as_fixed()
        assert fixed_field.value == fixed_value

    def test_register_errors_unfixable(self, make_and_check_field):
        """Check correct registering of unfixable errors."""
        field = make_and_check_field(UnknownField)
        reason = 'Unfixable'
        # pylint: disable-next=protected-access           # OK in tests
        with pytest.raises(EntrySyntaxError), field._register_errors():
            raise EntrySyntaxError(reason)
        assert not field.was_understood
        # pylint: disable-next=protected-access           # OK in tests
        assert field._not_understood == reason

    def test_register_errors_multiple(self, make_and_check_field):
        """Check that multiple different exceptions are recorded."""
        fixable, unfixable = 'Fixable', 'Unfixable'
        field = make_and_check_field(UnknownField, value='OK')
        assert field.was_understood
        assert not field.needs_fixing
        # pylint: disable-next=protected-access           # OK in tests
        with pytest.raises(EntrySyntaxError), field._register_errors():
            raise EntrySyntaxError(unfixable)
        with pytest.raises(FixableSyntaxError), field._register_errors():
            raise FixableSyntaxError(fixable)
        assert not field.was_understood
        assert field.needs_fixing
        # pylint: disable-next=protected-access           # OK in tests
        assert field._not_understood == unfixable
        # pylint: disable-next=protected-access           # OK in tests
        assert field._needs_fix.reason == fixable

    _comments = {
        'no comments': (FieldBase, 'no comment', 'no comment', ''),
        'comments': (FieldBase, 'Value # comment', 'Value ', '# comment'),
        'extra': (UnknownField, 'Value # comment', 'Value # comment', ''),
        'CommentLessField': (CommentLessField, 'Value # comment',
                             'Value # comment', ''),
        }

    @parametrize('field_cls,value_str,expect_value,expect_comment',
                 _comments.values(), ids=_comments)
    def test_strip_comments(self, field_cls, value_str,
                            expect_value, expect_comment):
        """Check expected splitting of comments."""
        # pylint: disable-next=protected-access           # OK in tests
        value, comments = field_cls._strip_comments(value_str)
        assert value == expect_value
        assert comments == expect_comment

    def test_update_from_field_raises(self):
        """Check complaints when calling _update_from_field on a non-field."""
        field = UnknownField('123')
        with pytest.raises(TypeError):
            # pylint: disable-next=protected-access       # OK in tests
            field._update_from_field('abcd')


class TestFieldBaseFromString:
    """Tests for the .from_string method of FieldBase."""

    def test_from_string_valid(self):
        """Check correct parsing from string."""
        value_str = '  test '
        field = FieldBase.from_string(value_str)
        field.check_value()
        # pylint: disable-next=unidiomatic-typecheck
        assert type(field) is UnknownField   # We want exact type match
        assert field.value == value_str.rstrip()
        assert not field.has_comments
        assert not field.is_missing
        assert not field.is_empty
        assert not field.needs_fixing
        assert field.was_understood

    _cant_parse = {
        'multiline': (FieldTag.NOTES, 'Some text\n On multiple lines'),
        'on line': (FieldTag.UNKNOWN, 'something without a label'),
        'invalid tag': (FieldTag.UNKNOWN, '# INVALID label'),
        }

    @parametrize('tag,value_str', _cant_parse.values(), ids=_cant_parse)
    def test_from_string_cant_parse(self, tag, value_str, remove_field_tag):
        """Check complaints when we can't parse a string."""
        remove_field_tag(tag)
        with pytest.raises(HistoryInfoError, match='No field can parse'):
            FieldBase.from_string(value_str)

    def test_from_string_cant_parse_wrong_regex(self, monkeypatch):
        """Check complaints when we can't parse a string for regex issues."""
        reason = 'No field can parse'
        monkeypatch.setattr(UnknownField, 'rgx_pattern', 'a__')
        with pytest.raises(HistoryInfoError, match=reason):
            FieldBase.from_string('does not match')

    def test_from_string_cant_parse_case(self, make_concrete_field,
                                         remove_field_tag):
        """Check handling of case in tags."""
        make_concrete_field(FieldBase, MockFieldTag.TAG_1)
        upper_str = f'{MockFieldTag.TAG_1.value} abcde'
        remove_field_tag(FieldTag.UNKNOWN, raising=True)
        with pytest.raises(HistoryInfoError, match='No field can parse'):
            FieldBase.from_string(upper_str.lower())

    @parametrize(tag=(*iter(FieldTag), *iter(MockFieldTag)))
    def test_from_string_subclass(self, tag, make_concrete_field, monkeypatch):
        """Check correct string parsing via a FieldBase subclass."""
        try:
            field_cls = FieldBase.for_tag(tag)
        except ValueError:  # A mock
            field_cls = make_concrete_field(FieldBase, tag)
        value = 'test_value'
        value_str = f'{tag.value}{value}'
        field = FieldBase.from_string(value_str)
        # pylint: disable-next=unidiomatic-typecheck
        assert type(field) is field_cls      # We want exact type match
        assert field.tag is tag
        assert field.value == value


class TestFieldBaseSubclasses:
    """Tests for custom subclasses of FieldBase."""

    @fixture(name='custom_fmt')
    def fixture_custom_fmt(self, make_concrete_field, make_and_check_field):
        """Return a subclass that modifies the standard string format."""
        @frozen
        class _CustomFormat(FieldBase):
            custom_: str = 'is customized:'
            def _format_string_value(self, value_str):
                return f'{self.tag.value} {self.custom_} {value_str}'
        field_cls = make_concrete_field(_CustomFormat, tag=MockFieldTag.TAG_3)
        return functools.partial(make_and_check_field, field_cls)

    @fixture(name='digit_str')
    def fixture_digit_str(self, make_concrete_field, make_and_check_field):
        """Return a subclass that accepts only string with digit values."""
        @frozen
        class _DigitOnly(FieldBase):
            def _check_str_value(self):
                super()._check_str_value()
                # pylint: disable-next=no-member  # Can't infer
                if not self.value.isdigit():
                    raise EntrySyntaxError('Not a valid digit')
        field_cls = make_concrete_field(_DigitOnly, tag=MockFieldTag.TAG_2)
        return functools.partial(make_and_check_field, field_cls)

    @fixture(name='store_str')
    def fixture_store_str(self, make_concrete_field, make_and_check_field,
                          monkeypatch):
        """Return a subclass that remembers its string value."""
        @frozen
        class _StoresStr(FieldBase):
            def _check_str_value(self):
                super()._check_str_value()
                set_frozen_attr(self, '_value_str',
                                f'custom+{self.value}+custom')
        field_cls = make_concrete_field(_StoresStr, tag=MockFieldTag.TAG_1)
        return functools.partial(make_and_check_field, field_cls)

    @fixture(name='wrong_check')
    def fixture_wrong_check(self, make_concrete_field, make_and_check_field):
        """Return a subclass that may call _check_str_value on non-strings."""
        @frozen
        class _WrongCheck(FieldBase):
            value: Union[int, str] = MissingField
            def _check_value(self):
                super()._check_value()
                self._check_str_value()
            def _check_int_value(self):
                pass
        field_cls = make_concrete_field(_WrongCheck, tag=MockFieldTag.TAG_1)
        return functools.partial(make_and_check_field, field_cls)

    def test_abstract_instance(self):
        """Check complaints when instantiating a tag-less subclass."""
        TagLess = type('TagLess', (FieldBase,), {})
        with pytest.raises(TypeError):
            TagLess()

    def test_already_handled(self):
        """Check complaints when subclassing FieldBase for the same tag."""
        with pytest.raises(ValueError):
            _ = type('AltreadyHandled', (FieldBase,), {}, tag=FieldTag.UNKNOWN)

    def test_custom_format_field(self, custom_fmt):
        """Check correct behavior of a subclass with custom format."""
        value = 'test value'
        field = custom_fmt(value)
        assert field.value == value
        _str = str(field)
        assert _str.endswith(value)
        assert _str.startswith(field.tag.value)
        sep = _str.replace(value, '').replace(field.tag.value, '')
        assert sep.strip() == field.custom_

    def test_digit_only_field(self, digit_str):
        """Check correct behavior of a digit-only subclass."""
        field = digit_str(value='1234')
        with not_raises((EntrySyntaxError, FixableSyntaxError)):
            field.check_value()
        field = digit_str(value='some text with 123 digits and words')
        with pytest.raises(EntrySyntaxError):
            field.check_value()

    def test_store_str(self, store_str):
        """Check correct behavior of a subclass that edits _value_str."""
        value, comment = 'test value', '## some comment chars   '
        from_value = store_str(value)
        assert from_value.value == value
        assert not from_value.has_comments
        _str = str(from_value)
        assert value in _str
        assert not _str.endswith(value)
        # pylint: disable-next=magic-value-comparison
        assert _str.count('custom') == 2

        value_str = f'{from_value.tag.value}{value}{comment}'
        from_string = FieldBase.from_string(value_str)
        assert type(from_string) is type(from_value)
        assert from_string.value == value
        assert from_string.has_comments
        assert str(value_str) != value_str.rstrip()  # Is customized

    _subclass = {
        'mock 1': (type('MockField_1', (FieldBase,), {}), MockFieldTag.TAG_1),
        'mock 2': (type('MockField_2', (FieldBase,), {}), MockFieldTag.TAG_2),
        'unknown': (UnknownField, FieldTag.UNKNOWN),
        }

    @parametrize('field_cls, tag', _subclass.values(), ids=_subclass)
    def test_subclass(self, field_cls, tag,
                      make_concrete_field,
                      make_and_check_field):
        """Check subclassing of FieldBase."""
        concrete_cls = make_concrete_field(field_cls, tag)
        value = 'value'
        field = make_and_check_field(concrete_cls, value)
        assert field.tag is tag
        assert field.value == value

    def test_subclass_from_string(self, remove_field_tag):
        """Check correct recognition of a subclass from a string value."""
        remove_field_tag(FieldTag.TENSOR_NUMS)
        @frozen
        class _MockTensorsField(FieldBase, tag=FieldTag.TENSOR_NUMS):
            pass
        value = ' \t test value'
        value_str = f'{FieldTag.TENSOR_NUMS.value}{value}  '
        field = FieldBase.from_string(value_str)
        field.check_value()
        # pylint: disable-next=unidiomatic-typecheck
        assert type(field) is _MockTensorsField  # Exact type match
        assert field.tag == FieldTag.TENSOR_NUMS
        assert field.value == value

    def test_wrong_str_check(self, wrong_check):
        """Check complaints when calling _check_str_value on non-strings."""
        str_field = wrong_check('test')
        with not_raises((EntrySyntaxError, FixableSyntaxError)):
            str_field.check_value()

        int_field = wrong_check(123)
        reason = DefaultMessage.NOT_STRING.value
        with pytest.raises(EntrySyntaxError, match=reason):
            int_field.check_value()


class TestMultilineField:
    """Tests for (concrete subclasses of) MultiLineField."""

    @fixture(name='multiline')
    def fixture_multiline(self, make_concrete_field, make_and_check_field):
        """Return a concrete subclass of MultiLineField."""
        field_cls = make_concrete_field(MultiLineField,
                                        tag=MockFieldTag.TAG_1)
        return functools.partial(make_and_check_field, field_cls)

    _init = {
        'two lines': 'line1\nline2',
        'special chars': 'line1\nline2\twith\tescape\\sequence',
        'three lines, comment': 'line1\n line2\n line3',
        'many lines': 'line\n' * 1000,
        'one line': 'single line without any line breaks',
        'empty': '',
        }

    @parametrize(value=_init.values(), ids=_init)
    def test_init(self, value, multiline):
        """Check correct initialization of a multi-line field."""
        field = multiline(value)
        _empty = not value
        _str = str(field)
        sep = _str.replace(value.rstrip(), '').replace(field.tag.value, '')
        assert field.value == (value or EmptyField)
        assert _str.endswith(value.rstrip())
        assert _str.startswith(field.tag.value)
        assert not sep.strip()
        assert not field.is_missing
        assert field.is_empty is _empty

    def test_format_faulty(self, multiline):
        """Check formatting of faulty field."""
        value = 'line one \n line two'
        field = multiline(value)
        fmt = field.format_faulty()
        assert all(line.rstrip() in fmt for line in value.splitlines())
        assert fmt.count('\n') == value.count('\n')

    def test_single_line(self, multiline):
        """Check there are no line breaks for a value without breaks."""
        field = multiline('no line breaks')
        # pylint: disable-next=magic-value-comparison
        assert '\n' not in str(field)


class TestUnknownField:
    """Tests for the UnknownField subclass of FieldBase."""

    def test_format_string_value(self, make_and_check_field):
        """Check expected outcome of _format_string_value."""
        value = '"unknown value" # and has ! comments'
        field = make_and_check_field(UnknownField, value=value)
        # pylint: disable-next=protected-access           # OK in tests
        formatted_value = field._format_string_value(field.value)
        assert formatted_value == value

    def test_minimal_input(self, make_and_check_field):
        """Check correct handling of a short input value."""
        minimal_input = '?'
        field = make_and_check_field(UnknownField, value=minimal_input)
        assert str(field) == minimal_input
