"""Tests for module field of viperleed.calc.bookkeeper.history.entry."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-08-30'
__license__ = 'GPLv3+'

from dataclasses import FrozenInstanceError
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
from viperleed.calc.bookkeeper.history.entry.field import NoneIsEmptyField
from viperleed.calc.bookkeeper.history.entry.field import UnknownField
from viperleed.calc.lib.dataclass_utils import frozen
from viperleed.calc.lib.dataclass_utils import set_frozen_attr

from .....helpers import not_raises
from .conftest import MockFieldTag


@fixture(name='unknown')
def make_unknown(make_field_factory):
    """Return a factory of UnknownField instances."""
    return make_field_factory(UnknownField)


class _TestFieldUtils:
    """Collection of utility functions for testing FieldBase."""

    test_cls = None

    def check_attrs(self, field_factory, attrs, *args, **kwargs):
        """Check that `field_factory(value)` has the given `attrs`.

        Parameters
        ----------
        field_factory : callable
            Takes *args, **kwargs and returns a FieldBase instance
            on which check_value was run.
        attrs : dict
            Keys are attribute names, values their expected value.
        *args : object
            Positional arguments passed to `field_factory`.
        **kwargs : object
            Keyword arguments passed to `field_factory`.

        Returns
        -------
        field : FieldBase
            The result of field_factory(*args, **kwargs)
        """
        field = field_factory(*args, **kwargs)
        for attr, expect in attrs.items():
            assert getattr(field, attr) == expect
        return field


class TestCommentlessField(_TestFieldUtils):
    """Tests for (concrete subclasses of) CommentLessField."""

    test_cls = CommentLessField

    @fixture(name='commentless')
    def fixture_commentless(self, make_concrete_field_instance):
        """Return a concrete subclass of CommentLessField."""
        return make_concrete_field_instance(self.test_cls)

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
        field = self.check_attrs(commentless, expect, value)
        assert not field.has_comments
        assert not field.needs_fixing

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


class TestFieldBase(_TestFieldUtils):
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
    def test_init(self, args, expect, make_field_factory):
        """Check initialization of a field."""
        field_cls, *args = args
        self.check_attrs(make_field_factory(field_cls), expect, *args)

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
    def test_big_value(self, value, unknown):
        """Check handling of large values."""
        field = unknown(value)
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
    def test_format_faulty_non_mandatory(self, value, expect, unknown):
        """Check formatting of faulty field."""
        field = unknown(value=value)
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

    def test_invalid_type_hint(self, make_concrete_field_instance):
        """Check complaints when an invalid type is given."""
        @frozen
        class _Dummy(FieldBase):
            value: str = MissingField
        factory = make_concrete_field_instance(_Dummy)
        field = factory(123)
        assert not field.was_understood

    def test_register_errors_fixable(self, unknown):
        """Check correct registering of fixable errors."""
        field = unknown()
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

    def test_register_errors_unfixable(self, unknown):
        """Check correct registering of unfixable errors."""
        field = unknown()
        reason = 'Unfixable'
        # pylint: disable-next=protected-access           # OK in tests
        with pytest.raises(EntrySyntaxError), field._register_errors():
            raise EntrySyntaxError(reason)
        assert not field.was_understood
        # pylint: disable-next=protected-access           # OK in tests
        assert field._not_understood == reason

    def test_register_errors_multiple(self, unknown):
        """Check that multiple different exceptions are recorded."""
        fixable, unfixable = 'Fixable', 'Unfixable'
        field = unknown(value='OK')
        assert field.was_understood
        assert not field.needs_fixing
        # pylint: disable-next=protected-access           # OK in tests
        with pytest.raises(EntrySyntaxError), field._register_errors():
            raise EntrySyntaxError(unfixable)
        # pylint: disable-next=protected-access           # OK in tests
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

    def test_from_string_cant_parse_lowercase(self, make_concrete_field,
                                              remove_field_tag):
        """Check handling of case in tags."""
        make_concrete_field(FieldBase, MockFieldTag.TAG_1)
        upper_str = f'{MockFieldTag.TAG_1.value} abcde'
        remove_field_tag(FieldTag.UNKNOWN, raising=True)
        with pytest.raises(HistoryInfoError, match='No field can parse'):
            FieldBase.from_string(upper_str.lower())

    @parametrize(tag=(*iter(FieldTag), *iter(MockFieldTag)))
    def test_from_string_subclass(self, tag, make_concrete_field):
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
        field_value = field.value
        if issubclass(field_cls, MultiLineField):  # Value is a tuple
            field_value, *rest = field_value
            assert not rest
        assert field_value == value


class TestFieldBaseSubclasses:
    """Tests for custom subclasses of FieldBase."""

    @fixture(name='custom_fmt')
    def fixture_custom_fmt(self, make_concrete_field_instance):
        """Return a subclass that modifies the standard string format."""
        @frozen
        class _CustomFormat(FieldBase):
            custom_: str = 'is customized:'
            def _format_string_value(self, value_str):
                return f'{self.tag.value} {self.custom_} {value_str}'
        return make_concrete_field_instance(_CustomFormat,
                                            tag=MockFieldTag.TAG_3)

    @fixture(name='digit_str')
    def fixture_digit_str(self, make_concrete_field_instance):
        """Return a subclass that accepts only string with digit values."""
        @frozen
        class _DigitOnly(FieldBase):
            def _check_str_value(self):
                super()._check_str_value()
                # pylint: disable-next=no-member  # Can't infer
                if not self.value.isdigit():
                    raise EntrySyntaxError('Not a valid digit')
        return make_concrete_field_instance(_DigitOnly, tag=MockFieldTag.TAG_2)

    @fixture(name='empty_ok')
    def fixture_empty_ok(self, make_concrete_field_instance):
        """Return a subclass that may call _check_str_value on non-strings."""
        @frozen
        class _EmptyOK(FieldBase):
            def _check_not_empty(self):
                try:
                    super()._check_not_empty()
                except EntrySyntaxError:
                    pass
        return make_concrete_field_instance(_EmptyOK)

    @fixture(name='store_str')
    def fixture_store_str(self, make_concrete_field_instance):
        """Return a subclass that remembers its string value."""
        @frozen
        class _StoresStr(FieldBase):
            def _check_str_value(self):
                super()._check_str_value()
                set_frozen_attr(self, '_value_str',
                                f'custom+{self.value}+custom')
        return make_concrete_field_instance(_StoresStr, tag=MockFieldTag.TAG_1)

    @fixture(name='wrong_check')
    def fixture_wrong_check(self, make_concrete_field_instance):
        """Return a subclass that may call _check_str_value on non-strings."""
        @frozen
        class _WrongCheck(FieldBase):
            value: Union[int, str] = MissingField
            def _check_value(self):
                super()._check_value()
                self._check_str_value()
            def _check_int_value(self):
                pass
        return make_concrete_field_instance(_WrongCheck)

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

    def test_empty_ok_field(self, empty_ok):
        """Check that there are no complaints for an empty acceptable field."""
        field = empty_ok('')
        assert field.is_empty
        with not_raises(EntrySyntaxError):
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
    def test_subclass(self, field_cls, tag, make_concrete_field_instance):
        """Check subclassing of FieldBase."""
        value = 'value'
        factory = make_concrete_field_instance(field_cls, tag=tag)
        field = factory(value)
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


class TestMultilineField(_TestFieldUtils):
    """Tests for (concrete subclasses of) MultiLineField."""

    test_cls = MultiLineField

    @fixture(name='multiline')
    def fixture_multiline(self, make_concrete_field_instance):
        """Return a concrete subclass of MultiLineField."""
        return make_concrete_field_instance(self.test_cls)

    _init = {
        'two lines': 'line1\nline2',
        'special chars': 'line1\nline2\twith\tescape\\sequence',
        'three lines, comment': 'line1   \n line2\n line3 # comment',
        'many lines': 'line\n' * 1000,
        'one line': 'single line without any line breaks',
        'empty': '',
        'tuple multi lines': (
            'line1', 'line2', '', 'line4 after an empty one',
            '', '', 'line # with comments', '', '', ''
            ),
        'tuple empty': (),
        }

    @parametrize(value=_init.values(), ids=_init)
    def test_init(self, value, multiline):
        """Check correct initialization of a multi-line field."""
        field = multiline(value)
        _empty = not value
        _value = EmptyField if _empty else value
        if isinstance(_value, str):
            assert not _empty
            # pylint: disable-next=redefined-variable-type
            _value = tuple(line.rstrip() for line in value.splitlines())
            _value += ('',) if value.endswith('\n') else ()
        _value_as_str = '' if _empty else '\n'.join(_value)
        assert field.value == _value
        assert not field.is_missing
        assert field.is_empty is _empty
        self._check_string_value(field, _value_as_str)

    _init_invalid = {
        'wrong type, dict': ({}, {'was_understood': False}),
        'wrong type, number': (1.234, {'was_understood': False}),
        'tuple, not strings': ((1, 2, '3'), {'was_understood': False}),
        }

    @parametrize('value,attrs', _init_invalid.values(), ids=_init_invalid)
    def test_init_invalid(self, value, attrs, multiline):
        """Check complaints when initializing with an invalid value."""
        self.check_attrs(multiline, attrs, value)

    @parametrize('value,_', _init_invalid.values(), ids=_init_invalid)
    def test_str_with_invalid_value(self, value, _, multiline):
        """Check the string version of a MultiLineField with invalid value."""
        field = multiline(value)
        # pylint: disable-next=protected-access           # OK in tests
        lines = field._prepare_lines_for_str()
        assert lines is None

    def test_check_tuple_too_early(self, multiline):
        """Check complaints when _check_tuple_value is called on a string."""
        value = 'abc'
        field = multiline(value)
        # Since multiline calls check_value, field is already OK.
        # Simulate a missing call by resetting the value to a string
        set_frozen_attr(field, 'value', value)
        with pytest.raises(TypeError):
            # pylint: disable-next=protected-access       # OK in tests
            field._check_tuple_value()

    def _check_string_value(self, field, value_as_str):
        """Check expected form of the string version of a multi-line field."""
        _str = str(field)
        sep = _str.replace(value_as_str, '').replace(field.tag.value, '')
        assert _str.endswith(value_as_str)
        assert _str.startswith(field.tag.value)
        assert not sep.strip()  # Only whitespace between tag and value

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


class TestNoneIsEmptyField(_TestFieldUtils):
    """Tests for the NoneIsEmptyField abstract subclass of FieldBase."""

    test_cls = NoneIsEmptyField

    @fixture(name='none_empty_field')
    def fixture_none_empty_field(self, make_concrete_field_instance):
        """Return an instance of a concrete subclass of NoneIsEmptyField."""
        return make_concrete_field_instance(self.test_cls)

    _init = {
        'None': (
            None,
            {'value': EmptyField, 'is_empty': True, 'was_understood': False},
            ),
        'non-None string': (
            'valid value',
            {'value': 'valid value', 'is_empty': False,
             'was_understood': True},
            ),
        'empty': (
            '',
            {'value': EmptyField, 'is_empty': True, 'was_understood': False},
            ),
        'empty dict': (
            {},
            {'value': {}, 'is_empty': False, 'was_understood': False},
            ),
        'empty list': (
            [],
            {'value': [], 'is_empty': False, 'was_understood': False},
            ),
        'empty tuple': (
            (),
            {'value': (), 'is_empty': False, 'was_understood': False},
            ),
        'false': (
            False,
            {'value': False, 'is_empty': False, 'was_understood': False},
            ),
        'zero': (
            0,
            {'value': 0, 'is_empty': False, 'was_understood': False},
            ),
        }

    @parametrize('value,attrs', _init.values(), ids=_init)
    def test_init(self, value, attrs, none_empty_field):
        """Check attributes after initialization and checking."""
        self.check_attrs(none_empty_field, attrs, value)


class TestUnknownField:
    """Tests for the UnknownField subclass of FieldBase."""

    def test_format_string_value(self, unknown):
        """Check expected outcome of _format_string_value."""
        value = '"unknown value" # and has ! comments'
        field = unknown(value=value)
        # pylint: disable-next=protected-access           # OK in tests
        formatted_value = field._format_string_value(field.value)
        assert formatted_value == value

    def test_minimal_input(self, unknown):
        """Check correct handling of a short input value."""
        minimal_input = '?'
        field = unknown(value=minimal_input)
        assert str(field) == minimal_input
