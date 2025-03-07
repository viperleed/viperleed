"""Tests for module time_field of viperleed.calc.bookkeeper.history.entry."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2024-09-40'
__license__ = 'GPLv3+'

from datetime import datetime
from enum import Enum
from enum import auto

import pytest
from pytest_cases import case
from pytest_cases import fixture
from pytest_cases import parametrize
from pytest_cases import parametrize_with_cases

from viperleed.calc.bookkeeper.history.entry.field import MissingField
from viperleed.calc.bookkeeper.history.entry.time_field import TimestampField
from viperleed.calc.bookkeeper.history.entry.time_field import TimestampFormat
from viperleed.calc.bookkeeper.history.errors import EntrySyntaxError

from .....helpers import exclude_tags


@fixture(name='make_time')
def fixture_make_time(make_field_factory):
    """Return a factory of TimestampField instances."""
    return make_field_factory(TimestampField)


class CasesTimestampField:
    """Test cases for the TimestampField class."""

    _formats = set(TimestampFormat)
    _formats.add(TimestampFormat.DEFAULT)
    all_formats = parametrize(fmt=_formats)

    class Tag(Enum):
        """Tags for test cases."""
        NO_FMT = auto()
        STR = auto()
        INVALID = auto()

    _dates = {
        'future': datetime(2050, 1, 1),
        'long ago': datetime(1900, 1, 1),
        'very far future': datetime.max,
        # The next one causes problems with format spec %y on Windows,
        # where, according to the exception, year must be >= 1900. The
        # datetime.min is (1, 1, 1, 0, 0).
        # 'very long ago': datetime.min,
        }

    @all_formats
    @parametrize(date_time=_dates.values(), ids=_dates)
    def case_date(self, make_time, date_time, fmt):
        """Return a TimestampField with a given date-time."""
        return make_time(date_time, time_format=fmt)

    _empty = ('', '    ')

    @parametrize(str_value=_empty)
    @case(tags=(Tag.STR, Tag.NO_FMT, Tag.INVALID))
    def case_empty(self, make_time, str_value):
        """Return an empty TimestampField."""
        return make_time(str_value)

    _invalid_value = {
        'set': set(),
        'float': 1.234,
        'int': 8,
        'dict': {},
        'list': [],
        'none': None,
        }

    @parametrize(value=_invalid_value.values(), ids=_invalid_value)
    @case(tags=(Tag.NO_FMT, Tag.INVALID))
    def case_invalid_value(self, make_time, value):
        """Return a TimestampField with an invalid value."""
        return make_time(value)

    @all_formats
    @parametrize(value=_invalid_value.values(), ids=_invalid_value)
    @case(tags=Tag.INVALID)
    def case_invalid_value_with_fmt(self, make_time, value, fmt):
        """Return a TimestampField with an invalid value."""
        return make_time(value, time_format=fmt)

    @case(tags=(Tag.NO_FMT, Tag.INVALID))
    def case_missing(self, make_time):
        """Return a missing TimestampField."""
        return make_time()

    @all_formats
    def case_now(self, make_time, fmt):
        """Return a TimestampField with value of now."""
        return make_time(datetime.now(), time_format=fmt)

    @case(tags=(Tag.NO_FMT, Tag.INVALID))
    def case_now_no_fmt(self, make_time):
        """Return a format-less TimestampField with value fo now."""
        return self.case_now(make_time, None)

    _string_valid = {  # With recognizable format
        '28.07.24 15:32:45': TimestampFormat.GERMAN,
        '2024-09-04 15:32:45': TimestampFormat.ISO,
        # pylint: disable-next=protected-access           # OK in tests
        'moved-240904-153245': TimestampFormat._CALC,
        # pylint: disable-next=protected-access           # OK in tests
        '  moved-240904-153245   ': TimestampFormat._CALC,
        '2024-02-29 12:00:00': TimestampFormat.ISO,   # Leap year 2024
        '0001-01-01 00:00:00': TimestampFormat.ISO,   # Far behind
        '9999-12-31 23:59:59': TimestampFormat.ISO,   # Far future
        }

    @parametrize(value=_string_valid)
    @case(tags=Tag.STR)
    def case_string(self, make_time, value):
        """Return a TimestampField from a string value."""
        return make_time(value)

    _string_invalid = {
        'spaces': 'invalid string with spaces',
        'no timestamp': 'invalid-timestamp',
        'extra text': '2024-08-26 14:55:00 extra',
        'other format': '31/12/2024 15:32:45',
        'no time': '2024-13-45',
        'no time, white space': '2024- 02-29   ',
        'white spaces': '2024-02 -29 15:32:45',  # Valid otherwise
        'invalid month': '2024-13-01 12:00:00',
        'invalid day': '2024-02-30 12:00:00',
        'invalid day, leap year': '2023-02-29 12:00:00',
        'invalid hours': '2024-02-29 25:00:00',
        'invalid minutes': '2024-02-29 22:61:12',
        'timezone GMT': '2024-02-29 12:00:00 +0100',
        'timezone Z': '2024-02-29 12:00:00 Z',
        'ISO swapped': '01-01-2024 12:00:00',
        'ISO wrong separators': '2024/01/01 12:00:00',
        'milliseconds': '2024-02-29 12:00:00.123456',
        'microseconds': '2024-02-29 12:00:00,123',
        'missing seconds': '2024-02-29 12:60',
        }

    @parametrize(value=_string_invalid.values(), ids=_string_invalid)
    @case(tags=(Tag.STR, Tag.NO_FMT, Tag.INVALID))
    def case_string_invalid(self, make_time, value):
        """Return a TimestampField from a string value."""
        return make_time(value)

    _inconsistent = {
        'empty': '',
        'missing': MissingField,
        **_string_invalid,
        }

    @all_formats
    @parametrize(value=_inconsistent.values(), ids=_inconsistent)
    @case(tags=Tag.INVALID)
    def case_value_fmt_inconsistent(self, make_time, value, fmt):
        """Return a field with inconsistent value and format."""
        return make_time(value, time_format=fmt)


all_cases = parametrize_with_cases('field', CasesTimestampField)
case_strings = parametrize_with_cases(
    'field',
    CasesTimestampField,
    has_tag=CasesTimestampField.Tag.STR,
    )
cases_valid = parametrize_with_cases(
    'field',
    CasesTimestampField,
    filter=exclude_tags(CasesTimestampField.Tag.INVALID)
    )
cases_with_fmt = parametrize_with_cases(
    'field',
    CasesTimestampField,
    filter=exclude_tags(CasesTimestampField.Tag.NO_FMT),
    )
cases_without_fmt = parametrize_with_cases(
    'field',
    CasesTimestampField,
    has_tag=CasesTimestampField.Tag.NO_FMT,
    )


class TestTimestampFormat:
    """Tests for the TimestampFormat enumeration."""

    def test_invalid(self):
        """Check complaints when accessing an invalid TimestampFormat."""
        invalid = 'INVALID_FORMAT'
        with pytest.raises(ValueError):
            TimestampFormat(invalid)
        with pytest.raises(KeyError):
            _ = TimestampFormat[invalid]

    _writable = {
        'GERMAN': True,
        'ISO': True,
        '_CALC': False,
        }

    @parametrize('fmt_name,expect', _writable.items(), ids=_writable)
    def test_writable(self, fmt_name, expect):
        """Test the writable property of TimestampFormat."""
        fmt = TimestampFormat[fmt_name]
        assert fmt.writable is expect


class TestTimestampField:
    """Tests for the TimestampField class."""

    @parametrize_with_cases('field', CasesTimestampField.case_empty)
    def test_init_empty(self, field):
        """Check initialization of an empty TimestampField."""
        assert field.is_empty

    @cases_without_fmt
    def test_init_invalid(self, field):
        """Check initialization with invalid arguments."""
        # pylint: disable-next=protected-access           # OK in tests
        _ = field._get_string_value()
        # pylint: disable-next=protected-access           # OK in tests
        assert field._value_str is None
        assert not field.was_understood
        assert not field.needs_fixing

    @case_strings
    def test_init_str_value(self, field):
        """Check format of a time_field created from a string value."""
        if field.needs_fixing:
            assert field.time_format is not TimestampFormat.DEFAULT
            self.test_as_fixed(field)
        elif field.time_format:
            assert field.time_format is TimestampFormat.DEFAULT
        else:
            with pytest.raises(EntrySyntaxError):
                field.check_value()

    @cases_valid
    def test_init_valid(self, field):
        """Check initialization of a TimestampField."""
        fmt = field.time_format
        # pylint: disable-next=protected-access           # OK in tests
        assert field._get_string_value() == field.value.strftime(fmt.value)

    @parametrize_with_cases(
        'field',
        CasesTimestampField.case_value_fmt_inconsistent
        )
    def test_init_value_and_format_mismatch(self, field):
        """Check complaints when time_format and value don't agree."""
        assert not field.was_understood
        assert not field.needs_fixing
        assert field.time_format is None

    @cases_valid
    def test_as_fixed(self, field):
        """Check the result of fixing a TimestampField."""
        fixed_field = field.as_fixed()
        # datetime.datetime objects are immutable. We never copy.
        assert fixed_field.value is field.value
        assert fixed_field.time_format is TimestampFormat.DEFAULT
        assert not fixed_field.needs_fixing

    _new_format = {
        # pylint: disable-next=protected-access           # OK in tests
        'not writable': (TimestampFormat._CALC, TimestampFormat._CALC),
        'from name': ('GERMAN', TimestampFormat.GERMAN),
        'from name, lowercase': ('iso', TimestampFormat.ISO),
        }

    @cases_valid
    @parametrize('new_fmt,expect_fmt', _new_format.values(), ids=_new_format)
    def test_with_format_other(self, new_fmt, expect_fmt, field):
        """Check reformatting of a TimestampField."""
        new_field = field.with_format(new_fmt)
        assert new_field.value is field.value
        assert new_field.time_format is expect_fmt

    _invalid_fmt = {
        None: TypeError,
        'invalid': ValueError,
        123: TypeError,
        }

    @all_cases
    @parametrize('invalid_fmt,exc', _invalid_fmt.items())
    def test_with_format_raises(self, invalid_fmt, exc, field):
        """Check complaints when requesting an invalid TimestampFormat."""
        with pytest.raises((EntrySyntaxError, exc)) as excinfo:
            field.with_format(invalid_fmt)
        if excinfo.type is EntrySyntaxError:  # field is problematic
            assert not field.was_understood

    @cases_valid
    def test_with_format_same(self, field):
        """Check with_format with the same format."""
        # Test with the same format
        same_field = field.with_format(field.time_format)
        assert same_field is field

    def test_sanitize_string_value_invalid(self, make_time):
        """Check no cleanup happens for a non-string value."""
        invalid = object()
        field = make_time(invalid)
        # pylint: disable-next=protected-access           # OK in tests
        field._sanitize_string_value()
        assert field.value is invalid
