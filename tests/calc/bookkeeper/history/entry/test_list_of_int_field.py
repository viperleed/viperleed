"""Test list_of_int_field of viperleed.calc.bookkeeper.history.entry."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2024-09-04'
__license__ = 'GPLv3+'

import ast
from inspect import isclass

import pytest
from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.calc.bookkeeper.history.entry.enums import FieldTag
from viperleed.calc.bookkeeper.history.entry.field import FieldBase
from viperleed.calc.bookkeeper.history.entry.field import DefaultMessage
from viperleed.calc.bookkeeper.history.entry.field import MissingField
from viperleed.calc.bookkeeper.history.entry.list_of_int_field import (
    CommaSeparatedIntsField,
    JobIdsField,
    ListOfIntsField,
    PositiveIntsField,
    RunInfoField,
    SpaceSeparatedIntsField,
    TensorNumsField,
    )
from viperleed.calc.bookkeeper.history.errors import EntrySyntaxError
from viperleed.calc.sections.calc_section import CalcSection

from .test_field import _TestFieldUtils


class _TestAbstractBase(_TestFieldUtils):
    """Base class for tests of abstract ListOfIntsField classes."""

    test_cls = None

    @fixture(name='list_of_int')
    def fixture_list_of_int(self, make_concrete_field_instance):
        """Return an instance of a concrete subclass of ListOfIntsField."""
        subclass = self.make_concrete_subclass()
        return make_concrete_field_instance(subclass)

    @classmethod
    def make_concrete_subclass(cls):
        """Return a concrete, tag-less subclass of self.test_cls."""
        concrete = cls.test_cls
        method = getattr(concrete, '_clean_up_string')
        is_abstract = getattr(method, '__isabstractmethod__', False)
        if not is_abstract:
            # test_cls is already concrete
            return concrete

        abstract = cls.test_cls
        if not isclass(abstract) or not issubclass(abstract, FieldBase):
            raise RuntimeError(f'{cls.__name__}: Must define a non-None '
                               'class attribute test_cls that must be a'
                               'class.')
        # pylint: disable-next=inherit-non-class,too-few-public-methods
        class _Subclass(abstract):
            def _clean_up_string(self):
                try:
                    return ast.literal_eval(self.value + ',')
                except (TypeError, ValueError,
                        SyntaxError, MemoryError) as exc:
                    raise EntrySyntaxError(str(exc)) from exc
        return _Subclass

    def test_is_sequence_of_int_raises(self, list_of_int):
        """Ensure complaints when calling _check_is_sequence_of_int."""
        field = list_of_int(123.456)
        with pytest.raises(EntrySyntaxError):
            # pylint: disable-next=protected-access       # OK in tests
            field._check_is_sequence_of_int()


class TestListOfIntsField(_TestAbstractBase):
    """Test for the abstract ListOfIntsField base class."""

    test_cls = ListOfIntsField

    init = {
        'empty set': (
            set(),
            {'was_understood': False, 'is_empty': True},
            ),
        'missing': (
            MissingField,
            {'is_missing': True, 'was_understood': True, 'is_empty': False},
            ),
        'none': (
            None,
            {'is_empty': True},
            ),
        'non-handled sequence': (
            {1, 2, 3},
            {'value': {1, 2, 3}, 'was_understood': False},
            ),
        'not a sequence': (
            123.456,
            {'value': 123.456, 'was_understood': False, '_value_str': None},
            ),
        'one item': (
            1,
            {'value': (1,), 'was_understood': True, '_value_str': '1'},
            ),
        'sequence of int strings': (
            ('1', '2', '3'),
            {'value': ('1', '2', '3'), 'was_understood': False,
             '_value_str': None},
            ),
        'valid string value': (
            '1, 2, 3',
            {'value': (1, 2, 3), 'was_understood': True,
             '_value_str': '1, 2, 3'},
            ),
        'valid': (
            [1, 2, 3],
            {'value': (1, 2, 3), 'was_understood': True,
             'needs_fixing': False, '_value_str': '1, 2, 3'},
            ),
        }

    @parametrize('value,attrs', init.values(), ids=init)
    def test_init(self, value, attrs, list_of_int):
        """Check attributes after initialization."""
        self.check_attrs(list_of_int, attrs, value)

    def test_no_loose_regex(self, list_of_int):
        """Check complaints when trying to clean up a string without regex."""
        field = list_of_int('abcd')
        with pytest.raises(NotImplementedError):
            # pylint: disable-next=protected-access       # OK in tests
            field._clean_and_validate_string_loose('separator')


class TestCommaSeparatedIntsField(_TestAbstractBase):
    """Tests for the tag-less CommaSeparatedIntsField subclass."""

    test_cls = CommaSeparatedIntsField

    _clean = {
        '1, 2, 3': {'value': (1, 2, 3), 'needs_fixing': False,
                    'was_understood': True},
        '1': {'value': (1,), '_value_str': '1', 'needs_fixing': False,
              'was_understood': True},
        12: {'value': (12,), '_value_str': '12', 'needs_fixing': False,
             'was_understood': True},
        # The next ones are invalid
        '1 2 3': {'value': '1 2 3', 'needs_fixing': True,
                  'was_understood': True},
        '1, 2 3': {'value': '1, 2 3', 'needs_fixing': True,
                   'was_understood': True},
        '1, two, 3': {'value': '1, two, 3', 'needs_fixing': False,
                      'was_understood': False},
        }

    @parametrize('value,attrs', _clean.items(), ids=_clean)
    def test_clean_up_string(self, value, attrs, list_of_int):
        """Check the result of cleaning up a string value."""
        self.check_attrs(list_of_int, attrs, value)


class TestPositiveIntsField(_TestAbstractBase):
    """Tests for the abstract PositiveIntsField subclass of ListOfIntsField."""

    test_cls = PositiveIntsField

    init = {
        **TestListOfIntsField.init,
        'all positive': (
            [10, 11, 12],
            {'value': (10, 11, 12), '_value_str': '10, 11, 12',
             'was_understood': True},
            ),
        'negative': (
            (-1, 2, 3),
            {'value': (-1, 2, 3), '_value_str': None, 'was_understood': False},
            ),
        }

    @parametrize('value,attrs', init.values(), ids=init)
    def test_init(self, value, attrs, list_of_int):
        """Check attributes after initialization."""
        self.check_attrs(list_of_int, attrs, value)


class TestSpaceSeparatedIntsField(_TestAbstractBase):
    """Tests for the tag-less SpaceSeparatedIntsField subclass."""

    test_cls = SpaceSeparatedIntsField
    init = {
        'commas': (
            '1, 2, 3',
            {'value': '1, 2, 3', '_value_str': None, 'was_understood': True,
             'needs_fixing': True},
            ),
        'invalid': (
            '1 two 3',
            {'value': '1 two 3', '_value_str': None, 'was_understood': False,
             'needs_fixing': False},
            ),
        'mixed': (
            '1 2, 3',
            {'value': '1 2, 3', '_value_str': None, 'was_understood': True,
             'needs_fixing': True},
            ),
        'one value': (
            159,
            {'value': (159,), '_value_str': '159', 'was_understood': True,
             'needs_fixing': False},
            ),
        'valid string': (
            '1 2 3',
            {'value': (1, 2, 3), '_value_str': '1 2 3', 'was_understood': True,
             'needs_fixing': False},
            ),
        'valid sequence': (
            [1, 2, 3],
            {'value': (1, 2, 3), '_value_str': '1 2 3', 'was_understood': True,
             'needs_fixing': False},
            ),
        }

    @parametrize('value,attrs', init.values(), ids=init)
    def test_clean_up_string(self, value, attrs, list_of_int):
        """Check the result of cleaning up a string value."""
        self.check_attrs(list_of_int, attrs, value)


class TestJobIdsField(_TestAbstractBase):
    """Tests for the concrete, tag-bearing JobIdsField subclass."""

    test_cls = JobIdsField
    init = {
        **TestPositiveIntsField.init,
        'missing': (
            MissingField,
            {'is_missing': True, 'was_understood': False, 'is_empty': False},
            ),
        }

    def test_class_attrs(self):
        """Check the expected class-level attributes of JobIdsField."""
        field = self.test_cls
        assert field.is_mandatory
        assert field.tag is FieldTag.JOB_NUMS

    @parametrize('value,attrs', init.values(), ids=init)
    def test_init(self, value, attrs, list_of_int):
        """Check attributes after initialization."""
        self.check_attrs(list_of_int, attrs, value)


class TestRunInfoField(_TestAbstractBase):
    """Tests for the concrete, tag-bearing RunInfoField subclass."""

    test_cls = RunInfoField

    def test_class_attrs(self):
        """Check the expected class-level attributes of JobIdsField."""
        field = self.test_cls
        assert not field.is_mandatory
        assert field.tag is FieldTag.RUN_INFO

    init = {
        **TestSpaceSeparatedIntsField.init,
        'invalid section': (
            (1, 999),
            {'value': (1, 999), 'was_understood': False},
            ),
        'one value': (
            CalcSection.REFCALC.value,
            {'value': (CalcSection.REFCALC.value,),
             '_value_str': f'{CalcSection.REFCALC.value}',
             'was_understood': True, 'needs_fixing': False},
            ),
        'section value': (
            [CalcSection.INITIALIZATION.value],
            {'value': (CalcSection.INITIALIZATION.value,),
             'was_understood': True, 'needs_fixing': False},
            ),
        }

    @parametrize('value,attrs', init.values(), ids=init)
    def test_init(self, value, attrs, list_of_int):
        """Check attributes after initialization."""
        self.check_attrs(list_of_int, attrs, value)


class TestTensorNumsField(TestCommaSeparatedIntsField):
    """Tests for the concrete, tag-bearing TensorNumsField subclass."""

    test_cls = TensorNumsField
    _none_attrs = {'value': 'None', '_value_str': 'None',
                   'was_understood': True, 'needs_fixing': False,
                   'no_tensors': True}
    init = {
        **TestPositiveIntsField.init,
        'init list': ([CalcSection.INITIALIZATION.value], _none_attrs),
        'init tuple': ((CalcSection.INITIALIZATION.value,), _none_attrs),
        'init value': (CalcSection.INITIALIZATION.value, _none_attrs),
        'missing': (
            MissingField,
            {'is_missing': True, 'was_understood': False, 'is_empty': False},
            ),
        'none': (None, _none_attrs,),
        'none string': ('None', _none_attrs),
        'none list': ([None], _none_attrs),
        'none tuple': ((None,), _none_attrs),
        }

    @parametrize('value,attrs', init.values(), ids=init)
    def test_init(self, value, attrs, list_of_int):
        """Check attributes after initialization."""
        self.check_attrs(list_of_int, attrs, value)

    def test_check_not_empty_with_empty_value(self):
        """Check complaints when an empty value is given."""
        field = self.test_cls(value=[])
        with pytest.raises(EntrySyntaxError, match=DefaultMessage.EMPTY.value):
            # pylint: disable-next=protected-access       # OK in tests
            field._check_not_empty()
