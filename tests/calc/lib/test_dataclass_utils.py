"""Tests for dataclasses_utils module of viperleed/calc/lib."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-08-04'
__license__ = 'GPLv3+'

from dataclasses import FrozenInstanceError
from dataclasses import InitVar
from dataclasses import dataclass
from dataclasses import field as data_field
from dataclasses import fields as data_fields
from typing import Any
from typing import ClassVar
from typing import Dict
from typing import List
from typing import NoReturn
from typing import Optional
from typing import Set
from typing import Tuple
from typing import Union

import pytest
from pytest_cases import parametrize

from viperleed.calc.lib.dataclass_utils import _find_type_origin
from viperleed.calc.lib.dataclass_utils import check_types
from viperleed.calc.lib.dataclass_utils import frozen
from viperleed.calc.lib.dataclass_utils import is_optional_field
from viperleed.calc.lib.dataclass_utils import non_init_field
from viperleed.calc.lib.dataclass_utils import replace_values
from viperleed.calc.lib.dataclass_utils import set_frozen_attr

from ...helpers import not_raises


@frozen
class SampleFrozenClass:
    """A frozen dataclass."""

    attr: int
    optional_attr: Optional[int] = None
    non_init: List[int] = non_init_field(
        default_factory=lambda: [1, 'a', {}]
        )


@dataclass
class SampleClass:
    """A non-frozen dataclass."""

    attr: int
    optional_attr: Optional[int] = None
    mutable_optional: List[Optional[int]] = data_field(default_factory=list)
    non_init_with_default: int = non_init_field(default=10)
    mutable_non_init: List[int] = non_init_field(default_factory=list)
    any_attr: Any = 3
    class_attr: ClassVar[int] = 42


@dataclass
class NestedClass:
    """A dataclass to be used as field factory for another one."""

    inner: int = 1


@dataclass
class ComplexClass:
    """A dataclass with non-trivial type hints."""

    attr: int
    optional_attr: Optional[int] = None
    list_non_init: List[Optional[int]] = non_init_field(default_factory=list)
    dict_non_init: Dict[str, List[int]] = non_init_field(default_factory=dict)
    init_only_var: InitVar[int] = 0
    class_attr: ClassVar[int] = 42
    dataclass_attr: NestedClass = non_init_field(default_factory=NestedClass)


def test_frozen_decorator():
    """Check correct behavior of the @frozen decorator."""
    instance = SampleFrozenClass(attr=1, optional_attr=2)
    with pytest.raises(FrozenInstanceError):
        instance.attr = 3


class TestFindTypeOrigin:
    """Tests for the _find_type_origin function."""

    _type_origins = {  # {hint: expected outcome}
        Optional[int]: (int, type(None)),
        List[int]: (list,),
        Union[List[int], Dict[str, Tuple[int, float]]]: (list, dict),
        Any: (object,),
        ...: (object,),
        ClassVar[int]: (int,),
        NoReturn: (),  # a _SpecialForm
        }
    _type_origins_str = {str(h): e for h, e in _type_origins.items()}

    @parametrize('hint,expect', _type_origins.items())
    def test_from_hint(self, hint, expect):
        """Check correct outcome of finding the origin of a type hint."""
        origin = tuple(_find_type_origin(hint))
        assert origin == expect

    @parametrize('hint,expect', _type_origins_str.items())
    @pytest.mark.xfail(reason='string annotations are still unsupported')
    def test_from_string(self, hint, expect):
        """Check correct outcome of finding the origin of a type hint."""
        origin = tuple(_find_type_origin(hint))
        assert origin == expect


class TestCheckTypes:
    """Collection of tests for the check_types function."""

    def test_valid_type(self):
        """Check that there are no complaints for a valid initialization."""
        instance = SampleClass(attr=1)
        with not_raises(TypeError):
            check_types(instance)

    _many_types = ({}, set(), tuple(), 3.5, 1+2j, object(), NestedClass())

    @parametrize(value=_many_types)
    def test_any_hint(self, value):
        """Check that hinting with Any never complains."""
        instance = SampleClass(attr=1, any_attr=value)
        with not_raises(TypeError):
            check_types(instance)

    def test_invalid_type(self):
        """Check complaints when an invalid type is given."""
        instance = SampleClass(attr='string', optional_attr=2)
        with pytest.raises(TypeError, match='attr'):
            check_types(instance)

    def test_dont_check_items(self):
        """Check complaints for an invalid type of list items."""
        instance = SampleClass(attr=1, mutable_optional=[None, 'string'])
        with not_raises(TypeError):  # We don't check items
            check_types(instance)

    def test_dont_check_initvar(self):
        """Check complaints for an invalid type of list items."""
        instance = ComplexClass(attr=1, init_only_var='something wrong')
        with not_raises(TypeError):
            check_types(instance)

    def test_init_only(self):
        """Check no complaints if only types of init fields are verified."""
        instance = SampleClass(1)
        instance.non_init_with_default = []  # should be int
        with pytest.raises(TypeError):
            check_types(instance)
        with not_raises(TypeError):
            check_types(instance, init_only=True)


def test_is_optional_field():
    """Check correctness of the result of is_optional_field."""
    fields = data_fields(SampleClass)
    assert not is_optional_field(fields[0])  # attr
    assert is_optional_field(fields[1])      # optional_attr


class TestNonInitField:
    """Tests for the non_init_field function."""

    def test_defaul_value(self):
        """Check the value of a non-initialization field with default."""
        instance = SampleClass(attr=5)
        # pylint: disable=magic-value-comparison  # OK for hard-coded
        assert instance.attr == 5
        assert instance.non_init_with_default == 10

    def test_repr_false(self):
        """Check that a non-initialization field does not show in repr."""
        instance = SampleClass(attr=5)
        # pylint: disable-next=magic-value-comparison  # OK for attr
        assert 'non_init_with_default' not in str(instance)

    def test_repr_true(self):
        """Check that a non-initialization field shows in repr if specified."""
        @dataclass
        class _TestClass:
            non_init: int = non_init_field(default=10, repr=True)

        instance = _TestClass()
        # pylint: disable-next=magic-value-comparison  # OK for attr
        assert 'non_init' in repr(instance)


class TestReplaceValues:
    """Tests for the replace_values function."""

    @dataclass
    class ReplaceInfo:  # pylint: disable=too-many-instance-attributes
        """Information about replacements."""

        cls_to_replace: ...
        init_args: Dict[str, Any]    # How to initialize cls_to_replace
        equal: Tuple[str]        # Attributes that should compare equal
        not_identical: Set[str]  # Equal attributes, but not "is"
        identical: List[str]     # Attributes that compare with "is"
        new_attr: str = None     # The attribute to be replaced
        new_value: Any = None    # The new value of new_attr
        replaced: Dict[str, Any] = non_init_field(repr=True,
                                                  default_factory=dict)

        def __post_init__(self):
            """Make replaced."""
            if self.new_attr is not None:
                self.replaced[self.new_attr] = self.new_value

    @staticmethod
    def _check_instances(new_instance, old_instance):
        """Test correct production of a replaced instance."""
        assert new_instance is not old_instance
        assert new_instance != old_instance

    @staticmethod
    def _check_identical_attrs(new_instance, old_instance, info):
        """Check that all identical attributes are really identical."""
        for identical_attr in info.identical:
            old_ = getattr(old_instance, identical_attr)
            assert getattr(new_instance, identical_attr) is old_

    @staticmethod
    def _check_equal_attrs(new_instance, old_instance, info):
        """Check that attributes are equal, and sometimes not identical."""
        for equal_attr in info.equal:
            old_ = getattr(old_instance, equal_attr)
            if equal_attr in info.not_identical:
                assert getattr(new_instance, equal_attr) is not old_
            else:
                assert getattr(new_instance, equal_attr) is old_
            assert getattr(new_instance, equal_attr) == old_

    @staticmethod
    def _prepare_replaced(info):
        """Return an instance of cls with init_args, and a replaced one."""
        init_args = info.init_args
        instance = info.cls_to_replace(**init_args)
        for arg, value in init_args.items():
            assert getattr(instance, arg) is value
        new_instance = replace_values(instance, **info.replaced)
        return instance, new_instance

    def _test_base(self, info):
        """Make an instance of class, its replaced version, and check them."""
        info.identical.remove(info.new_attr)
        instance, new_instance = self._prepare_replaced(info)
        self._check_instances(new_instance, instance)
        self._check_identical_attrs(new_instance, instance, info)
        self._check_equal_attrs(new_instance, instance, info)

        new_attr = info.new_attr
        assert getattr(new_instance, new_attr) != getattr(instance, new_attr)

    _simple_repl = {
        'attr': 15,
        'optional_attr': 'replaced',
        'any_attr': object(),
        }

    @parametrize('new_attr,new_value', _simple_repl.items(), ids=_simple_repl)
    def test_simple_replace(self, new_attr, new_value):
        """Check replacements on a simple class."""
        info = self.ReplaceInfo(
            cls_to_replace=SampleClass,
            init_args={
                'attr': object(),
                'optional_attr': object(),
                'mutable_optional': [1, 2, 3]
                },
            equal=('mutable_optional', 'mutable_non_init'),
            not_identical={'mutable_non_init',},
            identical=['attr', 'optional_attr', 'non_init_with_default',
                       'any_attr', 'class_attr'],
            new_attr=new_attr,
            new_value=new_value,
            )
        self._test_base(info)

    _frozen_repl = {
        'attr': 15,
        'optional_attr': 'replaced',
        }

    @parametrize('new_attr,new_value', _frozen_repl.items(), ids=_frozen_repl)
    def test_replace_frozen(self, new_attr, new_value):
        """Check correct replacement of attributes from a frozen class."""
        info = self.ReplaceInfo(
            cls_to_replace=SampleFrozenClass,
            init_args={'attr': object(),
                       'optional_attr': object()},
            equal=('non_init',),
            not_identical={'non_init',},
            identical=['attr', 'optional_attr'],
            new_attr=new_attr,
            new_value=new_value,
            )
        self._test_base(info)

    def test_replace_skip(self):
        """Check correct skipping of specified non-init fields."""
        info = self.ReplaceInfo(
            cls_to_replace=ComplexClass,
            init_args={
                'attr': object(),
                'optional_attr': object(),
                'init_only_var': object(),
                },
            equal=('list_non_init',
                   # 'dict_non_init',  # This is the one we skip
                   'dataclass_attr'),
            not_identical=('list_non_init',
                           'dataclass_attr'),
            identical=('attr',
                       # 'optional_attr',  # This is the one we replace
                       'init_only_var',
                       'class_attr'),
            new_attr='optional_attr',
            )
        replaced = {'optional_attr': NestedClass(),
                    'init_only_var': info.init_args['init_only_var']}
        skip = 'dict_non_init'
        new_attr = info.new_attr

        instance = ComplexClass(**info.init_args)
        instance.dict_non_init = {i: i for i in range(5)}
        new_instance = replace_values(instance, skip=skip, **replaced)
        self._check_instances(new_instance, instance)
        self._check_identical_attrs(new_instance, instance, info)
        self._check_equal_attrs(new_instance, instance, info)
        assert getattr(new_instance, new_attr) != getattr(instance, new_attr)
        assert getattr(new_instance, skip) != getattr(instance, skip)

    def test_invalid_attribute(self):
        """Check complaints when trying to replace non-existing attributes."""
        instance = SampleFrozenClass(1)
        with pytest.raises(TypeError):
            replace_values(instance, non_existent_attr=3)


class TestSetFrozenAttr:
    """Tests for the set_frozen_attr function."""

    def test_frozen(self):
        """Check correct setting of an attribute of a frozen dataclass."""
        instance = SampleFrozenClass(attr=1)
        set_frozen_attr(instance, 'attr', 10)
        # pylint: disable-next=magic-value-comparison  # OK for attr
        assert instance.attr == 10

    def test_not_frozen(self):
        """Check that setting a non-frozen attr works as expected."""
        @dataclass
        class _TestClass:
            attr : ...
            dependent_attr : ... = non_init_field()

            def __setattr__(self, attr, value):
                super().__setattr__(attr, value)
                super().__setattr__('dependent_attr', value)

        value = object()
        instance = _TestClass(1)
        set_frozen_attr(instance, 'attr', value)
        assert instance.attr is value
        assert instance.dependent_attr is value

    def test_raises_not_dataclass(self):
        """Check complaints when using set_frozen_attr on a non-dataclass."""
        with pytest.raises(TypeError):
            set_frozen_attr('not a dataclass', 'attr', 'value')

    def test_raises_attribute_error(self):
        """Check complaints when trying to set a non-existing attribute."""
        instance = SampleFrozenClass(attr=1)
        with pytest.raises(AttributeError):
            set_frozen_attr(instance, 'does_not_exist', 3)
