"""Tests for module _base of viperleed.tleedmlib.classes.rparams.special.

Created on 2023-10-27

@author: Michele Riva (@michele-riva)
"""

from dataclasses import dataclass

import pytest
from pytest_cases import fixture

from viperleed.tleedmlib.classes.rparams.special._base import (
    SpecialParameter, NotASpecialParameterError
)


class TestSpecialParameter:
    """Tests for SpecialParameter base class."""

    @fixture(name='make_special_class')
    def fixture_make_special_class(self):
        """Return a subclass of SpecialParameter for `param`."""
        def _make(param):
            @dataclass
            class TestParam(SpecialParameter, param=param):
                """A subclass of SpecialParameter."""
                value : object
            return TestParam
        return _make

    def test_get_subclass_valid_param(self, make_special_class):
        """Check correctness of the result of get_subclass."""
        new_cls = make_special_class('TEST_PARAM')
        assert SpecialParameter.get_subclass('TEST_PARAM') is new_cls

    def test_get_subclass_invalid_param(self):
        """Check complaints when trying to access an unknown parameter."""
        with pytest.raises(NotASpecialParameterError):
            SpecialParameter.get_subclass('INVALID_PARAM')

    def test_registers_subclass(self, make_special_class):
        """Check correct registration of a subclass."""
        new_cls = make_special_class('test_param')
        # pylint: disable-next=protected-access
        assert new_cls._subclasses['TEST_PARAM'] is new_cls

    def test_from_value(self, make_special_class):
        """Check the correct return value of the from_value method."""
        new_cls = make_special_class('TEST_PARAM')
        test_value = 42
        instance = new_cls.from_value(test_value)
        assert isinstance(instance, new_cls)
        assert instance.value == test_value
