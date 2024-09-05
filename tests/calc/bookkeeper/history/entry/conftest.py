"""Test configuration for tests/calc/bookkeeper/history/entry."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-08-30'
__license__ = 'GPLv3+'

from enum import Enum
import functools

from pytest_cases import fixture

from viperleed.calc.bookkeeper.history.errors import EntrySyntaxError
from viperleed.calc.bookkeeper.history.errors import FixableSyntaxError
from viperleed.calc.bookkeeper.history.entry.field import FieldBase
from viperleed.calc.lib.dataclass_utils import frozen

no_value = object()


class MockFieldTag(Enum):
    """A fake field tag, useful for tests only."""

    TAG_1 = '# TAG_1'
    TAG_2 = 'TAG_2:'
    TAG_3 = '~ TAG_3 --'


@fixture(name='make_field')
def fixture_make_field():
    """Return a FieldBase instance with a given value."""
    def _make(cls, value=no_value, **kwargs):
        return cls() if value is no_value else cls(value, **kwargs)
    return _make


@fixture(name='make_and_check_field')
def factory_make_and_check_field(make_field):
    """Return a FieldBase instance with a value after silently checking."""
    def _make(*args, **kwargs):
        field = make_field(*args, **kwargs)
        try:
            field.check_value()
        except (EntrySyntaxError, FixableSyntaxError):
            pass
        return field
    return _make


@fixture(name='make_concrete_field')
def factory_make_concrete_field(monkeypatch, request):
    """Return a concrete subclass of an abstract field."""
    def _make(abstract, tag):
        if abstract.tag is tag:
            return abstract

        def cleanup():
            # pylint: disable-next=protected-access       # OK in tests
            del FieldBase._subclasses[tag]

        request.addfinalizer(cleanup)
        if isinstance(tag, MockFieldTag):
            monkeypatch.setattr(
                'viperleed.calc.bookkeeper.history.entry.field.FieldTag',
                MockFieldTag
                )
        return frozen(type('_Concrete', (abstract,), {}, tag=tag))
    return _make


@fixture(name='make_field_factory')
def factory_make_field_factory(make_and_check_field):
    """Return a factory that produces instances of a field."""
    def _make(field_cls):
        return functools.partial(make_and_check_field, field_cls)
    return _make


@fixture
def make_concrete_field_instance(make_concrete_field, make_field_factory):
    """Return a factory that produces instances of a concrete field."""
    def _make(base_cls, tag=None):
        if tag is None:
            tag = base_cls.tag or MockFieldTag.TAG_1
        field_cls = make_concrete_field(base_cls, tag=tag)
        return make_field_factory(field_cls)
    return _make


@fixture
def remove_field_tag(monkeypatch):
    """Temporarily remove a tag from FieldBase."""
    def _remove(tag, raising=False):
        # pylint: disable-next=protected-access           # OK in tests
        monkeypatch.delitem(FieldBase._subclasses, tag, raising=raising)
    return _remove
