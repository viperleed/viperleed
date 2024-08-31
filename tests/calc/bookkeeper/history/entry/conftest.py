"""Test configuration for tests/calc/bookkeeper/history/entry."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-08-30'
__license__ = 'GPLv3+'

from pytest_cases import fixture

from viperleed.calc.bookkeeper.history.errors import EntrySyntaxError
from viperleed.calc.bookkeeper.history.errors import FixableSyntaxError
from viperleed.calc.bookkeeper.history.entry.field import FieldBase
from viperleed.calc.lib.dataclass_utils import frozen

no_value = object()

@fixture(name='make_field')
def fixture_make_field():
    """Return a FieldBase instance with a given value."""
    def _make(cls, value=no_value):
        return cls() if value is no_value else cls(value)
    return _make


@fixture
def make_and_check_field(make_field):
    """Return a FieldBase instance with a value after silently checking."""
    def _make(*args, **kwargs):
        field = make_field(*args, **kwargs)
        try:
            field.check_value()
        except (EntrySyntaxError, FixableSyntaxError):
            pass
        return field
    return _make


@fixture
def make_concrete_field(request):
    """Return a concrete subclass of an abstract field."""
    def _make(abstract, tag):
        if abstract.tag is tag:
            return abstract

        def cleanup():
            # pylint: disable-next=protected-access       # OK in tests
            del FieldBase._subclasses[tag]

        request.addfinalizer(cleanup)
        return frozen(type('_Concrete', (abstract,), {}, tag=tag))
    return _make


@fixture
def remove_field_tag(monkeypatch):
    """Temporarily remove a tag from FieldBase."""
    def _remove(tag, raising=False):
        # pylint: disable-next=protected-access           # OK in tests
        monkeypatch.delitem(FieldBase._subclasses, tag, raising=raising)
    return _remove
