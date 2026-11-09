"""Tests for the Qt modules exposed by viperleed.gui.qt_compat.

These tests use the Qt binding that is actually installed. They
are skipped if there is none.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-09'
__license__ = 'GPLv3+'

import importlib

import pytest
from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.gui.qt_compat.binding import PYQT5
from viperleed.gui.qt_compat.binding import QtBindingNotFoundError
from viperleed.gui.qt_compat.binding import get_binding

_PACKAGE = 'viperleed.gui.qt_compat'
_MOVED_TO_QTGUI = (
    'QAction',
    'QActionGroup',
    'QFileSystemModel',
    'QShortcut',
    'QUndoCommand',
    'QUndoGroup',
    'QUndoStack',
    )


@fixture(name='qt_binding')
def fixture_qt_binding(missing_dependency):
    """Return the name of the installed Qt binding, or skip."""
    try:
        return get_binding()
    except QtBindingNotFoundError as exc:
        missing_dependency('PyQt5/PyQt6', str(exc))
    raise AssertionError('Unreachable')


def _import_both(module, binding):
    """Return the `module` from qt_compat and from `binding`."""
    return (importlib.import_module(f'{_PACKAGE}.{module}'),
            importlib.import_module(f'{binding}.{module}'))


_REPRESENTATIVE = {
    'QtCore': 'QObject',
    'QtGui': 'QColor',
    'QtNetwork': 'QNetworkRequest',
    'QtSerialPort': 'QSerialPort',
    'QtWidgets': 'QWidget',
    }


@parametrize('module,name', _REPRESENTATIVE.items(), ids=_REPRESENTATIVE)
def test_module_contents(module, name, qt_binding):
    """Check that a qt_compat module exposes the binding's objects."""
    compat, original = _import_both(module, qt_binding)
    assert getattr(compat, name) is getattr(original, name)


@parametrize(name=_MOVED_TO_QTGUI)
def test_moved_to_qtgui(name, qt_binding):
    """Check that QtGui exposes classes that were in QtWidgets."""
    origin = 'QtWidgets' if qt_binding == PYQT5 else 'QtGui'
    qtg = importlib.import_module(f'{_PACKAGE}.QtGui')
    original = importlib.import_module(f'{qt_binding}.{origin}')
    assert getattr(qtg, name) is getattr(original, name)


@pytest.mark.usefixtures('qt_binding')
def test_scoped_enums():
    """Check that fully scoped enumeration members are available."""
    qtc = importlib.import_module(f'{_PACKAGE}.QtCore')
    flags = qtc.Qt.AlignmentFlag.AlignLeft | qtc.Qt.AlignmentFlag.AlignTop
    assert flags & qtc.Qt.AlignmentFlag.AlignTop
