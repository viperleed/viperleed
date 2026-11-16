"""Fixtures for the tests of viperleed.gui.

Modules that import Qt are imported via the fixtures defined here,
rather than at the top of test modules, so that test modules can be
collected without a Qt binding. Tests that use these fixtures are
skipped (or fail with --require-all) if no Qt binding is available.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-16'
__license__ = 'GPLv3+'

import importlib

from pytest_cases import fixture

from viperleed.gui.qt_compat.binding import QtBindingNotFoundError

_QT_COMPAT = 'viperleed.gui.qt_compat'


@fixture(name='import_with_qt')
def fixture_import_with_qt(missing_dependency):
    """Return a function that imports a module that needs Qt.

    The returned function takes the name of a module (e.g.,
    'viperleed.gui.qt_compat.QtCore') and returns the module. It
    skips the test (or fails with --require-all) if no Qt binding
    is available.
    """
    def _import(module_name):
        try:
            return importlib.import_module(module_name)
        except QtBindingNotFoundError as exc:
            missing_dependency('PyQt5/PyQt6', str(exc))
            raise  # Unreachable, as missing_dependency raises
    return _import


@fixture(name='qtc')
def fixture_qtc(import_with_qt):
    """Return the QtCore module of qt_compat."""
    return import_with_qt(f'{_QT_COMPAT}.QtCore')


@fixture(name='qtg')
def fixture_qtg(import_with_qt):
    """Return the QtGui module of qt_compat."""
    return import_with_qt(f'{_QT_COMPAT}.QtGui')


@fixture(name='qtw')
def fixture_qtw(import_with_qt):
    """Return the QtWidgets module of qt_compat."""
    return import_with_qt(f'{_QT_COMPAT}.QtWidgets')
