"""Fixtures for the tests of viperleed.gui.

Modules that import Qt are imported via the fixtures defined here,
rather than at the top of test modules, so that test modules can be
collected without a Qt binding. Tests that use these fixtures are
skipped (or fail with --require-all) if no Qt binding is available.
They fail if the QT_API environment variable prevents selecting one.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-16'
__license__ = 'GPLv3+'

import importlib

from pytest_cases import fixture

from viperleed.gui.qt_compat.binding import QtApiError
from viperleed.gui.qt_compat.binding import QtBindingNotFoundError

_QT_COMPAT = 'viperleed.gui.qt_compat'


@fixture(name='import_with_qt')
def fixture_import_with_qt(missing_dependency):
    """Return a function that imports a module that needs Qt.

    The returned function takes the name of a module (e.g.,
    'viperleed.gui.qt_compat.QtCore') and returns the module. It
    skips the test (or fails with --require-all) if no Qt binding
    is available. It fails if QT_API prevents selecting one.
    """
    def _import(module_name):
        try:
            return importlib.import_module(module_name)
        except QtApiError:  # QT_API is wrong: fail, rather than skip
            raise
        except QtBindingNotFoundError as exc:
            missing_dependency('PyQt', str(exc))
            raise  # Unreachable, as missing_dependency raises
    return _import


@fixture(name='no_blocking_dialogs')
def fixture_no_blocking_dialogs(qtw, mocker):
    """Make modal pop-ups return immediately, without showing.

    QMessageBox.question answers Yes (the return value can be changed
    via the returned mock). File dialogs are cancelled.
    """
    for method in ('getOpenFileName', 'getOpenFileNames', 'getSaveFileName'):
        mocker.patch.object(qtw.QFileDialog, method, return_value=('', ''))
    return mocker.patch.object(
        qtw.QMessageBox, 'question',
        return_value=qtw.QMessageBox.StandardButton.Yes,
        )


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


@fixture(name='show')
def fixture_show(qtbot, qtc):
    """Return a function that shows a widget, and waits until visible.

    Widgets deleted when closed (e.g., ViPErLEED plug-ins) are not
    registered to qtbot, which would close them again at teardown.
    Close them in the tests.
    """
    delete_on_close = qtc.Qt.WidgetAttribute.WA_DeleteOnClose

    def _show(widget):
        if not widget.testAttribute(delete_on_close):
            qtbot.addWidget(widget)
        with qtbot.waitExposed(widget):
            widget.show()
    return _show
