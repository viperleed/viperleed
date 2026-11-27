"""Fixtures for the tests of viperleed.gui.

Modules that import Qt are imported via the fixtures defined here,
rather than at the top of test modules, so that test modules can be
collected without PyQt5. Tests that use these fixtures are skipped
(or fail with --require-all) if PyQt5 is not installed. They fail if
PyQt5 is installed, but broken.

Fixtures
--------
import_with_qt (factory)
    Import a module that needs Qt, or skip the test without PyQt5.
no_blocking_dialogs
    Make modal pop-ups return immediately, without showing.
qtc, qtg, qtw
    The QtCore, QtGui and QtWidgets modules of PyQt5.
show (factory)
    Show a widget, and wait until it is visible.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-16'
__license__ = 'GPLv3+'

import importlib

from pytest_cases import fixture

from ..harness import qt


@fixture(name='import_with_qt')
def fixture_import_with_qt(missing_dependency):
    """Return a function that imports a module that needs Qt.

    The returned function takes the name of a module (e.g.,
    'viperleed.gui.widgets.spinboxes') and returns the module. It
    skips the test (or fails with --require-all) if PyQt5 is not
    installed. Any other ImportError is raised.
    """
    def _import(module_name):
        try:
            return importlib.import_module(module_name)
        except ImportError as exc:
            if not qt.is_not_installed(exc):  # E.g., PyQt5 is broken
                raise
            missing_dependency(qt.DEPENDENCY, str(exc))
            raise  # Unreachable, as missing_dependency raises
    return _import


@fixture(name='no_blocking_dialogs')
def fixture_no_blocking_dialogs(qtw, mocker):
    """Make modal pop-ups return immediately, without showing.

    QMessageBox.question answers Yes (the return value can be changed
    via the returned mock). File dialogs are cancelled, and so are
    other dialogs run via QDialog.exec (e.g., opened by a button that
    a test clicks by mistake), which would otherwise wait for a user.
    """
    for method in ('getOpenFileName', 'getOpenFileNames', 'getSaveFileName'):
        mocker.patch.object(qtw.QFileDialog, method, return_value=('', ''))
    for method in ('exec', 'exec_'):
        mocker.patch.object(qtw.QDialog, method,
                            return_value=qtw.QDialog.DialogCode.Rejected)
    return mocker.patch.object(
        qtw.QMessageBox, 'question',
        return_value=qtw.QMessageBox.StandardButton.Yes,
        )


@fixture(name='qtc')
def fixture_qtc(import_with_qt):
    """Return the QtCore module of PyQt5."""
    return import_with_qt(f'{qt.BINDING}.QtCore')


@fixture(name='qtg')
def fixture_qtg(import_with_qt):
    """Return the QtGui module of PyQt5."""
    return import_with_qt(f'{qt.BINDING}.QtGui')


@fixture(name='qtw')
def fixture_qtw(import_with_qt):
    """Return the QtWidgets module of PyQt5."""
    return import_with_qt(f'{qt.BINDING}.QtWidgets')


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
