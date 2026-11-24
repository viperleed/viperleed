"""Smoke tests for the GUI: import modules, open windows and widgets.

These tests do not check what the GUI shows. They make sure that
modules can be imported, and that widgets (including windows and
dialogs) can be created, shown, and closed without raising.
Exceptions raised in Qt slots and virtual methods (e.g., paintEvent)
make the tests fail, as pytest-qt captures them.

What is tested is collected dynamically, rather than listed: all the
modules of viperleed.gui, and all the QWidget subclasses that they
define. Widgets that need constructor arguments are skipped, unless
a factory is given in _FACTORIES. Tests of specific widgets belong
to the test modules of the modules that define them. Widgets with
known bugs that make these tests fail are listed in _KNOWN_BUGS, slow
ones in _SLOW_WIDGETS.

All tests need pytest-qt and PyQt5: they are skipped (or fail with
--require-all) otherwise (see tests/harness/qt.py and the
import_with_qt fixture in tests/gui/conftest.py). Modules are imported
at collection only if PyQt5 is installed, so that this file can be
collected without it.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-10'
__license__ = 'GPLv3+'

import importlib
import inspect
import pkgutil
from types import ModuleType

import numpy as np
import pytest
from pytest_cases import parametrize

import viperleed.gui

from ..harness import qt
from . import known_bugs

try:
    from PyQt5 import QtWidgets as qtw
except ImportError as exc:
    if not qt.is_not_installed(exc):  # Broken: fail at collection
        raise
    _HAS_QT = False  # All tests are skipped. See module docstring.
else:
    _HAS_QT = True

_NOT_IMPORTABLE = (
    'viperleed.gui.__main__',  # Runs the GUI
    )
# A parameter for when nothing can be collected (e.g., no PyQt5): the
# test is skipped via qtbot, or fails importing what is missing.
_PLACEHOLDER = 'nothing-collected'


def _import(dotted_name):
    """Return the object `dotted_name` (e.g., 'pkg.module:Class')."""
    module_name, _, obj_name = dotted_name.partition(':')
    module = importlib.import_module(module_name)
    return getattr(module, obj_name) if obj_name else module


def _import_all_gui_modules():
    """Return a dict of all the modules of viperleed.gui.

    Returns
    -------
    modules : dict
        Keys are module names. Values are the modules, the exceptions
        raised when importing them (reported by test_import), or None
        if PyQt5 is not installed.
    """
    names = (module.name
             for module in pkgutil.walk_packages(viperleed.gui.__path__,
                                                 prefix='viperleed.gui.')
             if module.name not in _NOT_IMPORTABLE)
    if not _HAS_QT:
        return dict.fromkeys(names)
    modules = {}
    for name in names:
        # Any exception is reported by test_import
        try:
            modules[name] = importlib.import_module(name)
        # pylint: disable-next=broad-exception-caught
        except Exception as exc:  # noqa: BLE001, PERF203
            modules[name] = exc
    return modules


_GUI_MODULES = _import_all_gui_modules()


def _iter_widget_classes():
    """Yield 'module:Class' for QWidget subclasses in viperleed.gui."""
    for module_name, module in _GUI_MODULES.items():
        if not isinstance(module, ModuleType):
            continue
        for obj in vars(module).values():
            if (isinstance(obj, type)
                    and issubclass(obj, qtw.QWidget)
                    and obj.__module__ == module_name):
                yield f'{module_name}:{obj.__qualname__}'


def _make_export_csv_dialog(cls):
    """Return an ExportCSVDialog for the default LEED pattern."""
    # ExportCSVDialog is written for the old LEEDPattern, used by
    # mainwindow_old. The new main window passes a new LEEDPattern.
    pattern = _import('viperleed.gui.leedsim.classes.oldleedpatterns:'
                      'LEEDPattern')
    default = _import('viperleed.gui.leedsim.dialogs.newfiledialog:'
                      'DEFAULT_STARTUP')
    return cls(pattern(default))


def _make_hover_annotation(cls):
    """Return a HoverAnnot with an arrow, as it cannot paint without."""
    annotation = cls()
    annotation.center = (50, 50)
    annotation.head_pos = (60, 60)
    return annotation


def _make_in_parent(cls):
    """Return a widget containing a `cls` (which needs a parent)."""
    parent = qtw.QWidget()
    layout = qtw.QVBoxLayout(parent)
    layout.addWidget(cls())
    return parent


def _make_surface_input(cls):
    """Return a SurfaceStructureInput for a square bulk lattice."""
    lattice = _import('viperleed.gui.classes.lattice2d:Lattice2D')
    return cls(bulk_lattice=lattice(np.eye(2)))


# How to create widgets that need constructor arguments or some set-up
# before they can be shown. Keys are 'module:Class', values callables
# that take the class and return a widget. Widgets that need arguments
# but are not listed here are skipped.
_FACTORIES = {
    'viperleed.gui.leedsim.dialogs.errorbox:ErrorBox':
        lambda cls: cls(error_while='testing'),
    'viperleed.gui.leedsim.dialogs.exportcsvdialog:ExportCSVDialog':
        _make_export_csv_dialog,
    'viperleed.gui.leedsim.widgets.hoverannot:HoverAnnot':
        _make_hover_annotation,
    'viperleed.gui.leedsim.widgets.leedcanvas:LEEDCanvas': _make_in_parent,
    'viperleed.gui.leedsim.widgets.matricespopup:MatricesPopup':
        lambda cls: cls([np.eye(2)], None),
    'viperleed.gui.leedsim.widgets.surfaceinput:SurfaceStructureInput':
        _make_surface_input,
    'viperleed.gui.widgets.basewidgets:PainterMatrix':
        lambda cls: cls(np.eye(2)),
    }

# Widgets with bugs that make test_open_and_close fail. Keys are
# 'module:Class', values the xfail marks of tests.gui.known_bugs.
_KNOWN_BUGS = {
    'viperleed.gui.widgets.spinboxes:TextBox':
        known_bugs.TEXTBOX_WITHOUT_ARGUMENTS,
    }


# Widgets whose test_open_and_close takes about SLOW_THRESHOLD or
# more (see tests/harness/tiers.py). They are marked here rather than
# in SLOW_CASES, as widgets are collected only if Qt can be imported:
# without Qt, entries in SLOW_CASES would be reported as stale.
_SLOW_WIDGETS = (
    'viperleed.gui.leedsim.widgets.matricespopup:MatricesPopup',
    )


def _widget_param(cls_name):
    """Return a pytest.param for widget `cls_name`, marked if needed."""
    marks = []
    if cls_name in _KNOWN_BUGS:
        marks.append(_KNOWN_BUGS[cls_name])
    if cls_name in _SLOW_WIDGETS:
        marks.append(pytest.mark.slow)
    return pytest.param(cls_name, marks=marks)


_WIDGETS = (tuple(_widget_param(name)
                  for name in sorted(_iter_widget_classes()))
            or (_PLACEHOLDER,))


def _create_widget(cls_name):
    """Return a new widget of class `cls_name`, or skip the test."""
    cls = _import(cls_name)
    if cls_name in _FACTORIES:
        return _FACTORIES[cls_name](cls)
    try:
        signature = inspect.signature(cls.__init__)
    except ValueError:  # Constructor of Qt: no mandatory arguments
        return cls()
    try:
        signature.bind(None)  # For self
    except TypeError:
        pytest.skip(f'{cls_name} needs constructor arguments. '
                    'Add a factory to _FACTORIES.')
    return cls()


def _skip_if_crashes_qt(widget):
    """Skip if showing `widget` would crash the Python interpreter.

    With Qt 5, showing any QMessageBox on a headless platform (at
    least 'offscreen' and 'minimal' on Windows) is an access violation.
    """
    headless = qtw.QApplication.platformName() in {'offscreen', 'minimal'}
    if headless and isinstance(widget, qtw.QMessageBox):
        pytest.skip('Qt 5 crashes when showing a QMessageBox on a '
                    'headless platform')


@parametrize(module=tuple(_GUI_MODULES))
def test_import(module, import_with_qt):
    """Check that a GUI module can be imported."""
    import_with_qt(module)


@pytest.mark.usefixtures('no_blocking_dialogs')
@parametrize(cls_name=_WIDGETS)
def test_open_and_close(cls_name, show):
    """Check that a widget can be created, shown, and closed."""
    widget = _create_widget(cls_name)
    _skip_if_crashes_qt(widget)
    show(widget)
    assert widget.close()
    assert not widget.isVisible()
