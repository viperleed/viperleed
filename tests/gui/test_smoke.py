"""Smoke tests for the GUI: import modules, open windows and widgets.

These tests do not check what the GUI shows. They make sure that
modules can be imported, and that widgets (including windows and
dialogs) can be created, shown, used along their main paths, and
closed with the installed Qt binding (PyQt5 or PyQt6) without raising.
Exceptions raised in Qt slots and virtual methods (e.g., paintEvent)
make the tests fail, as pytest-qt captures them.

What is tested is collected dynamically, rather than listed: all the
modules of viperleed.gui, all the QWidget subclasses that they define,
and all the plug-ins that the plug-in selector can open. Widgets that
need constructor arguments are skipped, unless their arguments are
given in _INIT_ARGS.

All tests need pytest-qt and a Qt binding: they are skipped (or fail
with --require-all) otherwise, via the qapp and qtbot fixtures (see
tests/qt_plugin.py). Modules are imported at collection only if a Qt
binding is available, so that this file can be collected without one.
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
from pytest_cases import fixture
from pytest_cases import parametrize

import viperleed.gui
from viperleed.gui.qt_compat.binding import PYQT6
from viperleed.gui.qt_compat.binding import QtBindingNotFoundError
from viperleed.gui.qt_compat.binding import get_binding

try:
    from viperleed.gui.qt_compat import QtCore as qtc
    from viperleed.gui.qt_compat import QtGui as qtg
    from viperleed.gui.qt_compat import QtWidgets as qtw
except QtBindingNotFoundError:
    pass  # All tests are skipped. See module docstring.

_ERROR_BOX = 'viperleed.gui.leedsim.dialogs.errorbox:ErrorBox'
_MAIN_WINDOW = 'viperleed.gui.leedsim.mainwindow'
_SELECTOR = 'viperleed.gui.selectplugin:ViPErLEEDSelectPlugin'
_NOT_IMPORTABLE = (
    'viperleed.gui.__main__',  # Runs the GUI
    )
# Message of the ImportError raised by mpl_graphics when PyQt5 is
# missing. TODO: remove once detect_graphics and mpl_graphics support
# PyQt6, together with NEEDS_MPL_CANVAS.
_MPL_CANVAS_NEEDS_PYQT5 = 'PyQt5 not found.'
NEEDS_MPL_CANVAS = pytest.mark.xfail(
    reason='matplotlib canvases need PyQt5 (detect_graphics.has_pyqt)',
    raises=ImportError,
    )
# A parameter for when nothing can be collected (e.g., no Qt binding):
# the test is skipped via qtbot, or fails importing what is missing.
_PLACEHOLDER = 'nothing-collected'


def _get_binding():
    """Return the installed Qt binding, or None."""
    try:
        return get_binding()
    except QtBindingNotFoundError:
        return None


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
        if no Qt binding is available.
    """
    names = (module.name
             for module in pkgutil.walk_packages(viperleed.gui.__path__,
                                                 prefix='viperleed.gui.')
             if module.name not in _NOT_IMPORTABLE)
    if _get_binding() is None:
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


def _needs_mpl_canvas(dotted_name):
    """Return whether the module of `dotted_name` misses PyQt5."""
    exc = _GUI_MODULES.get(dotted_name.partition(':')[0])
    return (_get_binding() == PYQT6
            and isinstance(exc, ImportError)
            and str(exc) == _MPL_CANVAS_NEEDS_PYQT5)


def _xfail_if_needs_mpl_canvas(dotted_name):
    """Return a decorator applying NEEDS_MPL_CANVAS if appropriate."""
    if _needs_mpl_canvas(dotted_name):
        return NEEDS_MPL_CANVAS
    return lambda decorated: decorated


def _param(value, dotted_name=None):
    """Return a pytest.param for `value`, marked if appropriate."""
    needs_mpl_canvas = _needs_mpl_canvas(dotted_name or value)
    marks = (NEEDS_MPL_CANVAS,) if needs_mpl_canvas else ()
    return pytest.param(value, marks=marks)


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


def _collect_plugins():
    """Return the names of the plug-ins that the selector can open."""
    selector_module = _GUI_MODULES.get(_SELECTOR.partition(':')[0])
    if not isinstance(selector_module, ModuleType):
        # No Qt (skipped), or the import fails (reported by test_import)
        return (_param(_PLACEHOLDER, _SELECTOR),)
    plugins = selector_module.ViPErLEEDSelectPlugin.modules
    return tuple(name for name, (_, cls) in plugins.items() if cls)


_PLUGINS = _collect_plugins()


def _default_parameters():
    """Return the LEED parameters used when starting from scratch."""
    params = _import('viperleed.gui.leedsim.classes.leedparameters:'
                     'LEEDParametersList')
    return params(_import('viperleed.gui.leedsim.dialogs.newfiledialog:'
                          'DEFAULT_STARTUP'))


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
    _ERROR_BOX: lambda cls: cls(error_while='testing'),
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

# Bugs found by the tests. Keys are 'module:Class', values the
# exception raised and the reason. Marked as strict xfail.
_KNOWN_BUGS = {
    'viperleed.gui.widgets.spinboxes:TextBox': (
        UnboundLocalError,
        'TextBox() without arguments uses contents before assignment',
        ),
    }


def _widget_param(cls_name):
    """Return a pytest.param for widget `cls_name`, marked if needed."""
    try:
        exc, reason = _KNOWN_BUGS[cls_name]
    except KeyError:
        return pytest.param(cls_name)
    return pytest.param(cls_name,
                        marks=pytest.mark.xfail(raises=exc, reason=reason))


_WIDGETS = (tuple(_widget_param(name)
                  for name in sorted(_iter_widget_classes()))
            or (_param(_PLACEHOLDER),))


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
    if (headless
            and _get_binding() != PYQT6
            and isinstance(widget, qtw.QMessageBox)):
        pytest.skip('Qt 5 crashes when showing a QMessageBox on a '
                    'headless platform')


@fixture(name='no_blocking_dialogs', autouse=True)
def fixture_no_blocking_dialogs(qapp, mocker):
    """Make modal pop-ups return immediately, without showing.

    QMessageBox.question answers Yes (the return value can be changed
    via the returned mock). File dialogs are cancelled.
    """
    _ = qapp  # Skips without Qt, before using qtw below
    for method in ('getOpenFileName', 'getOpenFileNames', 'getSaveFileName'):
        mocker.patch.object(qtw.QFileDialog, method, return_value=('', ''))
    return mocker.patch.object(
        qtw.QMessageBox, 'question',
        return_value=qtw.QMessageBox.StandardButton.Yes,
        )


@fixture(name='default_parameters')
def fixture_default_parameters():
    """Return the LEED parameters used when starting from scratch."""
    return _default_parameters()


def _is_deleted_on_close(widget):
    """Return whether closing `widget` deletes it."""
    return widget.testAttribute(qtc.Qt.WidgetAttribute.WA_DeleteOnClose)


def _show(qtbot, widget):
    """Show `widget`, and wait until it is visible.

    Widgets deleted when closed (e.g., ViPErLEED plug-ins) are not
    registered to qtbot, which would close them again at teardown.
    Close them in the tests.
    """
    if not _is_deleted_on_close(widget):
        qtbot.addWidget(widget)
    with qtbot.waitExposed(widget):
        widget.show()


@parametrize(module=tuple(_param(name) for name in _GUI_MODULES))
def test_import(module):
    """Check that a GUI module can be imported."""
    importlib.import_module(module)


@parametrize(cls_name=_WIDGETS)
def test_open_and_close(cls_name, qtbot):
    """Check that a widget can be created, shown, and closed."""
    widget = _create_widget(cls_name)
    _skip_if_crashes_qt(widget)
    _show(qtbot, widget)
    assert widget.close()
    assert not widget.isVisible()


@parametrize(plugin=_PLUGINS)
def test_open_plugin_from_selector(plugin, qtbot):
    """Check that the plug-in selector opens and closes a plug-in."""
    selector = _import(_SELECTOR)()
    _show(qtbot, selector)
    # pylint: disable=protected-access
    qtbot.mouseClick(selector._btns[plugin], qtc.Qt.MouseButton.LeftButton)
    module = selector._open_modules[plugin]
    assert module is not None
    qtbot.waitExposed(module)
    assert selector.close()
    assert selector._open_modules[plugin] is None


def _raise_value_error():
    """Raise a ValueError."""
    raise ValueError(_ERROR_MSG)


_ERROR_MSG = 'a test error'
_ERROR_TEXT = 'Some error'
_EXTRA_TEXT = ' Extra text.'


def test_error_box(qtbot, mocker):
    """Check that an ErrorBox is shown only when handling exceptions."""
    error_box = _import(_ERROR_BOX)
    box = error_box(error_while='testing', text=_ERROR_TEXT)
    qtbot.addWidget(box)
    show_modal = mocker.patch(  # Would block the test
        'viperleed.gui.qt_compat.QtWidgets.QMessageBox.exec',
        return_value=box.StandardButton.Ok,
        )
    assert box.exec() == box.StandardButton.Ok
    show_modal.assert_not_called()
    try:
        _raise_value_error()
    except ValueError:
        assert box.exec(_EXTRA_TEXT) == box.StandardButton.Ok
    show_modal.assert_called_once()
    assert box.text() == _ERROR_TEXT + _EXTRA_TEXT
    assert _ERROR_MSG in box.detailedText()


@_xfail_if_needs_mpl_canvas(_MAIN_WINDOW)
class TestLEEDPatternSimulator:
    """Tests for the main paths of the (new) LEED-pattern simulator."""

    @fixture(name='simulator')
    def fixture_simulator(self, qtbot, no_blocking_dialogs):
        """Yield a shown LEEDPatternSimulator, then close it."""
        cls = _import(f'{_MAIN_WINDOW}:LEEDPatternSimulator')
        simulator = cls()
        _show(qtbot, simulator)
        yield simulator
        # Do not save changes when closing
        no_blocking_dialogs.return_value = qtw.QMessageBox.StandardButton.No
        assert simulator.close()

    def test_new_pattern(self, simulator, default_parameters, qtbot):
        """Check that a LEED pattern can be created from parameters."""
        # pylint: disable=protected-access
        dialog = simulator._dialogs['file_new']
        with qtbot.waitSignal(dialog.editing_finished):
            dialog.editing_finished.emit(default_parameters)
        assert not simulator.saved

    def test_new_file_dialog_done(self, simulator, qtbot):
        """Check that the new-file dialog can be opened and accepted."""
        # pylint: disable=protected-access
        simulator._on_file_new_or_edit_pressed()
        dialog = simulator._dialogs['file_new']
        qtbot.waitExposed(dialog)
        with qtbot.waitSignal(dialog.editing_finished):
            qtbot.mouseClick(dialog._ctrls['done'],
                             qtc.Qt.MouseButton.LeftButton)
        assert not dialog.isVisible()


def test_drag_matrices_popup(qtbot):
    """Check that the pop-up with the matrices can be dragged."""
    popup = _create_widget(
        'viperleed.gui.leedsim.widgets.matricespopup:MatricesPopup'
        )
    _show(qtbot, popup)
    start = popup.pos()
    qtbot.mousePress(popup, qtc.Qt.MouseButton.LeftButton,
                     pos=qtc.QPoint(5, 5))
    # QTest.mouseMove does not report pressed buttons: send a move
    # event with the left button pressed via the event handler itself
    popup.mouseMoveEvent(_left_button_move_event(popup, qtc.QPoint(25, 35)))
    assert popup.pos() - start == qtc.QPoint(20, 30)


def _left_button_move_event(widget, local_pos):
    """Return a mouse-move event at `local_pos` with the left button."""
    left = qtc.Qt.MouseButton.LeftButton
    return qtg.QMouseEvent(qtc.QEvent.Type.MouseMove,
                           qtc.QPointF(local_pos),
                           qtc.QPointF(widget.mapToGlobal(local_pos)),
                           left, left, qtc.Qt.KeyboardModifier.NoModifier)


def test_wheel_on_text_box(qtbot):
    """Check that the mouse wheel changes the value of a TextBox."""
    text_box = _import('viperleed.gui.widgets.spinboxes:TextBox')('1.0')
    text_box.setStep(1, 'add')
    _show(qtbot, text_box)
    text_box.setFocus()
    qtbot.waitUntil(text_box.hasFocus)
    center = qtc.QPointF(text_box.rect().center())
    wheel = qtg.QWheelEvent(center, qtc.QPointF(text_box.mapToGlobal(
                                center.toPoint())),
                            qtc.QPoint(0, 0),
                            qtc.QPoint(0, 160),  # One step with a mouse
                            qtc.Qt.MouseButton.NoButton,
                            qtc.Qt.KeyboardModifier.NoModifier,
                            qtc.Qt.ScrollPhase.NoScrollPhase,
                            False)
    qtw.QApplication.sendEvent(text_box, wheel)
    assert float(text_box.text()) == pytest.approx(2.0)
