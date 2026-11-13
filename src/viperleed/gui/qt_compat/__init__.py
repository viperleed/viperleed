"""Package qt_compat of viperleed.gui.

======================================
  ViPErLEED Graphical User Interface
======================================

A thin compatibility layer that lets the graphical user interface run
with either PyQt5 or PyQt6. All Qt modules should be imported from
this package rather than from PyQt5/PyQt6 directly, e.g.,

    from viperleed.gui.qt_compat import QtCore as qtc
    from viperleed.gui.qt_compat import QtGui as qtg
    from viperleed.gui.qt_compat import QtWidgets as qtw

Code using this package must be written in a way that works with both
bindings. In particular:
- Use fully scoped enumeration members (e.g.,
  qtc.Qt.AlignmentFlag.AlignLeft rather than qtc.Qt.AlignLeft).
  PyQt5 supports them natively.
- Use .exec() rather than .exec_(). PyQt5 supports both.
- Take QAction, QActionGroup, QShortcut, QUndoCommand, QUndoGroup,
  QUndoStack, and QFileSystemModel from QtGui. They live in QtWidgets
  in PyQt5, but the QtGui module of this package also exposes them.
- Where no method works with both bindings, use the functions of the
  events module (e.g., get_global_position(mouse_event) instead of
  mouse_event.globalPos() or mouse_event.globalPosition()).

Which binding is used is decided once, the first time any of the
Qt modules is imported from this package. See the binding module.

Modules
-------
binding
    Selection of the Qt binding to be used. Does not import any Qt
    module unless a binding is requested.
events
    Functions for Qt events whose API differs between the bindings.
QtCore, QtGui, QtNetwork, QtSerialPort, QtWidgets
    The corresponding modules of the selected binding.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2025-11-10'
__license__ = 'GPLv3+'
