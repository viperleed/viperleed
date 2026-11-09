"""Module QtGui of viperleed.gui.qt_compat.

Exposes all the contents of the QtGui module of the Qt binding
selected by viperleed.gui.qt_compat.binding.get_binding.
In addition, it always exposes the classes that PyQt6 moved from
QtWidgets to QtGui, also when using PyQt5.
"""
# ruff: noqa: N999, F403
# Only one binding is installed, and pylint cannot introspect it.
# Hence, which of these suppressions is needed depends on the env.
# pylint: disable=useless-suppression
# pylint: disable=invalid-name,import-error,wildcard-import
# pylint: disable=unused-wildcard-import,unused-import

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-09'
__license__ = 'GPLv3+'

from viperleed.gui.qt_compat.binding import PYQT5 as _PYQT5
from viperleed.gui.qt_compat.binding import get_binding as _get_binding

_BINDING = _get_binding()

if _BINDING == _PYQT5:
    from PyQt5.QtGui import *

    # These classes moved from QtWidgets to QtGui in Qt6
    from PyQt5.QtWidgets import (  # noqa: F401
        QAction,
        QActionGroup,
        QFileSystemModel,
        QShortcut,
        QUndoCommand,
        QUndoGroup,
        QUndoStack,
        )
else:
    from PyQt6.QtGui import *
