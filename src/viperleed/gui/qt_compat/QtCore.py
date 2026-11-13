"""Module QtCore of viperleed.gui.qt_compat.

Exposes all the contents of the QtCore module of the Qt binding
selected by viperleed.gui.qt_compat.binding.get_binding.
"""
# ruff: noqa: N999, F403
# Only one binding is installed, and pylint cannot introspect it.
# Hence, which of these suppressions is needed depends on the env.
# pylint: disable=useless-suppression
# pylint: disable=invalid-name,import-error,wildcard-import
# pylint: disable=unused-wildcard-import

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
    from PyQt5.QtCore import *
else:
    from PyQt6.QtCore import *
