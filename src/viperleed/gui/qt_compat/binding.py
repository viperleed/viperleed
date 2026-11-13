"""Module binding of viperleed.gui.qt_compat.

Defines functionality for selecting which Python binding for Qt
(PyQt5 or PyQt6) is used by the graphical user interface.

The binding is selected once, the first time get_binding is called,
according to these rules:
- If the QT_API environment variable is set, its value decides.
  Accepted values are 'pyqt5' and 'pyqt6' (case insensitive). This
  is the same variable that matplotlib uses, which ensures that both
  use the same binding. An error is raised if the requested binding
  cannot be imported, or if QT_API has any other value.
- Otherwise, if one of the bindings was already imported (e.g., by
  user code or by another package), that one is used. Two different
  Qt bindings cannot coexist in the same process.
- Otherwise, the first binding that can be imported among PyQt6 and
  PyQt5 (in this order) is used.

Which bindings would be tried, in which order, is returned by
get_candidates, which imports none of them. This allows checking the
candidates in a subprocess first (see viperleed.gui.detect_graphics),
as a broken Qt installation may crash the interpreter when imported.
The binding found to work can then be selected with set_binding.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-09'
__license__ = 'GPLv3+'

import importlib
import os
import sys

PYQT5 = 'PyQt5'
PYQT6 = 'PyQt6'
QT_API_ENV = 'QT_API'
SUPPORTED_BINDINGS = (PYQT6, PYQT5)  # In order of preference
_QT_API_VALUES = {binding.lower(): binding for binding in SUPPORTED_BINDINGS}

_SELECTED_BINDING = None


def get_binding():
    """Return the name of the Qt binding in use.

    A binding is selected the first time this function is called.
    See the module docstring for the selection rules.

    Returns
    -------
    binding : str
        Either PYQT5 or PYQT6.

    Raises
    ------
    QtBindingNotFoundError
        If no usable binding is found, or if the QT_API environment
        variable has an unsupported value.
    """
    global _SELECTED_BINDING  # pylint: disable=global-statement
    if _SELECTED_BINDING is None:
        _SELECTED_BINDING = _select_binding()
    return _SELECTED_BINDING


def get_candidates():
    """Return the bindings to be tried, in order, without importing.

    Returns
    -------
    candidates : tuple
        The names of the bindings that get_binding would try to
        import, in the same order. A single binding if requested
        via QT_API or already imported.

    Raises
    ------
    QtBindingNotFoundError
        If the QT_API environment variable has an unsupported value.
    """
    requested = _get_requested_binding()
    if requested:
        return (requested,)
    already_imported = tuple(binding for binding in SUPPORTED_BINDINGS
                             if f'{binding}.QtCore' in sys.modules)
    return already_imported[:1] or SUPPORTED_BINDINGS


def set_binding(binding):
    """Select `binding`, without importing it.

    Use this only for a binding that is known to work, e.g., one that
    was successfully imported in a subprocess.

    Parameters
    ----------
    binding : str
        Either PYQT5 or PYQT6.

    Raises
    ------
    ValueError
        If `binding` is not supported, or if a different binding
        was already selected.
    """
    global _SELECTED_BINDING  # pylint: disable=global-statement
    if binding not in SUPPORTED_BINDINGS:
        raise ValueError(f'Unsupported Qt binding {binding!r}')
    if _SELECTED_BINDING not in (None, binding):
        raise ValueError(f'Cannot select {binding}: {_SELECTED_BINDING} '
                         'is already in use')
    _SELECTED_BINDING = binding


class QtBindingNotFoundError(ImportError):
    """No usable Qt binding could be found."""


def _get_requested_binding():
    """Return the binding requested via QT_API, or None."""
    requested = os.environ.get(QT_API_ENV, '').strip()
    if not requested:
        return None
    try:
        return _QT_API_VALUES[requested.lower()]
    except KeyError:
        supported = ', '.join(repr(v) for v in _QT_API_VALUES)
        raise QtBindingNotFoundError(
            f'Unsupported value {requested!r} for the {QT_API_ENV} '
            f'environment variable. ViPErLEED supports: {supported}.'
            ) from None


def _select_binding():
    """Return the name of the binding that should be used."""
    requested = _get_requested_binding()
    candidates = get_candidates()
    errors = []
    for binding in candidates:
        try:
            importlib.import_module(f'{binding}.QtCore')
        except ImportError as exc:
            errors.append(f'{binding}: {exc}')
            continue
        return binding
    msg = (f'Could not import {requested}, requested via '
           f'{QT_API_ENV}={os.environ[QT_API_ENV]!r}. ' if requested
           else f'Could not import any Qt binding among {candidates}. ')
    raise QtBindingNotFoundError(msg + 'Errors: ' + '; '.join(errors))
