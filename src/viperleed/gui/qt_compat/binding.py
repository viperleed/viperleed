"""Module binding of viperleed.gui.qt_compat.

Defines functionality for selecting which Python binding for Qt
(PyQt5 or PyQt6) is used by the graphical user interface.

The binding is selected once, the first time get_binding is called,
according to these rules:
- If the QT_API environment variable is set, its value decides.
  Accepted values are 'pyqt5' and 'pyqt6' (case insensitive). This
  is the same variable that matplotlib uses, which ensures that both
  use the same binding. A QtApiError is raised if the requested
  binding is not installed, if QT_API has any other value, or if
  another binding was already imported.
- Otherwise, if one of the bindings was already imported (e.g., by
  user code or by another package), that one is used. Two different
  Qt bindings cannot coexist in the same process.
- Otherwise, the first binding that can be imported among PyQt6 and
  PyQt5 (in this order) is used. A binding that is installed but
  broken (e.g., it misses some of its own modules), or older than
  MIN_VERSIONS, is skipped too. If no binding can be used and any
  of them is broken or too old, a plain ImportError is raised rather
  than a QtBindingNotFoundError, as the installation needs fixing
  rather than completing.

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
# The oldest (major, minor) versions of the bindings that work. PyQt6
# before 6.1 lacks names of enums used in viperleed.gui (e.g., it has
# Qt.MouseButtons rather than Qt.MouseButton), and some enums that
# Qt functions return as int (e.g., Qt.Key) are not int-based there.
# pylint: disable-next=consider-using-namedtuple-or-dataclass
MIN_VERSIONS = {PYQT5: (5, 15), PYQT6: (6, 1)}
_QT_API_VALUES = {binding.lower(): binding for binding in SUPPORTED_BINDINGS}

_SELECTED_BINDING = None


class QtBindingNotFoundError(ImportError):
    """No usable Qt binding could be found."""


class QtApiError(QtBindingNotFoundError):
    """The QT_API environment variable prevents selecting a binding."""


class QtBindingTooOldError(ImportError):
    """A Qt binding is older than its version in MIN_VERSIONS."""


def check_version(binding, qt_core):
    """Raise if `binding` is older than its MIN_VERSIONS.

    Parameters
    ----------
    binding : str
        Either PYQT5 or PYQT6.
    qt_core : module
        The (imported) QtCore module of `binding`.

    Raises
    ------
    QtBindingTooOldError
        If the version of `binding` is older than MIN_VERSIONS.
    """
    version = qt_core.PYQT_VERSION  # E.g., 0x060c00 for 6.12.0
    major_minor = (version >> 16, (version >> 8) & 0xFF)
    oldest = MIN_VERSIONS[binding]
    if major_minor < oldest:
        oldest_str = '.'.join(str(v) for v in oldest)
        raise QtBindingTooOldError(
            f'{binding} {qt_core.PYQT_VERSION_STR} is too old. ViPErLEED '
            f'needs {binding}>={oldest_str}.'
            )


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
    QtApiError
        If the QT_API environment variable has an unsupported value,
        conflicts with an already imported binding, or requests a
        binding that is not installed.
    QtBindingNotFoundError
        If no binding is installed.
    ImportError
        If no binding can be used, and any of them is installed but
        broken or too old (see MIN_VERSIONS).
    """
    if _SELECTED_BINDING is None:
        set_binding(_select_binding())
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
    QtApiError
        If the QT_API environment variable has an unsupported value,
        or if it requests a binding other than one already imported.
    """
    requested = _get_requested_binding()
    already_imported = tuple(binding for binding in SUPPORTED_BINDINGS
                             if f'{binding}.QtCore' in sys.modules)
    if requested and already_imported and requested not in already_imported:
        raise QtApiError(
            f'{QT_API_ENV}={os.environ[QT_API_ENV]!r} requests {requested}, '
            f'but {already_imported[0]} was already imported. Two Qt '
            'bindings cannot be used in the same process.'
            )
    if requested:
        return (requested,)
    return already_imported[:1] or SUPPORTED_BINDINGS


def import_binding(binding, modules=()):
    """Import `binding` and its `modules`, and check its version.

    This defines when a binding counts as not installed, broken, or
    too old, both for selecting a binding and for checking bindings
    in a subprocess (see viperleed.gui.detect_graphics). A binding
    counts as not installed only if its package, or its QtCore, is
    missing. The latter happens if the package is a namespace package
    without modules (e.g., if only PyQt5-sip and PyQt5-Qt5 are left
    after uninstalling PyQt5). Any other missing module (e.g.,
    PyQt5.sip, or QtWidgets) means that the binding is installed, but
    broken. QtCore is always imported, and the version is checked
    before importing other modules.

    Parameters
    ----------
    binding : str
        Either PYQT5 or PYQT6.
    modules : Sequence of str, optional
        Names of the Qt modules (e.g., 'QtWidgets') to import besides
        QtCore. Default is none.

    Raises
    ------
    QtBindingNotFoundError
        If `binding` is not installed.
    QtBindingTooOldError
        If `binding` is older than its MIN_VERSIONS.
    ImportError
        If `binding` is installed, but broken (e.g., it misses some
        of its own modules, or a library it needs fails to load).
    """
    qt_core = f'{binding}.QtCore'
    # Each module once, in this order (QtCore may be in modules)
    to_import = dict.fromkeys((binding, qt_core,
                               *(f'{binding}.{m}' for m in modules)))
    for module in to_import:
        try:
            imported = importlib.import_module(module)
        except ModuleNotFoundError as exc:
            if exc.name in (binding, qt_core):
                raise QtBindingNotFoundError(str(exc)) from exc
            raise
        if module == qt_core:
            check_version(binding, imported)


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


def _get_requested_binding():
    """Return the binding requested via QT_API, or None."""
    requested = os.environ.get(QT_API_ENV, '').strip()
    if not requested:
        return None
    try:
        return _QT_API_VALUES[requested.lower()]
    except KeyError:
        supported = ', '.join(repr(v) for v in _QT_API_VALUES)
        raise QtApiError(
            f'Unsupported value {requested!r} for the {QT_API_ENV} '
            f'environment variable. ViPErLEED supports: {supported}.'
            ) from None


def _select_binding():
    """Return the name of the binding that should be used."""
    requested = _get_requested_binding()
    candidates = get_candidates()
    errors = []
    broken = False
    for binding in candidates:
        try:
            import_binding(binding)
        except QtBindingNotFoundError as exc:
            errors.append(f'{binding}: {exc}')
        except ImportError as exc:  # Installed, but broken or too old
            errors.append(f'{binding}: {exc}')
            broken = True
        else:
            return binding
    msg = (f'Could not import {requested}, requested via '
           f'{QT_API_ENV}={os.environ[QT_API_ENV]!r}. ' if requested
           else f'Could not import any Qt binding among {candidates}. ')
    msg += 'Errors: ' + '; '.join(errors)
    if broken:
        raise ImportError(msg)
    if requested:
        raise QtApiError(msg)
    raise QtBindingNotFoundError(msg)
