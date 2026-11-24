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
  Qt bindings cannot be used together in the same process. If both
  were already imported, PyQt6 is used, as matplotlib does.
- Otherwise, the first binding that can be imported among PyQt6 and
  PyQt5 (in this order) is used. A binding that is installed but
  broken (e.g., it misses some of its own modules), or older than
  MIN_VERSIONS, is skipped too. If no binding can be used and any
  of them is broken or too old, a plain ImportError is raised rather
  than a QtBindingNotFoundError, as the installation needs fixing
  rather than completing.

Selecting a binding imports only its QtCore. Its other modules (e.g.,
QtWidgets) are not imported, as code that uses only QtCore (e.g.,
viperleed.calc) should not load the graphical ones. A binding that
lacks, e.g., QtWidgets is thus selected, and fails when QtWidgets is
imported. To check all the modules that are needed (e.g., by the
GUI) before selecting a binding, use import_binding, preferably in a
subprocess.

Which bindings would be tried, in which order, is returned by
get_candidates, which imports none of them. This allows checking the
candidates in a subprocess first, as a broken Qt installation may
crash the interpreter when imported. The binding found to work can
then be selected with set_binding, which accepts only candidates.
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

from viperleed.calc.lib.version import Version

PYQT5 = 'PyQt5'
PYQT6 = 'PyQt6'
QT_API_ENV = 'QT_API'
SUPPORTED_BINDINGS = (PYQT6, PYQT5)  # In order of preference
# The oldest (major, minor) versions of the bindings that work. PyQt6
# before 6.1 lacks names of enums used in viperleed.gui (e.g., it has
# Qt.MouseButtons rather than Qt.MouseButton), and some enums that
# Qt functions return as int (e.g., Qt.Key) are not int-based there.
MIN_VERSIONS = {PYQT5: Version('5.15'), PYQT6: Version('6.1')}
_QT_API_VALUES = {binding.lower(): binding for binding in SUPPORTED_BINDINGS}

_SELECTED_BINDING = None


class QtBindingNotFoundError(ImportError):
    """No usable Qt binding could be found."""


class QtApiError(QtBindingNotFoundError):
    """The QT_API environment variable prevents selecting a binding."""


class QtBindingTooOldError(ImportError):
    """A Qt binding is older than its version in MIN_VERSIONS."""

    def __init__(self, *args, qt_version=None):
        """Initialize an instance.

        Parameters
        ----------
        *args : object
            Passed on to ImportError, e.g., the error message.
        qt_version : str or None, optional
            The version of Qt used by the binding (as QT_VERSION_STR),
            if known. Default is None.
        """
        super().__init__(*args)
        self.qt_version = qt_version


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
    oldest = MIN_VERSIONS[binding]
    if Version(qt_core.PYQT_VERSION_STR) < oldest:
        raise QtBindingTooOldError(
            f'{binding} {qt_core.PYQT_VERSION_STR} is too old. ViPErLEED '
            f'needs {binding}>={oldest}.',
            qt_version=qt_core.QT_VERSION_STR,
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
    candidates, _ = _get_candidates()
    return candidates


def import_binding(binding, modules=()):
    """Import `binding` and its `modules`, and check its version.

    This defines when a binding counts as not installed, broken, or
    too old, both for selecting a binding and for checking bindings
    (e.g., in a subprocess, with all the modules the GUI needs). A
    binding counts as not installed only if its package, or its
    QtCore, is missing. The latter happens if the package is a
    namespace package without modules (e.g., if only PyQt5-sip and
    PyQt5-Qt5 are left after uninstalling PyQt5). Any other missing
    module (e.g., PyQt5.sip, or QtWidgets) means that the binding is
    installed, but broken. QtCore is always imported, and the version
    is checked before importing other modules.

    Parameters
    ----------
    binding : str
        Either PYQT5 or PYQT6.
    modules : Sequence of str, optional
        Names of the Qt modules (e.g., 'QtWidgets') to import besides
        QtCore. Default is none.

    Returns
    -------
    qt_core : module
        The QtCore module of `binding`.

    Raises
    ------
    QtBindingNotFoundError
        If `binding` is not installed.
    QtBindingTooOldError
        If `binding` is older than its MIN_VERSIONS.
    ImportError
        If `binding` is installed, but broken (e.g., it misses some
        of its own modules, a library it needs fails to load, or its
        sip module is incompatible with it).
    """
    qt_core = f'{binding}.QtCore'
    # Each module once, in this order (QtCore may be in modules)
    to_import = dict.fromkeys((binding, qt_core,
                               *(f'{binding}.{m}' for m in modules)))
    imported = {}
    for module in to_import:
        try:
            imported[module] = importlib.import_module(module)
        except ModuleNotFoundError as exc:
            if exc.name in (binding, qt_core):
                raise QtBindingNotFoundError(str(exc)) from exc
            raise
        except RuntimeError as exc:
            # E.g., 'the sip module implements API v12.x but the
            # PyQt6.QtCore module requires API v12.y'
            raise ImportError(f'{module}: {exc}') from exc
        if module == qt_core:
            check_version(binding, imported[module])
    return imported[qt_core]


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
    QtApiError
        If `binding` is not among the candidates (see get_candidates):
        if the QT_API environment variable has an unsupported value
        or requests another binding, or if another binding was
        already imported.
    """
    global _SELECTED_BINDING  # pylint: disable=global-statement
    if binding not in SUPPORTED_BINDINGS:
        raise ValueError(f'Unsupported Qt binding {binding!r}')
    if _SELECTED_BINDING == binding:
        return
    if _SELECTED_BINDING is not None:
        raise ValueError(f'Cannot select {binding}: {_SELECTED_BINDING} '
                         'is already in use')
    candidates, requested = _get_candidates()
    if binding not in candidates:  # Then it is a single one
        why = (f'{_describe_qt_api()} requests it' if requested
               else 'it was already imported')
        raise QtApiError(f'Cannot select {binding}: {candidates[0]} must '
                         f'be used, as {why}.')
    _SELECTED_BINDING = binding


def _describe_qt_api():
    """Return the value of QT_API, as 'QT_API=<value>', for messages."""
    return f'{QT_API_ENV}={os.environ[QT_API_ENV]!r}'


def _forget_binding(binding):
    """Remove the modules of `binding` from sys.modules.

    The modules stay loaded (e.g., QtCore of a binding that is too
    old), but nothing takes `binding` as in use any longer: neither
    get_candidates, nor matplotlib, which both prefer a binding whose
    QtCore was imported.
    """
    prefix = f'{binding}.'
    # Collect them first: sys.modules cannot change while iterating
    binding_modules = [name for name in sys.modules
                       if name == binding or name.startswith(prefix)]
    for name in binding_modules:
        del sys.modules[name]


def _get_candidates():
    """Return the candidate bindings, and the one requested via QT_API.

    Returns
    -------
    candidates : tuple
        The names of the bindings to try, in order (see get_candidates).
    requested : str or None
        The binding requested via QT_API, None if QT_API is not set.

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
            f'{_describe_qt_api()} requests {requested}, but '
            f'{already_imported[0]} was already imported. Two Qt '
            'bindings cannot be used in the same process.'
            )
    if requested:
        return (requested,), requested
    return already_imported[:1] or SUPPORTED_BINDINGS, None


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
    candidates, requested = _get_candidates()
    errors = []
    broken = False
    for binding in candidates:
        try:
            import_binding(binding)
        except QtBindingNotFoundError as exc:
            errors.append(f'{binding}: {exc}')
        except QtBindingTooOldError as exc:
            errors.append(f'{binding}: {exc}')
            broken = True
            _forget_binding(binding)  # Its QtCore was imported
        except ImportError as exc:  # Installed, but broken
            errors.append(f'{binding}: {exc}')
            broken = True
        else:
            return binding
    msg = (f'Could not import {requested}, requested via '
           f'{_describe_qt_api()}. ' if requested
           else f'Could not import any Qt binding among {candidates}. ')
    msg += 'Errors: ' + '; '.join(errors)
    if broken:
        raise ImportError(msg)
    if requested:
        raise QtApiError(msg)
    raise QtBindingNotFoundError(msg)
