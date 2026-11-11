"""Load pytest-qt only if a Qt binding is available.

pytest-qt aborts the whole test session if no Qt binding is installed,
even if no test needs Qt (e.g., for calc-only installations). Hence,
this plugin registers pytest-qt only when viperleed finds a usable Qt
binding, and blocks it otherwise. This happens when this plugin is
registered (in pytest_addoption, called right away for a new plugin),
that is, before pytest-qt is configured. pytest-qt is registered under
the same name as when pytest loads it automatically ('pytest-qt'), so
that it is never registered twice (e.g., in nested runs via pytester,
where automatic loading may come first or later).

pytest-qt is also told to use the same binding as viperleed.gui's
qt_compat (selected via the QT_API environment variable, see
viperleed.gui.qt_compat.binding): two bindings must never be imported
in the same process.

Nested pytest runs (via pytester) that do not load this plugin should
pass -p no:pytest-qt, as pytest-qt is otherwise loaded automatically.

When pytest-qt is not loaded (no Qt binding, or pytest-qt is not
installed), its qapp and qtbot fixtures are replaced by fixtures that
skip the test via missing_dependency (or fail with --require-all).
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-10'
__license__ = 'GPLv3+'

import importlib
import importlib.util
import os

import pytest

PYTEST_QT = 'pytest-qt'  # Name of the pytest11 entry point
PYTEST_QT_API_ENV = 'PYTEST_QT_API'
QT_COMPAT_BINDING = 'viperleed.gui.qt_compat.binding'


def _find_missing():
    """Return (name, detail) of what pytest-qt misses, or None."""
    try:
        # Not imported at the top, as the installed viperleed may not
        # have qt_compat (e.g., an editable install of another branch).
        # Then, the tests that need Qt are skipped, but all others run.
        binding_module = importlib.import_module(QT_COMPAT_BINDING)
    except ImportError as exc:
        return QT_COMPAT_BINDING, str(exc)
    try:
        binding = binding_module.get_binding()
    except binding_module.QtBindingNotFoundError as exc:
        return 'PyQt5/PyQt6', str(exc)
    if importlib.util.find_spec('pytestqt') is None:
        return PYTEST_QT, ''
    requested = os.environ.get(PYTEST_QT_API_ENV, '').strip().lower()
    if requested and requested != binding.lower():
        msg = (f'{PYTEST_QT_API_ENV}={requested!r} differs from the Qt '
               f'binding used by viperleed ({binding}). Select the binding '
               'via the QT_API environment variable instead.')
        raise pytest.UsageError(msg)
    os.environ[PYTEST_QT_API_ENV] = binding.lower()
    return None


_MISSING = _find_missing()


def pytest_addoption(parser, pluginmanager):
    """Register pytest-qt if it can be used, block it otherwise."""
    _ = parser  # Unused, but needed to get pluginmanager
    if _MISSING is not None:
        pluginmanager.set_blocked(PYTEST_QT)  # Also unregisters it
    elif pluginmanager.get_plugin(PYTEST_QT) is None:
        pluginmanager.register(importlib.import_module('pytestqt.plugin'),
                               PYTEST_QT)


if _MISSING is not None:
    @pytest.fixture(name='qapp')
    def fixture_qapp(missing_dependency):
        """Skip (or fail) the test, as pytest-qt cannot be used."""
        missing_dependency(*_MISSING)

    @pytest.fixture(name='qtbot')
    def fixture_qtbot(missing_dependency):
        """Skip (or fail) the test, as pytest-qt cannot be used."""
        missing_dependency(*_MISSING)
