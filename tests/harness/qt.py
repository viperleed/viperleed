"""Load pytest-qt only if PyQt5 is available.

pytest-qt aborts the whole test session if no Qt binding is installed,
even if no test needs Qt (e.g., for calc-only installations). Hence,
this plugin registers pytest-qt only when PyQt5, the Qt binding of
viperleed.gui, can be imported, and blocks it otherwise. This happens
when this plugin is registered (in pytest_addoption, called right away
for a new plugin), that is, before pytest-qt is configured. pytest-qt
is registered under the same name as when pytest loads it
automatically ('pytest-qt'), so that it is never registered twice
(e.g., in nested runs via pytester, where automatic loading may come
first or later).

pytest-qt is also told to use PyQt5 (via the PYTEST_QT_API environment
variable), as it would otherwise prefer other bindings, if installed
(e.g., PyQt6): two bindings must never be imported in the same process.

Nested pytest runs (via pytester) should load this plugin, as
pytest-qt is otherwise loaded automatically (see run_pytest in
tests/harness/tests/conftest.py).

When pytest-qt is not loaded (PyQt5 or pytest-qt is not installed),
its qapp and qtbot fixtures are replaced by fixtures that skip the
test via missing_dependency (or fail with --require-all). They make
the test fail instead if PyQt5 is installed but broken (e.g., a
library it needs cannot be loaded).

The tests run on Qt's 'offscreen' platform, so that no windows show
up while testing, and tests that create widgets also run without a
display. This plugin sets the QT_QPA_PLATFORM environment variable
accordingly, unless it is set already (e.g., in CI), or unless the
--show-windows option is given. As Qt takes the platform only once,
when the QApplication is created, the platform is the same for all
the tests of a session. Tests that need a real display can be marked
'needs_display': they are skipped on a headless platform. Tests that
need a headless platform can be marked 'needs_offscreen': they are
skipped otherwise.
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

BINDING = 'PyQt5'  # The Qt binding used by viperleed.gui
DEPENDENCY = 'PyQt'  # Its name among the optional dependencies
PYTEST_QT = 'pytest-qt'  # Name of the pytest11 entry point
PYTEST_QT_API_ENV = 'PYTEST_QT_API'
QT_QPA_PLATFORM_ENV = 'QT_QPA_PLATFORM'
SHOW_WINDOWS_OPTION = '--show-windows'

# Qt platforms that show no windows
_HEADLESS_PLATFORMS = (
    'minimal',
    'offscreen',
    )

# The binding is not installed if any of these is missing. QtCore
# is missing if only PyQt5-sip and PyQt5-Qt5 are left after
# uninstalling PyQt5: PyQt5 is then an empty namespace package.
_NOT_INSTALLED_IF_MISSING = (BINDING, f'{BINDING}.QtCore')


def is_not_installed(exc):
    """Return whether ImportError `exc` means that PyQt5 is missing.

    Any other ImportError raised when importing a module that needs
    PyQt5 (e.g., if one of the PyQt5 modules or a library that they
    need is missing) means that PyQt5 is installed, but broken.

    Parameters
    ----------
    exc : ImportError
        The exception raised when importing a module.

    Returns
    -------
    not_installed : bool
        Whether `exc` was raised as PyQt5 is not installed.
    """
    return (isinstance(exc, ModuleNotFoundError)
            and exc.name in _NOT_INSTALLED_IF_MISSING)


def pytest_addoption(parser, pluginmanager):
    """Add --show-windows. Register pytest-qt if it can be used."""
    parser.addoption(
        SHOW_WINDOWS_OPTION,
        action='store_true',
        help=('show the windows of the GUI tests, rather than running '
              f'them on the offscreen Qt platform ({QT_QPA_PLATFORM_ENV} '
              'is respected also without this option)'),
        )
    if _MISSING is not None:
        pluginmanager.set_blocked(PYTEST_QT)  # Also unregisters it
    elif pluginmanager.get_plugin(PYTEST_QT) is None:
        pluginmanager.register(importlib.import_module('pytestqt.plugin'),
                               PYTEST_QT)


def pytest_collection_modifyitems(config, items):
    """Skip the tests that cannot run on the current Qt platform."""
    _ = config  # Unused, but part of the hook
    headless = _is_headless()
    platform = os.environ.get(QT_QPA_PLATFORM_ENV)
    for item in items:
        if headless and item.get_closest_marker('needs_display'):
            item.add_marker(pytest.mark.skip(
                reason=(f'Needs a display, but Qt runs on {platform!r}. '
                        f'Use {SHOW_WINDOWS_OPTION} (and do not set '
                        f'{QT_QPA_PLATFORM_ENV} to a headless platform)')
                ))
        if not headless and item.get_closest_marker('needs_offscreen'):
            item.add_marker(pytest.mark.skip(
                reason='Needs a headless Qt platform, e.g., offscreen'
                ))


def pytest_configure(config):
    """Register the markers, and use the offscreen Qt platform."""
    config.addinivalue_line(
        'markers',
        'needs_display: needs a real display (run with '
        f'{SHOW_WINDOWS_OPTION})',
        )
    config.addinivalue_line(
        'markers',
        'needs_offscreen: needs a headless Qt platform (e.g., offscreen)',
        )
    # Before any QApplication is created (also by the subprocess
    # checks of viperleed.gui.detect_graphics, which inherit it)
    if not config.getoption(SHOW_WINDOWS_OPTION):
        os.environ.setdefault(QT_QPA_PLATFORM_ENV, 'offscreen')


def _find_missing():
    """Return what is missing to use pytest-qt with PyQt5, or None.

    If nothing is missing, also tell pytest-qt to use PyQt5 (via
    the PYTEST_QT_API environment variable).

    Returns
    -------
    missing : tuple or None
        None if pytest-qt can be used. Otherwise (name, detail, fail):
        the name of the optional dependency that is missing (either
        DEPENDENCY or PYTEST_QT), why it cannot be used, and whether
        tests that need it should fail rather than be skipped, i.e.,
        whether it is a problem to fix (e.g., a broken PyQt5) rather
        than an optional dependency that is not installed.

    Raises
    ------
    pytest.UsageError
        If PYTEST_QT_API requests a binding other than PyQt5.
    """
    for module in _NOT_INSTALLED_IF_MISSING:
        try:
            importlib.import_module(module)
        # RuntimeError: e.g., an incompatible sip module
        except (ImportError, RuntimeError) as exc:
            if is_not_installed(exc):
                return DEPENDENCY, f'{BINDING} is not installed', False
            broken = f'{BINDING} is installed, but broken. {exc}'
            return DEPENDENCY, broken, True
    if importlib.util.find_spec('pytestqt') is None:
        return PYTEST_QT, '', False
    requested = os.environ.get(PYTEST_QT_API_ENV, '').strip().lower()
    if requested and requested != BINDING.lower():
        msg = (f'{PYTEST_QT_API_ENV}={requested!r} is not supported: '
               f'viperleed.gui uses {BINDING}.')
        raise pytest.UsageError(msg)
    os.environ[PYTEST_QT_API_ENV] = BINDING.lower()
    return None


_MISSING = _find_missing()


if _MISSING is not None:
    @pytest.fixture(name='qapp')
    def fixture_qapp(missing_dependency):
        """Skip (or fail) the test, as pytest-qt cannot be used."""
        _report_missing(missing_dependency)

    @pytest.fixture(name='qtbot')
    def fixture_qtbot(missing_dependency):
        """Skip (or fail) the test, as pytest-qt cannot be used."""
        _report_missing(missing_dependency)


def _is_headless():
    """Return whether Qt runs on a platform that shows no windows."""
    platform = os.environ.get(QT_QPA_PLATFORM_ENV, '')
    # The name may be followed by options, e.g., 'offscreen:fontengine'
    return platform.split(':')[0].strip() in _HEADLESS_PLATFORMS


def _report_missing(missing_dependency):
    """Skip the current test, or fail it if Qt cannot be used."""
    name, detail, fail = _MISSING
    if fail:
        pytest.fail(f'Cannot use {name}: {detail}')
    missing_dependency(name, detail)
