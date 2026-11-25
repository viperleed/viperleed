"""Tests for module qt of tests.harness."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-23'
__license__ = 'GPLv3+'

import pytest
from pytest_cases import parametrize

from .. import optional_dependencies
from .. import qt
from .conftest import run_pytest

_BROKEN_FILES = {  # A PyQt5 package whose QtCore cannot be loaded
    'PyQt5/__init__': '',
    'PyQt5/QtCore': "raise ImportError('a library is missing')",
    }
# QT_QPA_PLATFORM before the run, options, and the test that runs
_DISPLAY_MARKS = {  # pylint: disable=R6101
    'default': (None, (), 'test_offscreen'),
    'minimal with options': ('minimal:a=1', (), 'test_offscreen'),
    'other platform': ('xcb', (), 'test_display'),
    'show': (None, (qt.SHOW_WINDOWS_OPTION,), 'test_display'),
    'show but offscreen': ('offscreen', (qt.SHOW_WINDOWS_OPTION,),
                           'test_offscreen'),
    }
_HIDE_PLUGIN = 'hide_pyqt5'
_HIDE_FILES = {  # A plugin, loaded before qt, that hides PyQt5
    _HIDE_PLUGIN: "import sys\nsys.modules['PyQt5'] = None\n",
    }
_MARKED_FILE = """
import pytest

@pytest.mark.needs_display
def test_display():
    pass

@pytest.mark.needs_offscreen
def test_offscreen():
    pass
"""
# QT_QPA_PLATFORM before the run, options, and its value in the tests
_PLATFORMS = {  # pylint: disable=R6101
    'default': (None, (), 'offscreen'),
    'set': ('xcb', (), 'xcb'),
    'show': (None, (qt.SHOW_WINDOWS_OPTION,), None),
    'show and set': ('offscreen', (qt.SHOW_WINDOWS_OPTION,), 'offscreen'),
    }
_REQUIRE = optional_dependencies.REQUIRE_ALL_OPTION
_SKIP_REASONS = {  # The test that runs: why the other one is skipped
    'test_display': '*Needs a headless Qt platform*',
    'test_offscreen': f'*Needs a display*{qt.SHOW_WINDOWS_OPTION}*',
    }
_TEST_FILE = """
def test_qtbot(qtbot):
    pass

def test_unrelated():
    pass
"""


def _run_without(pytester, files, *args):
    """Run the sample tests in a subprocess, with extra `files`."""
    pytester.makepyfile(test_sample=_TEST_FILE, **files)
    plugin = ('-p', _HIDE_PLUGIN) if _HIDE_PLUGIN in files else ()
    return run_pytest(pytester, *plugin, *args, '-rsf', own_process=True)


def _set_platform(monkeypatch, platform):
    """Set QT_QPA_PLATFORM to `platform`, or unset it if None."""
    if platform is None:
        monkeypatch.delenv(qt.QT_QPA_PLATFORM_ENV, raising=False)
    else:
        monkeypatch.setenv(qt.QT_QPA_PLATFORM_ENV, platform)


@pytest.mark.slow  # About 1 s (pytester, subprocess)
@parametrize(args=((), (_REQUIRE,)), ids=('default', 'require-all'))
def test_broken(pytester, args):
    """Check that tests needing Qt fail if PyQt5 is broken."""
    result = _run_without(pytester, _BROKEN_FILES, *args)
    result.assert_outcomes(passed=1, errors=1)  # Fail in qtbot
    result.stdout.fnmatch_lines([
        '*Cannot use PyQt: PyQt5 is installed, but broken.*library*',
        ])


@parametrize('platform,args,runs', _DISPLAY_MARKS.values(),
             ids=_DISPLAY_MARKS)
def test_display_marks(platform, args, runs, pytester, monkeypatch):
    """Check that only the tests suited to the Qt platform run."""
    _set_platform(monkeypatch, platform)
    pytester.makepyfile(_MARKED_FILE)
    result = run_pytest(pytester, *args, '-v', '-rs')
    result.assert_outcomes(passed=1, skipped=1)
    result.stdout.fnmatch_lines([f'*::{runs} PASSED*', _SKIP_REASONS[runs]])


@pytest.mark.slow  # About 1 s (pytester, subprocess)
@parametrize('args,outcome',
             (((), 'skipped'), ((_REQUIRE,), 'errors')),  # Fail in qtbot
             ids=('default', 'require-all'))
def test_not_installed(pytester, args, outcome):
    """Check that tests needing Qt skip (or fail) without PyQt5."""
    result = _run_without(pytester, _HIDE_FILES, *args)
    result.assert_outcomes(passed=1, **{outcome: 1})
    result.stdout.fnmatch_lines([
        '*optional dependency PyQt not available: PyQt5 is not installed*',
        ])


@parametrize('before,args,expect', _PLATFORMS.values(), ids=_PLATFORMS)
def test_platform(before, args, expect, pytester, monkeypatch):
    """Check that tests run offscreen, unless told otherwise."""
    _set_platform(monkeypatch, before)
    pytester.makepyfile(f"""
        import os

        def test_platform():
            platform = os.environ.get({qt.QT_QPA_PLATFORM_ENV!r})
            assert platform == {expect!r}
        """)
    result = run_pytest(pytester, *args)
    result.assert_outcomes(passed=1)


@pytest.mark.usefixtures('qapp')  # Skip if pytest-qt cannot be used
def test_pytest_qt_uses_pyqt5(pytester):
    """Check that pytest-qt is loaded once, and that it uses PyQt5."""
    pytester.makepyfile("""
        from pytestqt.qt_compat import qt_api

        def test_binding(qtbot, pytestconfig):
            assert qt_api.pytest_qt_api == 'pyqt5'
            manager = pytestconfig.pluginmanager
            plugins = [p for p in manager.get_plugins()
                       if getattr(p, '__name__', '') == 'pytestqt.plugin']
            assert len(plugins) == 1
        """)
    result = run_pytest(pytester)
    result.assert_outcomes(passed=1)


@pytest.mark.slow  # About 1 s (pytester, subprocess)
@pytest.mark.usefixtures('qapp')  # Skip if pytest-qt cannot be used
def test_pytest_qt_wrong_api(pytester, monkeypatch):
    """Check that pytest-qt cannot be told to use another binding."""
    monkeypatch.setenv(qt.PYTEST_QT_API_ENV, 'pyqt6')
    pytester.makepyfile(_TEST_FILE)
    result = run_pytest(pytester, own_process=True)
    assert result.ret == pytest.ExitCode.USAGE_ERROR
    result.stderr.fnmatch_lines([
        f"*{qt.PYTEST_QT_API_ENV}='pyqt6' is not supported*",
        ])
