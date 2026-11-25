"""Tests for module selectplugin of viperleed.gui.

The plug-ins that the selector can open are collected when the tests
are collected. The selector is imported then only if a Qt binding is
available, so that this file can be collected without one.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-16'
__license__ = 'GPLv3+'

import importlib

import pytest
from pytest_cases import fixture
from pytest_cases import parametrize

_SELECTOR = 'viperleed.gui.selectplugin'
# A parameter for when nothing can be collected (e.g., no Qt binding):
# the test is skipped without Qt, or fails importing what is missing.
_PLACEHOLDER = 'nothing-collected'


def _collect_plugins():
    """Return the names of the plug-ins that the selector can open."""
    try:
        selector = importlib.import_module(_SELECTOR)
    # Any exception is reported by test_smoke.py::test_import
    # pylint: disable-next=broad-exception-caught
    except Exception:  # noqa: BLE001
        return (_PLACEHOLDER,)
    plugins = selector.ViPErLEEDSelectPlugin.modules
    return tuple(name for name, (_, cls) in plugins.items() if cls)


@fixture(name='selector')
def fixture_selector(import_with_qt, show):
    """Return a shown ViPErLEEDSelectPlugin."""
    module = import_with_qt(_SELECTOR)
    selector = module.ViPErLEEDSelectPlugin()
    show(selector)
    return selector


@pytest.mark.usefixtures('no_blocking_dialogs')
@parametrize(plugin=_collect_plugins())
def test_open_plugin(plugin, selector, qtbot, qtc):
    """Check that the plug-in selector opens and closes a plug-in."""
    if plugin == _PLACEHOLDER:
        # Qt is available, but the selector failed to import when
        # collecting. test_smoke.py::test_import gives the reason.
        pytest.fail('No plug-in was collected: the plug-in selector '
                    'could not be imported when collecting the tests')
    # pylint: disable=protected-access
    qtbot.mouseClick(selector._btns[plugin], qtc.Qt.MouseButton.LeftButton)
    module = selector._open_modules[plugin]
    assert module is not None
    qtbot.waitExposed(module)
    assert selector.close()
    assert selector._open_modules[plugin] is None
