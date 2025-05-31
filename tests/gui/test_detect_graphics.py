"""Tests for module detect_graphics of viperleed.gui."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-05-22'
__license__ = 'GPLv3+'


import pytest
from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.gui.detect_graphics import _HAS_PYQT
from viperleed.gui.detect_graphics import _TIMEOUT
from viperleed.gui.detect_graphics import _init_dummy_qapp
from viperleed.gui.detect_graphics import has_graphics
from viperleed.gui.detect_graphics import has_pyqt

_MODULE = 'viperleed.gui.detect_graphics'


# TODO: would be better to find a way to mock at least _init_dummy_qapp
# instead of mocking multiprocessing.Process. My tests seem to reveal
# that, while has_graphics gets the correct target, that target is
# never actually called.


@fixture(name='mock_multiprocessing')
def fixture_mock_process(mocker):
    """Replace multiprocessing.Process with a mock.

    Notice that it is necessary to do so in order to correctly
    test the behavior of has_display. In fact, we can't mock
    detect_graphics.QApplication as the module is re-imported
    by multiprocessing. We could, in principle, mock the target
    of multiprocessing.Process, i.e., _init_dummy_qapp. However,
    I tried without success to do so.
    """
    def _mock(has_display):
        mock_process = mocker.MagicMock()
        mock_process_cls = mocker.MagicMock(return_value=mock_process)
        # Mimic the expected behavior when a display is not present
        mock_process.is_alive.return_value = not has_display
        mocker.patch(f'{_MODULE}.Process', mock_process_cls)
        return mock_process, mock_process_cls
    return _mock


def test_init_dummy_qapp(mocker):
    """Check calls in _init_dummy_qapp."""
    mock_qapp = mocker.patch(f'{_MODULE}.QApplication',
                             create=not _HAS_PYQT)
    with pytest.raises(SystemExit):
        _init_dummy_qapp()
    mock_qapp.assert_called_once_with([])


@parametrize(pyqt_present=(True, False))
@parametrize(has_display=(True, False))
def test_has_graphics(pyqt_present, has_display, mock_multiprocessing, mocker):
    """Ensure correct identification of graphics capabilities."""
    mocker.patch(f'{_MODULE}._HAS_PYQT', pyqt_present)
    proc, proc_cls = mock_multiprocessing(has_display)
    calls = {}
    if pyqt_present:
        calls = {
            'start': mocker.call(),
            'join': mocker.call(_TIMEOUT),
            'is_alive': mocker.call(),
            }
    if pyqt_present and not has_display:
        calls['terminate'] = mocker.call()

    result = has_graphics()
    assert result == (pyqt_present and has_display)

    # Check also internal calls
    if pyqt_present:
        proc_cls.assert_called_once_with(target=_init_dummy_qapp)
    for method_name, call in calls.items():
        method = getattr(proc, method_name)
        assert method.mock_calls == [call]


@parametrize(pyqt_present=(True, False))
def test_has_pyqt(pyqt_present, mocker):
    """Ensure correct detection of PyQt5 presence."""
    mocker.patch(f'{_MODULE}._HAS_PYQT', pyqt_present)
    assert has_pyqt() == pyqt_present
