"""Tests for PyQtSanity and PyQtSanityChecker of detect_graphics."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-17'
__license__ = 'GPLv3+'

from types import SimpleNamespace
from unittest.mock import ANY

from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.gui.detect_graphics import _DISPLAY_TIMEOUT
from viperleed.gui.detect_graphics import _EXIT_CODE_OFFSET
from viperleed.gui.detect_graphics import _QtProbe
from viperleed.gui.detect_graphics import _STARTUP_TIMEOUT
from viperleed.gui.detect_graphics import PyQtSanity
from viperleed.gui.detect_graphics import PyQtSanityChecker
from viperleed.gui.qt_compat.binding import PYQT5
from viperleed.gui.qt_compat.binding import QtBindingNotFoundError

from .conftest import BOTH
from .conftest import MODULE


# TODO: would be better to find a way to mock at least
# _QtProbe instead of mocking
# multiprocessing.Process. My tests seem to reveal
# that, while _check_pyqt_sanity gets the correct
# target, that target is never actually called.


@fixture(name='mock_multiprocessing')
def fixture_mock_multiprocessing(mocker):
    """Replace the multiprocessing tools of the checker with mocks.

    Notice that it is necessary to do so in order to correctly
    test the behavior of has_display. In fact, we can't mock
    detect_graphics.QApplication as the module is re-imported
    by multiprocessing. We could, in principle, mock the target
    of multiprocessing.Process, i.e., _init_dummy_qapp. However,
    I tried without success to do so.

    The returned factory takes the following arguments:
    ready : {'imports', 'exited', 'none'}
        Which of the waited-for objects is ready within the start-up
        time: the pipe (a binding was imported), the process sentinel
        (the subprocess exited before), or none of them.
    alive : bool
        Whether the subprocess is still alive after the QApplication
        time limit, i.e., whether creating the QApplication hangs.
    exitcode : object
        The exit code of the subprocess.
    """
    def _mock(ready='imports', alive=False, exitcode=0):
        imports_done = mocker.MagicMock(name='imports_done')
        imports_done.recv.return_value = PYQT5
        child_end = mocker.MagicMock(name='child_end')
        process = mocker.MagicMock(name='process', exitcode=exitcode)
        process.is_alive.return_value = alive
        ready_objects = {
            'imports': [imports_done],
            'exited': [process.sentinel],
            'none': [],
            }
        order = []  # Names of selected calls, in the order they're made
        process.start.side_effect = lambda: order.append('start')
        child_end.close.side_effect = lambda: order.append('close_child')

        def _wait(*_):
            order.append('wait')
            return ready_objects[ready]
        return SimpleNamespace(
            order=order,
            imports_done=imports_done,
            child_end=child_end,
            process=process,
            pipe=mocker.patch(f'{MODULE}.Pipe',
                              return_value=(imports_done, child_end)),
            process_cls=mocker.patch(f'{MODULE}.Process',
                                     return_value=process),
            wait=mocker.patch(f'{MODULE}.wait', side_effect=_wait),
            candidates=mocker.patch(f'{MODULE}.get_candidates',
                                    return_value=BOTH),
            )
    return _mock


@parametrize(sanity=PyQtSanity)
def test_sanity_importable(sanity):
    """Check which statuses allow importing Qt."""
    expect = sanity in {PyQtSanity.OK, PyQtSanity.NO_DISPLAY}
    assert sanity.importable == expect


class TestPyQtSanityChecker:
    """Tests for the PyQtSanityChecker class."""

    @fixture(autouse=True)
    def fresh_checker(self, mocker):
        """Forget previous checks."""
        mocker.patch.object(PyQtSanityChecker, '_status', None)
        mocker.patch.object(PyQtSanityChecker, '_binding', None)

    @fixture(name='set_binding', autouse=True)
    def fixture_set_binding(self, mocker):
        """Do not select bindings globally."""
        return mocker.patch(f'{MODULE}.set_binding')

    def test_checks_only_once(self, mocker):
        """Ensure that checking is skipped when already done."""
        mock_sanity = mocker.MagicMock()
        mock_do_check = mocker.patch.object(PyQtSanityChecker,
                                            '_check_pyqt_sanity',
                                            return_value=(mock_sanity, None))
        checker = PyQtSanityChecker()

        # Call it once: performs the check.
        assert checker() is mock_sanity
        # pylint: disable-next=protected-access
        assert PyQtSanityChecker._status is mock_sanity

        # Call it again: should skip the check
        assert checker() is mock_sanity
        mock_do_check.assert_called_once_with()

    @parametrize(binding=BOTH)
    def test_selects_binding(self, binding, set_binding, mocker):
        """Check that the binding found is selected for qt_compat."""
        mocker.patch.object(PyQtSanityChecker, '_check_pyqt_sanity',
                            return_value=(PyQtSanity.OK, binding))
        checker = PyQtSanityChecker()
        assert checker() is PyQtSanity.OK
        assert checker.binding == binding
        set_binding.assert_called_once_with(binding)

    def test_no_binding(self, set_binding, mocker):
        """Check that no binding is selected if none was found."""
        mocker.patch.object(PyQtSanityChecker, '_check_pyqt_sanity',
                            return_value=(PyQtSanity.NOT_FOUND, None))
        checker = PyQtSanityChecker()
        assert checker() is PyQtSanity.NOT_FOUND
        assert checker.binding is None
        set_binding.assert_not_called()

    def test_forget_status(self, mocker):
        """Ensure resetting re-checks."""
        mock_do_check = mocker.patch.object(PyQtSanityChecker,
                                            '_check_pyqt_sanity',
                                            return_value=(None, PYQT5))
        checker = PyQtSanityChecker()
        # Check once to ensure forget_status does the right thing
        checker()

        checker.forget_status()
        assert checker.binding is None
        checker()

        # pylint: disable-next=magic-value-comparison
        assert mock_do_check.call_count == 2

    def test_assume_ok(self, mocker):
        """Ensure that assuming a sane Qt binding skips checking."""
        mock_sanity = mocker.MagicMock()
        mock_do_check = mocker.patch.object(PyQtSanityChecker,
                                            '_check_pyqt_sanity',
                                            return_value=(mock_sanity, None))
        checker = PyQtSanityChecker()
        checker.assume_ok()
        assert checker() is PyQtSanity.OK
        assert checker.binding is None  # qt_compat selects one
        mock_do_check.assert_not_called()

        # forget_status undoes it
        checker.forget_status()
        assert checker() is mock_sanity
        mock_do_check.assert_called_once_with()

    _exit_codes = {
        # Offset exit codes are mapped to the PyQtSanity members
        **{_EXIT_CODE_OFFSET + s.value: s for s in PyQtSanity},
        # Any other exit code is a crash. In particular, the ones that
        # coincide with the values of PyQtSanity members, like 1 for
        # an unhandled exception, must not be taken as such members.
        **{s.value: PyQtSanity.RUNTIME_CRASH for s in PyQtSanity},
        -6: PyQtSanity.RUNTIME_CRASH,   # SIGABRT
        134: PyQtSanity.RUNTIME_CRASH,  # SIGABRT on Unix
        _EXIT_CODE_OFFSET - 1: PyQtSanity.RUNTIME_CRASH,
        _EXIT_CODE_OFFSET + len(PyQtSanity): PyQtSanity.RUNTIME_CRASH,
        }

    @staticmethod
    def check_common_calls(mocks):
        """Check the calls made whatever the outcome of the check."""
        mocks.pipe.assert_called_once_with(duplex=False)
        mocks.process_cls.assert_called_once_with(target=ANY,
                                                  args=(mocks.child_end,))
        probe = mocks.process_cls.call_args[1]['target']
        assert isinstance(probe, _QtProbe)
        assert probe._candidates == BOTH  # pylint: disable=W0212
        mocks.process.start.assert_called_once_with()
        mocks.child_end.close.assert_called_once_with()
        mocks.wait.assert_called_once_with(
            (mocks.imports_done, mocks.process.sentinel),
            _STARTUP_TIMEOUT,
            )
        mocks.imports_done.close.assert_called_once_with()
        # The parent closes its copy of the child's end of the
        # pipe right after starting, i.e., before waiting.
        assert mocks.order == ['start', 'close_child', 'wait']

    @parametrize(ready=('imports', 'exited'))
    @parametrize('exit_code,expect', _exit_codes.items())
    def test_check_sanity_exits(self, ready, exit_code, expect,
                                mock_multiprocessing):
        """Check the result when the subprocess exits by itself."""
        mocks = mock_multiprocessing(ready=ready, exitcode=exit_code)
        # pylint: disable-next=protected-access
        sanity, binding = PyQtSanityChecker._check_pyqt_sanity()
        assert sanity is expect
        # The binding is known only if it was sent before exiting,
        # and it is reported only if it can be used
        # pylint: disable-next=magic-value-comparison
        usable = ready == 'imports' and expect.importable
        expect_binding = PYQT5 if usable else None
        assert binding == expect_binding
        self.check_common_calls(mocks)
        mocks.process.join.assert_called_once_with(_DISPLAY_TIMEOUT)
        mocks.process.terminate.assert_not_called()

    def test_check_sanity_invalid_qt_api(self, mock_multiprocessing):
        """Check the result when QT_API has an invalid value."""
        mocks = mock_multiprocessing()
        mocks.candidates.side_effect = QtBindingNotFoundError
        # pylint: disable-next=protected-access
        result = PyQtSanityChecker._check_pyqt_sanity()
        assert result == (PyQtSanity.INVALID_QT_API, None)
        mocks.pipe.assert_not_called()
        mocks.process_cls.assert_not_called()

    def test_check_sanity_no_display(self, mock_multiprocessing):
        """Check the result when creating a QApplication hangs."""
        mocks = mock_multiprocessing(ready='imports', alive=True)
        # pylint: disable-next=protected-access
        result = PyQtSanityChecker._check_pyqt_sanity()
        assert result == (PyQtSanity.NO_DISPLAY, PYQT5)
        self.check_common_calls(mocks)
        mocks.process.join.assert_called_once_with(_DISPLAY_TIMEOUT)
        mocks.process.terminate.assert_called_once_with()

    def test_check_sanity_no_message(self, mock_multiprocessing):
        """Check the result when the subprocess dies without sending."""
        mocks = mock_multiprocessing(ready='imports',
                                     exitcode=-6)  # SIGABRT
        mocks.imports_done.recv.side_effect = EOFError
        # pylint: disable-next=protected-access
        result = PyQtSanityChecker._check_pyqt_sanity()
        assert result == (PyQtSanity.RUNTIME_CRASH, None)

    def test_check_sanity_unresponsive(self, mock_multiprocessing):
        """Check the result when importing does not finish in time."""
        mocks = mock_multiprocessing(ready='none', alive=True)
        # pylint: disable-next=protected-access
        result = PyQtSanityChecker._check_pyqt_sanity()
        assert result == (PyQtSanity.UNRESPONSIVE, None)
        self.check_common_calls(mocks)
        mocks.process.join.assert_not_called()
        mocks.process.terminate.assert_called_once_with()
