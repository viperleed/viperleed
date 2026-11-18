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
from viperleed.gui.qt_compat.binding import PYQT6
from viperleed.gui.qt_compat.binding import QtApiError

from .conftest import BOTH
from .conftest import MODULE

_FALLBACK = {PyQtSanity.NOT_FOUND,
             PyQtSanity.IMPORT_ERROR,
             PyQtSanity.RUNTIME_CRASH}
_STOP = set(PyQtSanity) - _FALLBACK - {PyQtSanity.INVALID_QT_API}


@fixture(name='mock_multiprocessing')
def fixture_mock_multiprocessing(mocker):
    """Replace the multiprocessing tools of the checker with mocks.

    The checker runs a _QtProbe in a subprocess. Mocks made here do
    not reach it (e.g., with the spawn start method, the subprocess
    imports detect_graphics anew), and the probe may hang or crash.
    Hence, Pipe, Process and wait are replaced, and their outcome is
    faked. The tests check that the subprocess gets a _QtProbe as its
    target. _QtProbe itself is tested by calling it in this process
    (see test_qt_probe.py).

    The returned factory takes the following arguments:
    ready : {'imports', 'exited', 'none'}
        Which of the waited-for objects is ready within the start-up
        time: the pipe (the binding was imported), the process sentinel
        (the subprocess exited before), or none of them.
    alive : bool
        Whether the subprocess is still alive after the QApplication
        time limit, i.e., whether creating the QApplication hangs.
    exitcode : object
        The exit code of the subprocess.
    """
    def _mock(ready='imports', alive=False, exitcode=_EXIT_CODE_OFFSET):
        imports_done = mocker.MagicMock(name='imports_done')
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
            )
    return _mock


@parametrize(sanity=PyQtSanity)
class TestPyQtSanity:
    """Tests for the PyQtSanity enumeration."""

    def test_allows_fallback(self, sanity):
        """Check after which statuses another binding is tried."""
        assert sanity.allows_fallback == (sanity in _FALLBACK)

    def test_importable(self, sanity):
        """Check which statuses allow importing Qt."""
        expect = sanity in {PyQtSanity.OK, PyQtSanity.NO_DISPLAY}
        assert sanity.importable == expect


class TestPyQtSanityChecker:
    """Tests for the public interface of PyQtSanityChecker."""

    @fixture(autouse=True)
    def fresh_checker(self, mocker):
        """Forget previous checks."""
        mocker.patch.object(PyQtSanityChecker, '_status', None)
        mocker.patch.object(PyQtSanityChecker, '_binding', None)
        mocker.patch.object(PyQtSanityChecker, '_results', ())

    @fixture(name='set_binding', autouse=True)
    def fixture_set_binding(self, mocker):
        """Do not select bindings globally."""
        return mocker.patch(f'{MODULE}.set_binding')

    @fixture(name='mock_check')
    def fixture_mock_check(self, mocker):
        """Return a function that fakes the results of the checks."""
        def _mock(*results):
            return mocker.patch.object(PyQtSanityChecker,
                                       '_check_bindings',
                                       return_value=results)
        return _mock

    def test_assume_ok(self, mock_check):
        """Ensure that assuming a sane Qt binding skips checking."""
        mock_do_check = mock_check((PYQT6, PyQtSanity.RUNTIME_CRASH))
        checker = PyQtSanityChecker()
        checker.assume_ok()
        assert checker() is PyQtSanity.OK
        assert checker.binding is None  # qt_compat selects one
        assert not checker.results
        mock_do_check.assert_not_called()

        # forget_status undoes it
        checker.forget_status()
        assert checker() is PyQtSanity.RUNTIME_CRASH
        mock_do_check.assert_called_once_with()

    def test_checks_only_once(self, mock_check):
        """Ensure that checking is skipped when already done."""
        mock_do_check = mock_check((PYQT6, PyQtSanity.OK))
        checker = PyQtSanityChecker()
        assert checker() is PyQtSanity.OK
        assert checker() is PyQtSanity.OK
        mock_do_check.assert_called_once_with()

    @parametrize(sanity=(PyQtSanity.IMPORT_ERROR,
                         PyQtSanity.RUNTIME_CRASH,
                         PyQtSanity.UNRESPONSIVE))
    def test_failed_binding(self, sanity, set_binding, mock_check):
        """Check that a binding that failed is reported, not selected."""
        mock_check((PYQT6, sanity))
        checker = PyQtSanityChecker()
        assert checker() is sanity
        assert checker.binding == PYQT6
        set_binding.assert_not_called()

    def test_forget_status(self, mock_check):
        """Ensure resetting re-checks."""
        mock_do_check = mock_check((PYQT5, PyQtSanity.OK))
        checker = PyQtSanityChecker()
        checker()
        checker.forget_status()
        assert checker.binding is None
        assert not checker.results
        checker()
        # pylint: disable-next=magic-value-comparison
        assert mock_do_check.call_count == 2

    def test_no_binding(self, set_binding, mock_check):
        """Check that no binding is selected if none was found."""
        mock_check((PYQT6, PyQtSanity.NOT_FOUND),
                   (PYQT5, PyQtSanity.NOT_FOUND))
        checker = PyQtSanityChecker()
        assert checker() is PyQtSanity.NOT_FOUND
        assert checker.binding is None
        set_binding.assert_not_called()

    def test_results(self, mock_check):
        """Check that the results of all bindings checked are kept."""
        results = ((PYQT6, PyQtSanity.RUNTIME_CRASH),
                   (PYQT5, PyQtSanity.OK))
        mock_check(*results)
        checker = PyQtSanityChecker()
        assert checker() is PyQtSanity.OK
        assert checker.results == results

    @parametrize(binding=BOTH)
    def test_selects_binding(self, binding, set_binding, mock_check):
        """Check that the binding found is selected for qt_compat."""
        mock_check((binding, PyQtSanity.OK))
        checker = PyQtSanityChecker()
        assert checker() is PyQtSanity.OK
        assert checker.binding == binding
        set_binding.assert_called_once_with(binding)


class TestCheckBinding:
    """Tests for checking one binding in a subprocess."""

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
    def check_common_calls(mocks, binding):
        """Check the calls made whatever the outcome of the check."""
        mocks.pipe.assert_called_once_with(duplex=False)
        mocks.process_cls.assert_called_once_with(target=ANY,
                                                  args=(mocks.child_end,))
        probe = mocks.process_cls.call_args[1]['target']
        assert isinstance(probe, _QtProbe)
        assert probe._binding == binding  # pylint: disable=W0212
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

    @parametrize(binding=BOTH)
    @parametrize(ready=('imports', 'exited'))
    @parametrize('exit_code,expect', _exit_codes.items())
    def test_exits(self, binding, ready, exit_code, expect,
                   mock_multiprocessing):
        """Check the result when the subprocess exits by itself."""
        mocks = mock_multiprocessing(ready=ready, exitcode=exit_code)
        # pylint: disable-next=protected-access
        assert PyQtSanityChecker._check_binding(binding) is expect
        self.check_common_calls(mocks, binding)
        mocks.process.join.assert_called_once_with(_DISPLAY_TIMEOUT)
        mocks.process.terminate.assert_not_called()

    def test_no_display(self, mock_multiprocessing, mocker):
        """Check the result when creating a QApplication hangs."""
        mocks = mock_multiprocessing(alive=True)
        # pylint: disable-next=protected-access
        result = PyQtSanityChecker._check_binding(PYQT6)
        assert result is PyQtSanity.NO_DISPLAY
        self.check_common_calls(mocks, PYQT6)
        mocks.process.terminate.assert_called_once_with()
        # Once for the QApplication, once after terminating it
        assert mocks.process.join.mock_calls == [
            mocker.call(_DISPLAY_TIMEOUT),
            mocker.call(_DISPLAY_TIMEOUT),
            ]

    def test_unresponsive(self, mock_multiprocessing):
        """Check the result when importing does not finish in time."""
        mocks = mock_multiprocessing(ready='none', alive=True)
        # pylint: disable-next=protected-access
        result = PyQtSanityChecker._check_binding(PYQT6)
        assert result is PyQtSanity.UNRESPONSIVE
        self.check_common_calls(mocks, PYQT6)
        mocks.process.terminate.assert_called_once_with()
        # Only after terminating, not to leave a zombie behind
        mocks.process.join.assert_called_once_with(_DISPLAY_TIMEOUT)


class TestCheckBindings:
    """Tests for checking the candidate bindings one at a time."""

    @fixture(name='mocks')
    def fixture_mocks(self, mocker):
        """Fake the candidates (all installed), and their checks."""
        return SimpleNamespace(
            candidates=mocker.patch(f'{MODULE}.get_candidates',
                                    return_value=BOTH),
            find_spec=mocker.patch('importlib.util.find_spec'),
            check=mocker.patch.object(PyQtSanityChecker, '_check_binding'),
            )

    def test_invalid_qt_api(self, mocks):
        """Check the results if QT_API prevents selecting a binding."""
        mocks.candidates.side_effect = QtApiError
        # pylint: disable-next=protected-access
        results = PyQtSanityChecker._check_bindings()
        assert results == ((None, PyQtSanity.INVALID_QT_API),)
        mocks.check.assert_not_called()

    def test_not_installed(self, mocks, mocker):
        """Check that no subprocess is started for a missing binding."""
        mocks.find_spec.side_effect = (
            lambda name: None if name == PYQT6 else mocker.MagicMock()
            )
        mocks.check.return_value = PyQtSanity.OK
        # pylint: disable-next=protected-access
        results = PyQtSanityChecker._check_bindings()
        assert results == ((PYQT6, PyQtSanity.NOT_FOUND),
                           (PYQT5, PyQtSanity.OK))
        mocks.check.assert_called_once_with(PYQT5)

    @parametrize(first=sorted(_STOP))
    def test_stops(self, first, mocks):
        """Check that no other binding is tried after some statuses."""
        mocks.check.return_value = first
        # pylint: disable-next=protected-access
        results = PyQtSanityChecker._check_bindings()
        assert results == ((PYQT6, first),)
        mocks.check.assert_called_once_with(PYQT6)

    @parametrize(first=sorted(_FALLBACK - {PyQtSanity.NOT_FOUND}))
    def test_tries_next(self, first, mocks):
        """Check that the next binding is tried after some failures."""
        mocks.check.side_effect = [first, PyQtSanity.OK]
        # pylint: disable-next=protected-access
        results = PyQtSanityChecker._check_bindings()
        assert results == ((PYQT6, first), (PYQT5, PyQtSanity.OK))


# The results of the checks, and the expected (status, binding)
# pylint: disable-next=consider-using-namedtuple-or-dataclass
_SUMMARIES = {
    'first OK': (((PYQT6, PyQtSanity.OK),),
                 (PyQtSanity.OK, PYQT6)),
    'fallback OK': (((PYQT6, PyQtSanity.RUNTIME_CRASH),
                     (PYQT5, PyQtSanity.OK)),
                    (PyQtSanity.OK, PYQT5)),
    'none installed': (((PYQT6, PyQtSanity.NOT_FOUND),
                        (PYQT5, PyQtSanity.NOT_FOUND)),
                       (PyQtSanity.NOT_FOUND, None)),
    'first unusable': (((PYQT6, PyQtSanity.IMPORT_ERROR),
                        (PYQT5, PyQtSanity.RUNTIME_CRASH)),
                       (PyQtSanity.IMPORT_ERROR, PYQT6)),
    'skip missing': (((PYQT6, PyQtSanity.NOT_FOUND),
                      (PYQT5, PyQtSanity.IMPORT_ERROR)),
                     (PyQtSanity.IMPORT_ERROR, PYQT5)),
    'last stops': (((PYQT6, PyQtSanity.RUNTIME_CRASH),
                    (PYQT5, PyQtSanity.NO_DISPLAY)),
                   (PyQtSanity.NO_DISPLAY, PYQT5)),
    'invalid QT_API': (((None, PyQtSanity.INVALID_QT_API),),
                       (PyQtSanity.INVALID_QT_API, None)),
    }


@parametrize('results,expect', _SUMMARIES.values(), ids=_SUMMARIES)
def test_summarize(results, expect):
    """Check the status and binding of all results."""
    # pylint: disable-next=protected-access
    assert PyQtSanityChecker._summarize(results) == expect
