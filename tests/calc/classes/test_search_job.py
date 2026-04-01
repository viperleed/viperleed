"""Test for viperleed.calc.classes.SearchJob."""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-07-17'
__license__ = 'GPLv3+'

from contextlib import contextmanager
import multiprocessing
import subprocess
import sys
import time

import psutil
import pytest
from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.calc.classes.search_job import COMMAND_NOT_FOUND_RETURN_CODE
from viperleed.calc.classes.search_job import IS_WINDOWS
from viperleed.calc.classes.search_job import KILLED_EXIT_CODE
from viperleed.calc.classes.search_job import SearchJob
from viperleed.calc.classes.search_job import SearchWorkerABC
from viperleed.calc.classes.search_job import _psutil_kill_process_tree


@pytest.mark.timeout(5)
class TestSearchJob:
    """Tests for the SearchJob class."""

    @fixture(name='pyjob')
    def fixture_pyjob(self):
        """Return a started SearchJob running python code."""
        def _make(code, **kwargs):
            job = SearchJob([sys.executable, '-c', code], **kwargs)
            job.start()
            return job
        return _make

    def test_termination(self, pyjob, subtests):
        """Check clean termination of an infinite SearchJob."""
        job = pyjob('while True: pass', log_path=None)
        time.sleep(0.5)  # give it time to start

        with subtests.test('is running'):
            assert job.is_running()

        job.terminate()
        job.wait()  # ensure cleanup

        with subtests.test('is not running after termination'):
            assert not job.is_running()
        assert job.returncode == KILLED_EXIT_CODE

    def test_invalid_command(self):
        """Check return code of a job with an invalid target."""
        job = SearchJob(['nonexistent_command_xyz'], '', log_path=None)
        job.start()
        job.wait()
        assert job.returncode == COMMAND_NOT_FOUND_RETURN_CODE

    def test_correct_return_code(self, pyjob):
        """Test that the job returns the correct return code."""
        exit_code = 42
        code = f'''
import sys
import time

# We have to wait a bit, otherwise the
# process dies before we can fetch the PID
time.sleep(0.5)
sys.exit({exit_code})
'''
        job = pyjob(code, log_path=None)
        job.wait()
        assert job.returncode == exit_code

    def test_stdout_and_stderr_to_log_file(self, pyjob, tmp_path, subtests):
        """Both stdout and stderr should be redirected into the log file.

        Note that we explicitly test for UTF-8 characters. To ensure
        compatibility on Windows, we reconfigure the sys streams to use
        UTF-8 encoding. This is necessary as the tests are run in a
        non-console environment, which, on Windows, defaults to the 'locale'
        (ANSI) encoding.
        """
        log_file = tmp_path / 'search_job.log'
        messages = {
            'stdout': 'test log output with utf-8 😄',
            'stderr': 'test error output with utf-8 😡',
            }
        code = f'''
import sys

# Reconfigure sys.stdout and sys.stderr and print
sys.stdout.reconfigure(encoding='utf-8')
sys.stderr.reconfigure(encoding='utf-8')

print({messages['stdout']!r})
print({messages['stderr']!r}, file=sys.stderr)
'''
        job = pyjob(code, log_path=log_file)
        job.wait()

        content = log_file.read_text(encoding='utf-8')

        for stream, text in messages.items():
            with subtests.test(stream=stream):
                assert text in content

    def test_job_kill_flag(self, pyjob, subtests):
        """Test that the job can be killed using the kill flag."""
        job = pyjob('while True: pass', log_path=None)
        time.sleep(0.5)  # give it time to start

        with subtests.test('is running'):
            assert job.is_running()
        job._kill_me_flag.value = True
        job.wait()  # ensure cleanup

        with subtests.test('is not running after kill'):
            assert not job.is_running()

    def test_notice_immediate_job_death(self, pyjob):
        """Test that the job is not running immediately after start."""
        job = pyjob('pass', log_path=None)
        job.wait()

        assert job.returncode == 0
        assert job._mp_proc.exitcode == 0

    @pytest.mark.timeout(10)  # .join waits 5 seconds
    def test_hard_to_kill(self, mocker):
        """Ensure terminate() kills a hard-to-die process."""
        # Make a process that blocks on an Event
        evt = multiprocessing.Event()
        job = SearchJob('dummy cmd')
        job._mp_proc = multiprocessing.Process(target=evt.wait)
        job._mp_proc.start()

        assert job.is_running()

        # Terminate should now enter the is_alive() branch
        job.terminate()
        job.wait()

        assert not job.is_running()



class TestSearchJobNotStarted:
    """Tests for SearchJob before .start is called."""

    def test_wait(self):
        """Check that waiting on a non-started job is valid."""
        job = SearchJob('')
        job.wait()
        assert not job.is_running()

    def test_terminate(self):
        """Check that waiting on a non-started job is valid."""
        job = SearchJob('')
        job.terminate()
        assert not job.is_running()
        assert job._kill_me_flag.value


class TestSearchWorkerABC:
    """Tests for the SearchWorkerABC class."""

    @fixture(name='worker')
    def fixture_worker(self, mocker):
        """Return a concrete SearchWorkerABC instance."""
        mock_proc = mocker.MagicMock(started=True, stopped=False)
        class _Concrete(SearchWorkerABC):
            def start_subprocess(self_, *_):
                return mock_proc
            def _terminate_subprocess(self_):
                mock_proc.stopped = True
        dummy_args = {
            'command': mocker.MagicMock(),
            'input_data': None,
            'log_path': None,
            'kill_flag': mocker.MagicMock(),
            'return_code': mocker.MagicMock(),
            }
        instance = _Concrete(**dummy_args)
        instance.mock_proc = mock_proc
        return instance

    @fixture(name='stdin_worker')
    def fixture_worker_with_stdin(self, worker):
        """Return a concrete SearchWorkerABC that needs piped STDIN."""
        worker.input_data = 'STDIN data'
        return worker

    _methods_raise = (
        'start_subprocess',
        '_feed_inputs',
        '_monitor_and_wait',
        )

    @parametrize(method=_methods_raise)
    def test_keyboard_interrupt(self, method, worker, mocker):
        """Check termination when running is interrupted."""
        mocker.patch.object(worker, method, side_effect=KeyboardInterrupt)
        worker.run()
        proc = worker._proc
        if proc:
            assert proc.stopped
        assert worker.return_code.value == KILLED_EXIT_CODE

    @parametrize(method=_methods_raise)
    def test_fails_to_run(self, method, worker, mocker):
        """Check termination when running is interrupted."""
        mocker.patch.object(worker, method, side_effect=Exception)
        with pytest.raises(Exception):
            worker.run()
        proc = worker._proc
        if proc:
            assert proc.stopped
        assert worker.return_code.value == KILLED_EXIT_CODE

    @parametrize(exc=(OSError, subprocess.SubprocessError))
    def test_feed_input_fails(self, exc, stdin_worker, caplog):
        """Ensure complaints when feeding inputs fails."""
        proc = stdin_worker.mock_proc.__enter__.return_value
        proc.communicate.side_effect = exc
        with pytest.raises(exc):
            stdin_worker.run()
        assert stdin_worker.return_code.value == KILLED_EXIT_CODE
        assert caplog.text

    def test_feed_input_success(self, stdin_worker):
        """Ensure successful feeding of data to STDIN."""
        stdin_worker.run()
        proc = stdin_worker._proc
        assert stdin_worker.return_code.value is proc.returncode

    def test_feed_input_timeout(self, stdin_worker):
        """Ensure successful feeding of data to STDIN."""
        proc = stdin_worker.mock_proc.__enter__.return_value
        # subprocess.TimeoutExpired expects two positional args
        exc = subprocess.TimeoutExpired('cmd', 'timeout')
        proc.communicate.side_effect = exc
        self.test_feed_input_success(stdin_worker)


@pytest.mark.timeout(5)
class TestPsutilKillProcessTree:
    """Tests for the _psutil_kill_process_tree function."""

    @fixture(name='pyproc')
    def fixture_pyproc(self):
        """Return a started subprocess that uses the current python."""
        @contextmanager
        def _make(code):
            cmd = [sys.executable, '-c', code]
            with subprocess.Popen(cmd) as proc:
                time.sleep(0.5)  # ensure it starts
                yield proc
        return _make

    @pytest.mark.skipif(IS_WINDOWS,
                        reason='Does .terminate() instead of .kill()')
    def test_hard_to_kill(self, pyproc):
        """Test that a process ignoring terminate() is killed via kill()."""
        code = '''
import signal, time
def ignore(*_):
    return None
signal.signal(signal.SIGTERM, ignore)
time.sleep(10)
'''
        with pyproc(code) as proc:
            ps_proc = psutil.Process(proc.pid)
            assert ps_proc.is_running()
            _psutil_kill_process_tree(proc.pid)
            assert not ps_proc.is_running()

    def test_parent_without_children(self, pyproc):
        """Test graceful termination of a parent process only, no children."""
        with pyproc('import time; time.sleep(10)') as proc:
            _psutil_kill_process_tree(proc.pid)
            time.sleep(0.2)  # Give OS a moment to terminate
            assert proc.poll() is not None  # process should be gone

    def test_parent_with_children(self, pyproc):
        """Test graceful termination of a parent process with children."""
        parent_code = '''
import subprocess, time, sys
sub = subprocess.Popen([sys.executable, '-c', 'import time; time.sleep(5)'])
time.sleep(5)
'''
        with pyproc(parent_code) as proc:
            parent_ps = psutil.Process(proc.pid)
            children_before = parent_ps.children(recursive=True)
            _psutil_kill_process_tree(proc.pid)
            assert not parent_ps.is_running()
            for child in children_before:
                assert not child.is_running()

    def test_parent_with_some_children(self, pyproc):
        """Test graceful termination when only some children survive."""
        parent_code = f'''
import subprocess, time, sys
sub = subprocess.Popen([sys.executable, '-c', 'import time; time.sleep(5)'])
sub_dies_early = subprocess.Popen([
    sys.executable, '-c', 'import time; time.sleep(1)'
    ])
time.sleep(5)
'''
        with pyproc(parent_code) as proc:
            parent_ps = psutil.Process(proc.pid)
            children_before = parent_ps.children(recursive=True)
            _psutil_kill_process_tree(proc.pid)
            assert not parent_ps.is_running()
            for child in children_before:
                assert not child.is_running()
