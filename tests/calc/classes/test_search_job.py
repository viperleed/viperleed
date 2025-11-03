"""Test for viperleed.calc.classes.SearchJob."""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-07-17'
__license__ = 'GPLv3+'

import sys
import time

import pytest
from pytest_cases import fixture

from viperleed.calc.classes.search_job import COMMAND_NOT_FOUND_RETURN_CODE
from viperleed.calc.classes.search_job import KILLED_EXIT_CODE
from viperleed.calc.classes.search_job import SearchJob



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
