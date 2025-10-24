"""Test for viperleed.calc.classes.SearchJob."""

__authors__ = ('Alexander M. Imre (@amimre)',)
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-07-17'
__license__ = 'GPLv3+'

import sys
import time

import pytest

from viperleed.calc.classes.search_job import (
    COMMAND_NOT_FOUND_RETURN_CODE,
    SearchJob,
)


@pytest.mark.timeout(5)
def test_termination(subtests):
    script = [sys.executable, '-c', 'while True: pass']
    job = SearchJob(script, '', log_path=None)
    job.start()
    time.sleep(0.5)  # give it time to start

    with subtests.test('is running'):
        assert job.is_running()

    job.terminate()
    job.wait()  # ensure cleanup

    with subtests.test('is not running after termination'):
        assert not job.is_running()
    assert job.returncode is not None


@pytest.mark.timeout(5)
def test_invalid_command():
    job = SearchJob(['nonexistent_command_xyz'], '', log_path=None)
    job.start()
    job.wait()
    assert job.returncode == COMMAND_NOT_FOUND_RETURN_CODE


@pytest.mark.timeout(5)
def test_correct_return_code():
    """Test that the job returns the correct return code."""
    script = [
        sys.executable,
        '-c',
        # we have to wait a bit, otherwise the process dies before we
        # can fetch the PID
        'import sys; import time; time.sleep(0.5); sys.exit(42)',
    ]
    job = SearchJob(script, '', log_path=None)
    job.start()
    job.wait()
    assert job.returncode == 42


@pytest.mark.timeout(5)
def test_writing_stdout_and_stderr_to_log_file(tmp_path, subtests):
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

    # reconfigure sys.stdout and sys.stderr and print
    code = (
        'import sys\n'
        'sys.stdout.reconfigure(encoding="utf-8")\n'
        'sys.stderr.reconfigure(encoding="utf-8")\n'
        f'print({messages["stdout"]!r})\n'
        f'print({messages["stderr"]!r}, file=sys.stderr)\n'
    )
    script = [sys.executable, '-c', code]

    job = SearchJob(script, '', log_path=log_file)
    job.start()
    job.wait()

    content = log_file.read_text(encoding='utf-8')

    for stream, text in messages.items():
        with subtests.test(stream=stream):
            assert text in content


@pytest.mark.timeout(5)
def test_job_kill_flag(subtests):
    """Test that the job can be killed using the kill flag."""
    script = [sys.executable, '-c', 'while True: pass']
    job = SearchJob(script, '', log_path=None)
    job.start()
    time.sleep(0.5)  # give it time to start

    with subtests.test('is running'):
        assert job.is_running()
    job._kill_me_flag.value = True
    job.wait()  # ensure cleanup

    with subtests.test('is not running after kill'):
        assert not job.is_running()


@pytest.mark.timeout(5)
def test_notice_immediate_job_death():
    """Test that the job is not running immediately after start."""
    script = [sys.executable, '-c', 'pass']
    job = SearchJob(script, '', log_path=None)
    job.start()
    job.wait()

    assert job.returncode == 0
    assert job._mp_proc.exitcode == 0
