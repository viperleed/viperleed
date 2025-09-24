"""Class SearchJob."""

__authors__ = ('Alexander M. Imre (@amimre)',)
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-07-17'
__license__ = 'GPLv3+'

from pathlib import Path
import sys
import time

import pytest

from viperleed.calc.classes.search_job import SearchJob


@pytest.mark.timeout(5)
def test_termination(subtests):
    script = [sys.executable, '-c', 'while True: pass']
    job = SearchJob(script, "", log_path=None)
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
    assert job.returncode == 127  # mimics "command not found"

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
def test_writing_to_log_file(tmp_path):
    """Test that the job writes to the log file."""
    log_file = tmp_path / 'search_job.log'
    script = [sys.executable, '-c', 'print("Test log output")']
    job = SearchJob(script, '', log_path=log_file)
    job.start()
    job.wait()

    with open(log_file, 'r') as f:
        content = f.read()

    assert 'Test log output' in content

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
def test_notice_immediate_job_death(capfd):
    """Test that the job is not running immediately after start."""
    # Note: we HAVE to use capfd here (rather than caplog, capsys or even
    # "with pytest.raises" because the error occurs in another multiprocessing
    # process and will not be captured)
    script = [sys.executable, '-c', 'pass']
    job = SearchJob(script, '', log_path=None)
    job.start()
    job.wait()

    captured = capfd.readouterr()
    assert 'psutil.NoSuchProcess: process PID not found' in captured.err
