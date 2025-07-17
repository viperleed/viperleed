"""Class SearchJob."""

__authors__ = ("Alexander M. Imre (@amimre)",)
__copyright__ = "Copyright (c) 2019-2025 ViPErLEED developers"
__created__ = "2025-07-17"
__license__ = "GPLv3+"

import pytest

import logging
import time

from viperleed.calc.classes.search_job import SearchJob

logger = logging.getLogger(__name__)


@pytest.mark.timeout(5)
def test_termination(subtests):
    script = ["python3", "-c", "while True: pass"]
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

def test_invalid_command():
    job = SearchJob(["nonexistent_command_xyz"], "", log_path=None)
    job.start()
    job.wait()
    assert job.returncode == 127  # mimics "command not found"
