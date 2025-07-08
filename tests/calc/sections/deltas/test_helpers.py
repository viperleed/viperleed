"""Tests for module deltas of viperleed.calc.sections.

This module collects tests for various helper functions
in the deltas module.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-07-08'
__license__ = 'GPLv3+'

import logging

from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.calc.lib.context import execute_in_dir
from viperleed.calc.sections.deltas import _prepare_log_file

LOG_HEADER = ('Logs from multiple delta calculations are collected here. '
              'Their order may not be preserved.\n')


class TestPrepareLogFile:
    """Tests for the _prepare_log_file helper."""

    def check_file_contents(self, file_path):
        """Ensure the delta log file has the expected contents."""
        contents = file_path.read_text(encoding='utf-8')
        assert contents == LOG_HEADER

    @fixture(name='mock_rpars')
    def fixture_rpars(self, mocker):
        """Return a fake Rparams."""
        return mocker.MagicMock(timestamp='20240625')

    @fixture(name='call')
    def fixture_call(self, mock_rpars, tmp_path):
        """Call _prepare_log_file in a temporary directory."""
        def _call(**kwargs):
            with execute_in_dir(tmp_path):
                result = _prepare_log_file(mock_rpars, **kwargs)
            assert result == f'delta-{mock_rpars.timestamp}.log'
            return tmp_path/result
        return _call

    @parametrize(subdomain=(True, False))
    def test_sucess(self, subdomain, call, caplog):
        """Check successful file creation and logging messages."""
        caplog.set_level(logging.INFO)
        log_file = call(subdomain=subdomain)
        self.check_file_contents(log_file)
        logging_msgs = (
            'Generating delta files...',
            'collected in delta-20240625.log',
            )
        if subdomain:
            assert not any(msg in caplog.text for msg in logging_msgs)
        else:
            assert all(msg in caplog.text for msg in logging_msgs)

    def test_oserror_warns(self, call, mocker, caplog):
        """Check that failure to create the log file does not raise."""
        caplog.set_level(logging.WARNING)
        mocker.patch('pathlib.Path.write_text', side_effect=OSError)
        call(subdomain=False)
        expect_log = (
            'Error creating delta log file',
            'will not affect execution',
            )
        assert all(msg in caplog.text for msg in expect_log)
