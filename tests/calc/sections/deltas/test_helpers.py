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
from viperleed.calc.sections.deltas import _remove_old_param_file


class TestPrepareLogFile:
    """Tests for the _prepare_log_file helper."""

    def check_file_contents(self, file_path):
        """Ensure the delta log file has the expected contents."""
        log_header = ('Logs from multiple delta calculations are collected '
                      'here. Their order may not be preserved.\n')
        contents = file_path.read_text(encoding='utf-8')
        assert contents == log_header

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


class TestRemoveOldParamFile:
    """Tests for the _remove_old_param_file helper function."""

    @fixture(name='call')
    def fixture_call(self, tmp_path):
        """Call _remove_old_param_file in a temporary directory."""
        def _call():
            old_param = tmp_path/'PARAM'
            renamed = old_param.with_name('PARAM-old')
            contents = 'test'
            old_param.write_text(contents)
            with execute_in_dir(tmp_path):
                _remove_old_param_file()
            return old_param, renamed, contents
        return _call

    def test_deletes_param_if_rename_fails(self, call, mocker):
        """Check deletion of an existing PARAM when renaming fails."""
        mocker.patch('os.rename', side_effect=OSError)
        param, renamed, _ = call()
        assert not param.exists()
        assert not renamed.exists()

    def test_does_nothing_if_no_param_file(self, tmp_path):
        """Check that no action is take if PARAM does not exist."""
        with execute_in_dir(tmp_path):
            _remove_old_param_file()
        assert not (tmp_path / 'PARAM').exists()
        assert not (tmp_path / 'PARAM-old').exists()

    def test_renames_existing_param_file(self, call):
        """Check successful renaming of a PARAM file."""
        param, renamed, contents = call()
        assert not param.exists()
        assert renamed.exists()
        assert renamed.read_text() == contents

    def test_warns_if_remove_fails(self, call, mocker, caplog):
        """Check warnings are emitted if removal of PARAM fails."""
        mocker.patch('os.rename', side_effect=OSError)
        mocker.patch('os.remove', side_effect=OSError)
        param, *_ = call()
        assert 'Cannot rename/remove old PARAM file' in caplog.text
        assert param.exists()
