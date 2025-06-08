"""Tests for module bookkeeper of viperleed.calc.bookkeeper.

Collects tests for simulating the default execution of bookkeeper
before and after viperleed.calc.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2024-11-15'
__license__ = 'GPLv3+'

import pytest
from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.calc.bookkeeper.bookkeeper import Bookkeeper
from viperleed.calc.bookkeeper.exit_code import BookkeeperExitCode
from viperleed.calc.bookkeeper.mode import BookkeeperMode
from viperleed.calc.constants import DEFAULT_HISTORY
from viperleed.calc.classes.rparams.rparams import Rparams
from viperleed.calc.cli import ViPErLEEDCalcCLI
from viperleed.calc.lib.context import execute_in_dir

from ..conftest import MOCK_TIMESTAMP
from .run_bookkeeper_base import _TestBookkeeperRunBase


class TestBookkeeperDuringCalc(_TestBookkeeperRunBase):
    """Test the same conditions as when bookkeeper runs around calc."""

    def test_run_around_calc(self,
                             mock_tree_before_calc_execution,
                             mock_tree_after_calc_execution,
                             mocker,
                             caplog):
        """Check reuse of bookkeeper in the default calls around calc."""
        tmp_path = mock_tree_before_calc_execution()
        bookkeeper = Bookkeeper(cwd=tmp_path)
        # Before calc, we run in CLEAR mode
        self.run_before_calc_exec_and_check(bookkeeper,
                                            caplog,
                                            check_archiving_required=False,
                                            mode='clear')
        # Then, we simulate a calc run that produces some output
        mock_tree_after_calc_execution()
        # After calc, we run in ARCHIVE mode
        after_calc = (
            bookkeeper,
            bookkeeper.cwd / DEFAULT_HISTORY / f't004.r001_{MOCK_TIMESTAMP}',
            mocker,
            )
        self.run_archive_after_calc_and_check(after_calc,
                                              caplog,
                                              check_archiving_required=False)

    def test_run_calc_no_input_files(self, tmp_path):
        """Check results when calc fails because of missing inputs."""
        cli = ViPErLEEDCalcCLI()
        with execute_in_dir(tmp_path):
            calc_error = cli([])
            assert calc_error
            bookkeeper = Bookkeeper()
            exit_code = bookkeeper.run('fix')
            assert exit_code is BookkeeperExitCode.NOTHING_TO_DO


class TestBookkeeperWhenCalcFails(_TestBookkeeperRunBase):
    """Test bookkeeper behavior when calc fails to execute (See #353)."""

    @fixture(name='run_calc_cli')
    def fixture_calc_cli(self, mock_tree_before_calc_execution):
        """Execute calc in a tree."""
        tmp_path = mock_tree_before_calc_execution()
        def _run():
            cli = ViPErLEEDCalcCLI()
            with execute_in_dir(tmp_path):
                cli([])
            return tmp_path
        return _run

    def check_has_unlabeled_and_ori(self, tmp_path, mocker):
        """Ensure that `tmp_path` contains both unlabeled and _ori files."""
        bookkeeper = mocker.MagicMock(cwd=tmp_path)
        self.check_root_inputs_untouched(bookkeeper)
        self.check_root_inputs_renamed_to_ori(bookkeeper)

    def test_exits_on_halting(self, run_calc_cli, mocker, caplog):
        """Test behavior when calc exits because of HALTING."""
        _module = 'viperleed.calc.run'
        rpars = Rparams()
        mocker.patch(f'{_module}._make_rpars_and_slab',
                     return_value=(rpars, mocker.MagicMock(name='slab')))
        mocker.patch(f'{_module}._set_tensorleed_source')
        mocker.patch(f'{_module}._set_system_name')
        rpars.halt = 5
        assert rpars.halt > rpars.HALTING
        tmp_path = run_calc_cli()
        # pylint: disable-next=magic-value-comparison
        assert 'Halting' in caplog.text
        self.check_has_unlabeled_and_ori(tmp_path, mocker)

    def test_fails_to_read_inputs(self, run_calc_cli, mocker):
        """Test behavior when calc exits early because of errors in inputs."""
        tmp_path = run_calc_cli()  # Mocked input files are invalid
        self.check_has_unlabeled_and_ori(tmp_path, mocker)

    _cli_raises = {
        'cant make work': ('_make_work_directory', OSError),
        'cant copy tensors': ('_copy_tensors_and_deltas_to_work', OSError),
        'cant copy inputs': ('_copy_input_files_to_work', OSError),
        }

    @parametrize('helper,exc', _cli_raises.values(), ids=_cli_raises)
    def test_no_call_on_cli_exception(self, helper, exc, run_calc_cli, mocker):
        """Ensure that Bookkeeper is not run if calc raises."""
        _module = 'viperleed.calc.cli'
        mock_bookkeeper = mocker.MagicMock(spec=Bookkeeper)
        mocker.patch(f'{_module}.Bookkeeper', return_value=mock_bookkeeper)
        mocker.patch(f'{_module}.{helper}', side_effect=exc)
        with pytest.raises(exc):
            run_calc_cli()
        clear_before = mocker.call(mode=BookkeeperMode.CLEAR,
                                   requires_user_confirmation=True)
        valid_calls = (
            [clear_before],  # Only clear before running
            [],      # No call at all if can't make work
            )
        assert mock_bookkeeper.run.mock_calls in valid_calls
