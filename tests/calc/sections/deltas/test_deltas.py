"""Tests for module deltas of viperleed.calc.sections.

This module defines tests for the deltas function, the main entry point
of a delta-amplitudes calculation.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-07-04'
__license__ = 'GPLv3+'

import logging
from pathlib import Path
import re

import pytest
from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.calc.constants import DEFAULT_DELTAS
from viperleed.calc.lib.context import execute_in_dir
from viperleed.calc.sections.deltas import DeltaCompileTask
from viperleed.calc.sections.deltas import DeltaRunTask
from viperleed.calc.sections.deltas import compile_delta
from viperleed.calc.sections.deltas import deltas
from viperleed.calc.sections.deltas import run_delta

from .conftest import _MOCK_COMPILER
from .conftest import _MODULE

use = pytest.mark.usefixtures


@fixture(name='call_in_tmp')
def fixture_call_in_tmp(rpars, mocks, tmp_path, mocker):
    """Call deltas in a temporary directory."""
    def _call(**kwargs):
        args = mocker.MagicMock(name='slab'), rpars
        with execute_in_dir(tmp_path):
            result = deltas(*args, **kwargs)
        assert result is None
        mocks['fetch_deltas'].assert_called_once_with(rpars.TENSOR_INDEX,
                                                      required=False)
        mocks['fetch_tensor'].assert_called_once_with(*args)
        mocks['find_varied_atoms'].assert_called_once_with(*args)
        todo, at_el_todo, vacancies = mocks['find_varied_atoms'].return_value
        if at_el_todo:
            mocks['remove_param'].assert_called_once_with()
            mocks['make_tasks'].assert_called_once_with(
                *args,
                at_el_todo,
                f'delta-{rpars.timestamp}.log',
                )
            mocks['sort_deltas'].assert_called_once_with(todo, vacancies)
            mocks['write_input'].assert_called_once_with(
                *mocks['make_tasks'].return_value,
                )
        else:
            mocks['remove_param'].assert_not_called()
            mocks['make_tasks'].assert_not_called()
            mocks['sort_deltas'].assert_not_called()
            mocks['pool'].assert_not_called()
            mocks['write_input'].assert_not_called()
    return _call


@fixture(name='mock_atoms_need_deltas')
def fixture_mock_atoms_need_deltas(mocks, mocker):
    """Replace the _find_atoms_that_need_deltas function with a mock."""
    mocks['find_varied_atoms'] = mocker.patch(
        f'{_MODULE}._find_atoms_that_need_deltas',
        return_value=tuple(mocker.MagicMock() for _ in range(3)),
        )


@fixture(name='mocks')
def fixture_mock_implementation(rpars, mocker):
    """Replace implementation details of the deltas function with mocks."""
    tl_path = rpars.get_tenserleed_directory.return_value.path
    return {
        'read_disp': mocker.patch(f'{_MODULE}.readDISPLACEMENTS_block'),
        'fetch_tensor': mocker.patch(f'{_MODULE}._ensure_tensors_loaded'),
        'fetch_deltas': mocker.patch(f'{_MODULE}.leedbase.getDeltas'),
        'rmtree': mocker.patch('shutil.rmtree'),
        'is_dir': mocker.patch('pathlib.Path.is_dir', return_value=True),
        'os.listdir': mocker.patch('os.listdir', return_value=[]),
        'collect_inputs': mocker.patch(
            'viperleed.calc.files.iodeltas.collect_static_input_files',
            return_value=(mocker.MagicMock(),)*3,
            ),
        'make_delta_input': mocker.patch(
            'viperleed.calc.files.iodeltas.generateDeltaInput',
            side_effect=(('din', 'short', f'param {i}') for i in range(99999)),
            ),
        'check delta': mocker.patch(f'{_MODULE}.iodeltas.checkDelta',
                                    return_value=False),
        'checksums': mocker.patch(f'{_MODULE}.validate_multiple_files'),
        'pool': mocker.patch(f'{_MODULE}.parallelization.monitoredPool'),
        'copy_log': mocker.patch(f'{_MODULE}.leedbase.copy_compile_log'),
        'remove_param': mocker.patch(f'{_MODULE}._remove_old_param_file'),
        'sort_deltas': mocker.patch(
            f'{_MODULE}._sort_current_deltas_by_element',
            ),
        'write_input': mocker.patch(
            'viperleed.calc.files.iodeltas.write_delta_input_file',
            ),
        'make_tasks': mocker.patch(
            f'{_MODULE}._assemble_tasks',
            return_value=(
                [DeltaCompileTask('param', tl_path, 'index')],
                [DeltaRunTask('comptask')],
                ),
            )
        }


@fixture(name='no_deltas_to_do')
def fixture_no_atoms_need_deltas(mocks, mocker):
    """Replace the _find_atoms_that_need_deltas function with a mock."""
    mocks['find_varied_atoms'] = mocker.patch(
        f'{_MODULE}._find_atoms_that_need_deltas',
        return_value=([], [], []),
        )


class TestDeltasCalls:
    """Tests for implementation calls by the deltas function."""

    @use('mock_atoms_need_deltas')
    @parametrize(suppress=(True, False))
    # pylint: disable-next=too-many-arguments  # 4/6 fixtures
    def test_called_for_subdomain(self,
                                  suppress,
                                  rpars,
                                  mocks,
                                  call_in_tmp,
                                  mocker):
        """Check differences between calls as "main" or as subdomain."""
        rpars.SUPPRESS_EXECUTION = suppress  # Only determines manifest
        rpars.FORTRAN_COMP[0] = ''
        log_info = mocker.patch(f'{_MODULE}.logger.info')
        call_in_tmp(subdomain=True)
        log_info.assert_not_called()
        rpars.getFortranComp.assert_not_called()
        rpars.setHaltingLevel.assert_not_called()
        rpars.updateCores.assert_not_called()
        assert (DEFAULT_DELTAS in rpars.manifest) != suppress
        not_called = (  # Because deltas_domains does so
            'checksums',
            'pool',
            'copy_log',
            'rmtree',
            )
        for mock_name in not_called:
            mocks[mock_name].assert_not_called()

    @use('mock_atoms_need_deltas')
    @parametrize(enable=(True, False))
    def test_checksum_validation(self, enable, rpars, mocks, call_in_tmp):
        """Ensure validation of Fortran sources is done if requested."""
        rpars.TL_IGNORE_CHECKSUM = enable
        rpars.FORTRAN_COMP = ['some compiler', '']
        call_in_tmp()
        check_call = mocks['checksums'].getattr('assert_called_once' if enable
                                                else 'assert_not_called')
        check_call()

    @use('mock_atoms_need_deltas')
    def test_compile_and_run(self, rpars, mocks, call_in_tmp, caplog):
        """Check a full, clean execution for a single domain."""
        caplog.set_level(logging.INFO)
        call_in_tmp(subdomain=False)
        rpars.updateCores.assert_called_once()
        assert DEFAULT_DELTAS in rpars.manifest
        n_parallel_calls = 2  # compile, then run
        assert mocks['pool'].call_count == n_parallel_calls
        mocks['copy_log'].assert_called()
        mocks['rmtree'].assert_called()
        expect_log = (
            r'Generating delta files\.\.\.',
            r'Delta log will be written to local subfolders, and collected in',
            r'Compiling fortran files\.\.\.',
            r'Running delta calculations\.\.\.',
            r'Delta calculations finished\.',
            )
        log = caplog.text
        assert all(re.search(pattern, log) for pattern in expect_log)

    @use('no_deltas_to_do')
    def test_displacements_already_read(self, rpars, call_in_tmp, mocks):
        """Check expected calls when deltas are calculated after refcalc."""
        rpars.disp_block_read = True
        call_in_tmp()
        mocks['read_disp'].assert_not_called()

        # Ensure no log file is present in the current directory
        assert not any(Path.cwd().glob('*.log'))

    def test_domains_main(self, rpars, mocker):
        """Check calls when deltas is called in a multi-domain calculation."""
        rpars.domainParams = ['Domain 1']
        domains_impl = mocker.patch(f'{_MODULE}.deltas_domains')
        result = deltas(None, rpars)
        domains_impl.assert_called_once_with(rpars)
        assert result is None

    @use('no_deltas_to_do')
    def test_reads_most_recent_displacements_block(self, rpars, mocks, mocker):
        """Check that the most-recent DISPLACEMENTS block is read."""
        rpars.disp_block_read = False
        mock_slab = mocker.MagicMock(name='slab')
        deltas(mock_slab, rpars)
        mocks['read_disp'].assert_called_once_with(
            rpars,
            mock_slab,
            rpars.disp_blocks[rpars.search_index],
            )
        assert rpars.disp_block_read

    @use('no_deltas_to_do')
    def test_run_without_refcalc(self, rpars, call_in_tmp):
        """Check calls when deltas execute as the first segment."""
        rpars.runHistory = [0, 2, 3, 12, 99, 57]  # No refcalc (== 1)
        call_in_tmp()

    @staticmethod
    def _check_stopped_exec_calls(rpars,
                                  mocks,
                                  n_parallel_calls,
                                  pool_run,
                                  expect_tasks):
        """Ensure calls during execution are as expected."""
        assert mocks['pool'].call_count == n_parallel_calls
        args = mocks['pool'].call_args[0]
        size = min(rpars.N_CORES, len(expect_tasks))
        assert args == (rpars, size, pool_run, expect_tasks)
        mocks['rmtree'].assert_not_called()
        mocks['copy_log'].assert_not_called()

    @use('mock_atoms_need_deltas')
    def test_stop_after_compile(self, rpars, call_in_tmp, mocks):
        """Check expected calls when deltas are only compiled."""
        def _set_stop(*_, **__):
            rpars.STOP = True

        mocks['pool'].side_effect = _set_stop
        call_in_tmp()
        self._check_stopped_exec_calls(
            rpars,
            mocks,
            n_parallel_calls=1,  # compile only
            pool_run=compile_delta,
            expect_tasks=mocks['make_tasks'].return_value[0],
            )

    @use('mock_atoms_need_deltas')
    def test_stop_after_run(self, rpars, call_in_tmp, mocks):
        """Check expected calls when deltas are only compiled."""
        def _set_stop_after_run(*args, **__):
            if args[2] is run_delta:
                rpars.STOP = True

        mocks['pool'].side_effect = _set_stop_after_run
        call_in_tmp()
        self._check_stopped_exec_calls(
            rpars,
            mocks,
            n_parallel_calls=2,  # compile, then run
            pool_run=run_delta,
            expect_tasks=mocks['make_tasks'].return_value[1],
            )

    @use('mock_atoms_need_deltas', 'mocks')
    def test_triggers_compiler_discovery_if_missing(self, rpars, call_in_tmp):
        """Ensure that a compiler is looked up if not available yet."""
        call_in_tmp()
        assert rpars.FORTRAN_COMP[0] is _MOCK_COMPILER


@use('mock_atoms_need_deltas')
class TestDeltasGenerateInput:
    """Tests for creation of all input files needed for delta calculations."""

    @fixture(autouse=True)
    def suppress_execution(self, rpars):
        """Ensure only inputs are created."""
        rpars.SUPPRESS_EXECUTION = True

    @use('mocks')
    def test_raises_halting_level(self, rpars, call_in_tmp):
        """Check setting of halting level upon SUPPRESS_EXECUTION."""
        call_in_tmp()
        rpars.setHaltingLevel.assert_called_once_with(3)

    @use('mocks')
    def test_writes_log_file(self, call_in_tmp, tmp_path):
        """Ensure that a delta log file is written."""
        call_in_tmp()
        assert any(tmp_path.glob('*.log'))


class TestDeltasRaises:
    """Tests for conditions that cause exceptions when deltas is called."""

    @use('mock_atoms_need_deltas')
    def test_compilation_fails(self,
                               rpars,
                               mocks,
                               call_in_tmp):
        """Check failure when compilation fails."""
        rpars.FORTRAN_COMP = ['failing compiler', '']
        mocks['pool'].side_effect = RuntimeError('compile fail')
        with pytest.raises(RuntimeError, match='compile fail'):
            call_in_tmp()
        n_tasks = len(mocks['make_tasks'].return_value[0])
        assert mocks['copy_log'].call_count == n_tasks

    @use('mock_atoms_need_deltas', 'mocks')
    def test_find_compiler_fails(self, rpars, call_in_tmp):
        """Check failure when finding a missing compiler fails."""
        rpars.FORTRAN_COMP[0] = ''
        rpars.getFortranComp.side_effect = Exception('no compiler')
        with pytest.raises(RuntimeError):
            call_in_tmp()


class TestExceptionsPropagated:
    """Test that exceptions in helpers are not caught."""

    with_exceptions = parametrize(exc=(Exception, BaseException))

    @staticmethod
    def check_propagates(exc, call_in_tmp, mocks, mock_name):
        """Ensure that exceptions in a given mock are propagated."""
        mocks[mock_name].side_effect = exc
        with pytest.raises(exc):
            call_in_tmp()

    @with_exceptions
    def test_deltas_domains(self, exc, rpars, mocker):
        """Check propagation of exceptions in deltas_domains."""
        mocker.patch(f'{_MODULE}.deltas_domains', side_effect=exc)
        rpars.domainParams = ['some domain']
        with pytest.raises(exc):
            deltas('mock_slab', rpars)

    @with_exceptions
    def test_read_displacements(self, exc, rpars, call_in_tmp, mocks):
        """Check propagation of exceptions in readDISPLACEMENTS_block."""
        rpars.disp_block_read = False
        self.check_propagates(exc, call_in_tmp, mocks, 'read_disp')

    _preliminary = (
        'fetch_tensor',
        'fetch_deltas',
        'find_varied_atoms',
        )
    _before_suppress = (
        *_preliminary,
        'remove_param',
        'make_tasks',
        'sort_deltas',
        'write_input',
        )
    _execution = (
        *_before_suppress,
        )

    @use('no_deltas_to_do')
    @with_exceptions
    @parametrize(mock_name=_preliminary)
    def test_preliminary_call(self, exc, mock_name, call_in_tmp, mocks):
        """Check propagation of exceptions in helpers called before tasks."""
        self.check_propagates(exc, call_in_tmp, mocks, mock_name)

    @use('mock_atoms_need_deltas')
    @with_exceptions
    @parametrize(mock_name=_before_suppress)
    # pylint: disable-next=too-many-arguments  # 3/6 fixtures
    def test_suppressed(self, exc, mock_name, rpars, call_in_tmp, mocks):
        """Ensure propagation of exceptions in helpers before early return."""
        rpars.SUPPRESS_EXECUTION = True
        self.check_propagates(exc, call_in_tmp, mocks, mock_name)

    @use('mock_atoms_need_deltas')
    @with_exceptions
    @parametrize(mock_name=_execution)
    # pylint: disable-next=too-many-arguments  # 3/6 fixtures
    def test_execution(self, exc, mock_name, rpars, call_in_tmp, mocks):
        """Ensure propagation of exceptions in helpers before early return."""
        rpars.SUPPRESS_EXECUTION = False
        self.check_propagates(exc, call_in_tmp, mocks, mock_name)
