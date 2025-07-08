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
def fixture_call_in_tmp(rpars, tmp_path):
    """Call deltas in a temporary directory."""
    def _call(**kwargs):
        with execute_in_dir(tmp_path):
            return deltas('mock_slab', rpars, **kwargs)
    return _call


@fixture(name='mock_atoms_need_deltas')
def fixture_mock_atoms_need_deltas(mocker):
    """Replace the _find_atoms_that_need_deltas function with a mock."""
    atoms_todo = [
        # Fake Atom objects with bare-bone attributes
        mocker.MagicMock(name=f'Fake Atom {i}', known_deltas=[], num=i)
        for i in range(6)
        ]
    atoms_todo[0].disp_occ = {'ElementOne': None,
                              'ElementTwo': None}
    atoms_todo[1].disp_occ = {'OnlyOneElement': None}
    atoms_todo[2].disp_occ = {'ElementOne': None,
                              'ElementTwo': None,
                              'ElementThree': None}
    atoms_todo[3].disp_occ = {'ElementOne': None}
    atoms_todo[4].disp_occ = {'OnlyOneElement': None}
    atoms_todo[5].disp_occ = {'OnlyOneElement': None}

    at_el_todo = [
        (atoms_todo[0], 'ElementOne'),
        (atoms_todo[0], 'ElementTwo'),
        (atoms_todo[1], 'OnlyOneElement'),
        (atoms_todo[2], 'ElementOne'),
        (atoms_todo[2], 'ElementTwo'),
        (atoms_todo[2], 'ElementThree'),
        (atoms_todo[2], 'Vac'),   # Auto-added
        (atoms_todo[3], 'Vac'),   # Auto-added
        (atoms_todo[3], 'ElementOne'),
        (atoms_todo[4], 'OnlyOneElement'),
        (atoms_todo[5], 'OnlyOneElement'),
        ]
    at_with_vacancies = [atoms_todo[2], atoms_todo[3]]
    return mocker.patch(
        f'{_MODULE}._find_atoms_that_need_deltas',
        return_value=(atoms_todo, at_el_todo, at_with_vacancies),
        )


@fixture(name='mocks')
def fixture_mock_implementation(mocker):
    """Replace implementation details of the deltas function with mocks."""
    return {
        'read_disp': mocker.patch(f'{_MODULE}.readDISPLACEMENTS_block'),
        'fetch tensor': mocker.patch(f'{_MODULE}._ensure_tensors_loaded'),
        'fetch deltas': mocker.patch(f'{_MODULE}.leedbase.getDeltas'),
        'rmtree': mocker.patch('shutil.rmtree'),
        'is_dir': mocker.patch('pathlib.Path.is_dir', return_value=True),
        'os.listdir': mocker.patch('os.listdir', return_value=[]),
        'collect_inputs': mocker.patch(
            f'{_MODULE}.iodeltas.collect_static_input_files',
            return_value=(mocker.MagicMock(),)*3,
            ),
        'delta input': mocker.patch(
            f'{_MODULE}.iodeltas.generateDeltaInput',
            side_effect=(('din', 'short', f'param {i}') for i in range(99999)),
            ),
        'check delta': mocker.patch(f'{_MODULE}.iodeltas.checkDelta',
                                    return_value=False),
        'checksums': mocker.patch(f'{_MODULE}.validate_multiple_files'),
        'pool': mocker.patch(f'{_MODULE}.parallelization.monitoredPool'),
        'copy_log': mocker.patch(f'{_MODULE}.leedbase.copy_compile_log'),
        'remove_param': mocker.patch(f'{_MODULE}._remove_old_param_file'),
        }


@fixture(name='no_deltas_to_do')
def fixture_no_atoms_need_deltas(mocker):
    """Replace the _find_atoms_that_need_deltas function with a mock."""
    return mocker.patch(f'{_MODULE}._find_atoms_that_need_deltas',
                        return_value=([], [], []))


class TestDeltasCalls:
    """Tests for implementation calls by the deltas function."""

    @use('mock_atoms_need_deltas')
    def test_called_for_subdomain(self, rpars, mocks, call_in_tmp, mocker):
        """Check differences between calls as "main" or as subdomain."""
        rpars.SUPPRESS_EXECUTION = True  # Skipped for subdomain
        rpars.FORTRAN_COMP[0] = ''
        log_info = mocker.patch(f'{_MODULE}.logger.info')
        result = call_in_tmp(subdomain=True)
        assert result is not None
        log_info.assert_not_called()
        rpars.getFortranComp.assert_not_called()
        rpars.setHaltingLevel.assert_not_called()
        rpars.updateCores.assert_not_called()
        assert DEFAULT_DELTAS in rpars.manifest
        not_called = (
            'checksums',
            'pool',
            'copy_log',
            )
        for mock_name in not_called:
            mocks[mock_name].assert_not_called()
        mocks['remove_param'].assert_called_once()
        mocks['fetch tensor'].assert_called_once()

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
        result = call_in_tmp(subdomain=False)
        assert result is None
        rpars.updateCores.assert_called_once()
        assert DEFAULT_DELTAS in rpars.manifest
        n_parallel_calls = 2  # compile, then run
        assert mocks['pool'].call_count == n_parallel_calls
        mocks['copy_log'].assert_called()
        mocks['rmtree'].assert_called()
        mocks['remove_param'].assert_called_once()
        mocks['fetch tensor'].assert_called_once()
        expect_log = (
            r'Generating delta files\.\.\.',
            r'Delta log will be written to local subfolders, and collected in',
            r'Compiling fortran files\.\.\.',
            r'Running delta calculations\.\.\.',
            r'Delta calculations finished\.',
            )
        log = caplog.text
        assert all(re.search(pattern, log) for pattern in expect_log)

    def test_domains_main(self, rpars, mocker):
        """Check calls when deltas is called in a multi-domain calculation."""
        rpars.domainParams = ['Domain 1']
        domains_impl = mocker.patch(f'{_MODULE}.deltas_domains')
        result = deltas(mocker.MagicMock(name='slab'), rpars)
        domains_impl.assert_called_once_with(rpars)
        assert result is None

    def test_prerun_calls_when_run_after_refcalc(self,
                                                 rpars,
                                                 no_deltas_to_do,
                                                 mocks,
                                                 mocker):
        """Check expected calls when deltas are calculated after refcalc."""
        rpars.disp_block_read = True
        mock_slab = mocker.MagicMock(name='slab')
        mocks['find_varied_atoms'] = no_deltas_to_do
        deltas(mock_slab, rpars)
        not_called = (
            'read_disp',
            'delta input',    # No deltas to calculate
            'remove_param',
            )
        called = {
            'fetch tensor': mocker.call(mock_slab, rpars),
            'fetch deltas': mocker.call(rpars.TENSOR_INDEX, required=False),
            'collect_inputs': mocker.call(mock_slab, rpars),
            'find_varied_atoms': mocker.call(mock_slab, rpars),
            }
        for mock_name in not_called:
            mocks[mock_name].assert_not_called()
        for mock_name, expect_call in called.items():
            assert mocks[mock_name].mock_calls == [expect_call]

        # Ensure no log file is present in the current directory
        assert not any(Path.cwd().glob('*.log'))

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
    def test_run_without_refcalc(self, rpars, mocks):
        """Check calls when deltas execute as the first segment."""
        rpars.runHistory = [0, 2, 3, 12, 99, 57]  # No refcalc (== 1)
        deltas('mock_slab', rpars)
        mocks['fetch tensor'].assert_called_once_with('mock_slab', rpars)

    def test_stop_after_compile(self,
                                rpars,
                                call_in_tmp,
                                mock_atoms_need_deltas,
                                mocks):
        """Check expected calls when deltas are only compiled."""
        def _set_stop(*_, **__):
            rpars.STOP = True

        mocks['pool'].side_effect = _set_stop
        call_in_tmp()
        mocks['pool'].assert_called_once()
        *static_args, tasks = mocks['pool'].call_args[0]
        n_tasks = len(mock_atoms_need_deltas.return_value[1])
        assert static_args == [rpars, rpars.N_CORES, compile_delta]
        assert len(tasks) == n_tasks
        mocks['remove_param'].assert_called_once()

    def test_stop_after_run(self,
                            rpars,
                            call_in_tmp,
                            mock_atoms_need_deltas,
                            mocks):
        """Check expected calls when deltas are only compiled."""
        def _set_stop_after_run(*args, **__):
            if args[2] is run_delta:
                rpars.STOP = True

        mocks['pool'].side_effect = _set_stop_after_run
        call_in_tmp()
        n_parallel_calls = 2  # compile, then run
        assert mocks['pool'].call_count == n_parallel_calls
        *static_args, tasks = mocks['pool'].call_args[0]
        n_tasks = len(mock_atoms_need_deltas.return_value[1])
        assert static_args == [rpars, rpars.N_CORES, run_delta]
        assert len(tasks) == n_tasks
        mocks['rmtree'].assert_not_called()
        mocks['copy_log'].assert_not_called()
        mocks['remove_param'].assert_called_once()

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
    def test_fails_to_write_delta_input(self, rpars, mocker, caplog):
        """Check that failure to write the input file is tolerated."""
        input_file_name = 'delta-input'
        def _fail_to_open_delta_input(fname, *_, **__):
            if fname == input_file_name:
                # pylint: disable-next=broad-exception-raised
                raise Exception('open failed')
        mocker.patch('builtins.open', side_effect=_fail_to_open_delta_input)
        deltas('mock_slab', rpars)
        expect_log = f'Failed to write file {input_file_name!r}'
        assert expect_log in caplog.text

    @use('mocks')
    def test_raises_halting_level(self, rpars, call_in_tmp):
        """Check setting of halting level upon SUPPRESS_EXECUTION."""
        call_in_tmp()
        rpars.setHaltingLevel.assert_called_once_with(3)

    @use('mock_atoms_need_deltas')
    def test_static_input_files_used(self, mocks, call_in_tmp):
        """Check that input files are passed on to the generation of inputs."""
        call_in_tmp()
        gen_delta_static_args = [
            args[0][-3:]  # The last three positional arguments
            for args in mocks['delta input'].call_args_list
            ]
        assert all(args == gen_delta_static_args[0]
                   for args in gen_delta_static_args)
        assert gen_delta_static_args[0] == mocks['collect_inputs'].return_value

    def test_generates_tasks_based_on_hash(self,
                                           mocks,
                                           call_in_tmp,
                                           tmp_path,
                                           mocker):
        """Check that compile tasks are created using the hash of PARAM."""
        ids = (0, 1, 2, 2, 3, 5, 7, 8, 7, 3, 2)
        mocks['delta input'].side_effect = tuple(
            ('din', f'short {i}', f'param {i}')
            for i in ids
            )

        instance_counts = {}
        def count_instance_creations(klass):
            old_init = klass.__init__
            instance_counts[klass] = 0
            def _mock_init(obj, *args, **kwargs):
                instance_counts[klass] += 1
                return old_init(obj, *args, **kwargs)
            return _mock_init

        mocker.patch(f'{_MODULE}.DeltaCompileTask.__init__',
                     count_instance_creations(DeltaCompileTask))
        mocker.patch(f'{_MODULE}.DeltaRunTask.__init__',
                     count_instance_creations(DeltaRunTask))
        call_in_tmp()
        assert mocks['delta input'].call_count == len(ids)
        assert instance_counts[DeltaCompileTask] == len(set(ids))
        assert instance_counts[DeltaRunTask] == len(ids)

        # Check the delta input file
        input_contents = (tmp_path/'delta-input').read_text(encoding='utf-8')
        expect_contents = '''\
# ABOUT THIS FILE:
# Input for the delta-calculations is collected here. The blocks of data are
# new 'PARAM' files, which are used to recompile the fortran code, and input
# for generation of specific DELTA files. Lines starting with '#' are comments
# on the function of the next block of data.
# In the DELTA file blocks, [AUXBEAMS] and [PHASESHIFTS] denote where the
# entire contents of the AUXBEAMS and PHASESHIFTS files should be inserted.

#### NEW 'PARAM' FILE: ####

param 0

#### INPUT for new DELTA file DEL_0_ElementOne_1: ####

short 0

#### NEW 'PARAM' FILE: ####

param 1

#### INPUT for new DELTA file DEL_0_ElementTwo_1: ####

short 1

#### NEW 'PARAM' FILE: ####

param 2

#### INPUT for new DELTA file DEL_1_OnlyOneElement_1: ####

short 2

#### INPUT for new DELTA file DEL_2_ElementOne_1: ####

short 2

#### INPUT for new DELTA file DEL_5_OnlyOneElement_1: ####

short 2

#### NEW 'PARAM' FILE: ####

param 3

#### INPUT for new DELTA file DEL_2_ElementTwo_1: ####

short 3

#### INPUT for new DELTA file DEL_4_OnlyOneElement_1: ####

short 3

#### NEW 'PARAM' FILE: ####

param 5

#### INPUT for new DELTA file DEL_2_ElementThree_1: ####

short 5

#### NEW 'PARAM' FILE: ####

param 7

#### INPUT for new DELTA file DEL_2_Vac_1: ####

short 7

#### INPUT for new DELTA file DEL_3_ElementOne_1: ####

short 7

#### NEW 'PARAM' FILE: ####

param 8

#### INPUT for new DELTA file DEL_3_Vac_1: ####

short 8
'''
        assert input_contents == expect_contents

    @use('mocks')
    def test_writes_log_file(self, call_in_tmp, tmp_path):
        """Ensure that a delta log file is written."""
        call_in_tmp()
        assert any(tmp_path.glob('*.log'))


class TestDeltasRaises:
    """Tests for conditions that cause exceptions when deltas is called."""

    def test_compilation_fails(self,
                               rpars,
                               mocks,
                               mock_atoms_need_deltas,
                               call_in_tmp):
        """Check failure when compilation fails."""
        rpars.FORTRAN_COMP = ['failing compiler', '']
        mocks['pool'].side_effect = RuntimeError('compile fail')
        with pytest.raises(RuntimeError, match='compile fail'):
            call_in_tmp()
        n_tasks = len(mock_atoms_need_deltas.return_value[1])
        assert mocks['copy_log'].call_count == n_tasks

    @use('mock_atoms_need_deltas', 'mocks')
    def test_find_compiler_fails(self, rpars, call_in_tmp):
        """Check failure when finding a missing compiler fails."""
        rpars.FORTRAN_COMP[0] = ''
        rpars.getFortranComp.side_effect = Exception('no compiler')
        with pytest.raises(RuntimeError):
            call_in_tmp()
