"""Tests for function search of viperleed.calc.sections.

This module contains tests for the implementation of the
search function of viperleed.calc.sections. They mock away
parts of the function that actually do execute a search.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-06-06'
__license__ = 'GPLv3+'

from pathlib import Path

import pytest
from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.calc.classes.rparams.rparams import Rparams
from viperleed.calc.lib.context import execute_in_dir
from viperleed.calc.lib.version import Version
from viperleed.calc.sections.search import search

from .conftest import _MODULE

use = pytest.mark.usefixtures


@fixture(name='mock_compile')
def fixture_mock_implementation_compilation(rpars,
                                            mock_no_run,
                                            tmp_path,
                                            mocker):
    """Replace details of function search till the compilation step."""
    mock_compile = {
        'tenserleed': mocker.patch.object(
            rpars,
            'get_tenserleed_directory',
            # Returns a files.tenserleed.TensErLEEDSource object
            return_value=mocker.MagicMock(path=tmp_path),
            ),
        'copy': mocker.patch('shutil.copy2'),
        'checksum': mocker.patch(f'{_MODULE}.validate_multiple_files'),
        'compile': mocker.patch(f'{_MODULE}.leedbase.fortran_compile_batch'),
        }
    return {**mock_no_run, **mock_compile}


@fixture(name='mock_no_run')
def fixture_mock_implementation_no_execution(mocker):
    """Replace details of function search for when no execution happens."""
    return {
        'read_disp': mocker.patch(f'{_MODULE}.readDISPLACEMENTS_block'),
        'get_deltas': mocker.patch(f'{_MODULE}.leedbase.getDeltas'),
        'write_rfinfo': mocker.patch(f'{_MODULE}.iosearch.writeRfInfo'),
        'write_input': mocker.patch(f'{_MODULE}.iosearch.generateSearchInput'),
        'write_output': mocker.patch(f'{_MODULE}.iosearch.writeSearchOutput'),
        }


@fixture(name='mock_run')
def fixture_mock_implementation_run(rpars, mock_compile, mocker):
    """Replace details of function search till the execution step."""
    def _make_sdtl(*_, **__):
        Path('SD.TL').touch()

    # Note: os calls are no longer necessary, as SearchJob
    # handles the process management via psutil and multiprocessing.

    # mock the timer to always expire immediately
    mock_eval_timer = mocker.MagicMock()
    mock_eval_timer.has_expired.return_value = True

    # mock the SearchJob, which handles the search process
    # create SD.TL when started, and then behaves as if it finished
    mock_search_job = mocker.MagicMock()
    mock_search_job.start.side_effect = _make_sdtl
    mock_search_job.is_running.side_effect = (True, False)

    # set required rpars attributes
    rpars.searchEvalTime = 0.05  # Seconds
    rpars.output_interval = 1    # Generations
    mock_run = {
        'sleep': mocker.patch('time.sleep'),  # Make tests faster
        'results': mocker.patch(f'{_MODULE}.processSearchResults'),
        'update_rpars': mocker.patch(f'{_MODULE}.parameters.update'),
        'eval_timer': mocker.patch(
            f'{_MODULE}.ExpiringTimerWithDeadline',
            return_value=mock_eval_timer,
        ),
        'search_job': mocker.patch(
            f'{_MODULE}.SearchJob',
            return_value=mock_search_job,
        ),
    }
    mocks = {**mock_compile, **mock_run}
    # When running, writeSearchOutput is called as part of
    # processSearchResults, which we mock away above. This
    # means it will never be called in the tests.
    del mocks['write_output']
    return mocks


@fixture(name='mock_src_files')
def fixture_mock_glob_src_files(tmp_path, mocker):
    """Fake globbing of Fortran source files."""
    fortran_src_files = {  # glob_pattern: files_to_return
        'search*': ('search_src_file', 'search.mpi_src_file'),
        'search.mpi*': ('search.mpi_src_file', ),
        'lib.search*': ('lib.search_src_file', 'lib.search.mpi_src_file'),
        'lib.search.mpi*': ('lib.search.mpi_src_file', ),
        'intarr_hashing*.f90': ('intarr_hashing_src_file.f90', ),
        '*': (),  # No files to clean up
        }
    fortran_src_files = {
        pattern: tuple(tmp_path/f for f in files)
        for pattern, files in fortran_src_files.items()
        }
    def _mock_glob(_, pattern):
        yield from fortran_src_files[pattern]
    mocker.patch('pathlib.Path.glob', _mock_glob)
    return fortran_src_files


@fixture(name='mpirun_available')
@parametrize(available=(True, False))
def fixture_mock_mpirun_available(available, mocker):
    """Fake the presence of mpirun on the system."""
    mocker.patch('shutil.which', return_value=available)
    return available


@fixture(name='n_cores')
@parametrize(n_cores=(None, 1, 999999))
def fixture_mock_n_cores(rpars, n_cores):
    """Assign rpars.N_CORES."""
    rpars.N_CORES = n_cores


@fixture(name='rpars')
def fixture_mock_rpars():
    """Return an Rparams with attributes set to dummy values."""
    rpars = Rparams()
    rpars.disp_blocks.append('avoid IndexError at readDISPLACEMENTS_block')
    return rpars


@fixture(name='tl_version')
@parametrize(version=('1.6.1', '1.7.3', '1.7.4', '2.0.0', '9999.99.99'))
def fixture_mock_tl_version(rpars, version):
    """Assign rpars.TL_VERSION."""
    rpars.TL_VERSION = Version(version)


@fixture(name='vary')
def fixture_mock_stuff_to_vary(rpars):
    """Fake the presence of some parameters to be varied."""
    rpars.indyPars = 5


class TestSearch:
    """Collection of tests for the search function."""

    @use('n_cores')
    def test_no_displacements(self, rpars, mock_no_run, tmp_path, mocker):
        """Check behavior when no displacements are defined."""
        with execute_in_dir(tmp_path):
            search(mocker.MagicMock(name='slab'), rpars)
        assert rpars.N_CORES is not None
        for mock in mock_no_run.values():
            mock.assert_called_once()

    @use('n_cores', 'vary', 'tl_version')
    def test_supress_exec(self, rpars, mock_no_run, tmp_path, mocker):
        """Check that nothing is executed when execution is suppressed."""
        rpars.SUPPRESS_EXECUTION = True
        not_called = mock_no_run.pop('write_output')
        self.test_no_displacements(rpars, mock_no_run, tmp_path, mocker)
        not_called.assert_not_called()
        assert rpars.halt >= 3

    @pytest.mark.timeout(2)
    @use('n_cores',
         'vary',
         'mpirun_available',
         'tl_version',
         'mock_src_files')
    def test_compile_and_run(self, rpars, mock_run, tmp_path, mocker):
        """Check that nothing is executed when execution is suppressed."""
        rpars.SUPPRESS_EXECUTION = False
        called_more_than_once = (
            # The following mocks are called multiple times,
            # while test_no_displacements uses assert_called_once.
            'copy',
        )
        mocks_called_multiple_times = [mock_run.pop(n)
                                       for n in called_more_than_once]
        self.test_no_displacements(rpars, mock_run, tmp_path, mocker)
        for mock in mocks_called_multiple_times:
            mock.assert_called()
        # TODO: check relevant calls
