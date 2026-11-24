"""Tests for function search of viperleed.calc.sections.

This module contains tests for the implementation of the
search function of viperleed.calc.sections. They mock away
parts of the function that actually do execute a search.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2025-06-06'
__license__ = 'GPLv3+'

from pathlib import Path
from shutil import which as sh_which  # Import it because we mock it

import pytest
from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.calc.classes.rparams.rparams import Rparams
from viperleed.calc.classes.rparams.special.energy_range import TheoEnergies
from viperleed.calc.classes.search_job import IS_WINDOWS
from viperleed.calc.files.beams import readOUTBEAMS
from viperleed.calc.lib.context import execute_in_dir
from viperleed.calc.lib.version import Version
from viperleed.calc.sections.search import search

from ....helpers import TEST_DATA
from ....helpers import write_delta_stub
from .conftest import _MODULE

_AG100_INPUTS = TEST_DATA / 'Ag(100)' / 'initialization'
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
    aliases = {'mpirun'}
    if IS_WINDOWS:
        aliases.add('mpiexec')

    def _mock_which(cmd, *others, **kwargs):
        if cmd in aliases:
            return 'path/to/mpirun' if available else None
        return sh_which(cmd, *others, **kwargs)

    mocker.patch('shutil.which', side_effect=_mock_which)
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


def _skip_on_missing_compiler(rpars, mocker, missing_dependency):
    """Skip (or fail) the test if search finds no Fortran compiler."""
    compilers = {  # Method of rpars: name of the dependency
        'getFortranComp': 'gfortran',
        'getFortranMpiComp': 'mpifort',
        }
    for method_name, dependency in compilers.items():
        find_or_skip = _find_compiler_or_skip(getattr(rpars, method_name),
                                              dependency,
                                              missing_dependency)
        mocker.patch.object(rpars, method_name, side_effect=find_or_skip)


def _find_compiler_or_skip(find_compiler, dependency, missing_dependency):
    """Return a wrapper of find_compiler that skips if none is found."""
    def _find_or_skip(*args, **kwargs):
        try:
            return find_compiler(*args, **kwargs)
        except FileNotFoundError:
            missing_dependency(dependency, 'compiler not found')
        return None  # Unreachable: missing_dependency raises
    return _find_or_skip


class TestSearchDeltaEnergies:
    """Tests for the check of delta energies at the start of search."""

    @fixture(name='run_search')
    def fixture_run_search(self, rpars, mock_no_run, tmp_path, mocker):
        """Return a function to run search with given delta files."""
        # As in Ag(100): EXPBEAMS from 48 eV, THEO_ENERGIES 50-350 eV
        rpars.expbeams = readOUTBEAMS(str(_AG100_INPUTS / 'EXPBEAMS.csv'))
        rpars.THEO_ENERGIES = TheoEnergies(50, 350, 3)

        def _run(*delta_ranges, old_format=False):
            # One delta file per atom, one atom per (e_min, e_max)
            rpars.search_atlist = []
            for i, (e_min, e_max) in enumerate(delta_ranges, start=1):
                delta_name = f'DEL_{i}_Ag_1'
                write_delta_stub(tmp_path / delta_name,
                                 range(e_min, e_max + 1, 3),
                                 old_format=old_format)
                rpars.search_atlist.append(
                    mocker.MagicMock(current_deltas=[delta_name])
                    )
            with execute_in_dir(tmp_path):
                search(mocker.MagicMock(name='slab'), rpars)
            return mock_no_run
        return _run

    @parametrize(e_min=(50, 47))
    @parametrize(e_max=(350, 353))
    def test_delta_covers_range(self, e_min, e_max, run_search):
        """Check that no complaint is raised for sufficient energies."""
        mocks = run_search((e_min, e_max))
        mocks['write_rfinfo'].assert_called_once()

    @parametrize('e_min, e_max', ((50, 347), (53, 350), (50, 101)))
    def test_delta_too_short(self, e_min, e_max, run_search):
        """Check complaints if delta energies are too short."""
        with pytest.raises(RuntimeError, match='DEL_1_Ag_1') as exc:
            run_search((e_min, e_max))
        assert f'{e_max:.2f}' in str(exc.value)
        assert 'from 50.00 to 350.00 eV' in str(exc.value)

    def test_second_file_too_short(self, run_search):
        """Check that all delta files are checked."""
        with pytest.raises(RuntimeError, match='DEL_2_Ag_1') as exc:
            run_search((50, 350), (50, 347))
        assert 'DEL_1_Ag_1' not in str(exc.value)

    def test_experiment_narrower(self, rpars, run_search):
        """Check that deltas need to cover the experiment range only."""
        for beam in rpars.expbeams:
            beam.intens = {e: intensity
                           for e, intensity in beam.intens.items()
                           if 100 <= e <= 300}
        mocks = run_search((89, 311))
        mocks['write_rfinfo'].assert_called_once()

    @parametrize(e_min=(50, 200, 302))
    def test_v16_covers_range(self, e_min, rpars, run_search):
        """Check delta files of TensErLEED v1.6 at high energies."""
        rpars.THEO_ENERGIES = TheoEnergies(e_min, 350, 3)
        mocks = run_search((e_min, 350), old_format=True)
        mocks['write_rfinfo'].assert_called_once()

    @parametrize(e_min=(50, 200, 302))
    def test_v16_too_short(self, e_min, rpars, run_search):
        """Check complaints for too-short v1.6 delta files."""
        rpars.THEO_ENERGIES = TheoEnergies(e_min, 350, 3)
        with pytest.raises(RuntimeError, match='DEL_1_Ag_1') as exc:
            run_search((e_min, 347), old_format=True)
        assert f'{e_min:.2f} to 347.00' in str(exc.value)
        assert f'from {e_min:.2f} to 350.00 eV' in str(exc.value)


class TestSearch:
    """Collection of tests for the search function."""

    @use('n_cores', 'tl_version')
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
         'tl_version',
         'mock_src_files',
         'mpirun_available')
    def test_compile_and_run(self, rpars, mock_run, tmp_path, mocker,
                             missing_dependency):
        """Check that nothing is executed when execution is suppressed."""
        rpars.SUPPRESS_EXECUTION = False
        _skip_on_missing_compiler(rpars, mocker, missing_dependency)
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


class TestSearchPdfOneGeneration:
    """Tests for plots of a search that ends after one generation."""

    @use('vary', 'mock_src_files', 'tl_version')
    def test_pdfs_written(self, rpars, mock_run, tmp_path, mocker,
                          missing_dependency):
        """Check PDFs are written if SD.TL has a single generation."""
        rpars.SUPPRESS_EXECUTION = False
        rpars.N_CORES = 1
        _skip_on_missing_compiler(rpars, mocker, missing_dependency)
        # The search finishes before SD.TL is ever evaluated
        mock_run['eval_timer'].return_value.has_expired.return_value = False
        sdtl_content = 'one generation'
        mocker.patch(f'{_MODULE}.iosearch.readSDTL_next',
                     return_value=(len(sdtl_content), sdtl_content))
        mocker.patch(f'{_MODULE}.iosearch.readSDTL_blocks',
                     return_value=[(100, [0.3, 0.4], [(100, (1, 2))])])
        write_progress = mocker.patch(
            f'{_MODULE}.searchpdf.writeSearchProgressPdf'
            )
        write_report = mocker.patch(
            f'{_MODULE}.searchpdf.writeSearchReportPdf'
            )
        with execute_in_dir(tmp_path):
            search(mocker.MagicMock(name='slab'), rpars)
        write_progress.assert_called_once()
        write_report.assert_called_once()
        gens = write_progress.call_args[0][1]
        assert gens == [100]

    @use('vary', 'mock_src_files', 'tl_version')
    def test_restart_after_one_generation(self, rpars, mock_run, tmp_path,
                                          mocker, missing_dependency):
        """Check the final read of SD.TL after a restarted search."""
        rpars.SUPPRESS_EXECUTION = False
        rpars.N_CORES = 1
        _skip_on_missing_compiler(rpars, mocker, missing_dependency)

        # One generation is read in the first run, which restarts as
        # GAUSSIAN_WIDTH changes. The second run ends before SD.TL is
        # evaluated, and has one more generation in the final read.
        def _change_gaussian_width(rpars_):
            rpars_.GAUSSIAN_WIDTH *= 0.5
        mock_run['update_rpars'].side_effect = _change_gaussian_width
        mock_run['search_job'].return_value.is_running.side_effect = (
            True, False,   # First run, SD.TL read in the loop
            False,         # Second run
            )
        sdtl_content = 'one generation'
        mocker.patch(f'{_MODULE}.iosearch.readSDTL_next',
                     return_value=(len(sdtl_content), sdtl_content))
        first_config, final_config = [(100, (1, 2))], [(100, (2, 1))]
        mocker.patch(f'{_MODULE}.iosearch.readSDTL_blocks',
                     side_effect=([(100, [0.3, 0.4], first_config)],
                                  [(7, [0.2, 0.3], final_config)]))
        mocker.patch(f'{_MODULE}.parabolaFit', return_value=(None, None))
        write_progress = mocker.patch(
            f'{_MODULE}.searchpdf.writeSearchProgressPdf'
            )
        mocker.patch(f'{_MODULE}.searchpdf.writeSearchReportPdf')
        with execute_in_dir(tmp_path):
            search(mocker.MagicMock(name='slab'), rpars)
        write_progress.assert_called_once()
        _, gens, _, lastconfig = write_progress.call_args[0]
        assert gens == [100, 107]
        assert lastconfig is final_config
