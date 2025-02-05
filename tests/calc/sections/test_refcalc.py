"""Tests for viperleed.calc section refcalc."""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-07-28'
__license__ = 'GPLv3+'

from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.calc.lib.context import execute_in_dir
from viperleed.calc.sections.refcalc import compile_refcalc


class TestRefCalc:
    """Test the successful outcome of a full-dynamical calculation."""

    def test_successful_run(self, refcalc_files):
        """Check that refcalc exits without errors."""
        assert not refcalc_files.failed
        assert refcalc_files.records is not None
        assert refcalc_files.records.get_last_state_for_section('refcalc')

    @parametrize(expected_file=('THEOBEAMS.csv',))
    def test_refcalc_files_present(self, refcalc_files, expected_file):
        """Ensure the expected output files were generated."""
        assert refcalc_files.expected_file_exists(expected_file)


class TestCompileRefcalc:
    """Tests for the compile_refcalc function."""

    @fixture(name='make_comptask')
    def factory_comptask(self, tmp_path, mocker):
        """Return a fake RefcalcCompileTask."""
        def _make(sources=None):
            if sources is None:
                sources = 'lib.f90', 'src.f90', None, 'muftin.f90'
            task = mocker.MagicMock(
                basedir=str(tmp_path),
                exename='test_exe',
                foldername='test_folder',
                fortran_comp=['gfortran', '-O2'],
                param='test_param',
                copy_source_files_to_local=mocker.MagicMock(),
                )
            task.__str__.return_value = f'RefcalcCompileTask {task.foldername}'
            task.get_source_files.return_value = [
                mocker.MagicMock(name=s) if isinstance(s, str) else s
                for s in sources
                ]
            return task
        return _make

    @fixture(name='run_compile')
    def fixture_run_compile(self, mocker, tmp_path):
        """Run compile_refcalc after replacing inner workings with mocks."""
        def _mock(open_raises=None,
                  mkdir_raises=None,
                  chdir_raises=None,
                  compile_raises=None):
            mocker.patch('builtins.open', create=True, side_effect=open_raises)
            mocker.patch('os.mkdir', side_effect=mkdir_raises)
            mocker.patch('os.chdir', side_effect=chdir_raises)
            mocker.patch('viperleed.calc.lib.leedbase.fortran_compile_batch',
                         side_effect=compile_raises)

        def _run(task, fails=False, **mocks):
            with execute_in_dir(tmp_path):
                _mock(**mocks)
                error = compile_refcalc(task)
            assert error if fails else not error
            return error
        return _run

    def test_success(self, make_comptask, run_compile, caplog):
        """Check the results of a successful execution of compile_refcalc."""
        comptask = make_comptask()
        run_compile(comptask)
        comptask.copy_source_files_to_local.assert_called()
        assert not caplog.text

    def test_work_exists_warning(self, make_comptask, run_compile,
                                 tmp_path, caplog):
        """Check warnings are emitted when the work directory exists."""
        (tmp_path / 'test_folder').mkdir()
        comptask = make_comptask(sources=(None, None, None, None))
        run_compile(comptask)
        assert 'Contents may get overwritten.' in caplog.text

    def test_fails_to_write_param(self, make_comptask, run_compile,
                                  caplog, mocker):
        """Check failure when PARAM can't be written to."""
        comptask = make_comptask()
        result = run_compile(comptask, fails=True, open_raises=OSError)

        expect_log = 'Error writing PARAM file'
        expect_results = (
            'Error encountered by RefcalcCompileTask test_folder',
            'while trying to write PARAM file.',
            )
        print(result)
        assert expect_log in caplog.text
        assert all(e in result for e in expect_results)

    def test_fails_to_copy_sources(self, make_comptask, run_compile, caplog):
        """Check failure when Fortran source files can't be fetched."""
        comptask = make_comptask()
        comptask.copy_source_files_to_local.side_effect = OSError
        result = run_compile(comptask, fails=True)

        expect_log = 'Error getting TensErLEED files for refcalc'
        expect_results = (
            'Error encountered by RefcalcCompileTask test_folder',
            'while trying to fetch fortran source files',
            )
        assert expect_log in caplog.text
        assert all(e in result for e in expect_results)

    def test_fails_to_compile(self, make_comptask, run_compile, caplog):
        """Check failure when compilation fails."""
        comptask = make_comptask(sources=('lib.f90', 'src.f90', None, None))
        result = run_compile(comptask, fails=True, compile_raises=Exception)

        expect_log =  'Error compiling fortran files'
        expect_results = (
            'Fortran compile error in RefcalcCompileTask test_folder',
            )
        assert expect_log in caplog.text
        assert all(e in result for e in expect_results)
