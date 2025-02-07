"""Tests for viperleed.calc section refcalc."""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-07-28'
__license__ = 'GPLv3+'

from pathlib import Path

import pytest
from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.calc.lib.context import execute_in_dir
from viperleed.calc.sections.refcalc import compile_refcalc
from viperleed.calc.sections.refcalc import run_refcalc

from ...helpers import filesystem_from_dict
from ...helpers import filesystem_to_dict


_MODULE = 'viperleed.calc.sections.refcalc'


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

    compile_func = compile_refcalc
    compiler_cls_name = 'RefcalcCompileTask'
    default_sources = 'lib.f90', 'src.f90', None, 'muftin.f90'
    section_name = 'refcalc'

    @fixture(name='make_comptask')
    def factory_comptask(self, mocker):
        """Return a fake RefcalcCompileTask."""
        def _make(sources=None):
            if sources is None:
                sources = self.default_sources
            task = mocker.MagicMock(
                exename='test_exe',
                foldername='test_folder',
                fortran_comp=['gfortran', '-O2'],
                param='test_param',
                copy_source_files_to_local=mocker.MagicMock(),
                )
            task.__str__.return_value = (
                f'{self.compiler_cls_name} {task.foldername}'
                )
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
            mocker.patch('pathlib.Path.open',
                         create=True,
                         side_effect=open_raises)
            mocker.patch('os.mkdir', side_effect=mkdir_raises)
            mocker.patch('os.chdir', side_effect=chdir_raises)
            mocker.patch('viperleed.calc.lib.leedbase.fortran_compile_batch',
                         side_effect=compile_raises)

        def _run(task, fails=False, **mocks):
            with execute_in_dir(tmp_path):
                _mock(**mocks)
                cls = type(self)
                error = cls.compile_func(task)  # Avoid to pass self
            assert error if fails else not error
            return error
        return _run

    def test_fails_to_compile(self, make_comptask, run_compile, caplog):
        """Check failure when compilation fails."""
        comptask = make_comptask(sources=('lib.f90', 'src.f90', None, None))
        result = run_compile(comptask, fails=True, compile_raises=Exception)

        expect_log =  'Error compiling fortran files'
        expect_results = (
            f'Fortran compile error in {self.compiler_cls_name} test_folder',
            )
        assert expect_log in caplog.text
        assert all(e in result for e in expect_results)

    def test_fails_to_copy_sources(self, make_comptask, run_compile, caplog):
        """Check failure when Fortran source files can't be fetched."""
        comptask = make_comptask()
        comptask.copy_source_files_to_local.side_effect = OSError
        result = run_compile(comptask, fails=True)

        expect_log = f'Error getting TensErLEED files for {self.section_name}'
        expect_results = (
            f'Error encountered by {self.compiler_cls_name} test_folder',
            'while trying to fetch fortran source files',
            )
        assert expect_log in caplog.text
        assert all(e in result for e in expect_results)

    def test_fails_to_write_param(self, make_comptask, run_compile, caplog):
        """Check failure when PARAM can't be written to."""
        comptask = make_comptask()
        result = run_compile(comptask, fails=True, open_raises=OSError)

        expect_log = 'Error writing PARAM file'
        expect_results = (
            f'Error encountered by {self.compiler_cls_name} test_folder',
            'while trying to write PARAM file.',
            )
        assert expect_log in caplog.text
        assert all(e in result for e in expect_results)

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
        expect_log = 'Contents may get overwritten.'
        assert expect_log in caplog.text



class TestRunRefcalc:
    """Tests for the run_refcalc function."""

    @fixture(name='runtask')
    def factory_runtask(self, mocker):
        """Return a fake RefcalcRunTask."""
        comptask = mocker.MagicMock(
            exename='test_exe',
            foldername='test_folder',
            )
        runtask = mocker.MagicMock(
            collect_at=None,
            comptask=comptask,
            energy=5.0,
            fin='test input',
            foldername='test_run_folder',
            logname='test.log',
            name='test_name',
            single_threaded=False,
            tl_version=mocker.MagicMock(return_value='1.7.3'),
            )
        return runtask

    @fixture(name='mock_implementation')
    def fixture_mock_implementation(self, mocker):
        """Replace implementation details of run_refcalc with mocks."""
        def _fake_edit_fin(task):
            return task.fin
        def _make(**kwargs):
            to_mock = {
                # NB: We purposely DO NOT mock os.chdir and os.mkdir
                # as we want to make sure to run things in a dedicated
                # folder! This would also mess around with the
                # execute_in_dir context that we use in the tests
                'pathlib.Path.open': ('open_raises', None),
                'pathlib.Path.read_text': ('path_read', None),
                'pathlib.Path.write_text': ('path_write', None),
                'shutil.copy2': ('copy_raises', None),
                'shutil.rmtree': ('rmtree_raises', None),
                'subprocess.run': ('subprocess_raises', None),
                f'{_MODULE}.edit_fin_energy_lmax': ('fin_raises',
                                                    _fake_edit_fin),
                }
            return {
                target: mocker.patch(target, side_effect=kwargs.get(*fails))
                for target, fails in to_mock.items()
                }
        return _make

    @fixture(name='run')
    def fixture_run(self, runtask, mock_implementation, tmp_path):
        """Execute run_refcalc, potentially mocking before."""
        def _run(fails=False, mock=True, **mocks):
            if mock or mocks:
                mock_implementation(**mocks)
            with execute_in_dir(tmp_path):
                error = run_refcalc(runtask)
            assert error if fails else not error
            return error
        return _run

    _amp_missing = {
        '2.0': 'Refcalc output file amp.out not found.',
        '1.7.3': 'Refcalc output file amp.out not found.',
        '1.7': None,
        }

    @parametrize('version,expect_log', _amp_missing.items(), ids=_amp_missing)
    # pylint: disable-next=too-many-arguments  # 3/6 fixtures
    def test_amp_out_missing(self, version, expect_log, runtask, run, caplog):
        """Check logging messages when failing to copy amp.out."""
        def _amp_not_found(src, _):
            # pylint: disable-next=magic-value-comparison
            if Path(src).name == 'amp.out':
                raise FileNotFoundError
        runtask.tl_version = version
        run(fails=False, copy_raises=_amp_not_found)
        assert expect_log in caplog.text if expect_log else not caplog.text

    @parametrize(dest=(None, 'another_place'))
    # pylint: disable-next=too-many-arguments  # 4/6 are fixtures
    def test_collects_files(self, dest, runtask, run, tmp_path, mocker):
        """Check that files are collected in the right place."""
        work_files = {  # Files that will be moved
            'T_1': 'tensor one',
            'T_2': 'tensor two',
            'fd.out': 'fd.out contents',
            'amp.out': 'amp.out contents',
            }
        root_files = {
            runtask.foldername: work_files,
            runtask.comptask.foldername: {  # A fake refcalc executable
                runtask.comptask.exename: '',
                },
            }
        # And the place where stuff is to be collected
        dest_path = tmp_path/dest if dest else tmp_path
        if dest:
            runtask.collect_at = dest_path
            root_files[dest] = {}
        filesystem_from_dict(root_files, tmp_path)

        mocker.patch(f'{_MODULE}.edit_fin_energy_lmax', lambda task: task.fin)
        mocker.patch('subprocess.run')
        with execute_in_dir(tmp_path):
            run(fails=False, mock=False)

        expect_files = {} if dest else root_files.copy()
        # All the work files should have been
        # moved to tmp_path and renamed there
        try:
            del expect_files[runtask.foldername]
        except KeyError:
            pass
        for fname, contents in work_files.items():
            file = Path(fname)
            expect_files[f'{file.stem}_5.00eV{file.suffix}'] = contents
        assert filesystem_to_dict(dest_path) == expect_files

    def test_fails_log_append(self, run, mock_implementation, caplog, mocker):
        """Check warnings when failing to extend the main log file."""
        def _open_log_fails(path, mode, *args, **kwargs):
            # pylint: disable-next=magic-value-comparison
            if path.suffix == '.log' and mode == 'a':
                raise OSError
            # pylint: disable-next=unspecified-encoding  # In **kwargs
            return open(path, mode, *args, **kwargs)
        mock_implementation()  # BEFORE replacing Path.open again
        mocker.patch('pathlib.Path.open', _open_log_fails)
        run(fails=False, mock=False)
        expect_log = 'Error writing refcalc log part'
        assert expect_log in caplog.text

    # pylint: disable-next=too-many-arguments  # All fixtures
    def test_fails_move_tensors(self, runtask, run, tmp_path, mocker, caplog):
        """Check failure when moving Tensor files fails."""
        # Make sure there are some Tensor files
        filesystem_from_dict({runtask.foldername: {'T_1': None}}, tmp_path)
        mocker.patch('shutil.move', side_effect=OSError)

        with execute_in_dir(tmp_path):
            result = run(fails=True, mock=True)
        expect_error = 'Failed to copy Tensor file out'
        expect_log = 'Failed to copy refcalc output file T_1'
        assert expect_error in result
        assert expect_log in caplog.text

    def test_fails_to_cleanup(self, run, caplog):
        """Check that failing to remove the work folder is not a failure."""
        run(fails=False, rmtree_raises=OSError)
        expect_log = 'Error deleting folder'
        assert expect_log in caplog.text

    def test_fails_to_copy_executable(self, run, caplog):
        """Check failure when copying the refcalc executable fails."""
        result = run(fails=True, copy_raises=OSError)
        expect_error = 'Failed to get refcalc executable'
        expect_log = 'Error getting refcalc executable'
        assert expect_error in result
        assert expect_log in caplog.text

    _copy_fails = {  # file: expect_error
        'fd.out': 'Failed to copy fd.out file out.',
        'amp.out': None,  # No error, only logging
        }

    @parametrize('file,expect_error', _copy_fails.items(), ids=_copy_fails)
    def test_fails_to_copy_output_file(self, file, expect_error, run, caplog):
        """Check complaints when failing to copy a refcalc output file."""
        def _copy_file_fails(src, _):
            if Path(src).name == file:
                raise OSError
        error = run(fails=expect_error, copy_raises=_copy_file_fails)
        expect_log = f'Failed to copy refcalc output file {file}'
        assert expect_error in error if expect_error else not error
        assert expect_log in caplog.text

    def test_fails_to_exec_refcalc(self, run, caplog):
        """Check failure when execution of the compiled refcalc fails."""
        result = run(fails=True, subprocess_raises=Exception)
        expect_error = 'Error during refcalc execution'
        expect_log = 'Error while executing reference calculation'
        assert expect_error in result
        assert expect_log in caplog.text

    def test_fails_to_mkdir(self, run, mocker):
        """Check complaints failing to create work."""
        mocker.patch('pathlib.Path.mkdir', side_effect=OSError)
        with pytest.raises(OSError):
            run()

    def test_fails_to_write_local_fin(self, run, mock_implementation,
                                      caplog, mocker):
        """Check that there are no complaints when writing FIN fails."""
        def _write_fin_fails(file, *_, **__):
            if file.name:
                raise OSError
        mock_implementation()
        mocker.patch('pathlib.Path.write_text', _write_fin_fails)
        run(fails=False, mock=False)
        assert not caplog.text

    def test_log_read_fails(self, run, caplog):
        """Check warnings when failing to read the local log file."""
        run(fails=False, path_read=OSError)
        expect_log = 'Could not read local refcalc log'
        assert expect_log in caplog.text

    def test_single_threaded_shortcut(self, runtask, run,
                                      mock_implementation, mocker):
        """Check that a single-threaded run executes fewer calls."""
        def _fail_copy_except_refcalc_exe(src, dst):
            if Path(src).name != runtask.comptask.exename:
                pytest.fail('copy2 was unexpectedly called. Arguments '
                            f'were src={src} and dst={dst}')
        mocks = mock_implementation(
            copy_raises=_fail_copy_except_refcalc_exe,
            path_read=OSError,
            path_write=OSError,
            )
        mocks['os.mkdir'] = mocker.patch('os.mkdir')
        mocks['os.chdir'] = mocker.patch('os.chdir')
        runtask.single_threaded = True
        run(fails=False, mock=False)
        not_called = (
            'os.mkdir',
            'shutil.rmtree',
            f'{_MODULE}.edit_fin_energy_lmax',
            )
        for mocked in not_called:
            mock = mocks[mocked]
            mock.assert_not_called()

    def test_success(self, run, mock_implementation, caplog):
        """Check expected result of a successful execution."""
        mocks = mock_implementation()
        run(mock=False)
        assert not caplog.text
        for mock in mocks.values():
            mock.assert_called()

    def test_work_exists_warning(self, run, runtask, caplog, tmp_path):
        """Check that warnings are emitted if work exists already."""
        work = tmp_path / runtask.foldername
        work.mkdir()
        run(fails=False, mock=True)
        expect_log = 'Contents may get overwritten.'
        assert expect_log in caplog.text
