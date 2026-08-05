"""Tests for viperleed.calc section refcalc.

This module contains tests for the implementation of the
compile_refcalc function. No reference calculation is
truly executed.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2023-07-28'
__license__ = 'GPLv3+'

from pytest_cases import fixture

from viperleed.calc.lib.context import execute_in_dir
from viperleed.calc.sections.refcalc import compile_refcalc


class TestCompileRefcalc:
    """Tests for the compile_refcalc function."""

    compile_func = compile_refcalc
    compiler_cls_name = 'RefcalcCompileTask'
    default_sources = 'lib.f90', 'src.f90', None, 'muftin.f90'
    section_name = 'refcalc'

    @fixture(name='make_comptask')
    def factory_comptask(self, mocker, tmp_path):
        """Return a fake RefcalcCompileTask."""
        def _mock_sources(sources):
            if sources is None:
                sources = self.default_sources
            mock_sources = []
            for src in sources:
                if isinstance(src, str):
                    mock_src = mocker.MagicMock()
                    mock_src.name = src
                else:
                    mock_src = src
                mock_sources.append(mock_src)
            return mock_sources

        def _make(sources=None):
            task = mocker.MagicMock(
                exename='test_exe',
                foldername='test_folder',
                fortran_comp=['gfortran', '-O2'],
                param='test_param',
                copy_source_files_to_local=mocker.MagicMock(),
                )
            task.exec_dir = tmp_path / task.foldername
            task.__str__.return_value = (
                f'{self.compiler_cls_name} {task.foldername}'
                )
            task.get_source_files.return_value = _mock_sources(sources)
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

        expect_log = 'Error compiling fortran files'
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
        work = tmp_path / 'test_folder'
        work.mkdir()
        comptask = make_comptask(sources=(None, None, None, None))
        run_compile(comptask)
        expect_log = 'Contents may get overwritten.'
        assert expect_log in caplog.text
