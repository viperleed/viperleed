"""Tests for viperleed.calc section deltas."""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-07-28'
__license__ = 'GPLv3+'

from pathlib import Path

from pytest_cases import fixture

from viperleed.calc.constants import DEFAULT_TENSORS
from viperleed.calc.lib.context import execute_in_dir
from viperleed.calc.sections.deltas import compile_delta
from viperleed.calc.sections.deltas import run_delta

from ...helpers import filesystem_from_dict
from ...helpers import filesystem_to_dict
from .test_refcalc import TestCompileRefcalc
from .test_refcalc import TestRunRefcalc


class TestDeltasAg100:
    """Test the successful outcome of a delta-amplitude calculation."""

    def test_successful_run(self, delta_files_ag100):
        """Check that delta-amplitude calculation exits without errors."""
        assert not delta_files_ag100.failed
        assert delta_files_ag100.records is not None
        assert delta_files_ag100.records.get_last_state_for_section('delta')

    def test_delta_input_written(self, delta_files_ag100):
        """Check that an input file was correctly written."""
        assert delta_files_ag100.expected_file_exists('delta-input')

    def test_deltas_zip_created(self, delta_files_ag100):
        """Check that an archive with delta-amplitude files was created."""
        assert delta_files_ag100.expected_file_exists('Deltas/Deltas_001.zip')


# Notice that the compile_delta and compile_refcalc functions are
# virtually identical, except for a few error messages and a few
# variable names. There is no need to add more tests. When doing
# #43, the tests can be given to the base-class method!
class TestCompileDelta(TestCompileRefcalc):
    """Tests for the compile_delta function."""

    compile_func = compile_delta
    compiler_cls_name = 'DeltaCompileTask'
    default_sources = 'delta.f', 'lib.tleed.f', 'lib.delta.f', 'GLOBAL'
    section_name = 'delta-amplitudes'


class TestRunDelta:
    """Tests for the run_delta function."""

    @fixture(name='runtask')
    def factory_runtask(self, mocker):
        """Return a fake DeltaRunTask."""
        comptask = mocker.MagicMock(
            exename='test_exe',
            foldername='test_folder',
            )
        runtask = mocker.MagicMock(
            comptask=comptask,
            deltaname='test_delta',
            tensorname='test_tensor',
            deltalogname='test_delta.log',
            din='test input',
            name='test_name',
            foldername='calculating_test_delta',
            )
        return runtask

    @fixture(name='mock_implementation')
    def fixture_mock_implementation(self, mocker):
        """Replace implementation details of run_delta with mocks."""
        def _make(**kwargs):
            to_mock = {
                # NB: We purposely DO NOT mock os.chdir and os.mkdir
                # as we want to make sure to run things in a dedicated
                # folder! This would also mess around with the
                # execute_in_dir context that we use in the tests
                'pathlib.Path.open': ('open_raises', None),
                'pathlib.Path.read_text': ('path_read', None),
                'shutil.copy2': ('copy_raises', None),
                'shutil.rmtree': ('rmtree_raises', None),
                'subprocess.run': ('subprocess_raises', None),
                }
            return {
                target: mocker.patch(target, side_effect=kwargs.get(*fails))
                for target, fails in to_mock.items()
                }
        return _make

    @fixture(name='run')
    def fixture_run(self, runtask, mock_implementation, tmp_path):
        """Execute run_delta, potentially mocking before."""
        def _run(fails=False, mock=True, **mocks):
            with execute_in_dir(tmp_path):
                if mock or mocks:
                    mock_implementation(**mocks)
                error = run_delta(runtask)
            assert error if fails else not error
            return error
        return _run

    # Here are tests that are common to run_delta and run_refcalc.
    # They may help highlighting how to generalize the functions
    # in #43. Notice that many of the tests below are implemented
    # only because of slight differences in the error messages
    # returned/logged (typically 'delta' and variations instead
    # of 'refcalc' and such).
    test_fails_to_mkdir = TestRunRefcalc.test_fails_to_mkdir
    test_fails_to_cleanup = TestRunRefcalc.test_fails_to_cleanup
    test_success = TestRunRefcalc.test_success
    test_work_exists_warning = TestRunRefcalc.test_work_exists_warning

    def test_collects_files(self, runtask, run, tmp_path, mocker):
        """Check that files are collected in the right place."""
        root_files = {
            DEFAULT_TENSORS : {runtask.tensorname: 'tensor contents'},
            runtask.comptask.foldername: {  # A fake delta executable
                runtask.comptask.exename: '',
                },
            runtask.foldername: {'DELWV': 'delta contents'},  # OUTPUT
            }
        filesystem_from_dict(root_files, tmp_path)

        mocker.patch('subprocess.run')
        with execute_in_dir(tmp_path):
            run(fails=False, mock=False)

        expect_files = root_files.copy()
        # All the work files should have been
        # moved to tmp_path and renamed there
        outputs = expect_files.pop(runtask.foldername)
        expect_files[runtask.deltaname] = outputs['DELWV']
        assert filesystem_to_dict(tmp_path) == expect_files

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
        expect_log = 'Error writing delta log part'
        assert expect_log in caplog.text

    def test_fails_on_tensor_not_found(self, run, caplog):
        """Check failure when no tensor file is available."""
        result = run(fails=True, copy_raises=FileNotFoundError)
        expect_error = 'Tensors not found'
        expect_log = 'Tensors file not found'
        assert expect_error in result
        assert expect_log in caplog.text

    def test_fails_to_collect_tensor(self, run, caplog):
        """Check failure when failing to collect Tensors."""
        result = run(fails=True, copy_raises=OSError)
        expect_error = 'Error copying Tensors file'
        expect_log = expect_error
        assert expect_log in caplog.text
        assert expect_error in result

    def test_fails_to_copy_executable(self, run, runtask, caplog):
        """Check failure when copying the delta executable fails."""
        def _copy_exe_fails(src, _):
            if Path(src).name == runtask.comptask.exename:
                raise OSError
        result = run(fails=True, copy_raises=_copy_exe_fails)
        expect_error = 'Failed to get delta executable'
        expect_log = 'Error getting delta executable'
        assert expect_error in result
        assert expect_log in caplog.text

    def test_fails_to_copy_output_file(self, run, caplog):
        """Check complaints when failing to copy a delta output file."""
        file = 'DELWV'
        def _copy_file_fails(src, _):
            if Path(src).name == file:
                raise OSError
        result = run(fails=True, copy_raises=_copy_file_fails)
        expect_log = f'Failed to copy delta output file {file}'
        expect_error = 'Failed to copy result file out'
        assert expect_error in result
        assert expect_log in caplog.text

    def test_fails_to_exec_delta(self, run, caplog):
        """Check failure when execution of the compiled delta fails."""
        result = run(fails=True, subprocess_raises=Exception)
        expect_error = 'Error during delta execution'
        expect_log = 'Error while executing delta-amplitudes calculation'
        assert expect_error in result
        assert expect_log in caplog.text

    def test_log_read_fails(self, run, caplog):
        """Check warnings when failing to read the local log file."""
        run(fails=False, path_read=OSError)
        expect_log = 'Could not read local delta log'
        assert expect_log in caplog.text
