"""Tests for module cleanup of viperleed.calc.section."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-02-03'
__license__ = 'GPLv3+'

from pathlib import Path
import re
import shutil
from zipfile import ZipFile

import pytest
from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.calc.classes.rparams import Rparams
from viperleed.calc.constants import DEFAULT_DELTAS
from viperleed.calc.constants import DEFAULT_OUT
from viperleed.calc.constants import DEFAULT_SUPP
from viperleed.calc.constants import DEFAULT_TENSORS
from viperleed.calc.constants import DEFAULT_WORK_HISTORY
from viperleed.calc.constants import LOG_PREFIX
from viperleed.calc.constants import ORIGINAL_INPUTS_DIR_NAME
from viperleed.calc.lib.context import execute_in_dir
from viperleed.calc.sections.calc_section import ALL_INPUT_FILES
from viperleed.calc.sections.cleanup import _OUT_FILES
from viperleed.calc.sections.cleanup import _SUPP_DIRS
from viperleed.calc.sections.cleanup import _SUPP_FILES
from viperleed.calc.sections.cleanup import OPTIONAL_INPUT_FILES
from viperleed.calc.sections.cleanup import cleanup
from viperleed.calc.sections.cleanup import organize_workdir
from viperleed.calc.sections.cleanup import prerun_clean
from viperleed.calc.sections.cleanup import preserve_original_inputs
from viperleed.calc.sections.cleanup import _collect_delta_files
from viperleed.calc.sections.cleanup import _collect_out_contents
from viperleed.calc.sections.cleanup import _collect_supp_contents
from viperleed.calc.sections.cleanup import _copy_files_and_directories
from viperleed.calc.sections.cleanup import _organize_all_work_directories
from viperleed.calc.sections.cleanup import _write_final_log_messages
from viperleed.calc.sections.cleanup import _write_manifest_file
from viperleed.calc.sections.cleanup import _zip_folder
from viperleed.calc.sections.cleanup import _zip_subfolders

from ...helpers import CustomTestException
from ...helpers import filesystem_from_dict
from ...helpers import filesystem_to_dict
from ...helpers import raises_test_exception

_MODULE = 'viperleed.calc.sections.cleanup'


@fixture(name='rpars')
def fixture_rpars():
    """Return an empty Rparams."""
    return Rparams()


@fixture(name='workdir')
def fixture_workdir(tmp_path):
    """Create a temporary work directory."""
    workdir = tmp_path / 'work_directory'
    workdir.mkdir()
    return workdir


class TestCleanup:
    """Tests for the cleanup function."""

    mocked = (  # All the stuff that is called inside cleanup
        'logger.info',
        '_organize_all_work_directories',
        '_write_manifest_file',
        '_write_final_log_messages',
        'close_all_handlers',
        'logging.shutdown',
        )

    @fixture(name='mock_implementation')
    def fixture_mock_implementation(self, mocker):
        """Replace implementation details of cleanup with mocks."""
        return {f: mocker.patch(f'{_MODULE}.{f}') for f in self.mocked}

    @fixture(name='manifest')
    def mock_manifest(self):
        """Return the contents of a sample manifest file."""
        return {'file1.txt', 'file2.txt', 'dir1'}

    @parametrize(mock_name=mocked)
    def test_raises(self, mock_name, rpars, manifest, mock_implementation):
        """Check complaints when the implementation raises exceptions."""
        mock = mock_implementation[mock_name]
        mock.side_effect = CustomTestException
        with pytest.raises(CustomTestException):
            cleanup(manifest, rpars)

    def test_success(self, rpars, manifest, mock_implementation):
        """Check a successful execution of cleanup."""
        cleanup(manifest, rpars)
        calls = {
            'logger.info': '\nStarting cleanup...',
            '_organize_all_work_directories': rpars,
            '_write_manifest_file': manifest,
            '_write_final_log_messages': rpars,
            'close_all_handlers': None,
            'logging.shutdown': None,
            }
        for func, arg in calls.items():
            mock = mock_implementation[func]
            if arg is None:
                mock.assert_called_once()
            else:
                mock.assert_called_once_with(arg)

    def test_no_rpars(self, manifest, mock_implementation, mocker):
        """Check calls when no Rparams is passed."""
        # Create a "singleton" that we can use to check that
        # cleanup creates an empty Rparams in this case
        fake_rpars = mocker.MagicMock()
        mocker.patch(f'{_MODULE}.Rparams', return_value=fake_rpars)

        cleanup(manifest)
        calls = {
            'logger.info': '\nStarting cleanup...',
            '_organize_all_work_directories': fake_rpars,
            '_write_manifest_file': fake_rpars.manifest,
            '_write_final_log_messages': fake_rpars,
            'close_all_handlers': None,
            'logging.shutdown': None,
            }
        for func, arg in calls.items():
            mock = mock_implementation[func]
            if arg is None:
                mock.assert_called_once()
            else:
                mock.assert_called_once_with(arg)


class TestPreserveOriginalInputs:
    """Tests for the preserve_original_inputs function."""

    @fixture(name='mock_implementation')
    def fixture_mock_implementation(self, mocker):
        """Return a fake shutil.copy2 and a logger."""
        def _mock(**kwargs):
            mock_copy = mocker.patch('shutil.copy2', **kwargs)
            mock_logger = mocker.patch(f'{_MODULE}.logger')
            return mock_copy, mock_logger
        return _mock

    @fixture(name='make_inputs')
    def fixture_make_inputs(self, tmp_path):
        """Create the given input files in a temporary directory."""
        def _make(files):
            for file in files:
                (tmp_path / file).touch()
            return tmp_path
        return _make

    @fixture(name='run_preserve')
    def fixture_run_preserve(self, rpars, make_inputs, mock_implementation):
        """Execute preserve_original_inputs."""
        def _run(input_files, **kwargs):
            tmp_path = make_inputs(input_files)
            mocked = mock_implementation(**kwargs)
            with execute_in_dir(tmp_path):
                preserve_original_inputs(rpars)
            assert (tmp_path/ORIGINAL_INPUTS_DIR_NAME).is_dir()
            return (*mocked, rpars)
        return _run

    cwd_inputs = {
        'simple': ALL_INPUT_FILES,
        'missing optional': ALL_INPUT_FILES - set(OPTIONAL_INPUT_FILES),
        }

    @parametrize(input_files=cwd_inputs.values(), ids=cwd_inputs)
    def test_success(self, input_files, run_preserve):
        """Check execution of a successful storage."""
        # There's one less because we have two acceptable EXPBEAMS
        nr_files_copied = len(input_files) - 1
        copy, logger, *_ = run_preserve(input_files)

        assert copy.call_count == nr_files_copied
        logger.warning.assert_not_called()

    def test_required_missing(self, run_preserve):
        """Check that a missing mandatory input file causes warnings."""
        missing_required = ('PARAMETERS',)
        copy, logger, rpars = run_preserve(missing_required)

        logger.warning.assert_called()
        assert copy.call_count == 1
        assert rpars.halt == 1

    def test_copy_fails(self, run_preserve):
        """Check expected complaints when failing to copy files."""
        _, logger, rpars = run_preserve(ALL_INPUT_FILES, side_effect=OSError)
        assert logger.warning.call_count == len(ALL_INPUT_FILES) - 1
        assert rpars.halt == 1

    def test_mkdir_fails(self, rpars, tmp_path, mocker):
        """Check complaints when failing to create original_inputs."""
        mocker.patch('pathlib.Path.mkdir', side_effect=OSError)
        with pytest.raises(OSError):
            preserve_original_inputs(rpars)
        assert not (tmp_path / ORIGINAL_INPUTS_DIR_NAME).exists()

    _skip_expbeams = {  # skip: keep
        'EXPBEAMS': 'EXPBEAMS.csv',
        'EXPBEAMS.csv': 'EXPBEAMS',
        }

    @parametrize('skip,expect', _skip_expbeams.items())
    def test_expbeams_selection(self, skip, expect, run_preserve, tmp_path):
        """Test that the first existing EXPBEAMS file is selected."""
        input_files = ALL_INPUT_FILES - {skip}
        nr_files_copied = len(input_files)
        copy, logger, *_ = run_preserve(input_files)

        assert copy.call_count == nr_files_copied
        logger.warning.assert_not_called()
        assert (tmp_path/expect).is_file()

        # Check that we copied the expected one to original_inputs
        # Notice that we can't just use assert_any_call, as the
        # arguments may be OS-dependent with Paths (e.g., on WSL)
        calls = {
            (Path(src).name, Path(dst).name)
            for (src, dst), *_ in copy.call_args_list
            }
        assert (expect, ORIGINAL_INPUTS_DIR_NAME) in calls


class TestWriteManifest:
    """Tests for the _write_manifest_file function."""

    _success = {
        'unique': ('file1.txt', 'file2.txt', 'file3.txt'),
        'duplicates': ('file1.txt', 'file2.txt', 'file1.txt'),
        }

    @parametrize(contents=_success.values(), ids=_success)
    def test_success(self, contents, tmp_path, caplog):
        """Check successful writing to the manifest file."""
        caplog.set_level(0)  # All messages
        with execute_in_dir(tmp_path):
            _write_manifest_file(contents)
        manifest = tmp_path/'manifest'
        assert manifest.is_file()

        # Check contents
        written = sorted(manifest.read_text().splitlines())
        assert written == sorted(set(contents))

        # Check logging
        expect_log = 'Wrote manifest file successfully.'
        assert expect_log in caplog.text

    def test_fails(self, tmp_path, mocker, caplog):
        """Test logging when opening manifest fails."""
        mocker.patch('pathlib.Path.write_text', side_effect=OSError)
        with execute_in_dir(tmp_path):
            _write_manifest_file(('should not write',))
        expect_log = 'Failed to write manifest file.'
        assert not (tmp_path/'manifest').is_file()
        assert expect_log in caplog.text

    def test_raises(self):
        """Check that only OSError is caught cleanly."""
        with raises_test_exception('pathlib.Path.write_text'):
            _write_manifest_file(('should not write',))


class TestPrerunClean:
    """Tests for the prerun_clean function."""

    executables = {  # Some compiled executables
        'refcalc-010203-040506': 'refcalc exe contents',
        'rfactor-999999-888888': 'rfactor exe contents',
        'search-000001-000002': 'search exe contents',
        'superpos-123456-123456': 'superpos exe contents',
        # And also their windows versions
        'refcalc-010203-040506.exe': 'refcalc exe windows contents',
        'rfactor-999999-888888.exe': 'rfactor exe windows contents',
        'search-000001-000002.exe': 'search exe windows contents',
        'superpos-123456-123456.exe': 'superpos windows exe contents',
        }
    old_out_files = {
        # Some old-style _OUT files
        'POSCAR_OUT': 'POSCAR_OUT contents',
        'PARAMETERS_OUT': 'PARAMETERS_OUT contents',
        'VIBROCC_OUT': 'VIBROCC_OUT contents',
        'R_OUT_refcalc_R=0.123': 'R_OUT contents',
        }
    silent_fail = {  # Those for which we don't complain
        'fortran-compile.log': '',
        }
    calc_tree = {
        DEFAULT_WORK_HISTORY: {'subfolder':{}},
        DEFAULT_OUT: {'out_file': 'out file contents'},
        DEFAULT_SUPP: {'supp_file': 'supp file contents'},
        **old_out_files,
        **executables,
        **silent_fail,
        }

    @fixture(name='run')
    def fixture_run_prerun_clean(self, tmp_path, rpars, work_tree):
        """Execute prerun_clean at tmp_path."""
        def _run(nothing_to_move=False):
            expect, log_name, move_oldruns = work_tree
            prerun_clean(rpars, logname=log_name)
            clean = filesystem_to_dict(tmp_path)
            if not nothing_to_move:
                move_oldruns.assert_called_once_with(rpars, prerun=True)
            return expect, clean
        return _run

    @fixture(name='work_tree')
    def fixture_work_tree(self, tmp_path, mocker):
        """Prepare a directory with leftovers from a calc run."""
        # Mock move_oldruns to avoid having
        # to deal with its implementation
        mock_move = mocker.patch(f'{_MODULE}.move_oldruns')
        # Input files should be untouched, as is the current log
        surviving = {f: f'{f} contents' for f in ALL_INPUT_FILES}
        log_name = 'current_log.log'
        surviving[log_name] = 'current log contents'
        # Since we mock out move_oldruns, the old log
        # is not moved to a new workhistory subfolder
        surviving['old_log.log'] = 'old log contents'
        tree = {**self.calc_tree, **surviving}
        filesystem_from_dict(tree, tmp_path)
        with execute_in_dir(tmp_path):
            yield surviving, log_name, mock_move

    def check_fails_to_remove_file(self, file, run, mocker):
        """Check failure to remove file."""
        _path_unlink = Path.unlink
        def _fails_to_remove(path):
            if path.name == file:
                raise OSError
            _path_unlink(path)
        mocker.patch('pathlib.Path.unlink', new=_fails_to_remove)
        expect, clean = run()
        # The failing file should still be there
        expect[file] = self.calc_tree[file]
        assert clean == expect

    @parametrize(folder=(DEFAULT_WORK_HISTORY, DEFAULT_OUT, DEFAULT_SUPP))
    def test_fails_to_rmdir(self, folder, run, caplog, mocker):
        """Check complaints when deleting a folder fails."""
        _rmtree = shutil.rmtree
        def _rmtree_raises(path):
            if Path(path).name == folder:
                raise OSError
            _rmtree(path)
        mocker.patch('shutil.rmtree', new=_rmtree_raises)
        expect, clean = run()
        # The stuff that failed to be removed should stay
        expect[folder] = self.calc_tree[folder]
        assert clean == expect
        expect_log = ('Failed to clear', folder)
        assert all(m in caplog.text for m in expect_log)

    def test_fails_to_move_oldruns(self, run, work_tree, caplog):
        """Check warnings when failing to sort previous executions."""
        *_, mock_move = work_tree
        mock_move.side_effect = OSError
        expect, clean = run()
        expect_log = 'old files may be lost'
        assert clean == expect
        assert expect_log in caplog.text

    def test_dont_catch_move_oldruns_bugs(self, run, work_tree):
        """Check that no generic exception is masked."""
        *_, mock_move = work_tree
        mock_move.side_effect = Exception
        with pytest.raises(Exception):
            run()

    @parametrize(file=executables)
    def test_fails_to_remove_executable(self, file, run, caplog, mocker):
        """Check warnings when failing to remove a compiled executable."""
        caplog.set_level(0)  # We emit debug
        self.check_fails_to_remove_file(file, run, mocker)
        expect_log = f'Failed to delete file {file}'
        assert expect_log in caplog.text

    @parametrize(file=silent_fail)
    def test_fails_to_remove_log(self, file, run, caplog, mocker):
        """Check there no complaints when failing to remove some log files."""
        caplog.set_level(0)  # Collect all messages
        self.check_fails_to_remove_file(file, run, mocker)
        assert not caplog.text

    @parametrize(file=old_out_files)
    def test_fails_to_remove_old_out_file(self, file, run, caplog, mocker):
        """Check warnings when failing to remove an old-style _OUT file."""
        self.check_fails_to_remove_file(file, run, mocker)
        expect_log = f'Failed to delete previous {file}'
        assert expect_log in caplog.text

    def test_no_need_to_move_oldruns(self, run, work_tree, tmp_path, caplog):
        """Check no complaints if there's no old stuff to archive."""
        expect, *_, move_oldruns = work_tree
        old_logs = ('old_log.log', *self.silent_fail)
        for log in old_logs:
            (tmp_path/log).unlink()
            try:
                del expect[log]
            except KeyError:
                pass
        _, clean = run(nothing_to_move=True)
        assert clean == expect
        assert not caplog.text
        move_oldruns.assert_not_called()

    def test_not_an_executable(self, run, tmp_path, caplog):
        """Check that non-executable files are retained."""
        not_an_exe = tmp_path/'refcalc-not_a_time'
        not_an_exe.touch()
        expect, clean = run()
        expect[not_an_exe.name] = ''
        assert clean == expect
        assert not caplog.text

    @parametrize(folder=(DEFAULT_WORK_HISTORY, DEFAULT_OUT, DEFAULT_SUPP))
    # pylint: disable-next=too-many-arguments  # 4/6 are fixtures
    def test_skips_missing_folder(self, folder, work_tree,
                                  tmp_path, run, caplog):
        """Check no complaints when folder is not found."""
        expect, *_ = work_tree
        shutil.rmtree(tmp_path/folder)
        _, clean = run()
        assert expect == clean
        assert not caplog.text

    @parametrize(file=silent_fail)
    # pylint: disable-next=too-many-arguments  # 4/6 are fixtures
    def test_skips_missing_log(self, file, work_tree, tmp_path, run, caplog):
        """Check no complaints when folder is not found."""
        expect, *_ = work_tree
        (tmp_path/file).unlink()
        _, clean = run()
        assert expect == clean
        assert not caplog.text

    def test_success(self, run, caplog):
        """Check the successful cleanup of a work directory."""
        expect, clean = run()
        assert clean == expect
        assert not caplog.text


class TestOrganizeWorkdir:
    """Tests for the organize_workdir function."""

    @fixture(name='mock_implementation')
    def fixture_mock_implementation(self, mocker):
        """Replace implementation details of organize_workdir with mocks."""
        helpers = (
            '_collect_delta_files',
            '_zip_subfolders',
            '_collect_supp_contents',
            '_collect_out_contents',
            )
        return {helper: mocker.patch(f'{_MODULE}.{helper}')
                for helper in helpers}

    @fixture(name='run')
    def fixture_run(self, rpars, workdir):
        """Execute organize_workdir with rpars and path=workdir."""
        def _run(**kwargs):
            organize_workdir(rpars, workdir, **kwargs)
        return _run

    def test_archive_and_delete(self, rpars, run, mock_implementation, mocker):
        """Check that the implementation is executed as expected."""
        rpars.TENSOR_INDEX = mocker.MagicMock()
        rpars.ZIP_COMPRESSION_LEVEL = level = mocker.MagicMock()
        run(delete_unzipped=False, tensors=True, deltas=True)
        calls = {
            '_collect_delta_files': (mocker.call(rpars.TENSOR_INDEX),),
            '_zip_subfolders': (
                mocker.call(DEFAULT_TENSORS, True, False, level),
                mocker.call(DEFAULT_DELTAS, True, False, level),
                ),
            '_collect_supp_contents': (mocker.call(),),
            '_collect_out_contents': (mocker.call(),),
            }
        for helper, mock in mock_implementation.items():
            mock.assert_has_calls(calls[helper])

    def test_archive_tensors(self, run, mock_implementation, mocker):
        """Check that the implementation is executed as expected."""
        run(delete_unzipped=False, tensors=True, deltas=False)
        calls = {
            '_collect_delta_files': (
                mocker.call(0),  # TENSOR_INDEX is None by default
                ),
            '_zip_subfolders': (
                mocker.call(DEFAULT_TENSORS, True, False, 2),
                ),
            '_collect_supp_contents': (mocker.call(),),
            '_collect_out_contents': (mocker.call(),),
            }
        for helper, mock in mock_implementation.items():
            mock.assert_has_calls(calls[helper])

    def test_archive_deltas(self, run, mock_implementation, mocker):
        """Check that the implementation is executed as expected."""
        run(delete_unzipped=False, tensors=False, deltas=True)
        calls = {
            '_collect_delta_files': (
                mocker.call(0),  # TENSOR_INDEX is None by default
                ),
            '_zip_subfolders': (
                mocker.call(DEFAULT_DELTAS, True, False, 2),
                ),
            '_collect_supp_contents': (mocker.call(),),
            '_collect_out_contents': (mocker.call(),),
            }
        for helper, mock in mock_implementation.items():
            mock.assert_has_calls(calls[helper])

    def test_delete_only(self, run, mock_implementation, mocker):
        """Check that the implementation is executed as expected."""
        run(delete_unzipped=True, tensors=False, deltas=False)
        calls = {
            '_collect_delta_files': (
                mocker.call(0),  # TENSOR_INDEX is None by default
                ),
            '_zip_subfolders': (
                mocker.call(DEFAULT_TENSORS, False, True, 2),
                mocker.call(DEFAULT_DELTAS, False, True, 2),
                ),
            '_collect_supp_contents': (mocker.call(),),
            '_collect_out_contents': (mocker.call(),),
            }
        for helper, mock in mock_implementation.items():
            mock.assert_has_calls(calls[helper])


class TestOrganizeAllWorkDirectories:
    """Tests for the _organize_all_work_directories function."""

    @fixture(name='organize')
    def mock_organize_workdir(self, mocker):
        """Replace the organize_workdir function with a mock."""
        return mocker.patch(f'{_MODULE}.organize_workdir')

    def test_domains(self, rpars, organize, mocker):
        """Check the expected calls for a multi-domain calculation."""
        rpars.closePdfReportFigs = mocker.MagicMock()
        main_call = mocker.call(rpars=rpars,
                                path='',
                                delete_unzipped=True,
                                tensors=False,
                                deltas=False)
        # Add some fake domains
        domain_1 = mocker.MagicMock()
        domain_1.workdir = mocker.MagicMock()
        domain_1.rp = Rparams()  # No Tensors/Deltas in manifest
        call_1 = mocker.call(rpars=domain_1.rp,
                             path=domain_1.workdir,
                             delete_unzipped=True,
                             tensors=False,
                             deltas=False)
        domain_2 = mocker.MagicMock()
        domain_2.workdir = mocker.MagicMock()
        domain_2.rp = Rparams()
        domain_2.rp.manifest = {DEFAULT_TENSORS, DEFAULT_DELTAS}
        call_2 = mocker.call(rpars=domain_2.rp,
                             path=domain_2.workdir,
                             delete_unzipped=True,
                             tensors=True,
                             deltas=True)
        rpars.domainParams = domain_1, domain_2
        calls = main_call, call_1, call_2

        _organize_all_work_directories(rpars)
        rpars.closePdfReportFigs.assert_called_once()
        organize.call_count == len(calls)
        organize.assert_has_calls(calls)

    def test_no_domains(self, rpars, organize, mocker):
        """Check the expected calls for a single-domain calculation."""
        rpars.closePdfReportFigs = mocker.MagicMock()
        rpars.manifest.add(DEFAULT_TENSORS)  # But not deltas

        kwargs = {
            'rpars': rpars,
            'path': '',
            'delete_unzipped': True,
            'tensors': True,  # Added to manifest
            'deltas': False,
            }
        _organize_all_work_directories(rpars)
        organize.assert_called_once_with(**kwargs)

    def test_rpars_empty(self, rpars, organize):
        """Check the expected calls with an empty Rparams."""
        _organize_all_work_directories(rpars)
        kwargs = {
            'rpars': rpars,
            'path': '',
            'delete_unzipped': True,
            'tensors': False,
            'deltas': False,
            }
        organize.assert_called_once_with(**kwargs)


class TestWriteFinalLogMessages:
    """Tests for the _write_final_log_messages function."""

    @staticmethod
    def check_has_records(caplog, records):
        """Ensure caplog has only certain records."""
        logged = tuple(r.getMessage() for r in caplog.records)
        print(logged)
        print(records)
        assert len(logged) == len(records)
        for log, expect in zip(logged, records):
            if isinstance(expect, str):
                assert log == expect
            else:
                assert expect.match(log)

    @fixture(name='rpars_filled')
    def fixture_rpars_filled(self, rpars, mocker):
        """Return an Rparams with relevant attributes set."""
        rpars.timer = mocker.MagicMock()
        rpars.timer.how_long.return_value = '1h 30m'
        rpars.runHistory = [1, 2, 3]
        rpars.stored_R = {'refcalc': [0.5, 0.4, 0.3],
                          'superpos': None}
        rpars.checklist = ['Check convergence', 'Verify inputs']
        return rpars

    def test_crashed_early(self, rpars, caplog):
        """Check logging messages when cleanup is called early."""
        caplog.set_level(0)  # All messages
        rpars.timer = None   # Like cleanup if called with None
        _write_final_log_messages(rpars)
        expect = (
            re.compile(r'\nFinishing execution at .*'
                       r'\nTotal elapsed time: unknown\n'),
            '',
            )
        self.check_has_records(caplog, expect)

    def test_domains(self, rpars_filled, caplog):
        """Check logging messages for a multi-domain calculation."""
        caplog.set_level(0)  # All messages
        rpars_filled.domainParams.append(1.234)
        _write_final_log_messages(rpars_filled)
        expect = (
            re.compile(r'\nFinishing execution at .*'
                       r'\nTotal elapsed time: 1h 30m\n'),
            'Executed segments: 1 2 3',
            'Final R (refcalc): 0.5000 (0.4000 / 0.3000)',
            re.compile(r'Domain.*bookkeeper.*--archive.'),
            '',
            '# The following issues should be checked before starting again:',
            '- Check convergence',
            '- Verify inputs',
            '',
            )
        self.check_has_records(caplog, expect)

    def test_no_checklist(self, rpars_filled, caplog):
        """Check that no checklist-related messages are emitted."""
        caplog.set_level(0)  # All messages
        rpars_filled.checklist = []
        rpars_filled.stored_R['refcalc'] = [0.5, 0, 0.5]
        _write_final_log_messages(rpars_filled)
        expect = (
            re.compile(r'\nFinishing execution at .*'
                       r'\nTotal elapsed time: 1h 30m\n'),
            'Executed segments: 1 2 3',
            'Final R (refcalc): 0.5000',
            '',
            )

    def test_no_domains(self, rpars_filled, caplog):
        """Check logging messages for a single-domain calculation."""
        caplog.set_level(0)  # All messages
        _write_final_log_messages(rpars_filled)
        expect = (
            re.compile(r'\nFinishing execution at .*'
                       r'\nTotal elapsed time: 1h 30m\n'),
            'Executed segments: 1 2 3',
            'Final R (refcalc): 0.5000 (0.4000 / 0.3000)',
            '',
            '# The following issues should be checked before starting again:',
            '- Check convergence',
            '- Verify inputs',
            '',
            )
        self.check_has_records(caplog, expect)


class TestOrganizeSuppOut:
    """Tests for the _organize_supp_out function."""

    # Stuff that goes to SUPP
    to_supp = {f: f'{f} contents' for f in _SUPP_FILES}
    to_supp.update({d: {'test': f'test file in {d}'} for d in _SUPP_DIRS})
    to_supp[f'not-a-{LOG_PREFIX}-log.log'] = 'some log contents'

    # Stuff that goes to OUT
    to_out = {f: f'{f} contents' for f in _OUT_FILES}
    to_out.update({f'{f}_OUT_blah': f'{f} contents'
                   for f in ('PARAMETERS', 'POSCAR', 'VIBROCC', 'R')})
    to_out['PARAMETERS'] = 'was a PARAMETERS file in work'

    # Stuff that is untouched
    stays = {
        f'{LOG_PREFIX}-main-log.log': 'fake calc log contents',
        'a-log-for-a-compile-step.log': 'fake compilation log',
        'another_folder': {'with-contents.txt': ''},
        }
    tree = {**stays, **to_out, **to_supp}

    # PARAMETERS is renamed when copying
    to_out['PARAMETERS_OUT'] = to_out.pop('PARAMETERS')

    @fixture(name='work_tree')
    def fixture_work_tree(self, workdir):
        """Create files to be moved to SUPP/OUT at workdir."""
        filesystem_from_dict(self.tree, workdir)
        return workdir

    @fixture(name='run')
    def fixture_run(self, work_tree):
        """Execute _organize_supp_out."""
        def _run(*which):
            expect = {**self.tree,
                      DEFAULT_SUPP: self.to_supp,
                      DEFAULT_OUT: self.to_out}
            with execute_in_dir(work_tree):
                for func in which:
                    func()
            organized = filesystem_to_dict(work_tree)
            return expect, organized
        return _run

    _ordered = {
        'out,supp': (_collect_out_contents, _collect_supp_contents),
        'supp,out': (_collect_supp_contents, _collect_out_contents),
        }

    @parametrize(which=_ordered.values(), ids=_ordered)
    def test_success(self, which, run, caplog):
        """Check a successful copy of contents to SUPP/OUT."""
        expect, organized = run(*which)
        assert organized == expect
        assert not caplog.text

    def test_fails_to_rename_parameters(self, run, caplog, mocker):
        """Check warnings when renaming OUT/PARAMETERS fails."""
        _copy = shutil.copy2
        def _rename_fails(src, dst):
            # pylint: disable-next=magic-value-comparison
            if Path(src).name == 'PARAMETERS':
                raise OSError
            _copy(src, dst)
        mocker.patch('shutil.copy2', _rename_fails)
        expect, organized = run(_collect_out_contents)
        expect_log = 'Error renaming'
        del expect[DEFAULT_SUPP]  # We only move OUT files
        expect[DEFAULT_OUT].pop('PARAMETERS_OUT')  # Not renamed
        assert organized == expect
        assert expect_log in caplog.text


class TestCopyFilesAndDirectories:
    """Tests for the _copy_files_and_directories helper."""

    tree = {
        'test.txt': 'test file contents',
        'test_directory': {'test_subfolder': {'test_file': 'another file'}},
        }

    @fixture(name='run')
    def fixture_run(self, workdir):
        """Run _copy_files_and_directories in workdir."""
        def _run(*args, make_tree=True):
            if make_tree:
                filesystem_from_dict(self.tree, workdir)
            with execute_in_dir(workdir):
                _copy_files_and_directories(*args)
            return filesystem_to_dict(workdir)
        return _run

    def test_copy_directories(self, run, workdir):
        """Check a successful copy of directories."""
        target = workdir / 'test_target'
        target.mkdir()
        dirs = [workdir/d for d, c in self.tree.items() if isinstance(c, dict)]
        copied = run([], dirs, target)
        expect = self.tree.copy()
        expect[target.name] = {d.name: self.tree[d.name] for d in dirs}
        assert copied == expect

    def test_copy_files(self, run, workdir):
        """Check a successful copy of files."""
        target = workdir / 'test_target'
        files = [workdir/f for f, c in self.tree.items() if isinstance(c, str)]
        copied = run(files, [], target)
        expect = self.tree.copy()
        expect[target.name] = {f.name: self.tree[f.name] for f in files}
        assert copied == expect

    def test_copy_non_existing(self, run, workdir):
        """Check no complaints when copying non-existing files/folders."""
        target = workdir / 'test_target'
        expect = {target.name: {}}
        files = (workdir/'non_existing_file.txt',)
        dirs = (workdir/'non_existing_dir',)
        copied = run(files, dirs, target, make_tree=False)
        assert copied == expect

    def test_copy_fails(self, run, workdir, caplog, mocker):
        """Check logging messages when copying fails."""
        copy_file = mocker.patch('shutil.copy2', side_effect=OSError)
        copy_dir = mocker.patch(f'{_MODULE}.copytree_exists_ok',
                                side_effect=OSError)
        expect = self.tree.copy()
        target = workdir / 'test_target'
        expect[target.name] = {}
        files = [workdir/f for f, c in self.tree.items() if isinstance(c, str)]
        dirs = [workdir/d for d, c in self.tree.items() if isinstance(c, dict)]
        copied = run(files, dirs, target)
        assert copied == expect
        assert copy_file.call_count == len(files)
        assert copy_dir.call_count == len(dirs)
        log_msg_file = 'Error moving test_target file {}'
        log_msg_dir = 'Error moving test_target directory {}'
        log = caplog.text
        assert all(log_msg_file.format(f.name) in log for f in files)
        assert all(log_msg_dir.format(d.name) in log for d in dirs)

    def test_mkdir_fails(self, run, workdir, caplog, mocker):
        """Check warnings when making the target fails."""
        files = [workdir/f for f, c in self.tree.items() if isinstance(c, str)]
        target = workdir / 'test_target'
        filesystem_from_dict(self.tree, workdir)  # BEFORE patching
        mocker.patch('pathlib.Path.mkdir', side_effect=OSError)
        expect_log = 'Error creating test_target folder'
        unchanged = run(files, (), target, make_tree=False)
        assert unchanged == self.tree
        assert expect_log in caplog.text


class TestCollectDeltaFiles:
    """Tests for the _collect_delta_files function."""

    @fixture(name='run')
    def fixture_run(self, workdir):
        """Execute _collect_delta_files in workdir."""
        def _run(tensor_index, tree=None):
            if tree:
                filesystem_from_dict(tree, workdir)
            with execute_in_dir(workdir):
                _collect_delta_files(tensor_index)
            return filesystem_to_dict(workdir)
        return _run

    def test_success(self, run, caplog):
        """Check results of a successful execution."""
        tree = {f: f'{f} contents'
                for f in ('DEL_1', 'DEL_two', 'DEL_tafile')}
        collected = run(123, tree)
        expect = {'Deltas': {'Deltas_123': tree}}
        assert collected == expect
        assert not caplog.text

    def test_nothing_to_collect(self, run, workdir, caplog):
        """Check that no folder is created if no files should be collected."""
        run(9)
        assert not (workdir/'Deltas'/'Deltas_009').exists()
        assert not caplog.text

    def test_delta_folder_exists(self, run, workdir, caplog):
        """Check no complaints if the target folder is there already."""
        (workdir/'Deltas'/'Deltas_123').mkdir(parents=True)
        self.test_success(run, caplog)

    def test_mkdir_fails(self, run, caplog, mocker):
        """Check warnings when failing to create the Deltas directory."""
        mocker.patch('pathlib.Path.mkdir', side_effect=OSError)
        tree = {'DEL_1': 'contents'}
        unchanged = run(5, tree)
        expect_log = 'Failed to create'
        assert unchanged == tree  # No changes
        assert expect_log in caplog.text

    def test_file_copy_fails(self, run, caplog, mocker):
        """Check complaints when moving files fails."""
        mock_move = mocker.patch('pathlib.Path.replace', side_effect=OSError)
        tree = {'DEL_1': '', 'DEL_2': ''}
        collected = run(75, tree)
        expect = {'Deltas': {'Deltas_075': {}},  # Fails to move all
                  **tree}
        expect_log = 'Error moving Delta'
        assert collected == expect
        assert mock_move.call_count == len(tree)
        assert expect_log in caplog.text


class TestZipSubfolders:
    """Tests for the _zip_subfolders function."""

    folder = 'test_zipping'
    contents = {
        'test_zipping_001': {'file_1': 'file_1 contents'},
        'test_zipping_999': {'file_999': 'file_999 contents'},
        'not-a-directory': 'some stray file in test_zipping',
        }
    tree = {folder: contents}

    @property
    def packed_all(self):
        """Return the expected tree after packing self.tree."""
        files = {f: c for f, c in self.contents.items() if isinstance(c, str)}
        files.update({f'{f}.zip': f'{f}.zip is a binary file'
                      for f in self.contents
                      if f not in files})
        return {self.folder: files}

    @fixture(name='work_tree')
    def fixture_work_tree(self, workdir):
        """Create a temporary tree at workdir."""
        filesystem_from_dict(self.tree, workdir)
        return workdir

    @fixture(name='run')
    def fixture_run(self, work_tree):
        """Execute _zip_subfolders in workdir."""
        def _run(**kwargs):
            with execute_in_dir(work_tree):
                _zip_subfolders(at_path=self.folder, **kwargs)
            return filesystem_to_dict(work_tree)
        return _run

    def test_delete_only(self, run, caplog):
        """Check no deletion of unzipped directories."""
        clean = run(delete_unzipped=True,
                    archive=False,
                    compression_level=2)
        expect = {self.folder: {
            f: c for f, c in self.contents.items() if isinstance(c, str)
            }}
        assert clean == expect
        assert not caplog.text

    def test_delete_fails(self, run, caplog, mocker):
        """Check complaints when deleting an unzipped folder fails."""
        mocker.patch('shutil.rmtree', side_effect=OSError)
        clean = run(delete_unzipped=True,
                    archive=True,
                    compression_level=2)
        expect = self.packed_all
        expect[self.folder].update(self.contents)
        expect_log = 'Error deleting unzipped'
        assert clean == expect
        assert expect_log in caplog.text

    def test_do_nothing(self, run, caplog):
        """Check correct behavior for not packing nor deleting."""
        unchanged = run(delete_unzipped=False,
                        archive=False,
                        compression_level=2)
        assert unchanged == self.tree
        assert not caplog.text

    def test_folder_does_not_exist(self, work_tree, caplog):
        """Check that nothing happens if called with a non-existing folder."""
        _zip_subfolders('does-not-exist', True, True, 2)
        unchanged = filesystem_to_dict(work_tree)
        assert unchanged == self.tree
        assert not caplog.text

    def test_pack_and_delete(self, run, workdir, caplog):
        """Check packing and deleting of both Tensors and Deltas."""
        clean = run(delete_unzipped=True,
                    archive=True,
                    compression_level=2)
        assert clean == self.packed_all
        assert not caplog.text

        # Check the contents of an archive
        archived = next(f for f, c in self.contents.items()
                        if isinstance(c, dict))
        target = workdir/'extracted'
        target.mkdir()
        archive_p = (workdir/self.folder/archived).with_suffix('.zip')
        shutil.unpack_archive(archive_p, target)
        extracted = filesystem_to_dict(target)
        assert extracted == self.contents[archived]

    def test_pack_fails(self, run, caplog, mocker):
        """Check nothing is deleted if packing is not successful."""
        mocker.patch(f'{_MODULE}.ZipFile', side_effect=OSError)
        unchanged = run(delete_unzipped=True,  # not honored
                        archive=True,
                        compression_level=2)
        expect_log = 'Error packing'
        assert unchanged == self.tree
        assert expect_log in caplog.text

    def test_pack_only(self, run, caplog):
        """Check no deletion of unzipped directories."""
        clean = run(delete_unzipped=False,
                    archive=True,
                    compression_level=2)
        expect = self.packed_all
        expect[self.folder].update(self.contents)
        assert clean == expect
        assert not caplog.text

    @pytest.mark.xfail(
        reason=('This supposedly invalid folder name is collected anyway. It '
                'was meant to fail the match.end check but it is not really '
                'clear what that check guards')
        )
    def test_invalid_name_format(self, run, workdir, caplog):
        """Check that a folder with invalid format is not processed."""
        tree = {self.folder: {'test_zipping_001-invalid-fmt': {}}}
        filesystem_from_dict(tree, workdir)
        clean = run(delete_unzipped=True,
                    archive=True,
                    compression_level=2)
        expect = self.packed_all
        expect[self.folder].update(tree[self.folder])
        assert clean == expect
        assert not caplog.text

    def test_not_a_valid_folder(self, run, workdir, caplog):
        """Check that a stray subfolder is not packed."""
        tree = {self.folder: {'not-a-test_zipping_001': {}}}
        filesystem_from_dict(tree, workdir)
        clean = run(delete_unzipped=True,
                    archive=True,
                    compression_level=2)
        expect = self.packed_all
        expect[self.folder].update(tree[self.folder])
        assert clean == expect
        assert not caplog.text


class TestZipFolder:
    """Tests for the _zip_folder function."""

    def test_does_not_pack_itself(self, tmp_path):
        """Check that an existing archive is not included."""
        folder = 'tst_folder'
        contents = {'file': 'file contents'}
        filesystem_from_dict({folder: contents}, tmp_path)
        # Create also an existing zip file
        archive = (tmp_path/folder).with_suffix('.zip')
        with ZipFile(archive, 'w') as zip_:
            zip_.writestr('pre-existing', 'contents were already there')

        _zip_folder(tmp_path/folder, 2)

        # Check that the archive contains what it should
        extract_into = tmp_path/'extracted'
        shutil.unpack_archive(archive, extract_into)
        extracted = filesystem_to_dict(extract_into)
        expect = {**contents,
                  'pre-existing': 'contents were already there'}
        assert extracted == expect
