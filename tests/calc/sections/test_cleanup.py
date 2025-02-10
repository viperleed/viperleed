"""Tests for module cleanup of viperleed.calc.section."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-02-03'
__license__ = 'GPLv3+'

import os
from pathlib import Path
import shutil

from pytest_cases import fixture
from pytest_cases import parametrize
import pytest

from viperleed.calc.classes.rparams import Rparams
from viperleed.calc.constants import DEFAULT_OUT
from viperleed.calc.constants import DEFAULT_SUPP
from viperleed.calc.constants import DEFAULT_WORK_HISTORY
from viperleed.calc.constants import ORIGINAL_INPUTS_DIR_NAME
from viperleed.calc.lib.context import execute_in_dir
from viperleed.calc.sections.calc_section import ALL_INPUT_FILES
from viperleed.calc.sections.cleanup import OPTIONAL_INPUT_FILES
from viperleed.calc.sections.cleanup import prerun_clean
from viperleed.calc.sections.cleanup import preserve_original_inputs
from viperleed.calc.sections.cleanup import _write_manifest_file

from ...helpers import filesystem_from_dict
from ...helpers import filesystem_to_dict
from ...helpers import raises_test_exception

_MODULE = 'viperleed.calc.sections.cleanup'


@fixture(name='rpars')
def fixture_rpars():
    """Return an empty Rparams."""
    return Rparams()


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
        'refcalc-010203-040506': 'refcalc exe contents',                        # TODO: Not OK for Windows. Add tests.
        'rfactor-999999-888888': 'rfactor exe contents',
        'search-000001-000002': 'search exe contents',
        'superpos-123456-123456': 'superpos exe contents',
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
        mock_move.side_effect = Exception
        expect, clean = run()
        expect_log = 'old files may be lost'
        assert clean == expect
        assert expect_log in caplog.text

    @parametrize(file=executables)
    def test_fails_to_remove_executable(self, file, run, caplog, mocker):
        """Check warnings when failing to remove a compiled executable."""
        _os_remove = os.remove
        def _fails_to_remove(path):
            if Path(path).name == file:
                raise OSError
            _os_remove(path)
        mocker.patch('os.remove', new=_fails_to_remove)
        caplog.set_level(0)  # We emit debug
        expect, clean = run()
        # The failing file should still be there
        expect[file] = self.calc_tree[file]
        expect_log = f'Failed to delete file {file}'
        assert clean == expect
        assert expect_log in caplog.text

    @parametrize(file=silent_fail)
    def test_fails_to_remove_log(self, file, run, caplog, mocker):
        """Check there no complaints when failing to remove some log files."""
        _os_remove = os.remove
        def _fails_to_remove(path):
            if Path(path).name == file:
                raise OSError
            _os_remove(path)
        mocker.patch('os.remove', new=_fails_to_remove)
        caplog.set_level(0)  # Collect all messages
        expect, clean = run()
        # The failing file should still be there
        expect[file] = self.calc_tree[file]
        assert clean == expect
        assert not caplog.text

    @parametrize(file=old_out_files)
    def test_fails_to_remove_old_out_file(self, file, run, caplog, mocker):
        """Check warnings when failing to remove an old-style _OUT file."""
        _path_unlink = Path.unlink
        def _fails_to_remove(path):
            if path.name == file:
                raise OSError
            _path_unlink(path)
        mocker.patch('pathlib.Path.unlink', new=_fails_to_remove)
        expect, clean = run()
        # The failing file should still be there
        expect[file] = self.calc_tree[file]
        expect_log = f'Failed to delete previous {file}'
        assert clean == expect
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
