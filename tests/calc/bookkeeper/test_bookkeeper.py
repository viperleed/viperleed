"""Tests for module viperleed.calc.bookkeeper."""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-08-02'
__license__ = 'GPLv3+'

import ast
from enum import Enum
import functools
import logging
from pathlib import Path
import shutil
import time

import pytest
from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.calc import DEFAULT_HISTORY
from viperleed.calc import DEFAULT_WORK_HISTORY
from viperleed.calc import LOG_PREFIX
from viperleed.calc import ORIGINAL_INPUTS_DIR_NAME
from viperleed.calc.bookkeeper.bookkeeper import _FROM_ROOT
from viperleed.calc.bookkeeper.bookkeeper import Bookkeeper
from viperleed.calc.bookkeeper.bookkeeper import BookkeeperExitCode
from viperleed.calc.bookkeeper.history.constants import HISTORY_INFO_NAME
from viperleed.calc.bookkeeper.history.entry.entry import _DISCARDED
from viperleed.calc.bookkeeper.log import BOOKIE_LOGFILE
from viperleed.calc.bookkeeper.mode import BookkeeperMode
from viperleed.calc.sections.cleanup import DEFAULT_OUT
from viperleed.calc.sections.cleanup import DEFAULT_SUPP
from viperleed.calc.sections.cleanup import PREVIOUS_LABEL

from ...helpers import execute_in_dir
from ...helpers import not_raises
from ...helpers import make_obj_raise
from ...helpers import raises_exception
from ...helpers import raises_test_exception
from .conftest import MOCK_INPUT_CONTENT
from .conftest import MOCK_OUT_CONTENT
from .conftest import MOCK_TIMESTAMP
from .conftest import MOCK_STATE_FILES
from .conftest import MOCK_WORKHISTORY
from .conftest import NOTES_TEST_CONTENT


_UPDATE_METHOD = 'update_from_cwd'
raises_oserror = functools.partial(raises_exception, exc=OSError)
make_raise_oserror = functools.partial(make_obj_raise, exc=OSError)
not_raises_oserror = functools.partial(not_raises, exc=OSError)



@fixture(name='after_archive')
def fixture_after_archive(after_calc_run):
    """Yield a temporary directory for testing the bookkeeper."""
    bookkeeper, *_ = after_calc_run
    bookkeeper.run(mode=BookkeeperMode.ARCHIVE)
    bookkeeper.update_from_cwd(silent=True)
    return after_calc_run


@fixture(name='before_calc_run')
def fixture_before_calc_run(tmp_path):
    """Yield a temporary directory for testing the bookkeeper.

    This represents a new calculation, i.e., before any viperleed calc or
    bookkeeper run.

    Parameters
    ----------
    tmp_path : fixture
        Path to a temporary directory.

    Yields
    ------
    bookkeeper : Bookkeeper
        A bookkeeper instance ready for running in tmp_path.
    """
    # create mock input files
    for file in MOCK_STATE_FILES:
        (tmp_path / file).write_text(MOCK_INPUT_CONTENT)

    bookkeeper = Bookkeeper(cwd=tmp_path)
    bookkeeper.update_from_cwd(silent=True)
    with execute_in_dir(tmp_path):
        yield bookkeeper
    # It would be nice to clean up, but the following line causes
    # a PermissionError. Likely because of logging keeping a hold
    # of the bookkeeper.log file
    # shutil.rmtree(tmp_path)


@fixture(name='history_path')
def fixture_history_path(mock_tree_after_calc_run):
    """Return the path to a history subfolder of `mock_tree_after_calc_run`."""
    return mock_tree_after_calc_run / DEFAULT_HISTORY


@fixture(name='history_run_path')
def fixture_history_run_path(history_path):
    """Return the path to a history run subfolder of `history_path`."""
    return history_path / f't000.r001_{MOCK_TIMESTAMP}'


class _TestBookkeeperRunBase:
    """Base class for checking correct execution of bookkeeper."""

    mode = None

    def check_history_exists(self, bookkeeper, history_run_path):
        """Test that history_path and directory/history.info exist."""
        assert history_run_path.is_dir()
        assert (bookkeeper.cwd / HISTORY_INFO_NAME).exists()

    def check_history_folder_empty(self, bookkeeper, *_):
        """Test that an empty history folder exists, except bookkeeper.log."""
        cwd = bookkeeper.cwd
        assert (cwd / DEFAULT_HISTORY).exists()
        history_contents = (f for f in (cwd / DEFAULT_HISTORY).iterdir()
                            if f.name != BOOKIE_LOGFILE)
        assert not any(history_contents)

    def check_no_warnings(self, caplog, contain=None, exclude_msgs=()):
        """Check that there are no warnings or errors."""
        if contain is None:
            def _record_faulty(record):
                return record.levelno >= logging.WARNING
        else:
            def _record_faulty(record):
                return (record.levelno >= logging.WARNING
                        and contain in str(record))
        assert not any(
            _record_faulty(rec)
            for rec in caplog.records
            if not any(exclude in str(rec) for exclude in exclude_msgs)
            )

    def check_out_files_in_history(self, *run):
        """Check that the expected state files are stored in 'OUT'."""
        *_, history_run_path = run
        for file in MOCK_STATE_FILES:
            hist_file = history_run_path / DEFAULT_OUT / f'{file}_OUT'
            assert hist_file.is_file()
            assert MOCK_OUT_CONTENT in hist_file.read_text()

    def check_root_inputs_replaced_by_out(self, bookkeeper, *_):
        """Check that the input files in root come from OUT."""
        cwd = bookkeeper.cwd
        for file in MOCK_STATE_FILES:
            out_content = (cwd / file).read_text()
            assert MOCK_OUT_CONTENT in out_content

    def check_root_inputs_untouched(self, bookkeeper, *_):
        """Check the the original state files have not been moved."""
        cwd = bookkeeper.cwd
        for file in MOCK_STATE_FILES:
            assert (cwd / file).is_file()
            input_content = (cwd / file).read_text()
            assert MOCK_INPUT_CONTENT in input_content

    def check_root_is_clean(self, bookkeeper, *_):
        """Check that no calc output is present in the main directory."""
        cwd = bookkeeper.cwd
        for file in MOCK_STATE_FILES:
            assert not (cwd / f'{file}_ori').is_file()
        assert not (cwd / DEFAULT_SUPP).exists()
        assert not (cwd / DEFAULT_OUT).exists()
        assert not any(cwd.glob('*.log'))

    def _test_after_archive_base(self, after_archive):
        """Check that running bookkeeper after ARCHIVE does basic stuff."""
        bookkeeper, *_ = after_archive
        # bookkeeper should not think that it needs archiving
        assert not bookkeeper.archiving_required
        bookkeeper.run(mode=self.mode)
        bookkeeper.update_from_cwd(silent=True)
        self.check_history_exists(*after_archive)
        self.check_out_files_in_history(*after_archive)
        self.check_root_is_clean(*after_archive)
        # Check that the workhistory directories are
        # also where they should be
        history = bookkeeper.top_level_history_path
        workhistory = bookkeeper.cwd / DEFAULT_WORK_HISTORY
        assert not workhistory.is_dir()
        for ori_name, hist_name in MOCK_WORKHISTORY.items():
            if hist_name is None:  # File should be deleted
                continue
            moved_dir = history/hist_name
            moved_file = moved_dir/'file'
            assert moved_dir.is_dir()
            assert moved_file.is_file()
            assert ori_name in moved_file.read_text()

    def _test_after_calc_run_base(self, after_calc_run):
        """Check that running bookkeeper after calc does some basic stuff."""
        bookkeeper, *_ = after_calc_run
        # bookkeeper should think that it needs archiving
        assert bookkeeper.archiving_required
        bookkeeper.run(mode=self.mode)
        bookkeeper.update_from_cwd(silent=True)
        self.check_history_exists(*after_calc_run)
        self.check_out_files_in_history(*after_calc_run)

    def test_before_calc_run(self, before_calc_run):
        """Check that running bookkeeper before calc does almost nothing."""
        bookkeeper = before_calc_run
        # bookkeeper should not think that it needs archiving
        assert not bookkeeper.archiving_required
        bookkeeper.run(mode=self.mode)
        # Bookkeeper should not do anything (except for logging)
        self.check_history_folder_empty(before_calc_run)
        self.check_root_inputs_untouched(before_calc_run)


class TestBookkeeperArchive(_TestBookkeeperRunBase):
    """Tests for correct behavior of ARCHIVE bookkeeper runs."""

    mode = BookkeeperMode.ARCHIVE

    def _check_root_inputs_renamed_to_ori(self, bookkeeper, *_):
        """Check that the input files have now a _ori suffix."""
        cwd = bookkeeper.cwd
        for file in MOCK_STATE_FILES:
            ori_file = cwd / f'{file}_ori'
            assert ori_file.is_file()
            input_content = (cwd / file).read_text()
            assert MOCK_OUT_CONTENT in input_content

    def test_archive_after_calc_run(self, after_calc_run, caplog):
        """Check correct storage of history files in ARCHIVE mode."""
        self._test_after_calc_run_base(after_calc_run)
        self._check_root_inputs_renamed_to_ori(*after_calc_run)
        self.check_root_inputs_replaced_by_out(*after_calc_run)
        self.check_no_warnings(caplog)

    def test_archive_again(self, after_archive, caplog):
        """Bookkeeper ARCHIVE after ARCHIVE should not do anything."""
        bookkeeper, *_ = after_archive
        # Write stuff to files to check they are not overwritten
        cwd = bookkeeper.cwd
        sentinel_text = 'something else'
        for file in MOCK_STATE_FILES:
            (cwd / file).write_text(sentinel_text)
        super()._test_after_archive_base(after_archive)
        for file in MOCK_STATE_FILES:
            assert (cwd / file).read_text() == sentinel_text
        self.check_no_warnings(caplog)

    # pylint: disable-next=arguments-differ
    def test_before_calc_run(self, before_calc_run, caplog):
        """Check no archiving happens before calc runs."""
        super().test_before_calc_run(before_calc_run)
        self.check_no_warnings(caplog)


class TestBookkeeperClear(_TestBookkeeperRunBase):
    """Tests for correct behavior of CLEAR bookkeeper runs."""

    mode = BookkeeperMode.CLEAR

    # pylint: disable-next=arguments-differ
    def test_before_calc_run(self, before_calc_run, caplog):
        """Check correct overwriting of input files in CLEAR mode."""
        super().test_before_calc_run(before_calc_run)
        self.check_no_warnings(caplog)

    def test_clear_after_archive(self, after_archive, caplog):
        """Check behavior of CLEAR after ARCHIVE (e.g., manual call)."""
        self._test_after_archive_base(after_archive)
        self.check_root_inputs_replaced_by_out(*after_archive)
        self.check_no_warnings(caplog)

    def test_clear_after_calc_run(self, after_calc_run, caplog):
        """Check behavior of CLEAR after a non-ARCHIVEd calc run.

        This may happen, for example, if the previous (calc or
        bookkeeper) execution crashed.

        Parameters
        ----------
        after_calc_run: fixture
        caplog: fixture

        Returns
        -------
        None.
        """
        self._test_after_calc_run_base(after_calc_run)
        self.check_no_warnings(caplog)
        self.check_root_is_clean(*after_calc_run)

        # Original SHOULD NOT be replaced by output:
        # ARCHIVE does not run only if the run crashed,
        # in which case we don't want to overwrite
        self.check_root_inputs_untouched(*after_calc_run)


class TestBookkeeperDiscard(_TestBookkeeperRunBase):
    """Tests for correct behavior of DISCARD bookkeeper runs."""

    mode = BookkeeperMode.DISCARD

    # pylint: disable-next=arguments-differ
    def test_before_calc_run(self, before_calc_run, caplog):
        """Check correct overwriting of input files in CLEAR mode."""
        super().test_before_calc_run(before_calc_run)
        self.check_no_warnings(
            caplog,
            exclude_msgs=('Failed to mark last entry as discarded',),
            )

    def test_discard_after_archive(self, after_archive, caplog):
        """Check reverting of state when DISCARDing an ARCHIVEd calc run."""
        self._test_after_archive_base(after_archive)

        # Original be replaced by output                                        # TODO: this does something else!
        bookkeeper, *_ = after_archive
        cwd = bookkeeper.cwd
        for file in MOCK_STATE_FILES:
            out_content = (cwd / file).read_text()
            assert MOCK_INPUT_CONTENT in out_content
        # A 'DISCARDED' note should be in history.info
        assert _DISCARDED in bookkeeper.history_info.path.read_text()
        # Some fields are knowingly faulty,
        # but we can still DISCARD them.
        faulty_entry_logs = (
            'Found entry with',
            'Could not understand',
            'Faulty entry is',
            )
        self.check_no_warnings(caplog, exclude_msgs=faulty_entry_logs)

    def test_discard_after_calc_run(self, after_calc_run, caplog):
        """Check behavior of DISCARD after a non-ARCHIVEd calc run.

        This may happen, for example, if the previous (calc or
        bookkeeper) execution crashed.

        Parameters
        ----------
        after_calc_run: fixture
        caplog: fixture

        Returns
        -------
        None.
        """
        self._test_after_calc_run_base(after_calc_run)
        self.check_no_warnings(caplog, exclude_msgs=('discarded',))
        self.check_root_is_clean(*after_calc_run)

        # Original SHOULD NOT be replaced by output:
        # ARCHIVE does not run only if the run crashed,
        # in which case we don't want to overwrite
        self.check_root_inputs_untouched(*after_calc_run)

        # A 'DISCARDED' note should be in history.info
        bookkeeper, *_ = after_calc_run
        assert bookkeeper.history_info.last_entry_was_discarded
        assert _DISCARDED in bookkeeper.history_info.path.read_text()


class TestBookkeeperDiscardFull(_TestBookkeeperRunBase):
    """Tests for correct behavior of DISCARD_FULL bookkeeper runs."""

    mode = BookkeeperMode.DISCARD_FULL

    # pylint: disable-next=arguments-differ
    def test_before_calc_run(self, before_calc_run, caplog):
        """Check correct overwriting of input files in DISCARD_FULL mode.

        This should do the same as a normal DISCARD.

        Parameters
        ----------
        before_calc_run: fixture
        caplog: fixture

        Returns
        -------
        None.
        """
        super().test_before_calc_run(before_calc_run)
        self.check_no_warnings(
            caplog,
            exclude_msgs=('remove',)   # Catches two expected warnings
            )

    def test_discard_full_after_archive(self, after_archive, caplog):
        """Check correct behavior of DISCARD_FULL on an ARCHIVEd calc run.

        Should revert state to before the run and purge the history.

        Parameters
        ----------
        after_archive: fixture
        caplog: fixture

        Returns
        -------
        None.
        """
        bookkeeper, *_, history_run_path = after_archive
        last_entry = bookkeeper.history_info.last_entry
        bookkeeper.run(mode=self.mode)
        if not last_entry.can_be_removed:
            # Should prevent the removal of the history
            assert history_run_path.is_dir()
            return
        assert not history_run_path.is_dir()
        self.check_root_is_clean(*after_archive)
        self.check_root_inputs_untouched(*after_archive)
        self.check_no_warnings(caplog)

    def test_discard_full_after_calc_run(self, after_calc_run, caplog):
        """Check behavior of DISCARD_FULL on a non-ARCHIVEd calc run."""
        bookkeeper, *_, history_run_path = after_calc_run
        cwd = bookkeeper.cwd
        notes_in_history_info = (
            NOTES_TEST_CONTENT in (cwd / HISTORY_INFO_NAME).read_text()
            )
        last_entry = bookkeeper.history_info.last_entry
        bookkeeper.run(mode=self.mode)
        # Since the run was not archived, the history should be empty
        assert not history_run_path.is_dir()
        if notes_in_history_info:
            # pylint: disable-next=magic-value-comparison
            assert 'last entry in history.info has user notes' in caplog.text
        elif last_entry and not last_entry.can_be_removed:
            expected_logs = (
                'contains invalid fields that could not be interpreted',
                'contains fields with non-standard format',
                'some expected fields were deleted',
                )
            assert any(msg in caplog.text for msg in expected_logs)
        else:
            expected_logs = (
                'could not identify directory to remove',
                'No entries to remove',
                )
            assert any(msg in caplog.text for msg in expected_logs)


class TestBookkeeperOthers:
    """Collections of various tests for bits not covered by other tests."""

    # Note about the disable: It is more convenient to modify the key
    # for the glob-all-logs pattern later, as we can then simply edit
    # if we ever change again the naming of log files.
    # pylint: disable-next=dict-init-mutate
    _logs = {
        'tleedm-231110-103910.log': {'r_super': '0.1582',
                                     'run_info': '0 3 31 12'},
        'viperleed-calc-231110-103910.log': {
            'r_super': '0.1582',
            'r_ref': '0.1582 (0.1354 / 0.1827)',
            'run_info': '0 1 11 2 3 31 12',
            },
        }
    _logs['*.log'] = _logs['viperleed-calc-231110-103910.log']

    @parametrize('pattern,expect', _logs.items(), ids=_logs)
    def test_infer_from_log(self, pattern, expect, data_path, tmp_path):
        """Check correct detection of information from a log file."""
        for logfile in (data_path/'bookkeeper').glob(pattern):
            shutil.copy2(logfile, tmp_path)
        bookkeeper = Bookkeeper(cwd=tmp_path)
        bookkeeper.update_from_cwd(silent=True)
        # pylint: disable-next=protected-access   # OK in tests
        log_info = bookkeeper._infer_run_info_from_log()
        assert bookkeeper.all_cwd_logs
        assert log_info == expect

    def test_no_state_files_in_out(self):
        """Check correct behavior when there is no file in OUT to be used."""
        bookkeeper = Bookkeeper()
        # pylint: disable-next=protected-access   # OK in tests
        bookkeeper._update_state_files_from_out()  # There's no out
        assert not any(Path(f).exists() for f in MOCK_STATE_FILES)

    def test_remove_tensors_deltas(self, tmp_path):
        """Check removal of tensor and delta files."""
        bookkeeper = Bookkeeper(cwd=tmp_path)
        removed_files = (
            tmp_path/'Tensors/Tensors_003.zip',
            tmp_path/'Deltas/Deltas_003.zip',
            )
        surviving_files= (
            tmp_path/'Tensors/Tensors_002.zip',
            tmp_path/'Tensors/Tensors_001.zip',
            tmp_path/'Deltas/Deltas_001.zip',
            )
        folders = (tmp_path/'Tensors',
                   tmp_path/'Deltas')
        for folder in folders:
            folder.mkdir(parents=True)
        for file in (*removed_files, *surviving_files):
            file.touch()
        bookkeeper.update_from_cwd()
        # pylint: disable-next=protected-access   # OK in tests
        bookkeeper._remove_tensors_and_deltas()
        assert not any(f.exists() for f in removed_files)
        assert all(f.exists() for f in surviving_files)
        assert all(f.exists() for f in folders)

    @fixture(name='funky_files')
    def fixture_funky_files(self, tmp_path):
        """Prepare a bunch of files/folders that will not be considered."""
        # Stuff that should not be copied over:
        not_collected_log = tmp_path/'not_a_log.log'
        not_collected_log.mkdir()

        (tmp_path/DEFAULT_OUT).mkdir()  # Otherwise 'nothing to do'

        # Stuff that should not be considered for the state
        workhistory = tmp_path/DEFAULT_WORK_HISTORY
        not_collected_dirs = (
            workhistory/'some_folder',
            workhistory/f't001.r002_{PREVIOUS_LABEL}',
            )
        for directory in not_collected_dirs:
            directory.mkdir(parents=True)

        # Some stuff in history that should not be considered
        history = tmp_path/DEFAULT_HISTORY
        history.mkdir()
        tensor_num_unused = 999
        invalid_history_stuff = (
            history/f't{tensor_num_unused}.r999_some_file',
            history/'some_stray_directory',
            )
        for path in invalid_history_stuff:
            # pylint: disable-next=magic-value-comparison  # 'file'
            make = getattr(path, 'touch' if 'file' in path.name else 'mkdir')
            make()
        return (tensor_num_unused, not_collected_log,
                not_collected_dirs, invalid_history_stuff)

    def test_funky_files(self, funky_files, tmp_path):
        """Check that funny files and directories are not considered."""
        bookkeeper = Bookkeeper(cwd=tmp_path)
        (tensor_num_unused,
         not_collected_log,
         not_collected_dirs,
         invalid_history_stuff) = funky_files

        bookkeeper.update_from_cwd(silent=True)
        history_dir = bookkeeper.history_dir
        history_info = bookkeeper.history_info
        assert tensor_num_unused not in bookkeeper.max_job_for_tensor
        assert not_collected_log not in bookkeeper.all_cwd_logs

        # About the disables: W0212 (protected-access) is OK in a test;
        # E1135 (unsupported-membership-test) is a pylint problem, as
        # it cannot infer that by the time we reach this one we have
        # a tuple in _paths['to_be_archived'].
        # pylint: disable-next=E1135,W0212
        assert not any(f in bookkeeper._paths['to_be_archived']
                       for f in not_collected_dirs)

        bookkeeper.run('archive')
        assert not_collected_log.exists()
        assert not (history_dir/not_collected_log.name).exists()
        for file in not_collected_dirs:
            assert (not file.exists() if PREVIOUS_LABEL in file.name
                    else file.is_dir())
        assert all(f.exists() for f in invalid_history_stuff)
        assert history_dir.is_dir()
        assert history_info.path.read_bytes()

        # Now discard should remove workhistory and all the stuff
        bookkeeper.run('discard_full')
        assert not_collected_log.exists()
        assert not (tmp_path/DEFAULT_WORK_HISTORY).exists()
        assert not history_dir.exists()
        assert not history_info.path.read_bytes()


class TestBookkeeperRaises:
    """"Collection of tests for various bookkeeper complaints."""

    def test_cant_make_history(self):
        """Check complaints when we fail to make the history directory."""
        bookkeeper = Bookkeeper()
        with raises_oserror('pathlib.Path.mkdir'):
            # pylint: disable-next=protected-access    # OK in tests
            bookkeeper._make_history_and_prepare_logger()

        with raises_test_exception('pathlib.Path.mkdir'):
            # pylint: disable-next=protected-access    # OK in tests
            bookkeeper._make_history_and_prepare_logger()

    def test_discard_full_cant_remove_folder(self, after_archive, caplog):
        """Check complaints when it's not possible to remove folders."""
        bookkeeper, *_ = after_archive
        # Purge notes otherwise we can't remove anything.
        info = bookkeeper.history_info.raw_contents
        info, _ = info.rsplit('Notes:', maxsplit=1)
        if info:
            info += 'Notes:\n'
        info = bookkeeper.history_info.path.write_text(info)
        bookkeeper.history_info.read()

        with raises_test_exception('shutil.rmtree'):
            # pylint: disable-next=protected-access    # OK in tests
            bookkeeper._run_discard_full_mode()
        with make_raise_oserror('shutil.rmtree'):
            # pylint: disable-next=protected-access    # OK in tests
            exit_code = bookkeeper._run_discard_full_mode()
            # pylint: disable-next=magic-value-comparison
            assert 'Failed to delete' in caplog.text
            assert exit_code is BookkeeperExitCode.FAIL

    def test_invalid_mode(self):
        """Check complaints when an invalid mode is used."""
        bookkeeper = Bookkeeper()
        with pytest.raises(ValueError):
            bookkeeper.run('invalid')

    def test_no_runner_implemented(self, monkeypatch):
        """Check complaints if a runner is not available for a mode."""
        class _MockMode(Enum):
            INVALID = 'invalid'
        monkeypatch.setattr(
            'viperleed.calc.bookkeeper.bookkeeper.BookkeeperMode',
            _MockMode
            )
        bookkeeper = Bookkeeper()
        with pytest.raises(NotImplementedError):
            bookkeeper.run('invalid')

    raises = object()
    logs = object()
    skips = object()
    _os_error = {
        '_copy_out_and_supp': ('shutil.copytree', logs),
        '_discard_workhistory_previous': ('shutil.rmtree', logs),
        '_make_and_copy_to_history': ('pathlib.Path.mkdir', raises),
        '_move_and_cleanup_workhistory(True)': ('shutil.rmtree', logs),
        '_read_and_clear_notes_file-read': ('pathlib.Path.read_text', logs),
        '_read_and_clear_notes_file-write': ('pathlib.Path.write_text', logs),
        '_read_most_recent_log': ('pathlib.Path.open', skips),
        '_remove_log_files': ('pathlib.Path.unlink', logs),
        '_remove_out_and_supp': ('shutil.rmtree', logs),
        '_replace_state_files_from_ori': ('pathlib.Path.replace', raises),
        }

    @parametrize(method_name=_os_error)
    def test_oserror(self, method_name, tmp_path, monkeypatch, caplog):
        """Check expected outcome of calling a method on a bookkeeper."""

        workhistory = tmp_path/DEFAULT_WORK_HISTORY
        workhistory.mkdir()
        (workhistory/f't001.r002_some_{PREVIOUS_LABEL}_folder').mkdir()
        (tmp_path/DEFAULT_OUT).mkdir()
        (tmp_path/'notes').touch()
        (tmp_path/f'{LOG_PREFIX}-20xxxx-xxxxxx.log').touch()
        (tmp_path/'PARAMETERS_ori').touch()

        bookkeeper = Bookkeeper(cwd=tmp_path)
        bookkeeper.update_from_cwd(silent=True)

        to_patch, action = self._os_error[method_name]
        method_name, *_ = method_name.split('-')
        if '(' in method_name:  # pylint: disable=magic-value-comparison
            method_name, args = method_name.split('(')
            args = args.replace(')', '') + ','
            args = ast.literal_eval(args)
        else:
            args = tuple()
        method = getattr(bookkeeper, method_name)
        with raises_test_exception(to_patch):
            method(*args)

        if action is self.raises:
            with raises_oserror(to_patch):
                method(*args)
        else:
            with make_raise_oserror(to_patch), not_raises_oserror():
                method(*args)
        # Check logging results, unless we silently skip
        if action is not self.skips:
            assert any(r for r in caplog.records
                       if r.levelno >= logging.WARNING)

    _attr_needs_update = (
        'all_cwd_logs',
        'archiving_required',
        'cwd_logs',
        'files_need_archiving',
        'history_dir',
        'history_dir_base_name',
        'history_info',
        'history_with_same_base_name_exists',
        'max_job_for_tensor',
        'tensor_number',
        'timestamp',
        )
    _method_needs_update = (  # Only those without args
        '_collect_files_to_archive',
        '_copy_input_files_from_original_inputs_or_cwd',
        '_copy_log_files',
        '_copy_out_and_supp',
        '_deal_with_workhistory_and_history_info',
        '_discard_workhistory_previous',
        '_find_base_name_for_history_subfolder',
        '_find_new_history_directory_name',
        '_get_conflicting_history_subfolders',
        '_infer_run_info_from_log',
        '_make_and_copy_to_history',
        '_move_workhistory_folders',
        '_remove_log_files',
        '_remove_tensors_and_deltas',
        '_run_archive_mode',
        '_run_clear_mode',
        '_run_discard_common',
        '_run_discard_full_mode',
        '_run_discard_mode',
        )

    @parametrize(attr=_attr_needs_update)
    def test_too_early_attribute_access(self, attr):
        """check that accessing attributes before update_from_cwd fails."""
        bookkeeper = Bookkeeper()
        with pytest.raises(AttributeError) as exc:
            getattr(bookkeeper, attr)
        assert _UPDATE_METHOD in str(exc)

    @parametrize(method_name=_attr_needs_update)
    def test_too_early_method_call(self, method_name):
        """check that accessing attributes before update_from_cwd fails."""
        bookkeeper = Bookkeeper()
        with pytest.raises(AttributeError) as exc:
            # Some methods raise already at getattr, some at call
            method = getattr(bookkeeper, method_name)
            method()
        assert _UPDATE_METHOD in str(exc)


class TestBookkeeperComplaints:
    """Tests for situations that do not raise but issue log warnings/errors."""

    def test_cwd_file_newer(self, caplog, tmp_path):
        """Check warnings when a cwd file is newer than an original_input."""
        bookkeeper = Bookkeeper(cwd=tmp_path)
        ori_inputs = tmp_path/DEFAULT_SUPP/ORIGINAL_INPUTS_DIR_NAME
        ori_inputs.mkdir(parents=True)
        (ori_inputs/'POSCAR').touch()
        time.sleep(0.05)
        (tmp_path/'POSCAR').touch()

        bookkeeper.update_from_cwd()
        # pylint: disable-next=protected-access  # OK in tests
        bookkeeper._copy_input_files_from_original_inputs_or_cwd()
        # pylint: disable-next=magic-value-comparison
        assert 'is newer' in caplog.text

    def test_copy_from_root(self, caplog, tmp_path):
        """Check warnings when no original_input is present."""
        bookkeeper = Bookkeeper(cwd=tmp_path)
        (tmp_path/'POSCAR').touch()
        (tmp_path/DEFAULT_OUT).mkdir()  # Otherwise 'nothing to do'
        bookkeeper.update_from_cwd()
        target = bookkeeper.history_dir
        bookkeeper.run('archive')
        assert _FROM_ROOT in caplog.text
        assert (tmp_path/'POSCAR').is_file()
        assert (target/'POSCAR_from_root').is_file()
