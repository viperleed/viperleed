"""Tests for module viperleed.calc.bookkeeper."""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-08-02'
__license__ = 'GPLv3+'

import logging

from pytest_cases import fixture

from viperleed.calc import DEFAULT_HISTORY
from viperleed.calc import ORIGINAL_INPUTS_DIR_NAME
from viperleed.calc.bookkeeper.bookkeeper import Bookkeeper
from viperleed.calc.bookkeeper.bookkeeper import store_input_files_to_history
from viperleed.calc.bookkeeper.constants import HISTORY_INFO_NAME
from viperleed.calc.bookkeeper.history import _DISCARDED
from viperleed.calc.bookkeeper.mode import BookkeeperMode
from viperleed.calc.sections.cleanup import DEFAULT_OUT
from viperleed.calc.sections.cleanup import DEFAULT_SUPP

from ...helpers import execute_in_dir
from .conftest import MOCK_INPUT_CONTENT
from .conftest import MOCK_ORIG_CONTENT
from .conftest import MOCK_OUT_CONTENT
from .conftest import MOCK_TIMESTAMP
from .conftest import MOCK_STATE_FILES
from .conftest import NOTES_TEST_CONTENT


@fixture(name='after_archive')
def fixture_bookkeeper_mock_dir_after_archive(after_run):
    """Yield a temporary directory for testing the bookkeeper."""
    bookkeeper, bookkeeper_mock_dir_after_run, history_path, history_run_path = after_run
    bookkeeper.run(mode=BookkeeperMode.ARCHIVE)
    return bookkeeper, bookkeeper_mock_dir_after_run, history_path, history_run_path


@fixture(name='before_run')
def fixture_bookkeeper_mock_dir_new(tmp_path):
    """Yield a temporary directory for testing the bookkeeper.
    This represents a new calculation, i.e., before any viperleed calc or
    bookkeeper run."""
    # create mock input files
    for file in MOCK_STATE_FILES:
        (tmp_path / file).write_text(MOCK_INPUT_CONTENT)

    bookkeeper = Bookkeeper(cwd=tmp_path)
    with execute_in_dir(tmp_path):
        yield bookkeeper, tmp_path
    # It would be nice to clean up, but the following line causes
    # a PermissionError. Likely because of logging keeping a hold
    # of the bookkeeper.log file
    # shutil.rmtree(tmp_path)


@fixture(name='history_path')
def fixture_history_path(bookkeeper_mock_dir_after_run):
    """Return the path to a history subfolder of `bookkeeper_mock_dir`."""
    return bookkeeper_mock_dir_after_run / DEFAULT_HISTORY


@fixture(name='history_path_run')
def fixture_history_path_run(history_path):
    """Return the path to a history run subfolder of `history_path`."""
    return history_path / f't000.r001_{MOCK_TIMESTAMP}'


def test_bookkeeper_mode_enum():
    """Check values of bookkeeper mode enum."""
    assert BookkeeperMode.ARCHIVE is BookkeeperMode('archive')
    assert BookkeeperMode.CLEAR is BookkeeperMode('clear')
    assert BookkeeperMode.DISCARD is BookkeeperMode('discard')
    assert BookkeeperMode.DISCARD_FULL is BookkeeperMode('discard_full')


def test_store_input_files_to_history(tmp_path, bookkeeper_mock_dir_after_run):
    """Check correct storage of original input files to history."""
    inputs_path = bookkeeper_mock_dir_after_run / DEFAULT_SUPP / ORIGINAL_INPUTS_DIR_NAME
    history = tmp_path / DEFAULT_HISTORY
    history.mkdir(parents=True, exist_ok=True)
    store_input_files_to_history(inputs_path, history)
    for file in MOCK_STATE_FILES:
        hist_file = history / file
        hist_file_content = hist_file.read_text()
        assert MOCK_ORIG_CONTENT in hist_file_content


def check_no_warnings(caplog, contain=None, exclude_msgs=()):
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


def check_history_exists(run):
    """Test that history_path and directory/history.info exist."""
    _, mock_dir, history_path, history_path_run = run
    assert history_path.exists()
    assert history_path_run.is_dir()
    assert (mock_dir / HISTORY_INFO_NAME).exists()


def check_out_files_in_history(run):
    """Check that the expected state files are stored in 'OUT'."""
    *_, history_path_run = run
    for file in MOCK_STATE_FILES:
        hist_file = history_path_run / DEFAULT_OUT / f'{file}_OUT'
        assert hist_file.is_file()
        hist_content = hist_file.read_text()
        assert MOCK_OUT_CONTENT in hist_content


def check_original_inputs_untouched(run):
    """Check the the original state files have not been moved."""
    _, mock_dir, *_ = run
    for file in MOCK_STATE_FILES:
        assert (mock_dir / file).is_file()
        input_content = (mock_dir / file).read_text()
        assert MOCK_INPUT_CONTENT in input_content


def check_main_directory_clean(run):
    """Check that no calc output is present in the main directory."""
    _, mock_dir, *_ = run
    for file in MOCK_STATE_FILES:
        ori_file = mock_dir / f'{file}_ori'
        assert not ori_file.is_file()
    assert not (mock_dir / DEFAULT_SUPP).exists()
    assert not (mock_dir / DEFAULT_OUT).exists()
    assert not any(mock_dir.glob('*.log'))


class TestBookkeeperArchive:
    def test_archive_after_run(self, after_run, caplog):
        """Check correct storage of history files in ARCHIVE mode."""
        bookkeeper, mock_dir, *_ = after_run
        # bookkeeper should think that it needs archiving
        assert bookkeeper.archiving_required
        bookkeeper.run(mode=BookkeeperMode.ARCHIVE)
        check_history_exists(after_run)
        check_out_files_in_history(after_run)
        # Original moved to _ori
        for file in MOCK_STATE_FILES:
            ori_file = mock_dir / f'{file}_ori'
            assert ori_file.is_file()
            input_content = (mock_dir / file).read_text()
            assert MOCK_OUT_CONTENT in input_content
        # Original replaced by output
        for file in MOCK_STATE_FILES:
            out_content = (mock_dir / file).read_text()
            assert MOCK_OUT_CONTENT in out_content
        check_no_warnings(caplog)

    def test_archive_new(self, before_run, caplog):
        bookkeeper, mock_dir = before_run
        assert bookkeeper.archiving_required is False
        bookkeeper.run(mode=BookkeeperMode.ARCHIVE)
        # Bookkeeper should not do anything (except for logging)
        assert (mock_dir / 'history').exists()
        history_contents = (f for f in (mock_dir / 'history').iterdir()
                            if f.name != 'bookkeeper.log')
        assert not any(history_contents)
        check_original_inputs_untouched(before_run)
        check_no_warnings(caplog)

    def test_archive_again(self, after_archive, caplog):
        """Bookkeeper ARCHIVE after ARCHIVE should not do anything."""
        bookkeeper, mock_dir, *_ = after_archive
        assert bookkeeper.archiving_required is False
        # write stuff to files to check they are not overwritten
        for file in MOCK_STATE_FILES:
            (mock_dir / file).write_text('something else')
        bookkeeper.run(mode=BookkeeperMode.ARCHIVE)
        for file in MOCK_STATE_FILES:
            assert (mock_dir / file).read_text() == 'something else'
        check_no_warnings(caplog)


class TestBookkeeperClear:
    def test_clear_new(self, before_run, caplog):
        """Check correct overwriting of input files in CLEAR mode."""
        bookkeeper, mock_dir = before_run
        # bookkeeper should not think that it needs archiving
        assert not bookkeeper.archiving_required
        bookkeeper.run(mode=BookkeeperMode.CLEAR)
        # Bookkeeper should not do anything (except for logging)
        history_contents = (f for f in (mock_dir / 'history').iterdir()
                            if f.name != 'bookkeeper.log')
        assert not any(history_contents)
        check_original_inputs_untouched(before_run)
        check_no_warnings(caplog)

    def test_clear_after_run(self, after_run, caplog):
        """Check correct behaviour of calling CLEAR after a run without ARCHiVE.
        E.g. if the run crashed."""
        bookkeeper, *_ = after_run
        # bookkeeper should think that it needs archiving
        assert bookkeeper.archiving_required
        bookkeeper.run(mode=BookkeeperMode.CLEAR)
        check_history_exists(after_run)
        check_out_files_in_history(after_run)
        check_main_directory_clean(after_run)
        # Original SHOULD NOT be replaced by output:
        # ARCHIVE does not run only if the run crashed,
        # in which case we don't want to overwrite
        check_original_inputs_untouched(after_run)
        check_no_warnings(caplog)

    def test_clear_after_archive(self, after_archive, caplog):
        """Check correct behaviour of calling CLEAR after a run and ARCHIVE.
        E.g. if the run crashed."""
        bookkeeper, mock_dir, *_ = after_archive
        # bookkeeper should not think that it needs archiving
        assert not bookkeeper.archiving_required
        bookkeeper.run(mode=BookkeeperMode.CLEAR)
        check_history_exists(after_archive)
        check_out_files_in_history(after_archive)
        check_main_directory_clean(after_archive)
        # Original be replaced by output
        for file in MOCK_STATE_FILES:
            out_content = (mock_dir / file).read_text()
            assert MOCK_OUT_CONTENT in out_content
        check_no_warnings(caplog)


class TestBookkeeperDiscard:
    def test_discard_new(self, before_run, caplog):
        """Check correct overwriting of input files in CLEAR mode."""
        bookkeeper, mock_dir = before_run
        # bookkeeper should not think that it needs archiving
        assert not bookkeeper.archiving_required
        bookkeeper.run(mode=BookkeeperMode.DISCARD)
        # Bookkeeper should not do anything (except for logging)
        history_contents = (f for f in (mock_dir / 'history').iterdir()
                            if f.name != 'bookkeeper.log')
        assert not any(history_contents)
        check_original_inputs_untouched(before_run)
        check_no_warnings(
            caplog,
            exclude_msgs=('Failed to mark last entry as discarded',),
            )

    def test_discard_after_run(self, after_run, caplog):
        """Calling DISCARD after a run without ARCHiVE. E.g. if the run crashed."""
        bookkeeper, *_ = after_run
        # bookkeeper should think that it needs archiving
        assert bookkeeper.archiving_required
        bookkeeper.run(mode=BookkeeperMode.DISCARD)
        check_history_exists(after_run)
        check_out_files_in_history(after_run)
        check_main_directory_clean(after_run)
        # Original SHOULD NOT be replaced by output:
        # ARCHIVE does not run only if the run crashed,
        # in which case we don't want to overwrite
        check_original_inputs_untouched(after_run)
        # A 'DISCARDED' note should be in history.info
        assert bookkeeper.history_info.last_entry_was_discarded
        assert _DISCARDED in bookkeeper.history_info.path.read_text()
        check_no_warnings(caplog, exclude_msgs=('discarded',))

    def test_discard_after_archive(self, after_archive, caplog):
        """Bookkeeper DISCARD after a run and ARCHiVE.
        Should revert state to before the run."""
        bookkeeper, mock_dir, *_ = after_archive
        # bookkeeper should not think that it needs archiving
        assert not bookkeeper.archiving_required
        bookkeeper.run(mode=BookkeeperMode.DISCARD)
        check_history_exists(after_archive)
        check_out_files_in_history(after_archive)
        check_main_directory_clean(after_archive)
        # Original be replaced by output                                        # TODO: this does something else!
        for file in MOCK_STATE_FILES:
            out_content = (mock_dir / file).read_text()
            assert MOCK_INPUT_CONTENT in out_content
        # A 'DISCARDED' note should be in history.info
        with bookkeeper.history_info.path.open() as f:
            assert _DISCARDED in f.read()
        # Some fields are knowingly faulty,
        # but we can still DISCARD them.
        faulty_entry_logs = (
            'Found entry with',
            'Could not understand',
            'Faulty entry is',
            )
        check_no_warnings(caplog, exclude_msgs=faulty_entry_logs)


class TestBookkeeperDiscardFull:
    def test_discard_full_new(self, before_run, caplog):
        """Check correct overwriting of input files in DISCARD_FULL mode.
        Should be the same as normal DISCARD in this case."""
        bookkeeper, mock_dir = before_run
        # bookkeeper should not think that it needs archiving
        assert not bookkeeper.archiving_required
        bookkeeper.run(mode=BookkeeperMode.DISCARD_FULL)
        # Bookkeeper should not do anything (except for logging)
        history_contents = (f for f in (mock_dir / 'history').iterdir()
                            if f.name != 'bookkeeper.log')
        assert not any(history_contents)
        check_original_inputs_untouched(before_run)
        check_no_warnings(
            caplog,
            exclude_msgs=('remove',)   # Catches two expected warnings
            )

    def test_discard_full_after_run(self, after_run, caplog):
        """Calling DISCARD_FULL after a run without ARCHiVE. E.g. if the run
        crashed."""
        bookkeeper, mock_dir, _, history_path_run = after_run
        notes_in_history_info = (
            NOTES_TEST_CONTENT in (mock_dir / HISTORY_INFO_NAME).read_text()
            )
        last_entry = bookkeeper.history_info.last_entry
        bookkeeper.run(mode=BookkeeperMode.DISCARD_FULL)
        # Since the run was not archived, the history should be empty
        assert not history_path_run.is_dir()
        if notes_in_history_info:
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

    def test_discard_full_after_archive(self, after_archive, caplog):
        """Bookkeeper DISCARD_FULL after a run and ARCHIVE.
        Should revert state to before the run and purge the history.
        """
        bookkeeper, *_, history_path_run = after_archive
        last_entry = bookkeeper.history_info.last_entry
        bookkeeper.run(mode=BookkeeperMode.DISCARD_FULL)
        if not last_entry.can_be_removed:
            # Should prevent the removal of the history
            assert history_path_run.is_dir()
            return
        assert not history_path_run.is_dir()
        check_main_directory_clean(after_archive)
        check_original_inputs_untouched(after_archive)
        check_no_warnings(caplog)
