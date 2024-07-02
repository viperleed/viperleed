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
    inputs_path = bookkeeper_mock_dir_after_run / 'SUPP' / ORIGINAL_INPUTS_DIR_NAME
    history = tmp_path / DEFAULT_HISTORY
    history.mkdir(parents=True, exist_ok=True)
    store_input_files_to_history(inputs_path, history)
    for file in MOCK_STATE_FILES:
        hist_file = history / file
        hist_file_content = hist_file.read_text()
        assert MOCK_ORIG_CONTENT in hist_file_content


class TestBookkeeperArchive:
    def test_bookkeeper_archive_after_run(self,
                                          after_run,
                                          caplog):
        """Check correct storage of history files in ARCHIVE mode."""
        bookkeeper, mock_dir, history_path, history_path_run = after_run
        # bookkeeper should think that it needs archiving
        assert bookkeeper.archiving_required is True
        bookkeeper.run(mode=BookkeeperMode.ARCHIVE)
        assert history_path.exists()
        assert history_path_run.is_dir()
        assert (mock_dir / 'history.info').exists()
        # Out stored in history
        for file in MOCK_STATE_FILES:
            hist_file = history_path_run / 'OUT' / f'{file}_OUT'
            hist_content = hist_file.read_text()
            assert MOCK_OUT_CONTENT in hist_content
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
        # Check that there are no errors or warnings in log
        assert not any(rec.levelno >= logging.WARNING for rec in caplog.records)

    def test_bookkeeper_archive_new(self,
                                    before_run,
                                    caplog):
        bookkeeper, mock_dir = before_run
        assert bookkeeper.archiving_required is False
        bookkeeper.run(mode=BookkeeperMode.ARCHIVE)
        # Bookkeeper should not do anything (except for logging)
        assert (mock_dir / 'history').exists()
        history_contents = (f for f in (mock_dir / 'history').iterdir()
                            if f.name != 'bookkeeper.log')
        assert not any(history_contents)
        # Originals untouched
        for file in MOCK_STATE_FILES:
            assert (mock_dir / file).is_file()
            input_content = (mock_dir / file).read_text()
            assert MOCK_INPUT_CONTENT in input_content
        # Check that there are no errors or warnings in log
        assert not any(rec.levelno >= logging.WARNING for rec in caplog.records)

    def test_bookkeeper_archive_again(self,
                                    after_archive,
                                    caplog):
        """Bookkeeper ARCHIVE after ARCHIVE should not do anything."""
        bookkeeper, mock_dir, *_ = after_archive
        assert bookkeeper.archiving_required is False
        # write stuff to files to check they are not overwritten
        for file in MOCK_STATE_FILES:
            (mock_dir / file).write_text('something else')
        bookkeeper.run(mode=BookkeeperMode.ARCHIVE)
        for file in MOCK_STATE_FILES:
            assert (mock_dir / file).read_text() == 'something else'
        # Check that there are no errors or warnings in log
        assert not any(rec.levelno >= logging.WARNING for rec in caplog.records)


class TestBookkeeperClear:
    def test_bookkeeper_clear_new(self, before_run, caplog):
        """Check correct overwriting of input files in CLEAR mode."""
        bookkeeper, mock_dir = before_run
        # bookkeeper should not think that it needs archiving
        assert bookkeeper.archiving_required is False
        bookkeeper.run(mode=BookkeeperMode.CLEAR)
        # Bookkeeper should not do anything (except for logging)
        history_contents = (f for f in (mock_dir / 'history').iterdir()
                            if f.name != 'bookkeeper.log')
        assert not any(history_contents)
        # Originals untouched
        for file in MOCK_STATE_FILES:
            assert (mock_dir / file).is_file()
            input_content = (mock_dir / file).read_text()
            assert MOCK_INPUT_CONTENT in input_content
        # Check that there are no errors or warnings in log
        assert not any(rec.levelno >= logging.WARNING for rec in caplog.records)

    def test_bookkeeper_clear_after_run(self,
                                          after_run,
                                          caplog):
        """Check correct behaviour of calling CLEAR after a run without ARCHiVE.
        E.g. if the run crashed."""
        bookkeeper, mock_dir, history_path, history_path_run = after_run
        # bookkeeper should think that it needs archiving
        assert bookkeeper.archiving_required is True
        bookkeeper.run(mode=BookkeeperMode.CLEAR)
        assert history_path.exists()
        assert history_path_run.is_dir()
        assert (mock_dir / 'history.info').exists()
        # Out stored in history
        for file in MOCK_STATE_FILES:
            hist_file = history_path_run / 'OUT' / f'{file}_OUT'
            hist_content = hist_file.read_text()
            assert MOCK_OUT_CONTENT in hist_content
        # _ori files, SUPP, OUT and logs should be removed
        for file in MOCK_STATE_FILES:
            ori_file = mock_dir / f'{file}_ori'
            assert not ori_file.is_file()
        assert not (mock_dir / 'SUPP').exists()
        assert not (mock_dir / 'OUT').exists()
        assert len(tuple(mock_dir.glob('*.log'))) == 0
        # Original SHOULD NOT be replaced by output
        # (ARCHIVE only does not run if the run crashed, in which case we
        # don't want to overwrite)
        for file in MOCK_STATE_FILES:
            out_content = (mock_dir / file).read_text()
            assert MOCK_INPUT_CONTENT in out_content
        # Check that there are no errors or warnings in log
        assert not any(rec.levelno >= logging.WARNING for rec in caplog.records)

    def test_bookkeeper_clear_after_archive(self,
                                            after_archive,
                                            caplog):
        """Check correct behaviour of calling CLEAR after a run and ARCHIVE.
        E.g. if the run crashed."""
        bookkeeper, mock_dir, history_path, history_path_run = after_archive
        # bookkeeper should not think that it needs archiving
        assert bookkeeper.archiving_required is False
        bookkeeper.run(mode=BookkeeperMode.CLEAR)
        assert history_path.exists()
        assert history_path_run.is_dir()
        assert (mock_dir / 'history.info').exists()
        # Out stored in history
        for file in MOCK_STATE_FILES:
            hist_file = history_path_run / 'OUT' / f'{file}_OUT'
            hist_content = hist_file.read_text()
            assert MOCK_OUT_CONTENT in hist_content
        # _ori files, SUPP, OUT and logs should be removed
        for file in MOCK_STATE_FILES:
            ori_file = mock_dir / f'{file}_ori'
            assert not ori_file.is_file()
        assert not (mock_dir / 'SUPP').exists()
        assert not (mock_dir / 'OUT').exists()
        assert len(tuple(mock_dir.glob('*.log'))) == 0
        # Original be replaced by output
        for file in MOCK_STATE_FILES:
            out_content = (mock_dir / file).read_text()
            assert MOCK_OUT_CONTENT in out_content
        # Check that there are no errors or warnings in log
        assert not any(rec.levelno >= logging.WARNING for rec in caplog.records)


class TestBookkeeperDiscard:
    def test_bookkeeper_discard_new(self, before_run, caplog):
        """Check correct overwriting of input files in CLEAR mode."""
        bookkeeper, mock_dir = before_run
        # bookkeeper should not think that it needs archiving
        assert bookkeeper.archiving_required is False
        bookkeeper.run(mode=BookkeeperMode.DISCARD)
        # Bookkeeper should not do anything (except for logging)
        history_contents = (f for f in (mock_dir / 'history').iterdir()
                            if f.name != 'bookkeeper.log')
        assert not any(history_contents)
        # Originals untouched
        for file in MOCK_STATE_FILES:
            assert (mock_dir / file).is_file()
            input_content = (mock_dir / file).read_text()
            assert MOCK_INPUT_CONTENT in input_content
        # Check that there are no errors or warnings in log
        assert not any(
            rec.levelno >= logging.WARNING
            and not "Failed to mark last entry as discarded" in str(rec)
            for rec in caplog.records)

    def test_bookkeeper_discard_after_run(self,
                                          after_run,
                                          caplog):
        """Calling DISCARD after a run without ARCHiVE. E.g. if the run crashed."""
        bookkeeper, mock_dir, history_path, history_path_run = after_run
        # bookkeeper should think that it needs archiving
        assert bookkeeper.archiving_required is True
        bookkeeper.run(mode=BookkeeperMode.DISCARD)
        assert history_path.exists()
        assert history_path_run.is_dir()
        assert (mock_dir / 'history.info').exists()
        # Out stored in history
        for file in MOCK_STATE_FILES:
            hist_file = history_path_run / 'OUT' / f'{file}_OUT'
            hist_content = hist_file.read_text()
            assert MOCK_OUT_CONTENT in hist_content
        # _ori files, SUPP, OUT and logs should be removed
        for file in MOCK_STATE_FILES:
            ori_file = mock_dir / f'{file}_ori'
            assert not ori_file.is_file()
        assert not (mock_dir / 'SUPP').exists()
        assert not (mock_dir / 'OUT').exists()
        assert len(tuple(mock_dir.glob('*.log'))) == 0
        # Original SHOULD NOT be replaced by output
        # (ARCHIVE only does not run if the run crashed, in which case we
        # don't want to overwrite)
        for file in MOCK_STATE_FILES:
            out_content = (mock_dir / file).read_text()
            assert MOCK_INPUT_CONTENT in out_content
        # A 'DISCARDED' note should be in history.info
        assert bookkeeper.history_info.last_entry_was_discarded
        assert _DISCARDED in bookkeeper.history_info.path.read_text()
        # Check that there are no errors or warnings in log
        assert not any(
            rec.levelno >= logging.WARNING
            and not "discarded" in str(rec)
            for rec in caplog.records)

    def test_bookkeeper_discard_after_archive(self,
                                            after_archive,
                                            caplog):
        """Bookkeeper DISCARD after a run and ARCHiVE.
        Should revert state to before the run."""
        bookkeeper, mock_dir, history_path, history_path_run = after_archive
        # bookkeeper should not think that it needs archiving
        assert bookkeeper.archiving_required is False
        bookkeeper.run(mode=BookkeeperMode.DISCARD)
        assert history_path.exists()
        assert history_path_run.is_dir()
        assert (mock_dir / 'history.info').exists()
        # Out stored in history
        for file in MOCK_STATE_FILES:
            hist_file = history_path_run / 'OUT' / f'{file}_OUT'
            hist_content = hist_file.read_text()
            assert MOCK_OUT_CONTENT in hist_content
        # _ori files, SUPP, OUT and logs should be removed
        for file in MOCK_STATE_FILES:
            ori_file = mock_dir / f'{file}_ori'
            assert not ori_file.is_file()
        assert not (mock_dir / 'SUPP').exists()
        assert not (mock_dir / 'OUT').exists()
        assert len(tuple(mock_dir.glob('*.log'))) == 0
        # Original be replaced by output
        for file in MOCK_STATE_FILES:
            out_content = (mock_dir / file).read_text()
            assert MOCK_INPUT_CONTENT in out_content
        # A 'DISCARDED' note should be in history.info
        with bookkeeper.history_info.path.open() as f:
            assert _DISCARDED in f.read()
        # Check that there are no errors or warnings in log
        assert not any(rec.levelno >= logging.WARNING for rec in caplog.records)


class TestBookkeeperDiscardFull:
    def test_bookkeeper_discard_full_new(self, before_run, caplog):
        """Check correct overwriting of input files in DISCARD_FULL mode.
        Should be the same as normal DISCARD in this case."""
        bookkeeper, mock_dir = before_run
        # bookkeeper should not think that it needs archiving
        assert bookkeeper.archiving_required is False
        bookkeeper.run(mode=BookkeeperMode.DISCARD_FULL)
        # Bookkeeper should not do anything (except for logging)
        history_contents = (f for f in (mock_dir / 'history').iterdir()
                            if f.name != 'bookkeeper.log')
        assert not any(history_contents)
        # Originals untouched
        for file in MOCK_STATE_FILES:
            assert (mock_dir / file).is_file()
            input_content = (mock_dir / file).read_text()
            assert MOCK_INPUT_CONTENT in input_content
        # Check that there are no errors or warnings in log
        assert not any(rec.levelno >= logging.WARNING and
                       not "remove" in str(rec) # catch two expected warnings
                       for rec in caplog.records)

    def test_bookkeeper_discard_full_after_run(self,
                                               after_run,
                                               caplog):
        """Calling DISCARD_FULL after a run without ARCHiVE. E.g. if the run
        crashed."""
        bookkeeper, mock_dir, history_path, history_path_run = after_run
        notes_in_history_info = (
            NOTES_TEST_CONTENT in (mock_dir / HISTORY_INFO_NAME).read_text()
            )
        bookkeeper.run(mode=BookkeeperMode.DISCARD_FULL)
        assert not history_path_run.is_dir()
        # Since teh run was not archived, the history should be empty
        if notes_in_history_info:
            assert any("last entry in history.info has user notes" in str(rec)
                       for rec in caplog.records)
        else:
            assert any("could not identify directory to remove" in str(rec)
                    for rec in caplog.records)

    def test_bookkeeper_discard_full_after_archive(self,
                                                   after_archive,
                                                   caplog):
        """Bookkeeper DISCARD_FULL after a run and ARCHIVE.
        Should revert state to before the run and purge the history.
        """
        bookkeeper, mock_dir, history_path, history_path_run = after_archive
        bookkeeper.run(mode=BookkeeperMode.DISCARD_FULL)
        if bookkeeper.history_info.last_entry_has_notes:
            # should prevent the removal of the history
            assert history_path_run.is_dir()
            return
        assert not history_path_run.is_dir()
        # _ori files, SUPP, OUT and logs should be removed
        for file in MOCK_STATE_FILES:
            ori_file = mock_dir / f'{file}_ori'
            assert not ori_file.is_file()
        assert not (mock_dir / 'SUPP').exists()
        assert not (mock_dir / 'OUT').exists()
        assert len(tuple(mock_dir.glob('*.log'))) == 0
        for file in MOCK_STATE_FILES:
            out_content = (mock_dir / file).read_text()
            assert MOCK_INPUT_CONTENT in out_content
        # Check that there are no errors or warnings in log
        assert not any(rec.levelno >= logging.WARNING for rec in caplog.records)
