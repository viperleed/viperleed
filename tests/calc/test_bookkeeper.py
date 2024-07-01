"""Tests for module viperleed.calc.bookkeeper."""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-08-02'
__license__ = 'GPLv3+'

import shutil
import logging

import pytest
from pytest_cases import fixture, parametrize

from viperleed.calc import DEFAULT_HISTORY
from viperleed.calc import ORIGINAL_INPUTS_DIR_NAME
from viperleed.calc.bookkeeper import _CALC_LOG_PREFIXES, HISTORY_INFO_NAME
from viperleed.calc.bookkeeper import HISTORY_INFO_SEPARATOR
from viperleed.calc.bookkeeper import Bookkeeper
from viperleed.calc.bookkeeper import BookkeeperMode
from viperleed.calc.bookkeeper import HistoryInfoFile
from viperleed.calc.bookkeeper import store_input_files_to_history

from ..helpers import execute_in_dir


MOCK_TIMESTAMP = '010203-040506'
MOCK_LOG_FILES = [f'{pre}-{MOCK_TIMESTAMP}.log' for pre in _CALC_LOG_PREFIXES]
MOCK_HISTORY_INFO_FILES = {
    'no history.info': None,
    'empty history.info': "",
    'entry with job name': (
        '# TENSORS   \n# JOB ID    \n# JOB NAME  test_jobname\n'
        '# TIME      03.02.01 04:03:06\n# FOLDER    t003.r001_010203-040506\n'
        'Notes: This is a test note.\n'),
    'entry without job name': (
        '# TENSORS   \n# JOB ID    \n'
        '# TIME      03.02.01 04:03:06\n# FOLDER    t003.r001_010203-040506\n'
        'Notes: This is a test note.\n'),
    'entry without note': (
        '# TENSORS   \n# JOB ID    \n'
        '# TIME      03.02.01 04:03:06\n# FOLDER    t003.r001_010203-040506\n'
        'Notes:\n'),
    'with RUN': (
        '# TENSORS   \n# JOB ID    \n'
        '# RUN       1 2 3\n# TIME      03.02.01 04:03:06\n# FOLDER    t003.r001_010203-040506\n'
        'Notes:\n'),
    'with R REF': (
        '# TENSORS   \n# JOB ID    \n'
        '# TIME      03.02.01 04:03:06\n# R REF     0.1234\n# FOLDER    t003.r001_010203-040506\n'
        'Notes:\n'),
    'with R SUPER': (
        '# TENSORS   \n# JOB ID    \n'
        '# TIME      03.02.01 04:03:06\n# R SUPER   0.1234\n# FOLDER    t003.r001_010203-040506\n'
        'Notes:\n'),
    'entry discarded': (
        '# TENSORS   \n# JOB ID    \n# JOB NAME  test_jobname\n'
        '# TIME      03.02.01 04:03:06\n# FOLDER    t003.r001_010203-040506\n'
        'Notes: This is a test note.\n'
        'DISCARDED\n'),
    'two entries without job name': (
        '# TENSORS   \n# JOB ID    \n'
        '# TIME      03.02.01 04:03:06\n# FOLDER    t003.r001_010203-040506\n'
        'Notes: This is a test note.\n'
        f'{HISTORY_INFO_SEPARATOR}\n'
        '# TENSORS   \n# JOB ID    \n'
        '# TIME      03.02.01 04:05:06\n# FOLDER    t003.r001_010203-040506\n'
        'Notes: This is a test note.\n'),
}
MOCK_JOB_NAMES = (None, 'test_jobname')
NOTES_TEST_CONTENT = 'This is a test note.'
ALT_HISTORY_NAME = 'history_alt_name'
MOCK_STATE_FILES = ('POSCAR', 'VIBROCC', 'PARAMETERS')
MOCK_INPUT_CONTENT = 'This is a test input file.'
MOCK_ORIG_CONTENT = 'This is a test original input file.'
MOCK_OUT_CONTENT = 'This is a test output file.'


@fixture(name='bookkeeper_mock_dir_after_run')
@parametrize(log_file_name=MOCK_LOG_FILES, ids=MOCK_LOG_FILES)
@parametrize(history_info_file=MOCK_HISTORY_INFO_FILES.values(),
             ids=MOCK_HISTORY_INFO_FILES.keys())
def fixture_bookkeeper_mock_dir_after_run(tmp_path, log_file_name,
                                          history_info_file):
    """Yield a temporary directory for testing the bookkeeper."""
    out_path = tmp_path / 'OUT'
    supp_path = tmp_path / 'SUPP'
    tensors_path = tmp_path / 'Tensors'
    deltas_path = tmp_path / 'Deltas'
    directories = (out_path, supp_path, tensors_path, deltas_path)
    for directory in directories:
        directory.mkdir(parents=True, exist_ok=True)
    # create mock log files
    (tmp_path / log_file_name).touch()
    # create mock notes file
    notes_file = tmp_path / 'notes.txt'
    notes_file.write_text(NOTES_TEST_CONTENT)
    # mock history.info file
    if history_info_file is not None:
        hist_info_path = tmp_path / HISTORY_INFO_NAME
        with open(hist_info_path, 'w') as f:
            f.write(history_info_file)
    # create mock Tensor and Delta files
    (tensors_path / 'Tensors_003.zip').touch()
    (deltas_path / 'Deltas_003.zip').touch()

    # create non-empty mock alt_history folder
    alt_history = tmp_path / ALT_HISTORY_NAME
    for run_name in ('t001.r001_20xxxx-xxxxxx', 't002.r002_20xxxx-xxxxxx'):
        (alt_history / run_name).mkdir(parents=True, exist_ok=True)

    # create mock input files
    for file in MOCK_STATE_FILES:
        (tmp_path / file).write_text(MOCK_INPUT_CONTENT)

    # create mock original_inputs folder
    original_inputs_path = supp_path / ORIGINAL_INPUTS_DIR_NAME
    original_inputs_path.mkdir(parents=True, exist_ok=True)
    for file in MOCK_STATE_FILES:
        (original_inputs_path / file).write_text(MOCK_ORIG_CONTENT)

    # create mock OUT folder
    for file in MOCK_STATE_FILES:
        out_file = out_path / f'{file}_OUT'
        out_file.write_text(MOCK_OUT_CONTENT)

    with execute_in_dir(tmp_path):
        yield tmp_path
    # It would be nice to clean up, but the following line causes
    # a PermissionError. Likely because of logging keeping a hold
    # of the bookkeeper.log file.
    # shutil.rmtree(tmp_path)


@fixture(name='after_run')
@parametrize(job_name=MOCK_JOB_NAMES)
@parametrize(history_name=(DEFAULT_HISTORY, ALT_HISTORY_NAME))
def fixture_after_run(bookkeeper_mock_dir_after_run, job_name, history_name):
    """Return the path to the temporary directory after the run."""
    bookkeeper = Bookkeeper(cwd=bookkeeper_mock_dir_after_run,
                            job_name=job_name,
                            history_name=history_name)
    history_path = bookkeeper_mock_dir_after_run / history_name
    dir_name = f't003.r001_{MOCK_TIMESTAMP}'
    if job_name is not None:
        dir_name += f'_{job_name}'
    history_run_path = history_path / dir_name
    return bookkeeper, bookkeeper_mock_dir_after_run, history_path, history_run_path


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


@fixture(name='history_info_file')
def fixture_history_info(after_run):
    bookkeeper, *_ = after_run
    history_info = bookkeeper.history_info
    return history_info, bookkeeper.cwd / HISTORY_INFO_NAME

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
        assert "DISCARDED" in bookkeeper.history_info.path.read_text()
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
            assert 'DISCARDED' in f.read()
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
            'This is a test note.' in (mock_dir / HISTORY_INFO_NAME).read_text())
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


class TestHistoryInfoFile:

    def test_history_info_read_contents(self,history_info_file):
        """Check that the history.info file is read correctly."""
        history_info, actual_file = history_info_file
        assert actual_file.exists()
        assert history_info.path == actual_file

    def test_history_info_entry_parsing(self,history_info_file):
        history_info, actual_file = history_info_file
        assert bool(history_info.last_entry) == bool(actual_file.read_text())

    def test_history_info_has_notes(self,history_info_file):
        """Check that the history.info file is read correctly."""
        history_info, actual_file = history_info_file
        if 'This is a test note.' in actual_file.read_text():
            assert history_info.last_entry.notes == 'This is a test note.'

    def test_history_info_discarded(self,history_info_file):
        """Check that the history.info file is read correctly."""
        history_info, actual_file = history_info_file
        discarded_in_entry = 'DISCARDED' in actual_file.read_text()
        assert history_info.last_entry_was_discarded == discarded_in_entry

    def test_history_info_regenerate_from_entries(self,history_info_file):
        """Check that the history.info file can be regenerated after parsing."""
        history_info, actual_file = history_info_file
        actual_text = actual_file.read_text()
        assert actual_text == history_info.raw_contents
        last_entry = history_info.last_entry
        if last_entry is None:
            return
        history_info.remove_last_entry()
        history_info.append_entry(last_entry)
        assert actual_text.strip() == history_info.raw_contents.strip()

    def test_history_info_discard_last_entry(self, history_info_file):
        """Check we can discard the last history.info entry."""
        history_info, actual_file = history_info_file
        last_entry = history_info.last_entry
        if last_entry is None or last_entry.discarded:
            return
        history_info.discard_last_entry()
        assert history_info.last_entry_was_discarded

    def test_history_info_remove_last_entry(self, history_info_file):
        """Check we can remove the last history.info entry."""
        history_info, *_ = history_info_file
        # check number of entries before and run checks accordingly
        n_entries = history_info.path.read_text().count('# TENSORS')
        if n_entries == 0:
            assert history_info.last_entry is None
            with pytest.raises(ValueError):
                history_info.remove_last_entry()
        else:
            assert history_info.last_entry is not None
            history_info.remove_last_entry()
            assert history_info.path.read_text().count('# TENSORS') == n_entries - 1
