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

from pytest_cases import fixture, parametrize

from viperleed.calc import DEFAULT_HISTORY
from viperleed.calc import DEFAULT_WORK
from viperleed.calc import ORIGINAL_INPUTS_DIR_NAME
from viperleed.calc.bookkeeper import BookkeeperMode
from viperleed.calc.bookkeeper import _CALC_LOG_PREFIXES
from viperleed.calc.bookkeeper import Bookkeeper
from viperleed.calc.bookkeeper import store_input_files_to_history

from ..helpers import execute_in_dir


MOCK_TIMESTAMP = '010203-040506'
MOCK_LOG_FILES = [f'{pre}-{MOCK_TIMESTAMP}.log' for pre in _CALC_LOG_PREFIXES]
MOCK_JOB_NAMES = (None, 'test_jobname')
NOTES_TEST_CONTENT = 'This is a test note.'
ALT_HISTORY_NAME = 'history_alt_name'
MOCK_STATE_FILES = ('POSCAR', 'VIBROCC', 'PARAMETERS')
MOCK_INPUT_CONTENT = 'This is a test input file.'
MOCK_ORIG_CONTENT = 'This is a test original input file.'
MOCK_OUT_CONTENT = 'This is a test output file.'


@fixture(name='bookkeeper_mock_dir_after_run')
@parametrize(log_file_name=MOCK_LOG_FILES)
def fixture_bookkeeper_mock_dir_after_run(tmp_path, log_file_name):
    """Yield a temporary directory for testing the bookkeeper."""
    work_path = tmp_path / DEFAULT_WORK
    out_path = tmp_path / 'OUT'
    supp_path = tmp_path / 'SUPP'
    tensors_path = tmp_path / 'Tensors'
    deltas_path = tmp_path / 'Deltas'
    directories = (work_path, out_path, supp_path, tensors_path, deltas_path)
    for directory in directories:
        directory.mkdir(parents=True, exist_ok=True)
    # create mock log files
    (tmp_path / log_file_name).touch()
    (work_path / log_file_name).touch()
    # create mock notes file
    notes_file = tmp_path / 'notes.txt'
    notes_file.write_text(NOTES_TEST_CONTENT)
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
    original_inputs_path = work_path / ORIGINAL_INPUTS_DIR_NAME
    original_inputs_path.mkdir(parents=True, exist_ok=True)
    for file in MOCK_STATE_FILES:
        (original_inputs_path / file).write_text(MOCK_ORIG_CONTENT)

    # create mock OUT folder
    for file in MOCK_STATE_FILES:
        out_file = out_path / f'{file}_OUT'
        out_file.write_text(MOCK_OUT_CONTENT)

    with execute_in_dir(tmp_path):
        yield tmp_path
    shutil.rmtree(tmp_path)


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
    shutil.rmtree(tmp_path)


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
    inputs_path = bookkeeper_mock_dir_after_run / DEFAULT_WORK / ORIGINAL_INPUTS_DIR_NAME
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
        # Bookkeeper should not do anything
        assert (mock_dir / 'history').exists()
        assert not tuple((mock_dir / 'history').iterdir())
        assert not (mock_dir / 'history.info').exists()
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
        # Bookkeeper should not do anything
        assert not tuple((mock_dir / 'history').iterdir())
        assert not (mock_dir / 'history.info').exists()
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
        # Bookkeeper should not do anything
        assert not tuple((mock_dir / 'history').iterdir())
        assert not (mock_dir / 'history.info').exists()
        # Originals untouched
        for file in MOCK_STATE_FILES:
            assert (mock_dir / file).is_file()
            input_content = (mock_dir / file).read_text()
            assert MOCK_INPUT_CONTENT in input_content
        # Check that there are no errors or warnings in log
        assert not any(rec.levelno >= logging.WARNING for rec in caplog.records)

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
        # Check that there are no errors or warnings in log
        assert not any(rec.levelno >= logging.WARNING for rec in caplog.records)

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
        # Bookkeeper should not do anything
        assert not tuple((mock_dir / 'history').iterdir())
        assert not (mock_dir / 'history.info').exists()
        # Originals untouched
        for file in MOCK_STATE_FILES:
            assert (mock_dir / file).is_file()
            input_content = (mock_dir / file).read_text()
            assert MOCK_INPUT_CONTENT in input_content
        # Check that there are no errors or warnings in log
        assert not any(rec.levelno >= logging.WARNING for rec in caplog.records)

    def test_bookkeeper_discard_full_after_run(self,
                                               after_run,
                                               caplog):
        """Calling DISCARD_FULL after a run without ARCHiVE. E.g. if the run
        crashed."""
        bookkeeper, mock_dir, history_path, history_path_run = after_run
        bookkeeper.run(mode=BookkeeperMode.DISCARD_FULL)
        assert not history_path_run.is_dir()
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

    def test_bookkeeper_discard_full_after_archive(self,
                                                   after_archive,
                                                   caplog):
        """Bookkeeper DISCARD_FULL after a run and ARCHIVE.
        Should revert state to before the run and purge the history.
        """
        bookkeeper, mock_dir, history_path, history_path_run = after_archive
        bookkeeper.run(mode=BookkeeperMode.DISCARD_FULL)
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
