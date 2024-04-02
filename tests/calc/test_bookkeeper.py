"""Tests for module viperleed.calc.bookkeeper."""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-08-02'
__license__ = 'GPLv3+'

import shutil

from pytest_cases import fixture, parametrize

from viperleed.calc import DEFAULT_HISTORY
from viperleed.calc import DEFAULT_WORK
from viperleed.calc import ORIGINAL_INPUTS_DIR_NAME
from viperleed.calc.bookkeeper import BookkeeperMode
from viperleed.calc.bookkeeper import _CALC_LOG_PREFIXES
from viperleed.calc.bookkeeper import bookkeeper
from viperleed.calc.bookkeeper import store_input_files_to_history

from ..helpers import execute_in_dir


MOCK_TIMESTAMP = '010203-040506'
MOCK_LOG_FILES = [f'{pre}-{MOCK_TIMESTAMP}.log' for pre in _CALC_LOG_PREFIXES]
NOTES_TEST_CONTENT = 'This is a test note.'
MOCK_FILES = ('POSCAR', 'VIBROCC')
MOCK_INPUT_CONTENT = 'This is a test input file.'
MOCK_ORIG_CONTENT = 'This is a test original input file.'
MOCK_OUT_CONTENT = 'This is a test output file.'


@fixture(name='bookkeeper_mock_dir')
@parametrize(log_file_name=MOCK_LOG_FILES)
def fixture_bookkeeper_mock_dir(tmp_path, log_file_name):
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

    # create mock input files
    for file in MOCK_FILES:
        (tmp_path / file).write_text(MOCK_INPUT_CONTENT)

    # create mock original_inputs folder
    original_inputs_path = work_path / ORIGINAL_INPUTS_DIR_NAME
    original_inputs_path.mkdir(parents=True, exist_ok=True)
    for file in MOCK_FILES:
        (original_inputs_path / file).write_text(MOCK_ORIG_CONTENT)

    # create mock OUT folder
    for file in MOCK_FILES:
        (out_path / f'{file}_OUT_{MOCK_TIMESTAMP}').write_text(MOCK_OUT_CONTENT)

    with execute_in_dir(tmp_path):
        yield tmp_path
    shutil.rmtree(tmp_path)


@fixture(name='history_path')
def fixture_history_path(bookkeeper_mock_dir):
    """Return the path to a history subfolder of `bookkeeper_mock_dir`."""
    return bookkeeper_mock_dir / DEFAULT_HISTORY


@fixture(name='history_path_run')
def fixture_history_path_run(history_path):
    """Return the path to a history run subfolder of history_path."""
    return history_path / f't000.r001_{MOCK_TIMESTAMP}'


def test_bookkeeper_mode_enum():
    """Check values of bookkeeper mode enum."""
    assert BookkeeperMode.DEFAULT is BookkeeperMode('default')
    assert BookkeeperMode.CONT is BookkeeperMode('cont')
    assert BookkeeperMode.DISCARD is BookkeeperMode('discard')


def test_store_input_files_to_history(tmp_path, bookkeeper_mock_dir):
    """Check correct storage of original input files to history."""
    inputs_path = bookkeeper_mock_dir / DEFAULT_WORK / ORIGINAL_INPUTS_DIR_NAME
    history = tmp_path / DEFAULT_HISTORY
    history.mkdir(parents=True, exist_ok=True)
    store_input_files_to_history(inputs_path, history)
    for file in MOCK_FILES:
        hist_file = history / file
        hist_file_content = hist_file.read_text()
        assert MOCK_ORIG_CONTENT in hist_file_content


def test_bookkeeper_default_mode(bookkeeper_mock_dir,
                                 history_path,
                                 history_path_run):
    """Check correct storage of history files in DEFAULT mode."""
    # os.chdir(bookkeeper_mock_dir)
    bookkeeper(mode=BookkeeperMode.DEFAULT)
    assert history_path.exists()
    assert history_path_run.is_dir()
    assert (bookkeeper_mock_dir / 'history.info').exists()
    # out stored in history
    for file in MOCK_FILES:
        hist_file = history_path_run / 'OUT' / f'{file}_OUT_{MOCK_TIMESTAMP}'
        hist_content = hist_file.read_text()
        assert MOCK_OUT_CONTENT in hist_content
    # original not overwritten
    for file in MOCK_FILES:
        input_content = (bookkeeper_mock_dir / file).read_text()
        assert MOCK_INPUT_CONTENT in input_content


def test_bookkeeper_cont_mode(bookkeeper_mock_dir, history_path):
    """Check correct overwriting of input files in continuation mode."""
    # os.chdir(bookkeeper_mock_dir)
    bookkeeper(mode=BookkeeperMode.CONT)
    assert history_path.exists()
    # make sure input was overwritten
    for file in MOCK_FILES:
        input_content = (bookkeeper_mock_dir / file).read_text()
        assert MOCK_OUT_CONTENT in input_content


def test_bookkeeper_discard_mode(bookkeeper_mock_dir, history_path_run):
    """Check correct skipping of history storage in DISCARD mode."""
    # os.chdir(bookkeeper_mock_dir)
    bookkeeper(mode=BookkeeperMode.DISCARD)
    assert not history_path_run.is_dir()
    # original not overwritten
    for file in MOCK_FILES:
        input_content = (bookkeeper_mock_dir / file).read_text()
        assert MOCK_INPUT_CONTENT in input_content


def test_bookkeeper_with_job_name(history_path):
    """Check correct history storage when a specific job name is set."""
    # os.chdir(bookkeeper_mock_dir)
    bookkeeper(mode='default', job_name='test_job')
    assert (history_path / f't000.r001_{MOCK_TIMESTAMP}_test_job').exists()


def test_bookkeeper_with_existing_history_and_alt_name(bookkeeper_mock_dir):
    """Check correct storage with a non-empty, differently named history."""
    # os.chdir(bookkeeper_mock_dir)
    # Create some existing history folders
    hist_name = 'history_alt_name'
    alt_history = bookkeeper_mock_dir / hist_name
    for run_name in ('t001.r001_20xxxx-xxxxxx', 't002.r002_20xxxx-xxxxxx'):
        (alt_history / run_name).mkdir(parents=True, exist_ok=True)
    (bookkeeper_mock_dir / 'Tensors' / 'Tensors_003.zip').touch()
    bookkeeper(mode='default', history_name=hist_name)
    assert (alt_history / f't003.r001_{MOCK_TIMESTAMP}').exists()
