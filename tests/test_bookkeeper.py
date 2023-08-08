"""Tests for functionality of bookkeeper.py module.

Created on 2023-08-02

@author: Alex M. Imre
"""

import os
import shutil
import pytest
from enum import Enum
from pathlib import Path

from viperleed.calc.bookkeeper import bookkeeper, BookkeeperMode, store_input_files_to_history

MOCK_TIMESTAMP = "010203-040506"
MOCK_LOG_FILE = f"tleedm-{MOCK_TIMESTAMP}.log"
NOTES_TEST_CONTENT = "This is a test note."
MOCK_FILES = ['POSCAR', 'VIBROCC']
MOCK_INPUT_CONTENT = "This is a test input file."
MOCK_ORIG_CONTENT = "This is a test original input file."
MOCK_OUT_CONTENT = "This is a test output file."

# Define the fixtures needed for testing

@pytest.fixture(scope="function")
def bookkeeper_mock_dir(tmp_path):
    """Creates a temporary directory for testing."""
    work_path = tmp_path / "work"
    OUT_path = tmp_path / "OUT"
    SUPP_path = tmp_path / "SUPP"
    tensors_path = tmp_path / "Tensors"
    deltas_path = tmp_path / "Deltas"
    for dir in [work_path, OUT_path, SUPP_path, tensors_path, deltas_path]:
        dir.mkdir(parents=True, exist_ok=True)
    # create mock log files
    (tmp_path / MOCK_LOG_FILE).touch()
    (work_path / MOCK_LOG_FILE).touch()
    # create mock notes file
    notes_file = tmp_path / "notes.txt"
    notes_file.write_text(NOTES_TEST_CONTENT)
    home = Path.cwd()

    # create mock input files
    for file in MOCK_FILES:
        (tmp_path / file).write_text(MOCK_INPUT_CONTENT)

    # create mock original_inputs folder
    original_inputs_path = work_path / "original_inputs"
    original_inputs_path.mkdir(parents=True, exist_ok=True)
    for file in MOCK_FILES:
        (original_inputs_path / file).write_text(MOCK_ORIG_CONTENT)

    # create mock OUT folder
    for file in MOCK_FILES:
        (OUT_path / f"{file}_OUT_{MOCK_TIMESTAMP}").write_text(MOCK_OUT_CONTENT)

    os.chdir(tmp_path)
    yield tmp_path
    os.chdir(home)
    shutil.rmtree(tmp_path)


def test_bookkeeper_mode_enum():
    assert BookkeeperMode.DEFAULT is BookkeeperMode('default')
    assert BookkeeperMode.CONT is BookkeeperMode('cont')
    assert BookkeeperMode.DISCARD is BookkeeperMode('discard')


def test_store_input_files_to_history(tmp_path, bookkeeper_mock_dir):
    inputs_path = bookkeeper_mock_dir / "work" / "original_inputs"
    history_path = tmp_path / "history"
    history_path.mkdir(parents=True, exist_ok=True)
    store_input_files_to_history(inputs_path, history_path)
    for file in MOCK_FILES:
        hist_file = history_path / file
        hist_file_content = hist_file.read_text()
        assert MOCK_ORIG_CONTENT in hist_file_content


def test_bookkeeper_default_mode(bookkeeper_mock_dir):
    os.chdir(bookkeeper_mock_dir)
    bookkeeper(mode=BookkeeperMode.DEFAULT)
    hist_path = bookkeeper_mock_dir / "history"
    hist_path_run = hist_path / f't000.r001_{MOCK_TIMESTAMP}'
    assert (hist_path).exists()
    assert (bookkeeper_mock_dir / "history.info").exists()
    assert (hist_path_run).is_dir()
    # out stored in history
    for file in MOCK_FILES:
        mock_hist_file = hist_path_run / 'OUT' / f'{file}_OUT_{MOCK_TIMESTAMP}'
        hist_content = mock_hist_file.read_text()
        assert MOCK_OUT_CONTENT in hist_content
    # original not overwritten
    for file in MOCK_FILES:
        input_content = (bookkeeper_mock_dir / file).read_text()
        assert MOCK_INPUT_CONTENT in input_content


def test_bookkeeper_cont_mode(bookkeeper_mock_dir):
    os.chdir(bookkeeper_mock_dir)
    bookkeeper(mode=BookkeeperMode.CONT)
    hist_path = bookkeeper_mock_dir / "history"
    hist_path_run = hist_path / f't000.r001_{MOCK_TIMESTAMP}'
    assert (hist_path).exists()
    # make sure input was overwritten
    for file in MOCK_FILES:
        input_content = (bookkeeper_mock_dir / file).read_text()
        assert MOCK_OUT_CONTENT in input_content


def test_bookkeeper_discard_mode(bookkeeper_mock_dir):
    os.chdir(bookkeeper_mock_dir)
    bookkeeper(mode=BookkeeperMode.DISCARD)
    hist_path = bookkeeper_mock_dir / "history"
    hist_path_run = hist_path / f't000.r001_{MOCK_TIMESTAMP}'
    assert not (hist_path_run).is_dir()
    # original not overwritten
    for file in MOCK_FILES:
        input_content = (bookkeeper_mock_dir / file).read_text()
        assert MOCK_INPUT_CONTENT in input_content


def test_bookkeeper_with_job_name(bookkeeper_mock_dir):
    os.chdir(bookkeeper_mock_dir)
    bookkeeper(mode="default", job_name="test_job")
    assert (bookkeeper_mock_dir / "history" / f"t000.r001_{MOCK_TIMESTAMP}_test_job").exists()


def test_bookkeeper_with_existing_history_and_alt_name(bookkeeper_mock_dir):
    os.chdir(bookkeeper_mock_dir)
    # Create some existing history folders
    Path(bookkeeper_mock_dir / "history_alt_name" / "t001.r001_20xxxx-xxxxxx").mkdir(parents=True, exist_ok=True)
    Path(bookkeeper_mock_dir / "history_alt_name" / "t002.r002_20xxxx-xxxxxx").mkdir(parents=True, exist_ok=True)
    (bookkeeper_mock_dir / 'Tensors' / 'Tensors_003.zip').touch()
    bookkeeper(mode="default", history_name="history_alt_name")
    assert (bookkeeper_mock_dir / "history_alt_name" / f"t003.r001_{MOCK_TIMESTAMP}").exists()
