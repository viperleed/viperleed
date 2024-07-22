"""Test configuration for viperleed.tests.calc.bookkeeper.

Fixtures
--------
after_run
    Return the path to a temporary directory after a bookkeeper
    execution.
bookkeeper_mock_dir_after_run
    Yield a temporary directory for testing the bookkeeper.
"""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-08-02'
__license__ = 'GPLv3+'


from pytest_cases import fixture
from pytest_cases import parametrize
from pytest_cases import parametrize_with_cases

from viperleed.calc import DEFAULT_HISTORY
from viperleed.calc import ORIGINAL_INPUTS_DIR_NAME
from viperleed.calc.bookkeeper.bookkeeper import Bookkeeper
from viperleed.calc.bookkeeper.bookkeeper import CALC_LOG_PREFIXES
from viperleed.calc.bookkeeper.constants import HISTORY_INFO_NAME
from viperleed.calc.sections.cleanup import DEFAULT_OUT
from viperleed.calc.sections.cleanup import DEFAULT_SUPP

from ...helpers import execute_in_dir
from . import cases_bookkeeper
from .cases_bookkeeper import BookkeeperTag as Tag
from .cases_bookkeeper import NOTES_TEST_CONTENT

_QUICK = False  # Select only a subset of tests? (~factor 7 difference)


ALT_HISTORY_NAME = 'history_alt_name'
MOCK_INPUT_CONTENT = 'This is a test input file.'
MOCK_JOB_NAMES = (None, 'test_jobname')
MOCK_ORIG_CONTENT = 'This is a test original input file.'
MOCK_OUT_CONTENT = 'This is a test output file.'
MOCK_STATE_FILES = ('POSCAR', 'VIBROCC', 'PARAMETERS')
MOCK_TIMESTAMP = '010203-040506'
MOCK_LOG_FILES = [f'{pre}-{MOCK_TIMESTAMP}.log' for pre in CALC_LOG_PREFIXES]


if _QUICK:
    with_history_name = parametrize(history_name=(DEFAULT_HISTORY,))
    with_jobs = parametrize(job_name=(MOCK_JOB_NAMES[0],))
    with_logs = parametrize(log_file_name=(MOCK_LOG_FILES[0],),
                            ids=(MOCK_LOG_FILES[0],))
else:
    with_history_name = parametrize(
        history_name=(DEFAULT_HISTORY, ALT_HISTORY_NAME)
        )
    with_jobs = parametrize(job_name=MOCK_JOB_NAMES)
    with_logs = parametrize(log_file_name=MOCK_LOG_FILES, ids=MOCK_LOG_FILES)


@fixture(name='bookkeeper_mock_dir_after_run')
@with_logs
@parametrize_with_cases('history_info_contents',
                        cases=cases_bookkeeper,
                        has_tag=Tag.BOOKKEEPER)
def fixture_bookkeeper_mock_dir_after_run(tmp_path, log_file_name,
                                          history_info_contents):
    """Yield a temporary directory for testing the bookkeeper."""
    out_path = tmp_path / DEFAULT_OUT
    supp_path = tmp_path / DEFAULT_SUPP
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
    if history_info_contents is not None:
        hist_info_path = tmp_path / HISTORY_INFO_NAME
        with open(hist_info_path, 'w') as f:
            f.write(history_info_contents)
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


@fixture
@with_jobs
@with_history_name
def after_run(bookkeeper_mock_dir_after_run, job_name, history_name):
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
