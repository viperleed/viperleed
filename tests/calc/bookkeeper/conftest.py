"""Test configuration for viperleed.tests.calc.bookkeeper.

Fixtures
--------
after_calc_execution
    Return the path to a temporary directory after a bookkeeper
    execution.
mock_tree_after_calc_execution
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
from viperleed.calc import DEFAULT_WORK_HISTORY
from viperleed.calc import ORIGINAL_INPUTS_DIR_NAME
from viperleed.calc.bookkeeper.bookkeeper import Bookkeeper
from viperleed.calc.bookkeeper.bookkeeper import CALC_LOG_PREFIXES
from viperleed.calc.bookkeeper.history.constants import HISTORY_INFO_NAME
from viperleed.calc.sections.cleanup import DEFAULT_OUT
from viperleed.calc.sections.cleanup import DEFAULT_SUPP
from viperleed.calc.sections.cleanup import PREVIOUS_LABEL

from ...helpers import execute_in_dir
from . import cases_bookkeeper
from .cases_bookkeeper import BookkeeperTag as Tag
from .cases_bookkeeper import NOTES_TEST_CONTENT


ALT_HISTORY_NAME = 'history_alt_name'
MOCK_INPUT_CONTENT = 'This is a test input file.'
MOCK_JOB_NAMES = (None, 'test_jobname')
MOCK_ORIG_CONTENT = 'This is a test original input file.'
MOCK_OUT_CONTENT = 'This is a test output file.'
MOCK_STATE_FILES = ('POSCAR', 'VIBROCC', 'PARAMETERS')
MOCK_TIMESTAMP = '010203-040506'
MOCK_LOG_FILES = [f'{pre}-{MOCK_TIMESTAMP}.log' for pre in CALC_LOG_PREFIXES]
MOCK_WORKHISTORY = {  # name in workhistory: name in history
    f't003.r001_RDS_{MOCK_TIMESTAMP}': f't003.r001.001_RDS_{MOCK_TIMESTAMP}',
    f't004.r005_DS_{MOCK_TIMESTAMP}': f't004.r001.005_DS_{MOCK_TIMESTAMP}',
    f't003.r000_{PREVIOUS_LABEL}_xxx': None,  # deleted
    }


with_history_name = parametrize(
    history_name=(DEFAULT_HISTORY, ALT_HISTORY_NAME)
    )
with_jobs = parametrize(job_name=MOCK_JOB_NAMES)
with_logs = parametrize(log_file_name=MOCK_LOG_FILES)


@fixture(name='mock_tree_after_calc_execution')
@with_logs
@parametrize_with_cases('history_info_contents',
                        cases=cases_bookkeeper,
                        has_tag=Tag.BOOKKEEPER)
def fixture_mock_tree_after_calc_execution(tmp_path, log_file_name,
                                           history_info_contents):
    """Yield a temporary directory for testing the bookkeeper."""
    deltas_path = tmp_path / 'Deltas'
    tensors_path = tmp_path / 'Tensors'
    out_path = tmp_path / DEFAULT_OUT
    supp_path = tmp_path / DEFAULT_SUPP
    original_inputs_path = supp_path / ORIGINAL_INPUTS_DIR_NAME
    directories = [
        deltas_path,
        tensors_path,
        out_path,
        supp_path,
        original_inputs_path,
        tmp_path / ALT_HISTORY_NAME / 't001.r001_20xxxx-xxxxxx',
        tmp_path / ALT_HISTORY_NAME / 't002.r002_20xxxx-xxxxxx',
        ]
    directories.extend(tmp_path/DEFAULT_WORK_HISTORY/f
                       for f in MOCK_WORKHISTORY)

    for directory in directories:
        directory.mkdir(parents=True, exist_ok=True)

    files = {  # path: contents
        tmp_path / log_file_name: None,
        tmp_path / 'notes.txt': NOTES_TEST_CONTENT,
        deltas_path / 'Deltas_003.zip': None,
        tensors_path / 'Tensors_003.zip': None,
        }
    # Inputs in root
    files.update((tmp_path / f, MOCK_INPUT_CONTENT)
                 for f in MOCK_STATE_FILES)
    # Original inputs in SUPP
    files.update((original_inputs_path / f, MOCK_ORIG_CONTENT)
                 for f in MOCK_STATE_FILES)
    # OUT
    files.update((out_path / f'{f}_OUT', MOCK_OUT_CONTENT)
                 for f in MOCK_STATE_FILES)
    # history.info
    if history_info_contents is not None:
        files[tmp_path / HISTORY_INFO_NAME] = history_info_contents

    # workhistory
    files.update((tmp_path/DEFAULT_WORK_HISTORY/f/'file', f)
                 for f in MOCK_WORKHISTORY)

    for file, contents in files.items():
        if contents is None:
            file.touch()
        else:
            file.write_text(contents, encoding='utf-8')

    with execute_in_dir(tmp_path):
        yield tmp_path
    # It would be nice to clean up, but the following line causes
    # a PermissionError. Likely because of logging keeping a hold
    # of the bookkeeper.log file.
    # shutil.rmtree(tmp_path)


@fixture
@with_jobs
@with_history_name
def after_calc_execution(mock_tree_after_calc_execution,
                         job_name, history_name):
    """Return the bookkeeper, and the path to the history subfolder."""
    bookkeeper = Bookkeeper(cwd=mock_tree_after_calc_execution,
                            job_name=job_name,
                            history_name=history_name)
    bookkeeper.update_from_cwd(silent=True)
    history_path = bookkeeper.cwd / history_name
    dir_name = f't003.r001_{MOCK_TIMESTAMP}'
    if job_name is not None:
        dir_name += f'_{job_name}'
    history_run_path = history_path / dir_name
    return bookkeeper, history_run_path
