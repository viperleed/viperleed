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

from viperleed.calc.bookkeeper.bookkeeper import Bookkeeper
from viperleed.calc.bookkeeper.constants import CALC_LOG_PREFIXES
from viperleed.calc.bookkeeper.history.constants import HISTORY_INFO_NAME
from viperleed.calc.constants import DEFAULT_DELTAS
from viperleed.calc.constants import DEFAULT_HISTORY
from viperleed.calc.constants import DEFAULT_OUT
from viperleed.calc.constants import DEFAULT_SUPP
from viperleed.calc.constants import DEFAULT_TENSORS
from viperleed.calc.constants import DEFAULT_WORK_HISTORY
from viperleed.calc.constants import ORIGINAL_INPUTS_DIR_NAME
from viperleed.calc.sections.cleanup import PREVIOUS_LABEL

from ...helpers import execute_in_dir
from ...helpers import filesystem_from_dict
from .history import cases_history
from .history.entry import cases_entry
from .history.entry.cases_entry import NOTES_TEST_CONTENT
from .tag import BookkeeperTag as Tag


MOCK_INPUT_CONTENT = 'This is a test input file.'
MOCK_ORIG_CONTENT = 'This is a test original input file.'
MOCK_OUT_CONTENT = 'This is a test output file.'
MOCK_STATE_FILES = ('POSCAR', 'VIBROCC', 'PARAMETERS')
MOCK_TIMESTAMP = '210203-040506'
MOCK_LOG_FILES = [f'{pre}-{MOCK_TIMESTAMP}.log' for pre in CALC_LOG_PREFIXES]
MOCK_WORKHISTORY = {  # name in workhistory: name in history
    f't003.r005_DS_{MOCK_TIMESTAMP}': f't003.r001.005_DS_{MOCK_TIMESTAMP}',
    f't003.r000_{PREVIOUS_LABEL}_xxx': None,  # deleted
    f't004.r001_RDS_{MOCK_TIMESTAMP}': f't004.r001.001_RDS_{MOCK_TIMESTAMP}',
    f't004.r002_RDS_{MOCK_TIMESTAMP}': f't004.r001.002_RDS_{MOCK_TIMESTAMP}',
    }


@fixture(name='mock_tree_after_calc_execution')
@parametrize(log_file_name=MOCK_LOG_FILES)
@parametrize(with_notes=(True,False))
@parametrize_with_cases('history_info_contents',
                        cases=(cases_entry, cases_history),
                        has_tag=Tag.BOOKKEEPER)
def fixture_mock_tree_after_calc_execution(log_file_name, with_notes,
                                           history_info_contents, tmp_path):
    """Yield a temporary directory for testing the bookkeeper."""
    root_contents = {
        DEFAULT_DELTAS: {f'{DEFAULT_DELTAS}_004.zip': None},
        DEFAULT_TENSORS: {f'{DEFAULT_TENSORS}_004.zip': None},
        log_file_name: None,

        # Input files
        **{f: MOCK_INPUT_CONTENT for f in MOCK_STATE_FILES},

        # OUT files
        DEFAULT_OUT : {f'{f}_OUT': MOCK_OUT_CONTENT for f in MOCK_STATE_FILES},

        # Original inputs in SUPP
        f'{DEFAULT_SUPP}/{ORIGINAL_INPUTS_DIR_NAME}': {
            f: MOCK_ORIG_CONTENT for f in MOCK_STATE_FILES
            },

        # Pre-existing (empty) history directories
        f'{DEFAULT_HISTORY}/t001.r001_20xxxx-xxxxxx/some_directory': {},
        f'{DEFAULT_HISTORY}/t002.r002_20xxxx-xxxxxx/some_directory': {},

        # workhistory subfolders, with dummy files
        **{f'{DEFAULT_WORK_HISTORY}/{f}': {'file': f}
           for f in MOCK_WORKHISTORY},
        }
    # Files/folders that depend on the arguments
    if with_notes:
        root_contents['notes.txt'] = NOTES_TEST_CONTENT
    if history_info_contents is not None:  # history.info
        root_contents[HISTORY_INFO_NAME] = history_info_contents

    # Actually create files and folders
    filesystem_from_dict(root_contents, tmp_path)

    with execute_in_dir(tmp_path):
        yield tmp_path


@fixture
def after_calc_execution(mock_tree_after_calc_execution):
    """Return the bookkeeper, and the path to the history subfolder."""
    bookkeeper = Bookkeeper(cwd=mock_tree_after_calc_execution)
    bookkeeper.update_from_cwd(silent=True)
    history_path = bookkeeper.cwd / DEFAULT_HISTORY
    history_run_path = history_path / f't004.r001_{MOCK_TIMESTAMP}'
    return bookkeeper, history_run_path
