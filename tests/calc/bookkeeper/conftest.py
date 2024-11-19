"""Test configuration for viperleed.tests.calc.bookkeeper.

Fixtures
--------
after_archive
    Prepare a directory like the one after ARCHIVE was executed.
after_calc_execution
    Prepare a directory like the one after calc executes.
before_calc_execution
    Return a bookkeeper ready to run in a directory with calc inputs.
mock_tree_after_calc_execution
    Factory that produces and returns a temporary directory with
    contents like after a calc run.
mock_tree_before_calc_execution
    Factory that produces and returns a temporary directory with
    input files like those of calc.
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
from viperleed.calc.bookkeeper.mode import BookkeeperMode
from viperleed.calc.constants import DEFAULT_DELTAS
from viperleed.calc.constants import DEFAULT_HISTORY
from viperleed.calc.constants import DEFAULT_OUT
from viperleed.calc.constants import DEFAULT_SUPP
from viperleed.calc.constants import DEFAULT_TENSORS
from viperleed.calc.constants import DEFAULT_WORK_HISTORY
from viperleed.calc.constants import ORIGINAL_INPUTS_DIR_NAME
from viperleed.calc.sections.cleanup import PREVIOUS_LABEL

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
def factory_mock_tree_after_calc_execution(log_file_name, with_notes,
                                           history_info_contents, tmp_path):
    """Return a temporary directory with contents like after a calc run."""
    def _make():
        root_contents = {
            DEFAULT_DELTAS: {f'{DEFAULT_DELTAS}_004.zip': None},
            DEFAULT_TENSORS: {f'{DEFAULT_TENSORS}_004.zip': None},
            log_file_name: None,

            # Input files
            **{f: MOCK_INPUT_CONTENT for f in MOCK_STATE_FILES},

            # OUT files
            DEFAULT_OUT : {f'{f}_OUT': MOCK_OUT_CONTENT
                           for f in MOCK_STATE_FILES},

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

        return tmp_path
    return _make


@fixture(name='mock_tree_before_calc_execution')
def factory_mock_tree_before_calc_execution(tmp_path):
    """Return a directory with input files like before calc runs."""
    def make_():
        input_files = {f: MOCK_INPUT_CONTENT for f in MOCK_STATE_FILES}
        filesystem_from_dict(input_files, tmp_path)
        return tmp_path
    return make_


@fixture(name='after_archive')
def fixture_after_archive(after_calc_execution):
    """Prepare a directory like the one after ARCHIVE was executed.

    Parameters
    ----------
    after_calc_execution : fixture
        A bookkeeper ready to run in a directory with contents
        similar to those after calc has finished running, and
        a path to the main history directory that should be
        created by bookkeeper in ARCHIVE mode.

    Returns
    -------
    bookkeeper : Bookkeeper
        A Bookkeeper ready to run in a root directory that
        has been just archived. It is also up to date with
        the contents of such directory.
    history_run_path : Path
        Path to the main history subfolder created by `bookkeeper`
        as a result of the archiving.
    """
    bookkeeper, *_ = after_calc_execution
    bookkeeper.run(mode=BookkeeperMode.ARCHIVE)
    bookkeeper.update_from_cwd(silent=True)
    return after_calc_execution


@fixture(name='after_calc_execution')
def fixture_after_calc_execution(mock_tree_after_calc_execution):
    """Prepare a directory like the one after calc executes.

    Parameters
    ----------
    mock_tree_after_calc_execution : fixture
        Factory that produces a temporary directory containing
        files like those produced by run_calc when called.

    Returns
    -------
    bookkeeper : Bookkeeper
        A Bookkeeper ready to run in a the temporary directory.
        It is also up to date with the contents of such directory.
    history_run_path : Path
        Path to the main history subfolder that would be created
        by `bookkeeper` as a result of the archiving.
    """
    tmp_path = mock_tree_after_calc_execution()
    bookkeeper = Bookkeeper(cwd=tmp_path)
    bookkeeper.update_from_cwd(silent=True)
    history_path = bookkeeper.cwd / DEFAULT_HISTORY
    history_run_path = history_path / f't004.r001_{MOCK_TIMESTAMP}'
    return bookkeeper, history_run_path


@fixture(name='before_calc_execution')
def fixture_before_calc_execution(mock_tree_before_calc_execution):
    """Return a bookkeeper ready to run in a directory with calc inputs.

    This represents a new calculation, i.e., before any viperleed.calc
    or bookkeeper run. However, the bookkeeper is up to date about the
    contents of its root directory.

    Parameters
    ----------
    mock_tree_before_calc_execution : fixture
        Factory that produces a temporary directory containing
        calc input files when called.

    Returns
    ------
    bookkeeper : Bookkeeper
        A bookkeeper instance ready for running in the
        temporary directory.
    """
    # Create mock input files
    tmp_path = mock_tree_before_calc_execution()
    bookkeeper = Bookkeeper(cwd=tmp_path)
    bookkeeper.update_from_cwd(silent=True)
    return bookkeeper
