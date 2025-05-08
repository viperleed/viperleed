"""Test configuration for viperleed.tests.calc.bookkeeper.

Fixtures
--------
after_archive
    Prepare a directory like the one after ARCHIVE was executed.
after_bookkeper_run
    Factory that runs bookkeeper in a given mode on a given
    pre-run tree.
after_calc_execution
    Prepare a directory like the one after calc executes.
before_calc_execution
    Return a bookkeeper ready to run in a directory with calc inputs.
check_methods_called
    Patch selected methods of an object, and ensure that calling
    another method also calls the patched ones.
mock_attributes
    Replace multiple attributes of an object with MagicMock()s.
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
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2023-08-02'
__license__ = 'GPLv3+'

from copy import deepcopy
from operator import attrgetter

from pytest_cases import fixture
from pytest_cases import fixture_ref
from pytest_cases import parametrize
from pytest_cases import parametrize_with_cases

from viperleed.calc.bookkeeper.bookkeeper import Bookkeeper
from viperleed.calc.bookkeeper.constants import CALC_LOG_PREFIXES
from viperleed.calc.bookkeeper.constants import EDITED_SUFFIX
from viperleed.calc.bookkeeper.constants import ORI_SUFFIX
from viperleed.calc.bookkeeper.history.constants import HISTORY_INFO_NAME
from viperleed.calc.bookkeeper.history.meta import _METADATA_NAME
from viperleed.calc.bookkeeper.mode import BookkeeperMode
from viperleed.calc.constants import DEFAULT_DELTAS
from viperleed.calc.constants import DEFAULT_HISTORY
from viperleed.calc.constants import DEFAULT_OUT
from viperleed.calc.constants import DEFAULT_SUPP
from viperleed.calc.constants import DEFAULT_TENSORS
from viperleed.calc.constants import DEFAULT_WORK_HISTORY
from viperleed.calc.constants import LOG_PREFIX
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
MOCK_OUT_SUFFIXED_CONTENT = 'This is a test output file with _OUT in its name.'
MOCK_STATE_FILES = ('POSCAR', 'VIBROCC', 'PARAMETERS')
# Notice that the year for the timestamp is such that all files
# created in the tests are considered "not later than". This
# is to prevent labeling files as "_edited". As of 2024, "68"
# is translated to 2068, as per POSIX specification of stptime:
# https://pubs.opengroup.org/onlinepubs/9799919799/functions/strptime.html
MOCK_TIMESTAMP = '680203-040506'
MOCK_LOG_FILES = [f'{pre}-{MOCK_TIMESTAMP}.log' for pre in CALC_LOG_PREFIXES]
MOCK_WORKHISTORY = {  # name in workhistory: name in history
    # One run involving an already existing, previous tensor
    f't002.r009_DS_{MOCK_TIMESTAMP}': f't002.r003.009_DS_{MOCK_TIMESTAMP}',
    # One run involving a tensor that is not yet in history
    f't003.r005_DS_{MOCK_TIMESTAMP}': f't003.r001.005_DS_{MOCK_TIMESTAMP}',
    # A previous one, that should be deleted
    f't003.r000_{PREVIOUS_LABEL}_xxx': None,
    # And a couple runs from the current tensor
    f't004.r001_RDS_{MOCK_TIMESTAMP}': f't004.r001.001_RDS_{MOCK_TIMESTAMP}',
    f't004.r002_RDS_{MOCK_TIMESTAMP}': f't004.r001.002_RDS_{MOCK_TIMESTAMP}',
    }


@fixture(name='after_bookkeper_run')
def factory_after_bookkeeper_run():
    """Prepare a directory like the one after `mode` was executed."""
    def _make(before_this_run, mode, **kwargs):
        bookkeeper, *_ = before_this_run
        bookkeeper.run(mode=mode, **kwargs)
        bookkeeper.update_from_cwd(silent=True)
        return before_this_run
    return _make


@fixture(name='_make_root_tree')
def factory_make_root_tree(tmp_path):
    """Populate a temporary directory with files."""
    default_root_contents = {
        DEFAULT_DELTAS: {f'{DEFAULT_DELTAS}_004.zip': None},
        DEFAULT_TENSORS: {f'{DEFAULT_TENSORS}_004.zip': None},

        # Input files
        **{f: MOCK_INPUT_CONTENT for f in MOCK_STATE_FILES},

        # OUT files
        DEFAULT_OUT : {f: MOCK_OUT_CONTENT for f in MOCK_STATE_FILES},

        # Original inputs in SUPP
        f'{DEFAULT_SUPP}/{ORIGINAL_INPUTS_DIR_NAME}': {
            f: MOCK_ORIG_CONTENT for f in MOCK_STATE_FILES
            },

        # Pre-existing (empty) history directories
        DEFAULT_HISTORY: {'t001.r001_200101-111213/some_directory': {},
                          't001.r002.001_DS_200212-141516/some_directory': {},
                          't002.r002_200212-141516/some_directory': {}},

        # workhistory subfolders, with dummy files
        DEFAULT_WORK_HISTORY: {f: {'file': f}
                               for f in MOCK_WORKHISTORY},
        }

    def _make(skip=None, at_path=tmp_path, **kwargs):
        if skip is None:
            skip = set()
        def _populate_root():
            root_contents = {f: c for f, c in default_root_contents.items()
                             if f not in skip}
            root_contents.update(kwargs)
            filesystem_from_dict(root_contents, at_path)
            return at_path
        return _populate_root

    return _make


@fixture(name='mock_tree_after_calc_execution')
@parametrize(log_file_name=MOCK_LOG_FILES)
@parametrize(with_notes=(True, False))
@parametrize_with_cases('history_info_contents',
                        cases=(cases_entry, cases_history),
                        has_tag=Tag.BOOKKEEPER)
def factory_mock_tree_after_calc_execution(log_file_name,
                                           with_notes,
                                           history_info_contents,
                                           _make_root_tree):
    """Return a temporary directory with contents like after a calc run."""
    kwargs = {
        log_file_name: None,
        }
    # Files/folders that depend on the arguments
    if with_notes:
        kwargs['notes.txt'] = NOTES_TEST_CONTENT
    if history_info_contents is not None:  # history.info
        kwargs[HISTORY_INFO_NAME] = history_info_contents
    return _make_root_tree(**kwargs)


@fixture(name='mock_tree_after_calc_execution_with_out_suffix')
def factory_mock_tree_after_calc_execution_out_suffix(_make_root_tree):
    """Create files like those after calc HAD run before we dropped _OUT."""
    old_style_out = {
        f'{f}_OUT_010203-040506': MOCK_OUT_SUFFIXED_CONTENT
        for f in ('POSCAR', 'VIBROCC')
        }
    # In < v0.13.0, there was no PARAMETERS in OUT. It was in root.
    # old_style_out['PARAMETERS'] = MOCK_OUT_CONTENT
    kwargs = {
        f'{LOG_PREFIX}_{MOCK_TIMESTAMP}.log': None,

        # OUT files, with an old-style _OUT suffix
        DEFAULT_OUT : old_style_out,
        }
    return _make_root_tree(**kwargs)


@fixture(name='mock_tree_before_calc_execution')
def factory_mock_tree_before_calc_execution(tmp_path):
    """Return a directory with input files like before calc runs."""
    def make_():
        input_files = {f: MOCK_INPUT_CONTENT for f in MOCK_STATE_FILES}
        filesystem_from_dict(input_files, tmp_path)
        return tmp_path
    return make_


@fixture(name='after_archive')
def fixture_after_archive(after_calc_execution, after_bookkeper_run):
    """Prepare a directory like the one after ARCHIVE was executed.

    Parameters
    ----------
    after_calc_execution : fixture
        A bookkeeper ready to run in a directory with contents
        similar to those after calc has finished running, and
        a path to the main history directory that should be
        created by bookkeeper in ARCHIVE mode.
    after_bookkeper_run : fixture
        The after_bookkeper_run fixture factory, that will execute
        bookkeeper in ARCHIVE mode on the tree `after_calc_execution`.

    Returns
    -------
    bookkeeper : Bookkeeper
        A Bookkeeper ready to run in a root directory that
        has been just archived. It is also up to date with
        the contents of such directory.
    history_run_path : Path
        Path to the main history subfolder created by `bookkeeper`
        as a result of the archiving.
    mocker : fixture
        The pytest-mock fixture.
    """
    return after_bookkeper_run(after_calc_execution, BookkeeperMode.ARCHIVE)


_tree_after_calc_exec = (
    fixture_ref('mock_tree_after_calc_execution'),
    fixture_ref('mock_tree_after_calc_execution_with_out_suffix'),
    )

@fixture(name='after_calc_execution')
@parametrize(make_calc_tree=_tree_after_calc_exec)
def fixture_after_calc_execution(make_calc_tree, mocker):
    """Prepare a directory like the one after calc executes.

    Parameters
    ----------
    make_calc_tree : fixture
        Factory that produces a temporary directory containing
        files like those produced by run_calc when called.
    mocker : fixture
        The pytest-mock fixture.

    Returns
    -------
    bookkeeper : Bookkeeper
        A Bookkeeper ready to run in a the temporary directory.
        It is also up to date with the contents of such directory.
    history_run_path : Path
        Path to the main history subfolder that would be created
        by `bookkeeper` as a result of the archiving.
    mocker : fixture
        The pytest-mock fixture.
    """
    tmp_path = make_calc_tree()
    bookkeeper = Bookkeeper(cwd=tmp_path)
    bookkeeper.update_from_cwd(silent=True)
    history_path = bookkeeper.cwd / DEFAULT_HISTORY
    history_run_path = history_path / f't004.r001_{MOCK_TIMESTAMP}'
    return bookkeeper, history_run_path, mocker


@fixture(name='after_calc_with_edited_file')
def fixture_after_calc_with_edited(tmp_path):
    """Prepare a directory like the one after calc executes.

    However, modify one file as if the user had done so after
    calc started running. The tree is structured as in
    https://github.com/viperleed/viperleed/pull/198#issuecomment-2506549827.

    Parameters
    ----------
    tmp_path : fixture
        The pytest fixture for creating a temporary directory.

    Returns
    -------
    bookkeeper : Bookkeeper
        A bookkeeper instance ready to run in `tmp_path`.
    before : dict
        A dictionary representation of the current contents
        of `tmp_path`.
    archived : dict
        A dictionary representation of the expected contents
        of `tmp_path` after a successful ARCHIVE run.
    """
    # Use a time in the past. All files will be created AFTER this,
    # but we also check their contents to decide whether they were
    # edited since.
    calc_started = '250128-235532'
    before = {
        f'{LOG_PREFIX}-{calc_started}.log': '',
        'PARAMETERS': MOCK_ORIG_CONTENT,
        'VIBROCC': MOCK_ORIG_CONTENT,
        'POSCAR': MOCK_INPUT_CONTENT,  # This will be marked as _edited
        DEFAULT_SUPP : {
            ORIGINAL_INPUTS_DIR_NAME: {
                'PARAMETERS': MOCK_ORIG_CONTENT,  # Same as root file
                'VIBROCC': MOCK_ORIG_CONTENT,     # Same as root file
                'POSCAR': MOCK_ORIG_CONTENT,  # Differs from root one
                },
            },
        DEFAULT_OUT : {
            'PARAMETERS': MOCK_OUT_CONTENT,
            # Notice that there is no POSCAR nor VIBROCC in OUT
            },
        }

    # Now what we expect from an ARCHIVE run:
    # (1) History folder structure
    in_history = deepcopy(before)
    in_history['POSCAR'] = (
        # The only difference is the POSCAR file, which is pulled from
        # original_inputs rather than from the root directory
        before[DEFAULT_SUPP][ORIGINAL_INPUTS_DIR_NAME]['POSCAR']
        )
    archived = deepcopy(before)
    archived[DEFAULT_HISTORY] = {f't000.r001_{calc_started}': in_history}

    # (2) The rest of the root
    archived[f'PARAMETERS{ORI_SUFFIX}'] = before['PARAMETERS']
    archived['PARAMETERS'] = before[DEFAULT_OUT]['PARAMETERS']
    archived[f'VIBROCC{ORI_SUFFIX}'] = before['VIBROCC']
    # No VIBROCC in OUT, so new input stays the same
    archived[f'POSCAR{EDITED_SUFFIX}'] = before['POSCAR']
    archived['POSCAR'] = (  # No OUT. Pulled from original_inputs.
        before[DEFAULT_SUPP][ORIGINAL_INPUTS_DIR_NAME]['POSCAR']
        )

    filesystem_from_dict(before, tmp_path)
    bookkeeper = Bookkeeper(cwd=tmp_path)
    bookkeeper.update_from_cwd(silent=True)
    return bookkeeper, before, archived


@fixture(name='archived_domains')
def fixture_archived_domains(after_bookkeper_run,
                             domains_after_calc_execution):
    """Run ARCHIVE in a multi-domain root."""
    *_, domains = domains_after_calc_execution
    return after_bookkeper_run(domains_after_calc_execution,
                               'archive',
                               domains=domains)


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


@fixture(name='check_methods_called')
def fixture_check_methods_called(mock_attributes):
    """Check that calling `method_name` also calls other methods."""
    def check_methods_called(obj, method_name, **called):
        mock_attributes(obj, *called)
        method = attrgetter(method_name)(obj)
        method()

        for mocked_attr, called_attr in called.items():
            method = attrgetter(called_attr or mocked_attr)(obj)
            method.assert_called_once()
    return check_methods_called


@fixture(name='domains_after_calc_execution')
def fixture_domains_after_calc_execution(_make_root_tree, mocker):
    """Create a temporary tree with the results of a domain calculation."""
    older_timestamp = '010203-040506'
    root_kwargs = {
        # No Deltas/Tensors folders in the root of a domains
        # calculation. Skip workhistory because the default
        # folder names are inconsistent.
        'skip': {DEFAULT_DELTAS, DEFAULT_TENSORS,
                 DEFAULT_WORK_HISTORY},
        # Since no Tensors in the root, all history folders here have
        # a t000 index.
        DEFAULT_HISTORY: {f't000.r003_{older_timestamp}': {
            _METADATA_NAME: '''
[archived]
hash = c36f643af62014636ed9d7cd0905b486

[domains]
domains = (('Domain_1', '6bb858254252e34861f400ca03ec4866'),
           ('Domain_two', 'b9e631993d9da320902b52047547041d'),
           ('was_already_here', 'f37f2f0bc6f3c7b2196e211614a72c24'))
''',
            }},
        f'{LOG_PREFIX}_{MOCK_TIMESTAMP}.log': 'Executed segments: 0 1 2',
        }
    root_factory = _make_root_tree(**root_kwargs)
    tmp_path = root_factory()
    main_history_run = tmp_path/DEFAULT_HISTORY/f't000.r004_{MOCK_TIMESTAMP}'

    # Now domains
    domain_kwargs = {
        # One domain with an auto-generated name and no user flag
        'Domain_1': {
            DEFAULT_TENSORS: {'Tensors_004.zip': None,
                              'Tensors_005.zip': None},
            DEFAULT_HISTORY: {f't005.r004_{older_timestamp}': {
                _METADATA_NAME: f'''
[archived]
hash = 6bb858254252e34861f400ca03ec4866

[domains]
main = ('{tmp_path.as_posix()}', 'c36f643af62014636ed9d7cd0905b486')
'''}},
            'history_run': f't005.r005_{MOCK_TIMESTAMP}',
            },
        # One domain with an auto-generated name, and a user flag
        'Domain_two': {
            DEFAULT_TENSORS: {'Tensors_004.zip': None,
                              'Tensors_006.zip': None},
            DEFAULT_HISTORY: {f't006.r008_{older_timestamp}': {
                _METADATA_NAME: f'''
[archived]
hash = b9e631993d9da320902b52047547041d

[domains]
main = ('{tmp_path.as_posix()}', 'c36f643af62014636ed9d7cd0905b486')
'''}},
            'history_run': f't006.r009_{MOCK_TIMESTAMP}',
            },
        # And one domain that was not pulled from somewhere else
        'was_already_here': {
            DEFAULT_TENSORS: {'Tensors_099.zip': None,
                              'Tensors_105.zip': None},  # New
            DEFAULT_HISTORY: {f't099.r001_{older_timestamp}': {
            _METADATA_NAME: f'''
[archived]
hash = f37f2f0bc6f3c7b2196e211614a72c24

[domains]
main = ('{tmp_path.as_posix()}', 'c36f643af62014636ed9d7cd0905b486')
'''}},
            'history_run': f't105.r001_{MOCK_TIMESTAMP}',
            },
        }
    domain_info = {}
    for domain_folder, kwargs in domain_kwargs.items():
        kwargs['at_path'] = domain_root = tmp_path/domain_folder
        domain_info[domain_root] = {
            'history_run': (
                domain_root/DEFAULT_HISTORY/kwargs.pop('history_run')
                ),
            }
        make_domain = _make_root_tree(**kwargs)
        make_domain()

    bookkeeper = Bookkeeper(tmp_path)
    bookkeeper.update_from_cwd(silent=True)
    return bookkeeper, main_history_run, mocker, domain_info


@fixture(name='manual_run_one_domain')
def fixture_manual_run_one_domain(domains_after_calc_execution):
    """Produce a tree as after a user manually ran bookkeeper in one domain."""
    def _run_one_domain(mode):
        *_, domains = domains_after_calc_execution
        manual_run = next(iter(domains))
        bookkeeper = Bookkeeper(manual_run)
        bookkeeper.run(mode)
        return domains_after_calc_execution, manual_run
    return _run_one_domain


@fixture(name='mock_attributes')
def fixture_mock_attributes(mocker):
    """Replace attributes of `obj` with MagicMock(s)."""
    def _patch(obj, *attrs):
        return {attr_name: mocker.patch.object(obj, attr_name)
                for attr_name in attrs}
    return _patch
