"""Tests for module viperleed.calc.bookkeeper.

Collects tests for implementation details of Bookkeeper, i.e.,
all those that are not related to the overall behavior of
running bookkeeper in a specific mode.
"""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-08-02'
__license__ = 'GPLv3+'

import ast
from enum import Enum
import functools
import logging
from operator import attrgetter
from pathlib import Path
import shutil
import time

import pytest
from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.calc.bookkeeper.bookkeeper import _FROM_ROOT
from viperleed.calc.bookkeeper.bookkeeper import Bookkeeper
from viperleed.calc.bookkeeper.bookkeeper import BookkeeperExitCode
from viperleed.calc.bookkeeper.mode import BookkeeperMode as Mode
from viperleed.calc.constants import DEFAULT_DELTAS
from viperleed.calc.constants import DEFAULT_HISTORY
from viperleed.calc.constants import DEFAULT_OUT
from viperleed.calc.constants import DEFAULT_SUPP
from viperleed.calc.constants import DEFAULT_TENSORS
from viperleed.calc.constants import DEFAULT_WORK_HISTORY
from viperleed.calc.constants import LOG_PREFIX
from viperleed.calc.constants import ORIGINAL_INPUTS_DIR_NAME
from viperleed.calc.sections.cleanup import PREVIOUS_LABEL

from ....helpers import filesystem_from_dict
from ....helpers import not_raises
from ....helpers import make_obj_raise
from ....helpers import raises_exception
from ....helpers import raises_test_exception
from ..conftest import MOCK_TIMESTAMP
from ..conftest import MOCK_STATE_FILES


_UPDATE_METHODS = (
    'update_from_cwd',     # Bookkeeper
    'run',                 # Bookkeeper
    'collect',             # root_explorer.LogFiles
    'collect_info',        # root_explorer.RootExplorer
    'collect_subfolders',  # history.explorer.HistoryExplorer
    'prepare_info_file',   # history.explorer.HistoryExplorer
    )
raises_oserror = functools.partial(raises_exception, exc=OSError)
make_raise_oserror = functools.partial(make_obj_raise, exc=OSError)
not_raises_oserror = functools.partial(not_raises, exc=OSError)


def check_too_early():
    """Ensure an AttributeError is raised for a too-early getattr."""
    match_re = '|'.join(_UPDATE_METHODS)
    return pytest.raises(AttributeError, match=match_re)


class MockInput:  # pylint: disable=too-few-public-methods
    """Fake replacement for the input built-in function."""

    def __init__(self, *responses):
        """Initialize with some expected user responses."""
        self._responses = iter(responses)

    def __call__(self, *_):
        """Return a user response."""
        return next(self._responses, 'yes')


class TestBookkeeperComplaints:
    """Tests for situations that do not raise but issue log warnings/errors."""

    def test_cwd_file_newer(self, caplog, tmp_path):
        """Check warnings when a cwd file is newer than an original_input."""
        bookkeeper = Bookkeeper(cwd=tmp_path)
        ori_inputs = {
            f'{DEFAULT_SUPP}/{ORIGINAL_INPUTS_DIR_NAME}': {'POSCAR': None},
            }
        filesystem_from_dict(ori_inputs, tmp_path)
        time.sleep(0.05)
        (tmp_path/'POSCAR').touch()

        bookkeeper.update_from_cwd()
        # pylint: disable-next=protected-access  # OK in tests
        bookkeeper._mode = Mode.ARCHIVE
        # pylint: disable-next=protected-access  # OK in tests
        bookkeeper._copy_input_files_from_original_inputs_or_cwd()
        # pylint: disable-next=magic-value-comparison
        assert 'is newer' in caplog.text

    def test_copy_from_root(self, caplog, tmp_path):
        """Check warnings when no original_input is present."""
        bookkeeper = Bookkeeper(cwd=tmp_path)
        root_tree = {
            # A log file to ensure timestamps don't shift by 1 sec. See
            # also comments in fixture_funky_files. This is also enough
            # to ensure that there is something to archive.
            f'{LOG_PREFIX}-{MOCK_TIMESTAMP}.log': None,
            'POSCAR': None,
            }
        filesystem_from_dict(root_tree, tmp_path)
        bookkeeper.update_from_cwd()
        target = bookkeeper.history.new_folder.path
        bookkeeper.run('archive')
        assert _FROM_ROOT in caplog.text
        assert (tmp_path/'POSCAR').is_file()
        assert (target/'POSCAR_from_root').is_file()


class TestBookkeeperOthers:
    """Collections of various tests for bits not covered by other tests."""

    # Note about the disable: It is more convenient to modify the key
    # for the glob-all-logs pattern later, as we can then simply edit
    # if we ever change again the naming of log files.
    # pylint: disable-next=dict-init-mutate
    _logs = {
        'tleedm-231110-103910.log': {'r_super': '0.1582',
                                     'run_info': '0 3 31 12'},
        'viperleed-calc-231110-103910.log': {
            'r_super': '0.1582',
            'r_ref': '0.1582 (0.1354 / 0.1827)',
            'run_info': '0 1 11 2 3 31 12',
            },
        }
    _logs['*.log'] = _logs['viperleed-calc-231110-103910.log']

    @parametrize('pattern,expect', _logs.items(), ids=_logs)
    def test_infer_from_log(self, pattern, expect, data_path, tmp_path):
        """Check correct detection of information from a log file."""
        for logfile in (data_path/'bookkeeper').glob(pattern):
            shutil.copy2(logfile, tmp_path)
        bookkeeper = Bookkeeper(cwd=tmp_path)
        bookkeeper.update_from_cwd(silent=True)
        # pylint: disable-next=protected-access   # OK in tests
        logs = bookkeeper._root.logs
        log_info = logs.infer_run_info()
        assert logs.files
        assert log_info == expect

    def test_no_state_files_in_out(self):
        """Check correct behavior when there is no file in OUT to be used."""
        bookkeeper = Bookkeeper()
        # pylint: disable-next=protected-access           # OK in tests
        bookkeeper._root.update_state_files_from_out()    # No OUT dir
        assert not any(Path(f).exists() for f in MOCK_STATE_FILES)

    def test_remove_tensors_deltas(self, tmp_path):
        """Check removal of tensor and delta files."""
        bookkeeper = Bookkeeper(cwd=tmp_path)
        root_tree = {
            DEFAULT_TENSORS: {
                f'{DEFAULT_TENSORS}_003.zip': None,   # This is removed
                f'{DEFAULT_TENSORS}_002.zip': 'contents',    # This not
                f'{DEFAULT_TENSORS}_001.zip': 'contents',    # This not
                f'{DEFAULT_TENSORS}_003': {},  # Unzipped folder, stays
                },
            DEFAULT_DELTAS: {
                f'{DEFAULT_DELTAS}_003.zip': None,    # This is removed
                f'{DEFAULT_DELTAS}_001.zip': 'contents',     # This not
                f'{DEFAULT_DELTAS}_003': {},   # Unzipped folder, stays
                },
            # History information, needed for removal of Tensors
            DEFAULT_HISTORY: {   # Dummy empty files for hashing
                't003.r001_other_stuff': {'file': 'contents'},
                't002.r001_first_run': {'file': 'contents'},
                't002.r005_other_run': {'file': 'contents'},
                't001.r001_first_tensor': {'file': 'contents'},
                },
            }
        original_paths = filesystem_from_dict(root_tree, tmp_path)
        removed_files = {tmp_path/folder/file
                         for folder in (DEFAULT_TENSORS, DEFAULT_DELTAS)
                         for file, contents in root_tree[folder].items()
                         if contents is None}
        bookkeeper.update_from_cwd()
        # pylint: disable-next=protected-access   # OK in tests
        bookkeeper._root.remove_tensors_and_deltas()
        assert not any(f.exists() for f in removed_files)
        assert all(f.exists()
                   for f in original_paths
                   if f not in removed_files)

    @fixture(name='funky_files')
    def fixture_funky_files(self, tmp_path):
        """Prepare a bunch of files/folders that will not be considered."""
        tensor_num_unused = 999
        not_collected_log = tmp_path/'not_a_log.log'
        root_tree = {
            DEFAULT_OUT : {},  # Otherwise 'nothing to do'
            # A log file, just to make sure we don't need to infer
            # a timestamp, as this may happen to be different for
            # the runs we do in test_funky_files, and may lead to
            # unexpected failures.
            f'{LOG_PREFIX}-{MOCK_TIMESTAMP}.log': None,
            # Stuff that should not be copied over:
            # - a directory with the name of a log file
            not_collected_log.name: {},
            # - workhistory subfolders
            DEFAULT_WORK_HISTORY: {
                'some_folder': {},
                f't001.r002_{PREVIOUS_LABEL}': {},
                },
            # - history contents with invalid names
            DEFAULT_HISTORY: {
                f't{tensor_num_unused}.r999_some_file': None,
                'some_stray_directory': {},
                },
            }
        invalid_workhistory = tuple(
            tmp_path/DEFAULT_WORK_HISTORY/entry
            for entry in root_tree[DEFAULT_WORK_HISTORY]
            )
        invalid_history_stuff = tuple(
            tmp_path/DEFAULT_HISTORY/entry
            for entry in root_tree[DEFAULT_HISTORY]
            )
        filesystem_from_dict(root_tree, tmp_path)
        return (tensor_num_unused, not_collected_log,
                invalid_workhistory, invalid_history_stuff)

    def test_funky_files(self, funky_files, tmp_path):
        """Check that funny files and directories are not considered."""
        bookkeeper = Bookkeeper(cwd=tmp_path)
        (tensor_num_unused,
         not_collected_log,
         invalid_workhistory,
         invalid_history_stuff) = funky_files

        bookkeeper.update_from_cwd(silent=True)
        history_dir = bookkeeper.history.new_folder.path
        history_info = bookkeeper.history.info
        # pylint: disable-next=protected-access   # OK in tests
        logs = bookkeeper._root.logs.files
        assert tensor_num_unused not in bookkeeper.max_job_for_tensor
        assert not_collected_log not in logs

        # pylint: disable-next=protected-access           # OK in tests
        assert not any(f in bookkeeper._root._files_to_archive
                       for f in invalid_workhistory)

        bookkeeper.run('archive')
        # Funny workhistory stuff untouched, but '_previous' removed.
        for folder in invalid_workhistory:
            assert (not folder.exists() if PREVIOUS_LABEL in folder.name
                    else folder.is_dir())
        self._check_funky_files_untouched(history_dir,
                                          not_collected_log,
                                          invalid_history_stuff)
        assert history_dir.is_dir()
        assert history_info.path.read_bytes()

        # Now discard should remove workhistory, the archived
        # folder, and its corresponding history.info entry.
        bookkeeper.run('discard_full', requires_user_confirmation=False)
        self._check_funky_files_untouched(history_dir,
                                          not_collected_log,
                                          invalid_history_stuff)
        assert not (tmp_path/DEFAULT_WORK_HISTORY).exists()                     # TODO: this sometimes fails???
        assert not history_dir.exists()
        assert not history_info.path.read_bytes()

    @staticmethod
    def _check_funky_files_untouched(history_dir,
                                     not_collected_log,
                                     invalid_history_stuff):
        """Check that funny files and folders are not modified."""
        assert all(f.exists() for f in invalid_history_stuff)
        assert not_collected_log.exists()
        assert not (history_dir/not_collected_log.name).exists()

    _user_replies = {
        'no reply': ('', False),  # No by default
        'invalid reply, then no': (
            'please do not',
            'NoPe',  # This is the one that is used
            False,
            ),
        'confirmed': ('YES please', True),
        'invalid reply, then yes': (
            'maybe',
            'y',     # This is the one that is used
            True,
            ),
        }

    @parametrize(replies_and_expect=_user_replies.values(), ids=_user_replies)
    def test_user_confirmed(self, replies_and_expect, mock_path, mocker):
        """Check the result of asking user confirmation to proceed."""
        bookkeeper = Bookkeeper(mock_path)
        # pylint: disable-next=protected-access           # OK in tests
        bookkeeper._requires_user_confirmation = True
        *replies, expect = replies_and_expect
        mocker.patch('builtins.input', new=MockInput(*replies))
        # pylint: disable-next=protected-access           # OK in tests
        assert bookkeeper._user_confirmed() == expect


class TestBookkeeperRaises:
    """"Collection of tests for various bookkeeper complaints."""

    @staticmethod
    def _to_method_and_args(method_and_args):
        """Return a method name and arguments from a string."""
        # pylint: disable-next=magic-value-comparison
        if '(' not in method_and_args:
            return method_and_args, tuple()

        method_name, args = method_and_args.split('(')
        args = args.replace(')', '') + ','
        return method_name, ast.literal_eval(args)

    def test_cant_make_history(self):
        """Check complaints when we fail to make the history directory."""
        bookkeeper = Bookkeeper()
        with raises_oserror('pathlib.Path.mkdir'):
            # pylint: disable-next=protected-access    # OK in tests
            bookkeeper._make_history_and_prepare_logger()

        with raises_test_exception('pathlib.Path.mkdir'):
            # pylint: disable-next=protected-access    # OK in tests
            bookkeeper._make_history_and_prepare_logger()

    def test_discard_full_cant_remove_folder(self, after_archive, caplog):
        """Check complaints when it's not possible to remove folders."""
        bookkeeper, *_ = after_archive
        # Purge notes otherwise we can't remove anything.
        info = bookkeeper.history.info.raw_contents
        info, _ = info.rsplit('Notes:', maxsplit=1)
        if info:
            info += 'Notes:\n'
        info = bookkeeper.history.info.path.write_text(info)
        bookkeeper.history.info.read()

        with raises_test_exception('shutil.rmtree'):
            # pylint: disable-next=protected-access    # OK in tests
            bookkeeper._run_discard_full_mode()
        with make_raise_oserror('shutil.rmtree'):
            # pylint: disable-next=protected-access    # OK in tests
            exit_code = bookkeeper._run_discard_full_mode()
            # pylint: disable-next=magic-value-comparison
            assert 'Failed to delete' in caplog.text
            assert exit_code is BookkeeperExitCode.FAIL

    def test_invalid_mode(self):
        """Check complaints when an invalid mode is used."""
        bookkeeper = Bookkeeper()
        with pytest.raises(ValueError):
            bookkeeper.run('invalid')

    def test_no_runner_implemented(self, monkeypatch):
        """Check complaints if a runner is not available for a mode."""
        class _MockMode(Enum):
            INVALID = 'invalid'
        monkeypatch.setattr(
            'viperleed.calc.bookkeeper.bookkeeper.BookkeeperMode',
            _MockMode
            )
        bookkeeper = Bookkeeper()
        with pytest.raises(NotImplementedError):
            bookkeeper.run('invalid')

    raises = object()
    logs = object()
    skips = object()
    _os_error = {
        '_copy_out_and_supp': ('shutil.copytree', logs),
        '_make_and_copy_to_history': ('pathlib.Path.mkdir', raises),
        '_root.read_and_clear_notes_file-read': ('pathlib.Path.read_text',
                                                 logs),
        '_root.read_and_clear_notes_file-write': ('pathlib.Path.write_text',
                                                  logs),
        '_root.logs._read_most_recent': ('pathlib.Path.open', skips),
        '_root.logs.discard': ('pathlib.Path.unlink', logs),
        '_root._remove_out_and_supp': ('shutil.rmtree', logs),
        '_root._remove_ori_files': ('pathlib.Path.unlink', logs),
        '_root._replace_state_files_from_ori': ('pathlib.Path.replace',
                                                raises),
        '_workhistory._discard_previous': ('shutil.rmtree', logs),
        '_workhistory.move_current_and_cleanup(None)': ('shutil.rmtree', logs),
        }

    @parametrize(method_name=_os_error)
    def test_oserror(self, method_name, tmp_path, caplog):
        """Check complaints when functions raise OSError."""
        # Create some files and directories
        root_tree = {
            DEFAULT_WORK_HISTORY: {
                f't001.r002_some_{PREVIOUS_LABEL}_folder': {},
                },
            DEFAULT_OUT: {},
            'notes': None,
            f'{LOG_PREFIX}-20xxxx-xxxxxx.log': None,
            'PARAMETERS_ori': None,
            }
        filesystem_from_dict(root_tree, tmp_path)

        bookkeeper = Bookkeeper(cwd=tmp_path)
        bookkeeper.update_from_cwd(silent=True)

        to_patch, action = self._os_error[method_name]
        method_name, *_ = method_name.split('-')
        method_name, args = self._to_method_and_args(method_name)

        # Notice the use of operator.attrgetter instead of getattr,
        # as the latter does not handle dotted attributes
        method = attrgetter(method_name)(bookkeeper)
        with raises_test_exception(to_patch):
            method(*args)

        if action is self.raises:
            with raises_oserror(to_patch):
                method(*args)
        else:
            with make_raise_oserror(to_patch), not_raises_oserror():
                method(*args)
        # Check logging results, unless we silently skip
        if action is not self.skips:
            assert any(r for r in caplog.records
                       if r.levelno >= logging.WARNING)

    _attr_needs_update = (
        'archiving_required',
        'files_need_archiving',
        'history.info',
        'history.new_folder',
        'history.new_folder.name',
        'history.new_folder.exists',
        'max_job_for_tensor',
        'tensor_number',
        'timestamp',
        '_root.logs',
        )
    _method_needs_update = (
        '_archive_to_history_and_add_info_entry',
        '_copy_input_files_from_original_inputs_or_cwd',
        '_copy_log_files',
        '_copy_out_and_supp',
        '_make_and_copy_to_history',
        'history._find_name_for_new_history_subfolder(None, None)',
        'history.find_new_history_directory(None, None)',
        '_root.logs.discard',
        '_root.revert_to_previous_calc_run',
        '_root.remove_tensors_and_deltas',
        '_run_archive_mode',
        '_run_clear_mode',
        '_run_discard_full_mode',
        '_run_discard_mode',
        )

    @parametrize(attr=_attr_needs_update)
    def test_too_early_attribute_access(self, attr):
        """Check that accessing attributes before update_from_cwd fails."""
        bookkeeper = Bookkeeper()
        with check_too_early():
            attrgetter(attr)(bookkeeper)

    @parametrize(method_name=_method_needs_update)
    def test_too_early_method_call(self, method_name):
        """Check that accessing attributes before update_from_cwd fails."""
        bookkeeper = Bookkeeper()
        method_name, args = self._to_method_and_args(method_name)
        with check_too_early():
            # Some methods raise already at getattr, some at call
            method = attrgetter(method_name)(bookkeeper)
            method(*args)