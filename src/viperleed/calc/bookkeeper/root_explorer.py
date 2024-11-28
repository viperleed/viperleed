"""Module root_explorer of viperleed.calc.bookkeeper.

Defines objects useful for identifying, reading, and modifying files
and folders present in the directory in which Bookkeeper runs.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-10-14'
__license__ = 'GPLv3+'

from collections import namedtuple
from contextlib import nullcontext
from datetime import datetime
from functools import partial
from operator import attrgetter
from pathlib import Path
import shutil
import re

from viperleed.calc.constants import DEFAULT_DELTAS
from viperleed.calc.constants import DEFAULT_OUT
from viperleed.calc.constants import DEFAULT_SUPP
from viperleed.calc.constants import DEFAULT_TENSORS
from viperleed.calc.constants import ORIGINAL_INPUTS_DIR_NAME
from viperleed.calc.lib.leedbase import getMaxTensorIndex
from viperleed.calc.lib.log_utils import logging_silent
from viperleed.calc.lib.time_utils import DateTimeFormat
from viperleed.calc.sections.calc_section import ALL_INPUT_FILES

from .constants import CALC_LOG_PREFIXES
from .constants import EDITED_SUFFIX
from .constants import STATE_FILES
from .history.explorer import HistoryExplorer
from .history.workhistory import WorkhistoryHandler
from .log import LOGGER
from .utils import discard_files
from .utils import file_contents_identical
from .utils import make_property
from .utils import needs_update_for_attr


# Regular expressions for parsing the log file
_RFAC_RE = r'[\d.( )/]+'
_LOG_FILE_RE = {
    'run_info': re.compile(r'Executed segments:\s*(?P<run_info>[\d ]+)\s*'),
    'r_ref': re.compile(
        rf'Final R \(refcalc\)\s*:\s*(?P<r_ref>{_RFAC_RE})\s*'
        ),
    'r_super': re.compile(
        rf'Final R \(superpos\)\s*:\s*(?P<r_super>{_RFAC_RE})\s*'
        ),
    }


# TODO: read_and_clear_notes_file  -- support multiple files
class RootExplorer:
    """A class for handling files and folders in the bookkeeper's CWD."""

    _needs_collect = partial(needs_update_for_attr, updater='collect_info')

    def __init__(self, path, bookkeeper):
        """Initialize this explorer at `path` for a `bookkeeper`."""
        self._path = Path(path).resolve()
        self._logs = None              # LogFiles, set in collect
        self._files_to_archive = None  # See _collect_files_to_archive
        self.tensors = TensorAndDeltaInfo(self.path)
        self.history = HistoryExplorer(self.path)
        self.workhistory = WorkhistoryHandler(root=self.path,
                                              bookkeeper=bookkeeper)

    # Simple read-only properties
    logs = make_property('_logs', needs_update=True, updater='collect_info')
    path = make_property('_path')

    @property
    @_needs_collect('_logs')
    def calc_timestamp(self):
        """Return the timestamp of the most recent calc run, or None."""
        try:
            return self.logs.most_recent.timestamp
        except AttributeError:
            return None

    @property
    @_needs_collect('_files_to_archive')
    def needs_archiving(self):
        """Return whether there are files/folder that should go to history."""
        return any(self._files_to_archive)

    @property
    def orig_inputs_dir(self):
        """Return the path to the folder containing untouched input files."""
        return self.path / DEFAULT_SUPP / ORIGINAL_INPUTS_DIR_NAME

    def clear_for_next_calc_run(self):
        """Clean up the root directory for a clean calc execution.

        Notice that workhistory files and the workhistory directory
        itself are left untouched.

        Returns
        -------
        None.
        """
        self.logs.discard()
        self._remove_out_and_supp()
        self._remove_ori_files()
        # Let's not complain again about funny stuff. We have already
        # done so when collect_info was called the first time. Notice
        # that this must have happened, as this method relies on the
        # logs attribute to be up to date.
        self.collect_info(silent=True)

    def collect_info(self, silent=False):
        """Collect information from the root directory."""
        log_context = logging_silent if silent else nullcontext
        with log_context():
            self._logs = LogFiles(self.path)
            self._logs.collect()
            self.tensors.collect()
            self._collect_files_to_archive()
            self.history.collect_subfolders()

    def infer_run_info(self):
        """Return a dictionary of information read from the newest calc log."""
        return self.logs.infer_run_info()

    def list_files_to_replace(self):
        """Return a tuple of pairs of files that will be replaced.

        Returns
        -------
        to_replace : tuple of tuples
            Pairs of files that will be replaced. Each pair has format
            (<file_that_will_be_deleted>, <file_that will_be_renamed>).
            All items are Paths to the files.
        """
        ori_files = {self.path/file: self.path/f'{file}_ori'
                     for file in STATE_FILES}
        ori_files = {file: ori_file for file, ori_file in ori_files.items()
                     if ori_file.exists()}
        return tuple(ori_files.items())

    def list_paths_to_discard(self):
        """Return a list of files and folders that will be discarded."""
        to_discard = (
            self.path / DEFAULT_OUT,
            self.path / DEFAULT_SUPP,
            *self.logs.files,
            *self.tensors.list_paths_to_discard(self.history)
            )
        return tuple(p for p in to_discard if p.exists())

    def mark_edited_files(self):
        """Mark as _edited the files that were modified since calc started."""
        try:
            calc_started = datetime.strptime(self.calc_timestamp,
                                             DateTimeFormat.FILE_SUFFIX.value)
        except TypeError:  # No calc log in root
            assert self.calc_timestamp is None
            return

        ori_path = self.orig_inputs_dir
        for filename in ALL_INPUT_FILES:
            cwd_file = self.path / filename
            try:
                mod_time = datetime.utcfromtimestamp(cwd_file.stat().st_mtime)
            except FileNotFoundError:
                continue
            if mod_time < calc_started:
                continue
            # The input file was modified. See whether its contents
            # differ from the one we have copied in original_inputs.
            # If they differ (or there is no file in original_inputs),
            # label this file as having been user-edited.
            if not file_contents_identical(cwd_file, ori_path / filename):
                LOGGER.info(f"Renaming {filename} to _edited: mod={mod_time}, calc={calc_started}")
                cwd_file.replace(self.path / (filename + EDITED_SUFFIX))

    def read_and_clear_notes_file(self):
        """Return notes read from file. Clear the file contents."""
        notes_path = next(self.path.glob('notes*'), None)
        if notes_path is None:
            return ''
        try:
            notes = notes_path.read_text(encoding='utf-8')
        except OSError:
            LOGGER.error(f'Error: Failed to read {notes_path.name} file.',
                         exc_info=True)
            return ''
        try:
            notes_path.write_text('', encoding='utf-8')
        except OSError:
            LOGGER.error(f'Failed to clear the {notes_path.name} '
                         'file after reading.', exc_info=True)
        return notes

    def remove_tensors_and_deltas(self):
        """Delete relevant Tensors and Deltas archives."""
        self.tensors.discard(self.history)

    def revert_to_previous_calc_run(self):
        """Remove files generated by the previous calc run in self.path.

        Notice that files in the workhistory directory, as
        well as the workhistory directory are not touched.

        Raises
        ------
        OSError
            If replacing input files with their '_ori'-suffixed
            version fails.
        """
        self.logs.discard()
        self._remove_out_and_supp()
        self._replace_state_files_from_ori()
        # Let's not complain again about funny stuff. We have already
        # done so when collect_info was called the first time. Notice
        # that this must have happened, as this method relies on the
        # logs attribute to be up to date.
        self.collect_info(silent=True)

    def update_state_files_from_out(self):
        """Copy state files from OUT to root. Rename old to '_ori'."""
        out_path = self.path / DEFAULT_OUT
        if not out_path.is_dir():
            return
        failed = {}
        for file in STATE_FILES:
            out_file = out_path / f'{file}_OUT'
            cwd_file = self.path / file
            if not out_file.is_file():
                continue
            # NB: all the moving around is local to the same folder
            # tree, so we don't really need shutil. Path.replace always
            # works for the same file system. Notice also the use of
            # replace rather than rename, as the behavior of rename
            # is not identical for all platforms.
            try:
                cwd_file.replace(self.path / f'{file}_ori')
            except OSError as exc:
                failed[file] = exc
                continue
            # The OUT file is copied, however: .replace would move it
            shutil.copy2(out_file, cwd_file)
        if not failed:
            return
        # Report failures, so people can do it manually.
        LOGGER.error('Failed to copy OUT files to the current directory '
                     'and to rename the original inputs there.')
        LOGGER.error('The following files could not be processed:')
        for file in failed:
            LOGGER.error(f'{file} -x-> {file}_ori. '
                         f'{out_path.name}/{file}_OUT -x-> {file}')
        exc_msgs = (f'{file}: raised {type(exc).__name__} - {exc}'
                    for file, exc in failed.items())
        raise OSError('\n'.join(exc_msgs))

    @_needs_collect('_logs')
    def _collect_files_to_archive(self):
        """Return a tuple of files/folders to be archived in history."""
        to_archive = (
             # OUT and SUPP, if present
            self.path / DEFAULT_OUT,
            self.path / DEFAULT_SUPP,
            # Any calc logs
            *self.logs.calc,
            # And workhistory folders
            *self.workhistory.find_current_directories(contains='r'),
            )
        self._files_to_archive = tuple(p for p in to_archive
                                       if p.is_file() or p.is_dir())

    def _remove_ori_files(self):
        """Delete '_ori'-suffixed files from root."""
        ori_files = (self.path / f'{file}_ori' for file in STATE_FILES)
        discard_files(*ori_files)

    def _remove_out_and_supp(self):
        """Delete the SUPP and OUT directories from root."""
        discard_files(self.path / DEFAULT_OUT, self.path / DEFAULT_SUPP)

    def _replace_state_files_from_ori(self):
        """Replace input files with their '_ori'-suffixed version."""
        for state_file, ori_file in self.list_files_to_replace():
            try:
                ori_file.replace(state_file)
            except OSError:
                LOGGER.error(f'Failed to rename {ori_file.name} '
                             f'to {state_file.name}.')
                raise


LogInfo = namedtuple('LogInfo', ('timestamp', 'lines'))


class LogFiles:
    """Container to manage log files found in the bookkeeper's root."""

    _needs_collect = partial(needs_update_for_attr, updater='collect')

    def __init__(self, path):
        """Initialize an instance at path."""
        self._path = path
        self._calc = None    # Log files for viperleed.calc
        self._others = None  # Other log files found at path
        self.most_recent = None    # LogInfo from most recent calc log

    calc = make_property('_calc', needs_update=True)

    @property
    @_needs_collect('_calc')
    def files(self):
        """Return paths to all the log files in the root directory."""
        # tuple() is the start value. It would be nicer to specify
        # it as a keyword, but this was only introduced in py38.
        return self._calc + self._others

    def collect(self):
        """Collect and store internally information about log files."""
        self._collect_logs()
        self._read_most_recent()

    def discard(self):
        """Delete all log files at self._path."""
        discard_files(*self.files)
        self.collect()

    # NB: the next one needs _calc, not most_recent, as the latter may
    # be None also if there is no calc log file in the root directory.
    @_needs_collect('_calc')
    def infer_run_info(self):
        """Return a dictionary of information read from the newest calc log."""
        try:
            log_lines = self.most_recent.lines
        except AttributeError:  # No log in root
            log_lines = ()
        matched = {k: False for k in _LOG_FILE_RE}
        for line in reversed(log_lines):  # Info is at the end
            for info, already_matched in matched.items():
                if already_matched:
                    continue
                matched[info] = _LOG_FILE_RE[info].match(line)
            if all(matched.values()):
                break
        return {k: match[k] for k, match in matched.items() if match}

    def _collect_logs(self):
        """Find all the log files in the root folder. Store them internally."""
        calc_logs, other_logs = [], []
        for file in self._path.glob('*.log'):
            if not file.is_file():
                continue
            container = (calc_logs if file.name.startswith(CALC_LOG_PREFIXES)
                         else other_logs)
            container.append(file)
        self._calc = tuple(calc_logs)
        self._others = tuple(other_logs)

    @_needs_collect('_calc')
    def _read_most_recent(self):
        """Read information from the most recent log file, if available."""
        split_logs = {}  # Path to most recent log for each prefix
        for prefix in CALC_LOG_PREFIXES:  # newest to oldest
            calc_logs = (f for f in self._calc if f.name.startswith(prefix))
            try:
                split_logs[prefix] = max(calc_logs, key=attrgetter('name'))
            except ValueError:
                pass
        try:
            most_recent_log = next(iter(split_logs.values()))
        except StopIteration:  # No log files
            return

        timestamp = most_recent_log.name[-17:-4]
        last_log_lines = ()
        try:  # pylint: disable=too-many-try-statements
            with most_recent_log.open('r', encoding='utf-8') as log_file:
                last_log_lines = tuple(log_file.readlines())
        except OSError:
            pass
        self.most_recent = LogInfo(timestamp, last_log_lines)


class TensorAndDeltaInfo:
    """A container of information about zip Tensor (and Delta) files."""

    def __init__(self, root):
        """Initialize an instance that handles a `root` path."""
        self._root = root
        self._most_recent = None  # Index of most recent Tensor file

    most_recent = make_property('_most_recent',
                                needs_update=True,
                                updater='collect')

    def collect(self):
        """Collect Tensor/Deltas information from self.path."""
        self._most_recent = getMaxTensorIndex(home=self._root, zip_only=True)

    def discard(self, history):
        """Delete the most recent Tensor and Delta zip files."""
        discard_files(*self.list_paths_to_discard(history))

    def list_paths_to_discard(self, history):
        """Return a tuple of paths to Tensor/Delta files to be discarded."""
        # The tensors that may be removed are those of
        # the last run archived in the history folder
        tensor_nums = (f.tensor_num for f in history.last_folder_and_siblings)

        # However, we should only discard those that have only one run,
        # otherwise we may be loosing info from previous calculations.
        max_run = history.max_run_per_tensor
        tensor_nums = {t for t in tensor_nums if max_run[t] <= 1}

        tensor_fmt = f'{DEFAULT_TENSORS}/{DEFAULT_TENSORS}_{{ind:03d}}.zip'
        delta_fmt = f'{DEFAULT_DELTAS}/{DEFAULT_DELTAS}_{{ind:03d}}.zip'
        to_discard = [tensor_fmt.format(ind=t) for t in tensor_nums]
        to_discard.extend(delta_fmt.format(ind=t) for t in tensor_nums)

        # Make them into paths, keeps those that exist, sort by name
        to_discard_paths = (self._root / f for f in to_discard)
        to_discard_paths = (f for f in to_discard_paths if f.is_file())
        return tuple(sorted(to_discard_paths, key=attrgetter('name')))
