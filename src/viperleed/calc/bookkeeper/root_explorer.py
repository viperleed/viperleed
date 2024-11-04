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
from functools import partial
from operator import attrgetter
from pathlib import Path
import shutil
import re

from viperleed.calc.constants import DEFAULT_OUT
from viperleed.calc.constants import DEFAULT_SUPP
from viperleed.calc.constants import DEFAULT_WORK_HISTORY

from .constants import CALC_LOG_PREFIXES
from .constants import STATE_FILES
from .history.workhistory import WorkhistoryHandler
from .log import LOGGER
from .utils import discard_files
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
        # Tuples of paths to log files found in self.path
        self._logs = None              # LogFiles, set in collect
        self._files_to_archive = None  # See _collect_files_to_archive
        self.workhistory = WorkhistoryHandler(
            work_history_path=self.path / DEFAULT_WORK_HISTORY,
            bookkeeper=bookkeeper,
            )

    # Simple read-only properties
    path = make_property('_path')
    logs = make_property('_logs', needs_update=True, updater='collect_info')

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
        self.collect_info()

    def collect_info(self):
        """Collect information from the root directory."""
        self._logs = LogFiles(self.path)
        self._logs.collect()
        self._collect_files_to_archive()

    def infer_run_info(self):
        """Return a dictionary of information read from the newest calc log."""
        return self.logs.infer_run_info()

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
        self.collect_info()

    def update_state_files_from_out(self):
        """Move state files from OUT to root. Rename old to '_ori'."""
        out_path = self.path / DEFAULT_OUT
        if not out_path.is_dir():
            return
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
            cwd_file.replace(self.path / f'{file}_ori')
            # The OUT file is copied, however: .replace would move it
            shutil.copy2(out_file, cwd_file)

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
        for file in STATE_FILES:
            ori_file = self.path / f'{file}_ori'
            try:
                ori_file.replace(self.path / file)
            except FileNotFoundError:
                pass
            except OSError:
                LOGGER.error(f'Failed to move {ori_file} to {file}.')
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
