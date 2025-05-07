"""Module log of viperleed.calc.bookkeeper."""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2020-01-30'
__license__ = 'GPLv3+'

from collections import namedtuple
from functools import partial
import logging
from operator import attrgetter
import re

from viperleed.calc.bookkeeper.constants import CALC_LOG_PREFIXES
from viperleed.calc.bookkeeper.utils import discard_files
from viperleed.calc.bookkeeper.utils import make_property
from viperleed.calc.bookkeeper.utils import needs_update_for_attr
from viperleed.calc.lib.log_utils import get_handlers
from viperleed.calc.lib.string_utils import parent_name
from viperleed.calc.lib.version import Version


BOOKIE_LOGFILE = 'bookkeeper.log'  # Persistent among runs
LOGGER = logging.getLogger(parent_name(__name__))

# Regular expressions for parsing the log file
_CALC_VERSION_RE = re.compile(r'This is ViPErLEED version (?P<version>[\d.]+)')
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

# Container of information about one log file. Used by LogFiles below
LogInfo = namedtuple('LogInfo', ('timestamp', 'lines'))


def add_bookkeeper_logfile(at_path):
    """Attach a FileHandler to LOGGER for the bookkeeper log file."""
    bookkeeper_log = at_path / BOOKIE_LOGFILE
    has_log = get_handlers(LOGGER,
                           logging.FileHandler,
                           baseFilename=str(bookkeeper_log))
    if not any(has_log):
        LOGGER.addHandler(logging.FileHandler(bookkeeper_log, mode='a'))


def ensure_has_stream_handler():
    """Attach a stream handler to LOGGER unless one exists already."""
    if not any(get_handlers(LOGGER, logging.StreamHandler)):
        LOGGER.addHandler(logging.StreamHandler())
    LOGGER.setLevel(logging.INFO)
    LOGGER.propagate = True


def remove_bookkeeper_logfile(at_path):
    """Remove the FileHandler(s) `at_path` for LOGGER, if any."""
    handlers = get_handlers(LOGGER,
                            logging.FileHandler,
                            baseFilename=str(at_path/BOOKIE_LOGFILE))
    for handler in handlers:
        handler.flush()
        LOGGER.removeHandler(handler)
        handler.close()


class LogFiles:
    """Container to manage log files found in a folder."""

    _needs_collect = partial(needs_update_for_attr, updater='collect')

    def __init__(self, path):
        """Initialize an instance with a folder's `path`."""
        self._path = path
        self._calc = None    # Log files for viperleed.calc
        self._others = None  # Other log files found at path
        self._calc_version = None  # Set in _infer_calc_version
        self.most_recent = None    # LogInfo from most recent calc log

    calc = make_property('_calc', needs_update=True)
    version = make_property('_calc_version')

    @property
    @_needs_collect('_calc')
    def files(self):
        """Return paths to all the log files in the root directory."""
        return self._calc + self._others

    def collect(self):
        """Collect and store internally information about log files."""
        self._collect_logs()
        self._read_most_recent()
        self._infer_calc_version()

    def discard(self):
        """Delete all log files at self._path.

        Returns
        -------
        bool
            Whether any log file was deleted.
        """
        logs_discarded = discard_files(*self.files)
        self.collect()
        return logs_discarded

    # NB: the next one needs _calc, not most_recent, as the latter may
    # be None also if there is no calc log file in the root directory.
    @_needs_collect('_calc')
    def infer_run_info(self):
        """Return a dictionary of information read from the newest calc log."""
        try:
            log_lines = self.most_recent.lines
        except AttributeError:  # No log in root
            assert self.most_recent is None
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
    def _infer_calc_version(self):
        """Find out the version with which the calc log was created."""
        if not self.most_recent:
            return
        for line in self.most_recent.lines:
            match_ = _CALC_VERSION_RE.match(line)
            if not match_:
                continue
            try:
                self._calc_version = Version(match_['version'])
            except ValueError:  # Malformed version
                continue
            return

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
