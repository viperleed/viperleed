"""ViPErLEED bookkeeper module of package calc."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2020-01-30'
__license__ = 'GPLv3+'

from collections import defaultdict
from contextlib import nullcontext
from enum import IntEnum
from functools import wraps
import logging
from operator import attrgetter
from pathlib import Path
import re
import shutil
import time

from viperleed.calc import DEFAULT_HISTORY
from viperleed.calc import DEFAULT_WORK_HISTORY
from viperleed.calc import LOG_PREFIX
from viperleed.calc import ORIGINAL_INPUTS_DIR_NAME
from viperleed.calc.lib.base import logging_silent
from viperleed.calc.lib.leedbase import getMaxTensorIndex
from viperleed.calc.sections.calc_section import ALL_INPUT_FILES
from viperleed.calc.sections.cleanup import DEFAULT_OUT
from viperleed.calc.sections.cleanup import DEFAULT_SUPP
from viperleed.calc.sections.cleanup import PREVIOUS_LABEL

from .constants import HISTORY_INFO_NAME
from .constants import LOGGER
from .history import HistoryInfoEntry
from .history import HistoryInfoFile
from .history import NoHistoryEntryError
from .history import PureCommentEntry
from .mode import BookkeeperMode


BOOKIE_LOGFILE = 'bookkeeper.log'  # Persistent among runs
CALC_LOG_PREFIXES = (
    LOG_PREFIX,
    # The next ones are for backwards compatibility. Must
    # be in historical order of use: Most recent first!
    'tleedm',
    )
HIST_FOLDER_RE = re.compile(
    r't(?P<tensor_num>[0-9]{3}).r(?P<job_num>[0-9]{3})_'
    )
STATE_FILES = ('PARAMETERS', 'POSCAR', 'VIBROCC')

# Input files that may be generated at runtime - don't warn if missing
RUNTIME_GENERATED_INPUT_FILES = ('IVBEAMS', 'PHASESHIFTS')

# Regular expressions for parsing the log file
_LOG_FILE_RE = {
    'run_info': re.compile(r'Executed segments:\s*(<run_info>\w+)\s*'),
    'r_ref': re.compile(r'Final R (refcalc)\s*:\s*(<r_ref>\w+)\s*'),
    'r_super': re.compile(r'Final R (superpos)\s*:\s*(<r_super>\w+)\s*'),
    }


def _get_attr_or_dict_item(attr):
    """Return an attribute getter for 'attr_name' or 'dict_name[item_name]'.

    Parameters
    ----------
    attr : str
        The attribute of self to be looked up. If the value is None,
        the decorated function raises AttributeError suggesting to
        call update_from_cwd beforehand. `attr` may also have the form
        'dict_name[dict_item]'. In this case, self.dict_name[dict_item]
        is looked up instead.

    Returns
    -------
    getter : callable
        A function that returns `attr` from its first argument.
    cleaned_attr : str
        'attr_name' or 'item_name', depending on the form of `attr`.
    """
    if '[' in attr:  # pylint: disable=magic-value-comparison
        dict_, attr = attr.split('[')
        attr = attr.replace(']', '')
        def _getattr(self):
            container = getattr(self, dict_)
            return container[attr]
    else:
        def _getattr(self):
            return getattr(self, attr)
    return _getattr, attr


def _make_property(attr, needs_update=False):
    """Return a property getter that looks up `attr`."""
    getter, property_name = _get_attr_or_dict_item(attr)
    if needs_update:
        getter = _needs_update_for_attr(attr, attr_name=property_name)(getter)
    return property(getter)


def _needs_update_for_attr(attr, attr_name=None):
    """Return a decorator that complains if `attr` is None.

    Parameters
    ----------
    attr : str
        The attribute of self to be looked up. If the value is None,
        the decorated function raises AttributeError suggesting to
        call update_from_cwd beforehand. `attr` may also have the form
        'dict_name[dict_item]'. In this case, self.dict_name[dict_item]
        is looked up instead.
    attr_name : str, optional
        The name of the decorated function. Automatically fetched from
        the decorated function if not given. Default is None.

    Returns
    -------
    _decorator : callable
        A decorator to apply to a function that raises AttributeError
        before calling the function if `attr` is None.
    """
    _getattr, *_ = _get_attr_or_dict_item(attr)
    def _decorator(func):
        func_name = attr_name or func.__name__
        @wraps(func)
        def _wrapper(self, *args, **kwargs):
            if _getattr(self) is None:
                raise AttributeError(
                    f'{type(self).__name__} has no {func_name} yet. '
                    'Call update_from_cwd() beforehand.'
                    )
            return func(self, *args, **kwargs)
        return _wrapper
    return _decorator


class _FileNotOlderError(Exception):
    """Exception used internally for file age checks."""


class BookeeperExitCode(IntEnum):
    """Exit code of the bookkeeper."""
    SUCCESS = 0
    NOTHING_TO_DO = -1
    FAIL = 1


# TODO: catch errors from history.info, at .read()
class Bookkeeper:
    """Bookkeeper to archive or discard the most recent viperleed calc run."""

    def __init__(self,
                 job_name=None,
                 history_name=DEFAULT_HISTORY,
                 work_history_name=DEFAULT_WORK_HISTORY,
                 cwd=Path.cwd()):
        """Initialize the bookkeeper.

        Parameters
        ----------
        job_name : str or None, optional
            Custom name to append to the stored folder and to
            history.info. If not given or None, no extra name
            is added. Default is None.
        history_name : str, optional
            The name of the folder in the current directory where the
            most recent run should be archived. Default is 'history'.
        work_history_name : str, optional
            The name of the workhistory subfolder of the current
            directory where intermediate runs may have been stored.
            Results are also copied from here to ./history_name.
            Default is 'workhistory'.
        cwd : Path, optional
            The current working directory. Default is the current
            directory.
        """
        self._folder_names = {
            'top_level_history': history_name,
            'workhistory': work_history_name,
            'history_dir_base_name': None,   # Set in update_from_cwd
            'history_dir': None,             # Set in update_from_cwd
            }
        self._history_info = None   # Set in update_from_cwd
        self._paths = {
            'cwd': Path(cwd),
            # The next ones are set in update_from_cwd
            'calc_logs': None,      # tuple
            'other_logs': None,     # tuple
            'to_be_archived': None  # tuple
            }
        self._state_info = {
            'job_name': job_name,
            'logger_prepared': False,
            # The next ones are set in update_from_cwd
            'last_log_lines': None,       # tuple
            'max_job_for_tensor': None,   # defaultdict[int]
            'tensor_number': None,        # int
            'timestamp': None,            # str
            'history_with_same_base_name_exists': None,  # bool
            }

    # Simple dynamic @properties
    cwd = _make_property('_paths[cwd]')
    history_dir_base_name = _make_property(
        # Notice that this is not necessarily the name of the
        # history subfolder that is created. Use .history_dir
        # for that.
        '_folder_names[history_dir_base_name]',
        needs_update=True,
        )
    history_with_same_base_name_exists = _make_property(
        '_state_info[history_with_same_base_name_exists]',
        needs_update=True,
        )
    job_name = _make_property('_state_info[job_name]')
    max_job_for_tensor = _make_property('_state_info[max_job_for_tensor]',
                                        needs_update=True)
    tensor_number = _make_property('_state_info[tensor_number]',
                                   needs_update=True)
    timestamp = _make_property('_state_info[timestamp]',
                               needs_update=True)

    @property
    def all_cwd_logs(self):
        """Return paths to all the log files in the root directory."""
        return sum(self.cwd_logs, start=tuple())

    @property
    def archiving_required(self):
        """Check if archiving is required."""
        return (self.files_need_archiving
                and not self.history_with_same_base_name_exists)

    @property
    @_needs_update_for_attr('_paths[calc_logs]')
    def cwd_logs(self):
        """Return log files, split into those of calc and all the others."""
        return self._paths['calc_logs'], self._paths['other_logs']

    @property
    @_needs_update_for_attr('_paths[to_be_archived]')
    def files_need_archiving(self):
        """Check if there are any files that need archiving."""
        return bool(self._paths['to_be_archived'])

    @property
    @_needs_update_for_attr('_folder_names[history_dir]')
    def history_dir(self):
        """Return the path to the history subfolder that we deal with.

        This is the directory to which files will be moved upon run().

        Returns
        -------
        Path
        """
        return self.top_level_history_path / self._folder_names['history_dir']

    @property
    @_needs_update_for_attr('_history_info')
    def history_info(self):
        """Return the HistoryInfoFile handling history.info."""
        return self._history_info

    @property
    def orig_inputs_dir(self):
        """Return the path to the folder containing untouched input files."""
        return self.cwd / DEFAULT_SUPP / ORIGINAL_INPUTS_DIR_NAME

    @property
    def top_level_history_path(self):
        """Return the path to the 'history' folder."""
        return self.cwd / self._folder_names['top_level_history']

    @property
    def work_history_path(self):
        """Return the path to the 'workhistory' folder."""
        return self.cwd / self._folder_names['workhistory']

    def run(self, mode):
        """Run the bookkeeper in the given mode.

        Parameters
        ----------
        mode : str or BookkeeperMode
            Which bookkeeper mode to use. See help(BookkeeperMode).

        Returns
        -------
        exit_code : BookeeperExitCode

        Raises
        ------
        ValueError
            If `mode` is not a valid bookkeeper mode.
        NotImplementedError
            If `mode` is a valid bookkeeper mode, but
            there is no method named _run_<mode>_mode.
        OSError
            If creation of the history folder or any of the subfolders
            where results are to be stored fails.
        """
        try:
            mode = BookkeeperMode(mode)
        except ValueError as exc:                                               # TODO: untested raise
            raise ValueError(f'Unknown mode {mode}') from exc

        try:
            runner = getattr(self, f'_run_{mode.name.lower()}_mode')
        except AttributeError as exc:                                           # TODO: untested raise
            raise NotImplementedError from exc

        LOGGER.info(f'Running bookkeeper in {mode.name} mode.')
        self.update_from_cwd()
        try:
            return runner()
        finally:
            # Clean up all the state attributes so we
            # don't risk giving the wrong information.
            self._clean_state()

    def update_from_cwd(self, silent=False):
        """Update timestamp, tensor number, log lines, etc. from root.

        This method is always called automatically before every run
        of the bookkeeper.

        Parameters
        ----------
        silent : bool, optional
            Whether logging messages should be silenced or not.
            The default is False.

        Returns
        -------
        None.
        """
        context = logging_silent if silent else nullcontext
        with context():
            self._make_history_and_prepare_logger()

            # Figure out the number of the tensor (it's the most recent
            # one) and the highest run number currently stored for each
            # tensor in history.
            tensor_number = getMaxTensorIndex(home=self.cwd, zip_only=True)
            self._state_info['tensor_number'] = tensor_number
            self._collect_max_run_per_tensor()

            # Infer timestamp from log file, if possible
            self._collect_log_files()
            timestamp, log_lines = self._read_most_recent_log()
            self._state_info['timestamp'] = timestamp
            self._state_info['last_log_lines'] = log_lines

            # Infer which history subfolder we should be handling
            self._find_base_name_for_history_subfolder()
            self._find_new_history_directory_name()

            # history.info handler (may create history.info file)
            self._history_info = HistoryInfoFile(self.cwd / HISTORY_INFO_NAME,
                                                 create_new=True)
            self.history_info.read()

            # Collect all the necessary files
            self._collect_files_to_archive()

    def _clean_state(self):
        """Rebuild (almost) the same state as after __init__.

        The only surviving attributes are the ones given
        at __init__ and the state that concerns logging.

        Returns
        -------
        None.
        """
        kwargs = {
            'job_name': self.job_name,
            'history_name': self._folder_names['top_level_history'],
            'work_history_name': self._folder_names['workhistory'],
            'cwd': self.cwd
            }
        logger_prepared = self._state_info['logger_prepared']

        # Note on the disable: while it is indeed not very elegant to
        # call __init__ again, it is the simplest way to (i) make it
        # clear in __init__ which attributes we have, and (ii) avoid
        # code repetition here. The alternative would be to set all
        # attributes and keys to None in __init__, and call another
        # method both in __init__ and here. This seems a lot of
        # complication just to reset (almost) everything to None.
        # pylint: disable-next=unnecessary-dunder-call
        self.__init__(**kwargs)
        self._state_info['logger_prepared'] = logger_prepared

    @_needs_update_for_attr('_paths[calc_logs]')
    def _collect_files_to_archive(self):
        """Scan the root directory for files to be stored to history."""
        # Check for OUT and SUPP
        files_to_archive = [self.cwd / DEFAULT_OUT, self.cwd / DEFAULT_SUPP]
        # Any calc logs
        files_to_archive.extend(self._paths['calc_logs'])
        # Workhistory folders
        if self.work_history_path.is_dir():                                     # TODO: untested
            files_to_archive.extend(
                self._get_current_workhistory_directories(contains='r')
                )
        self._paths['to_be_archived'] = tuple(p for p in files_to_archive
                                              if p.is_file() or p.is_dir())

    def _collect_log_files(self):
        """Find the path to log files in root, store them internally."""
        calc_logs, other_logs = [], []
        for file in self.cwd.glob('*.log'):
            if not file.is_file():
                continue
            container = (calc_logs if file.name.startswith(CALC_LOG_PREFIXES)
                         else other_logs)
            container.append(file)
        self._paths['calc_logs'] = tuple(calc_logs)
        self._paths['other_logs'] = tuple(other_logs)

    def _collect_max_run_per_tensor(self):
        """Find the maximum run numbers for all directories in history."""
        history_path = self.top_level_history_path
        max_jobs = defaultdict(int)  # max. job number per tensor number
        for directory in history_path.iterdir():
            match = HIST_FOLDER_RE.match(directory.name)
            if not directory.is_dir() or not match:                             # TODO: untested
                continue
            tensor_num = int(match['tensor_num'])
            job_num = int(match['job_num'])
            max_jobs[tensor_num] = max(max_jobs[tensor_num], job_num)
        self._state_info['max_job_for_tensor'] = max_jobs

    def _copy_input_files_from_original_inputs_or_cwd(self, use_ori=False):
        """Copy input files to the history subfolder.

        The files are preferentially collected from the original_inputs
        subfolder of SUPP, only if they are not found, they are taken
        from the root directory (unless they are auto-generated ones).
        If the root directory contains a newer version of the file
        than in original_inputs, warn the user (but still copy the one
        from original_inputs).

        Parameters
        ----------
        use_ori : bool, optional
            Whether '_ori'-suffixed files should be copied rather
            than 'bare' ones. Default is False.

        Returns
        -------
        None.
        """
        for file in ALL_INPUT_FILES:
            original_file = self.orig_inputs_dir / file
            cwd_file = self.cwd / file
            if file in STATE_FILES and use_ori:
                cwd_file = self.cwd / f'{file}_ori'

            copy_file, with_name = None, None                                   # TODO: do we want to maintain the _ori suffix?? If not we can replace the second None with file.
            if original_file.is_file() and cwd_file.is_file():
                # Copy original, but warn if cwd is newer
                copy_file = original_file
                try:
                    _check_newer(older=cwd_file, newer=original_file)
                except _FileNotOlderError:
                    LOGGER.warning(
                        f'File {file} from {ORIGINAL_INPUTS_DIR_NAME} was '
                        'copied to history, but the file in the input '
                        'directory is newer.'
                        )
            elif original_file.is_file():  # Just copy original
                copy_file = original_file
            elif file in RUNTIME_GENERATED_INPUT_FILES:
                # Ignore optional inputs in root. If the user gave
                # these explicitly, they have already been copied
                # by calc to original_inputs
                continue
            elif cwd_file.is_file():  # Copy cwd and warn                       # TODO: untested
                copy_file, with_name = cwd_file, f'{cwd_file.name}_from_root'
                LOGGER.warning(
                    f'File {file} not found in {ORIGINAL_INPUTS_DIR_NAME}. '
                    'Using file from root directory instead and renaming to '
                    f'{with_name}.'
                    )
            if copy_file:
                self._copy_one_file_to_history(copy_file, with_name=with_name)

    def _copy_log_files(self):
        """Copy all the log files in root to history."""
        for file in self.all_cwd_logs:
            self._copy_one_file_to_history(file)

    def _copy_one_file_to_history(self, file_path, with_name=None):
        """Copy `file_path` to history, optionally renaming."""
        dest_name = with_name or file_path.name
        try:
            shutil.copy2(file_path, self.history_dir / dest_name)
        except OSError:                                                         # TODO: untested
            LOGGER.error(f'Failed to copy {file_path} to history.')

    def _copy_out_and_supp(self):
        """Copy OUT and SUPP directories to history."""
        for name in (DEFAULT_SUPP, DEFAULT_OUT):
            directory = self.cwd / name
            if not directory.is_dir():                                          # TODO: untested
                LOGGER.warning(f'Could not find {name} directory in '
                               f'{self.cwd}. It will not be copied '
                               'to history.')
                continue
            try:
                shutil.copytree(directory, self.history_dir / name)
            except OSError:                                                     # TODO: untested
                LOGGER.error(f'Failed to copy {name} directory to history.')

    def _deal_with_workhistory_and_history_info(self, discard=False):
        """Move work-history subfolders and update the history.info file."""
        tensor_nums = self._move_and_cleanup_workhistory(discard)
        tensor_nums.add(self.tensor_number)                                     # TODO: how about sorting on tensor numbers?

        self.history_info.append_entry(HistoryInfoEntry(
            tensor_nums=list(tensor_nums),
            job_nums=[self.max_job_for_tensor[tensor] + 1
                      for tensor in tensor_nums],
            timestamp=self.timestamp,
            discarded=discard,
            folder_name=self.history_dir.name,
            notes=self._read_and_clear_notes_file(),
            # Optional ones
            job_name=self.job_name,
            **self._infer_run_info_from_log()
            ))

    def _discard_workhistory_previous(self):                                    # TODO: untested
        """Remove 'previous'-labelled directories in work history."""
        work_hist_prev = self._get_workhistory_directories(
            contains=PREVIOUS_LABEL
            )
        for directory in work_hist_prev:
            try:
                shutil.rmtree(directory)
            except OSError:
                LOGGER.error(f'Failed to delete {directory} directory '
                             f'from {self.work_history_path}',
                             exc_info=True)

    def _find_base_name_for_history_subfolder(self):
        """Store internally the potential name of the history subfolder.

        This is not necessarily the name of the history directory that
        is actually created. Use self.history_dir.name instead. See
        also help(self._get_new_history_directory_name).

        Returns
        -------
        None.
        """
        suffix = (self.timestamp +
                  ('' if self.job_name is None else f'_{self.job_name}'))
        job_number = self.max_job_for_tensor[self.tensor_number]
        dir_name_fmt = f't{self.tensor_number:03d}.r{{job:03d}}_{suffix}'
        # If there is already a folder with the same name and correct
        # timestamp/suffix, we take that, otherwise, we increase the
        # job number: it's a new run.
        dir_name = dir_name_fmt.format(job=job_number)
        if not (self.top_level_history_path / dir_name).is_dir():
            job_number += 1
        base_name = dir_name_fmt.format(job=job_number)
        self._folder_names['history_dir_base_name'] = base_name

        # Now that we have a base name, we can also check if there
        # are already history folders with the same base name.
        self._state_info['history_with_same_base_name_exists'] = any(           # TODO: would it be enough to use the info above? We already know if there is a folder with the same base name.
            self._get_conflicting_history_subfolders()
            )

    def _find_new_history_directory_name(self):
        """Store internally the name of the history directory to be used.

        The path to the directory is then available via self.history_dir.

        Folder name has form 'tXXX.rYYY_<suffix>'. <suffix> can vary.
        It may be:
        - if there was a log file, and the folder is not
          present in history:
              <log_timestamp>[_<job_name>]
        - if there was a log file, but there is already a folder
          in history with the same timestamp and job_name (typically
          that means that bookkeeper was called multiple times with
          the same job_name on the same run):
              <log_timestamp>[_<job_name>]_moved-<bookie_timestamp>
        - if there was no log file:
              moved-<bookie_timestamp>[_<job_name>]
        - in the unlikely event that the folder in the previous point
          already exists (this would mean that the bookkeeper is called
          twice with the same job_name within one second):
              moved-<bookie_timestamp>[_<job_name>]_moved-<bookie_timestamp2>
          where <bookie_timestamp2> may be slightly later than
          <bookie_timestamp>, but is most likely the same.

        Returns
        -------
        None.
        """
        dir_name = self.history_dir_base_name
        if (self.top_level_history_path / dir_name).is_dir():                   # TODO: untested
            bookkeeper_timestamp = time.strftime('%y%m%d-%H%M%S',
                                                 time.localtime())
            dir_name = f'{dir_name}_moved-{bookkeeper_timestamp}'
        self._folder_names['history_dir'] = dir_name

    def _get_conflicting_history_subfolders(self):
        """Return an iterator of history subfolders with the same base name."""
        folders = self.top_level_history_path.glob(
            f'{self.history_dir_base_name}*'
            )
        return (f for f in folders if f.is_dir())

    def _get_current_workhistory_directories(self, contains=''):                # TODO: untested
        """Return a generator of subfolders in the current workhistory folder.

        As compared with _get_workhistory_directories, only those
        directories not marked as 'previous' are returned.

        Parameters
        ----------
        contains : str, optional
            Select only those subfolders whose name contains this
            string. Default is an empty string, corresponding to
            no filtering other than the one described in Returns.

        Returns
        -------
        subfolders : generator
            When iterated over, it yields paths to the immediate
            subfolders of the current workhistory folder whose name
            matches the HIST_FOLDER_RE regular expression, contains
            `contains`, but does not contain 'previous'.
        """
        directories = self._get_workhistory_directories(contains=contains)
        return (d for d in directories if PREVIOUS_LABEL not in d.name)

    def _get_workhistory_directories(self, contains=''):                        # TODO: untested
        """Return a generator of subfolders in the current workhistory folder.

        Parameters
        ----------
        contains : str, optional
            Select only those subdirectories whose name contains this
            string. Default is an empty string, corresponding to no
            filtering other than the one described in Returns.

        Returns
        -------
        subfolders : generator
            When iterated over, it yields paths to the immediate
            subdirectories of the current work-history folder whose
            name matches the HIST_FOLDER_RE regular expression, and
            whose name include `contains`.
        """
        work_history_path = self.work_history_path
        globbed = (work_history_path.glob(f'*{contains}*') if contains
                   else work_history_path.iterdir())
        return (d for d in globbed
                if d.is_dir() and HIST_FOLDER_RE.match(d.name))

    @_needs_update_for_attr('_state_info[last_log_lines]')
    def _infer_run_info_from_log(self):                                         # TODO: untested -- log_lines always empty?
        """Return a dictionary of information read from the calc log."""
        matched = {k: False for k in _LOG_FILE_RE}
        log_lines = self._state_info['last_log_lines'] or tuple()
        for line in reversed(log_lines):  # Info is at the end
            for info, already_matched in matched.items():
                if already_matched:
                    continue
                matched[info] = _LOG_FILE_RE[info].match(line)
            if all(matched.values()):
                break
        return {k: match[k] for k, match in matched.items() if match}

    def _make_and_copy_to_history(self, use_ori=False):
        """Create new history subfolder and copy all files there."""
        try:
            self.history_dir.mkdir()
        except OSError:
            LOGGER.error('Error: Could not create target directory '
                         f'{self.history_dir}\n Stopping...')
            raise
        self._copy_out_and_supp()
        self._copy_input_files_from_original_inputs_or_cwd(use_ori)
        self._copy_log_files()

    def _make_history_and_prepare_logger(self):
        """Make history folder and add handlers to the bookkeeper logger."""
        if self._state_info['logger_prepared']:
            return

        # Attach a stream handler to logger if not already present
        if not any(isinstance(h, logging.StreamHandler)
                   for h in LOGGER.handlers):
            LOGGER.addHandler(logging.StreamHandler())
        LOGGER.setLevel(logging.INFO)
        LOGGER.propagate = True

        try:  # Make top level history folder if not there yet
            self.top_level_history_path.mkdir(exist_ok=True)
        except OSError:                                                         # TODO: untested raise
            LOGGER.error('Error creating history folder.')
            raise

        # Attach file handler for history/bookkeeper.log
        bookkeeper_log = self.top_level_history_path / BOOKIE_LOGFILE
        LOGGER.addHandler(logging.FileHandler(bookkeeper_log, mode='a'))
        LOGGER.info(  # Log only once per instance
            '\n### Bookeeper running at '
            f'{time.strftime("%y%m%d-%H%M%S", time.localtime())} ###'
            )
        self._state_info['logger_prepared'] = True

    def _move_and_cleanup_workhistory(self, discard):
        """Move files from the current work-history folder, then clean up.

        If the current work-history folder is empty, it is always
        deleted. Otherwise, only if `discard` is True.

        Parameters
        ----------
        discard : bool
            Whether files should only be deleted, without copying them.

        Returns
        -------
        tensor_nums : set of int
            Indices of tensors found in the current work-history folder
            that have been moved to history as new history entries.
        """
        workhistory = self.work_history_path
        tensor_nums = set()
        if not workhistory.is_dir():
            return tensor_nums

        # Always remove any 'previous'-labelled folders
        self._discard_workhistory_previous()
        if not discard:
            tensor_nums = self._move_workhistory_folders()
        is_empty = not any(workhistory.iterdir())
        if is_empty or discard:
            try:
                shutil.rmtree(workhistory)
            except OSError:
                err_ = f'Failed to {{}} {workhistory.name} folder'
                err_ = err_.format('discard' if discard else 'delete empty')
                LOGGER.error(err_, exc_info=True)
        return tensor_nums

    def _move_workhistory_folders(self):
        """Move relevant folders from the current work history to history.

        Returns
        -------
        tensor_nums : set
            The tensor numbers of the folders that have been moved.
        """
        tensor_nums = set()
        directories = self._get_current_workhistory_directories(
            contains=self.timestamp
            )
        for directory in directories:
            match = HIST_FOLDER_RE.match(directory.name)
            if not match:  # Should always match
                continue
            tensor_num = int(match['tensor_num'])
            try:
                search_num = int(directory.name[6:9])                           # TODO: is this right? 6:9 should be job_num
            except (ValueError, IndexError):
                continue

            job_num = self.max_job_for_tensor[tensor_num] + 1
            newname = (
                f't{tensor_num:03d}.r{job_num:03d}.{search_num:03d}'
                + directory.name[9:]
                )
            target = self.top_level_history_path / newname
            if target.exists():
                LOGGER.error(f'Error: Failed to move {directory} to {target}: '
                             'Target path already exists. Stopping...')
                raise FileExistsError
            try:
                directory.replace(target)
            except OSError:
                LOGGER.error(f'Error: Failed to move {directory}.',
                             exc_info=True)
            tensor_nums.add(tensor_num)
        return tensor_nums

    def _read_and_clear_notes_file(self):
        """Return notes read from file. Clear the file contents."""
        notes_path = next(self.cwd.glob('notes*'), None)
        if notes_path is None:
            return ''                                                           # TODO: untested
        try:
            notes = notes_path.read_text(encoding='utf-8')
        except OSError:                                                         # TODO: untested
            LOGGER.error(f'Error: Failed to read {notes_path.name} file.',
                         exc_info=True)
            return ''
        try:
            notes_path.write_text('', encoding='utf-8')
        except OSError:                                                         # TODO: untested
            LOGGER.error(f'Failed to clear the {notes_path.name} '
                         'file after reading.', exc_info=True)
        return notes

    def _read_most_recent_log(self):
        """Read timestamp and lines from the most-recent calc log file."""
        last_log_lines = []
        split_logs = {}
        for prefix in CALC_LOG_PREFIXES:  # newest to oldest
            logs = (f for f in self._paths['calc_logs']
                    if f.name.startswith(prefix))
            try:
                split_logs[prefix] = max(logs, key=attrgetter('name'))
            except ValueError:
                pass
        try:
            most_recent_log = next(iter(split_logs.values()))
        except StopIteration:  # No log files
            timestamp = time.strftime('%y%m%d-%H%M%S', time.localtime())
            old_timestamp = f'moved-{timestamp}'
            return old_timestamp, last_log_lines

        old_timestamp = most_recent_log.name[-17:-4]
        try:  # pylint: disable=too-many-try-statements
            with most_recent_log.open('r', encoding='utf-8') as log_file:
                last_log_lines = log_file.readlines()
        except OSError:                                                         # TODO: untested
            pass
        return old_timestamp, last_log_lines

    def _remove_log_files(self):
        """Delete all log files in root."""
        _discard_files(*self.all_cwd_logs)

    def _remove_tensors_and_deltas(self):                                       # TODO: untested
        """Delete the most recent tensor and delta files."""
        tensor_file = f'Tensors/Tensors_{self.tensor_number:03d}.zip'
        delta_file = f'Deltas/Deltas_{self.tensor_number:03d}.zip'
        _discard_files(self.cwd / tensor_file, self.cwd / delta_file)

    def _remove_ori_files(self):
        """Delete '_ori'-suffixed files from root."""
        ori_files = (self.cwd / f'{file}_ori' for file in STATE_FILES)
        _discard_files(*ori_files)

    def _remove_out_and_supp(self):
        """Delete the SUPP and OUT directories from root."""
        _discard_files(self.cwd / DEFAULT_OUT, self.cwd / DEFAULT_SUPP)

    def _replace_state_files_from_ori(self):
        """Replace input files with their '_ori'-suffixed version."""
        for file in STATE_FILES:
            ori_file = self.cwd / f'{file}_ori'
            if ori_file.is_file():
                try:
                    ori_file.replace(self.cwd / file)
                except OSError:                                                 # TODO: untested
                    LOGGER.error(f'Failed to move {ori_file} to {file}.')
                    raise

    def _run_archive_mode(self):
        if self.history_with_same_base_name_exists:
            LOGGER.info(
                f'History directory for run {self.history_dir_base_name} '
                'exists. Exiting without doing anything.'
                )
            return BookeeperExitCode.NOTHING_TO_DO
        if not self.files_need_archiving:
            LOGGER.info('No files to be moved to history. Exiting '
                        'without doing anything.')
            return BookeeperExitCode.NOTHING_TO_DO
        self._make_and_copy_to_history(use_ori=False)

        # move old files to _ori and replace from OUT
        self._update_state_files_from_out()

        # workhistory and history.info
        self._deal_with_workhistory_and_history_info(discard=False)

        return BookeeperExitCode.SUCCESS

    def _run_clear_mode(self):
        if (not self.history_with_same_base_name_exists
            and self.files_need_archiving):
            LOGGER.info(f'History folder {self.history_dir} does not '
                        'yet exist. Running archive mode first.')
            self._make_and_copy_to_history(use_ori=True)

            # workhistory and history.info
            self._deal_with_workhistory_and_history_info(discard=False)

        # remove OUT, SUPP, logs and _ori files
        self._remove_log_files()
        self._remove_out_and_supp()
        self._remove_ori_files()

        return BookeeperExitCode.SUCCESS

    def _run_discard_common(self):
        """Removes files that get discarded for both DISCARD and DISCARD_FULL"""
        # replace input files from _ori
        self._replace_state_files_from_ori()

        if self.tensor_number not in self.max_job_for_tensor:
            self._remove_tensors_and_deltas()

        # remove OUT, SUPP and logs
        self._remove_log_files()
        self._remove_out_and_supp()

    def _run_discard_full_mode(self):
        # Check for notes in history.info
        if self.history_info.last_entry_has_notes:
            LOGGER.warning(f'The last entry in {HISTORY_INFO_NAME} has user '
                           'notes. If you really want to purge the last run, '
                           'remove the notes first.')
            return BookeeperExitCode.FAIL
        # Check whether there is an entry at all
        last_entry = self.history_info.last_entry
        if not last_entry:
            LOGGER.error('Error: Failed to remove last entry from '
                         f'{HISTORY_INFO_NAME}: No entries to remove.')
            return BookeeperExitCode.NOTHING_TO_DO
        # And check if there was some user edit in the last one
        if not last_entry.can_be_removed:
            if isinstance(last_entry, PureCommentEntry):
                err_ = 'is a comment-only field'
            elif not last_entry.was_understood:  # Can't be auto-fixed
                err_ = 'contains invalid fields that could not be interpreted'
            elif last_entry.misses_mandatory_fields:                            # TODO: untested
                err_ = 'some expected fields were deleted'
            else:    # Needs fixing, but can be done in --fixup mode            # TODO: implement
                assert not last_entry.has_notes  # Checked above
                err_ = ('contains fields with non-standard format (run '
                        'bookkeeper --fixup to automatically fix it)')
            LOGGER.error('Error: Failed to remove last entry from '
                         f'{HISTORY_INFO_NAME}: {err_}. Please '
                         'proceed manually.')
            return BookeeperExitCode.FAIL
                                                                                # TODO: !!! UNTESTED
        # The directory we want to remove is not self.history_dir (since that   # TODO: don't we also want to purge all the intermediate ones from workhistory?
        # would be _moved-<timestamp>), but the one with the same base name!
        dir_to_remove = (self.top_level_history_path
                         / self.history_dir_base_name)
        if not dir_to_remove.is_dir():
            LOGGER.error('FULL_DISCARD mode failed: could not identify '
                         'directory to remove. Please proceed manually.')
            return BookeeperExitCode.FAIL

        # Remove history folder                                                 # TODO: how about the workhistory ones?
        try:
            shutil.rmtree(dir_to_remove)
        except OSError:
            LOGGER.error(f'Error: Failed to delete {dir_to_remove}.')
            return BookeeperExitCode.FAIL
        self._run_discard_common()
        # And the history entry from history.info
        self.history_info.remove_last_entry()
        return BookeeperExitCode.SUCCESS

    def _run_discard_mode(self):
        if (not self.history_with_same_base_name_exists
                and self.files_need_archiving):
            LOGGER.info(f'History folder {self.history_dir} does not '
                        'yet exist. Running archive mode first.')
            self._make_and_copy_to_history(use_ori=False)

            # workhistory and history.info
            self._deal_with_workhistory_and_history_info(discard=True)

        self._run_discard_common()
        try:
            self.history_info.discard_last_entry()
        except NoHistoryEntryError as exc:
            LOGGER.error('Error: Failed to mark last entry as '
                         f'discarded in {HISTORY_INFO_NAME}: {exc}')
            return BookeeperExitCode.FAIL
        return BookeeperExitCode.SUCCESS

    def _update_state_files_from_out(self):
        """Move state files from OUT to root. Rename old to '_ori'."""
        out_path = self.cwd / 'OUT'
        if not out_path.is_dir():
            return                                                              # TODO: untested
        for file in STATE_FILES:
            out_file = out_path / f'{file}_OUT'
            cwd_file = self.cwd / file
            if not out_file.is_file():
                continue
            # NB: all the moving around is local to the same folder
            # tree, so we don't really need shutil. Path.replace always
            # works for the same file system. Notice also the use of
            # replace rather than rename, as the behavior of rename
            # is not identical for all platforms.
            cwd_file.replace(self.cwd / f'{file}_ori')
            out_file.replace(cwd_file)


def _check_newer(older, newer):
    """Raise if file `older` is newer than file `newer`."""
    newer_timestamp = newer.stat().st_mtime
    older_timestamp = older.stat().st_mtime
    if newer_timestamp < older_timestamp:                                       # TODO: untested
        raise _FileNotOlderError


def _discard_files(*file_paths):
    """Delete files at `file_paths`. Log if they can't be deleted."""
    for file in file_paths:
        _discard_one_file(file)


def _discard_one_file(file):
    """Delete a file or directory, if it exists."""
    if not file.exists():
        return
    if file.is_file():
        try:
            file.unlink()
        except OSError:                                                         # TODO: untested
            LOGGER.error(f'Failed to discard file {file.name}.')
        return
    assert file.is_dir()  # Should be a directory
    try:
        shutil.rmtree(file)
    except OSError:                                                             # TODO: untested
        LOGGER.error(f'Failed to discard directory {file.name}.')
