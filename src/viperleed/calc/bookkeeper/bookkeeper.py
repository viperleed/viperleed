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
import logging
from operator import attrgetter
import os
from pathlib import Path
import re
import shutil
import time

from viperleed.calc import DEFAULT_HISTORY
from viperleed.calc import DEFAULT_WORK_HISTORY
from viperleed.calc import LOG_PREFIX
from viperleed.calc import ORIGINAL_INPUTS_DIR_NAME
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
    'tleedm',   # For backwards compatibility
    )
HIST_FOLDER_RE = re.compile(
    r't(?P<tensor_num>[0-9]{3}).r(?P<job_num>[0-9]{3})_'
    )
STATE_FILES = ('PARAMETERS', 'POSCAR', 'VIBROCC')

# Input files that may be generated at runtime - don't warn if missing
RUNTIME_GENERATED_INPUT_FILES = ('IVBEAMS', 'PHASESHIFTS')


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
        self.cwd = Path(cwd)
        self.top_level_history_path = self.cwd / history_name
        self.work_history_path = self.cwd / work_history_name
        self.orig_inputs_dir = (self.cwd / DEFAULT_SUPP
                                / ORIGINAL_INPUTS_DIR_NAME)

        self.job_name = job_name
        self._state_info = {
            'logger_prepared': False,
            }
        self._make_history_and_prepare_logger()

        # history.info handler - creates file if not yet there
        self.history_info = HistoryInfoFile(self.cwd / HISTORY_INFO_NAME,
                                            create_new=True)
        self.update_from_cwd()                                                  # TODO: this is annoying. Logs stuff before running. Can we do it in run instead?

    @property
    def all_cwd_logs(self):
        """Return paths to all the log files in the root directory."""
        return sum(self.cwd_logs, start=[])

    @property
    def base_history_dir_name(self):
        """The name of the history directory based on timestamp, number and name.

        This is not necessarily the name of the history directory! Use
        self.history_dir.name instead."""
        suffix = (self.timestamp +
                  ('' if self.job_name is None else f'_{self.job_name}'))
        job_number = self.max_job_for_tensor[self.tensor_number]
        dir_name_fmt = f't{self.tensor_number:03d}.r{{job:03d}}_{suffix}'
        # If there is already a folder with the same name and correct
        # timestamp/suffix, we take that, otherwise, we increase the
        # job number
        dir_name = dir_name_fmt.format(job=job_number)
        if not (self.top_level_history_path / dir_name).is_dir():
            job_number += 1
        return dir_name_fmt.format(job=job_number)

    @property
    def cwd_logs(self):
        calc_logs, other_logs = [], []
        for file in self.cwd.glob('*.log'):
            if not file.is_file():
                continue
            container = (calc_logs if file.name.startswith(CALC_LOG_PREFIXES)
                         else other_logs)
            container.append(file)
        return calc_logs, other_logs





    @property
    def files_needs_archiving(self):
        """Check if there are any files that need archiving."""
        # check for OUT and SUPP
        files_to_archive = [self.cwd / DEFAULT_OUT, self.cwd / DEFAULT_SUPP]
        # any calc logs
        files_to_archive.extend(self.cwd_logs[0])
        # check workhistory
        if self.work_history_path.is_dir():                                     # TODO: untested
            files_to_archive.extend(_get_current_workhistory_directories(
                self.work_history_path, contains='r'))
        files_to_archive = [file for file in files_to_archive
                            if file.is_file() or file.is_dir()]
        return bool(files_to_archive)

    @property
    def history_with_same_base_name_exists(self):
        """Check if a history folder with the same base name exists."""
        all_history_entires = self.top_level_history_path.glob('*')
        return any(entry.is_dir() and
                   entry.name.startswith(self.base_history_dir_name)
                   for entry in all_history_entires)

    @property
    def archiving_required(self):
        """Check if archiving is required."""
        return (self.files_needs_archiving and
                not self.history_with_same_base_name_exists)

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
            If `mode` is not a valid bookkeeper mode, but
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
        return runner()

    def update_from_cwd(self):
        """Update timestamp, tensor number, log lines, etc. from root.

        This method is called in __init__, but can also be called manually if
        a new run happens during the lifetime of the bookkeeper.
        """
        # Figure out the number of the tensor (it's the most recent one)
        # and the highest run number currently stored for each tensor in
        # history_path
        self.tensor_number = getMaxTensorIndex(home=self.cwd, zip_only=True)
        self.max_job_for_tensor = self._find_max_run_per_tensor()

        # Infer timestamp from log file, if possible
        self.timestamp, self.last_log_lines = self._read_most_recent_log()

        # get history dir to deal with
        self.history_dir = (self.top_level_history_path /
                            self._get_new_history_directory_name())

        # Read the current state of the history.info file
        self.history_info.read()

    def _copy_input_files_from_original_inputs_and_cwd(self, use_ori=False):
        """Copy files from original_inputs_path to target_path."""
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
            elif cwd_file.is_file():  # copy cwd and warn                       # TODO: untested
                copy_file, with_name = cwd_file, f'{cwd_file.name}_from_root'
                LOGGER.warning(
                    f'File {file} not found in {ORIGINAL_INPUTS_DIR_NAME}. '
                    'Using file from root directory instead and renaming to '
                    f'{with_name}.'
                    )
            if copy_file:
                self._copy_one_file_to_history(copy_file, with_name=with_name)

    def _copy_log_files_to_history(self):
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

        run_info, r_ref, r_super = self._infer_run_info_from_log()
        tensor_nums.add(self.tensor_number)                                     # TODO: how about sorting on tensor numbers?

        self.history_info.append_entry(HistoryInfoEntry(
            tensor_nums=list(tensor_nums),
            job_nums=[self.max_job_for_tensor[tensor] + 1
                      for tensor in tensor_nums],
            timestamp=self.timestamp,
            discarded=discard,
            folder_name=self.history_dir.name,
            notes=self._read_and_clear_notes_file(),
            # optionals
            job_name=self.job_name,
            run_info=run_info,
            r_ref=r_ref,
            r_super=r_super,
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

    def _find_max_run_per_tensor(self):
        """Return maximum run numbers for all directories in 'history'.

        Returns
        -------
        max_nums : defaultdict(int)
            The maximum run number currently stored in the 'history'
            folder for each of the tensors in 'history'.
        """
        history_path = self.top_level_history_path
        max_nums = defaultdict(int)  # max. job number per tensor number
        for directory in history_path.iterdir():
            match = HIST_FOLDER_RE.match(directory.name)
            if not directory.is_dir() or not match:                             # TODO: untested
                continue
            tensor_num = int(match['tensor_num'])
            job_num = int(match['job_num'])
            max_nums[tensor_num] = max(max_nums[tensor_num], job_num)
        return max_nums

    def _get_current_workhistory_directories(self, contains=''):                # TODO: untested
        """Return a generator of subfolders in the current workhistory folder.

        As compared with _get_workhistory_directories, only those
        directories not marked as 'previous' are returned.

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
            subdirectories of the current workhistory folder whose
            name matches the HIST_FOLDER_RE regular expression, whose
            name include `contains`, but does not contain 'previous'.
        """
        directories = self._get_workhistory_directories(contains=contains)
        return (d for d in directories if PREVIOUS_LABEL not in d.name)

    def _get_new_history_directory_name(self):
        """Return the name of a history directory for a given run.

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
        """
        dir_name = self.base_history_dir_name
        if (self.top_level_history_path / dir_name).is_dir():                   # TODO: untested
            bookkeeper_timestamp = time.strftime('%y%m%d-%H%M%S',
                                                 time.localtime())
            dir_name = f'{dir_name}_moved-{bookkeeper_timestamp}'
        return dir_name

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

    def _move_and_cleanup_workhistory(self, discard):                           # TODO: untested
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

    def _move_workhistory_folders(self):                                        # TODO: untested
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

            tensor_nums.add(tensor_num)
            job_num = self.max_job_for_tensor[tensor_num] + 1
            newname = (
                f't{tensor_num:03d}.r{job_num:03d}.{search_num:03d}'
                + directory.name[9:]
                )
            try:
                shutil.move(directory, self.top_level_history_path / newname)
            except OSError:
                LOGGER.error(f'Error: Failed to move {directory}.')
        return tensor_nums

    def _infer_run_info_from_log(self):                                         # TODO: untested -- log_lines always empty?
        run_info, r_ref, r_super = None, None, None
        for line in self.last_log_lines:
            if line.startswith('Executed segments: '):
                run_info = line.split('Executed segments: ')[1].strip()
                break
        for line in self.last_log_lines:
            if 'Final R (refcalc)' in line:
                r_ref = float(line.split(':', maxsplit=1)[1].strip())
                break
        for line in self.last_log_lines:
            if 'Final R (superpos)' in line:
                r_super = float(line.split(':', maxsplit=1)[1].strip())
                break
        return run_info, r_ref, r_super

    def _make_and_copy_to_history(self, use_ori=False):
        """Create new history subfolder and copy all files there."""
        try:
            self.history_dir.mkdir()
        except OSError:                                                         # TODO: untested
            LOGGER.error('Error: Could not create target directory '
                         f'{self.history_dir}\n Stopping...')
            raise
        self._copy_out_and_supp()
        self._copy_input_files_from_original_inputs_and_cwd(use_ori)
        self._copy_log_files_to_history()

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
        try:
            most_recent_log = max(
                (log_file
                 for prefix in CALC_LOG_PREFIXES
                 for log_file in self.cwd.glob(f'{prefix}*.log')
                 if log_file.is_file()),
                key=attrgetter('name')
                )
        except ValueError:  # No log files
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
                f'History directory for run {self.base_history_dir_name} '
                'exists. Exiting without doing anything.'
                )
        if not self.files_needs_archiving:
            return BookeeperExitCode.NOTHING_TO_DO
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
            and self.files_needs_archiving):
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
            return BookeeperExitCode.FAIL
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
                         / self.base_history_dir_name)
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
                and self.files_needs_archiving):
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
