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
from dataclasses import dataclass
from enum import Enum
from operator import attrgetter
from pathlib import Path
from typing import List, Optional
import logging
import os
import re
import shutil
import time

from viperleed.calc import DEFAULT_HISTORY
from viperleed.calc import DEFAULT_WORK_HISTORY
from viperleed.calc import LOG_PREFIX
from viperleed.calc import ORIGINAL_INPUTS_DIR_NAME
from viperleed.calc.lib.leedbase import getMaxTensorIndex
from viperleed.calc.sections.calc_section import ALL_INPUT_FILES
from viperleed.calc.sections.cleanup import PREVIOUS_LABEL
from viperleed.cli_base import ViPErLEEDCLI
from viperleed.calc.sections.cleanup import DEFAULT_OUT, DEFAULT_SUPP

logger = logging.getLogger(__name__)

_CALC_LOG_PREFIXES = (
    LOG_PREFIX,
    'tleedm',   # For backwards compatibility
    )
HIST_FOLDER_RE = re.compile(
    r't(?P<tensor_num>[0-9]{3}).r(?P<job_num>[0-9]{3})_'
    )
HISTORY_INFO_NAME = 'history.info'
_HISTORY_INFO_SPACING = 12  # For the leftmost field in history.info
HISTORY_INFO_SEPARATOR = '\n###########\n'

STATE_FILES = ('PARAMETERS', 'POSCAR', 'VIBROCC')

class BookkeeperMode(Enum):
    """Enumeration of bookkeeper modes.

    Attributes
    ----------
    ARCHIVE
        Store last run in history. Overwrite PARAMETERS, POSCAR & VIBROCC from
        OUT. Runs after run_calc by default.
    CLEAR
        Clear the input directory of last run.
        Runs before run_calc by default.
    DISCARD
        Re-start from the same input as the previous run. The discarded run is
        kept in history. Has to be run manually after run_calc.
    DISCARD_FULL
        Discard previous run as if it never happened and removes it from 
        history. Has to be run manually after run_calc.
    """
    ARCHIVE = 'archive'
    CLEAR = 'clear'
    DISCARD = 'discard'
    DISCARD_FULL = 'discard_full'

    @property
    def discard(self):
        """Return whether this is mode DISCARD."""
        return self is BookkeeperMode.DISCARD


class Bookkeeper():
    """Bookkeeper to archive or discard the most recent viperleed calc run.
    """
    def __init__(self,
                 job_name=None,
                 history_name=DEFAULT_HISTORY,
                 work_history_name=DEFAULT_WORK_HISTORY,
                 cwd=Path.cwd()):
        """Initialize the bookkeeper.

        Parameters
        ----------
        job_name : str or None, optional
            Custom name to append to the stored folder and to history.info.
            If not given or None, no extra name is added. Default is None.
        history_name : str, optional
            The name of the folder in the current directory where the
            most recent run should be archived. Default is 'history'.
        work_history_name : str, optional
            The name of the workhistory subfolder of ./work where
            intermediate runs may have been stored. Results are also
            copied from here to ./history_name. Default is workhistory.
        cwd : Path, optional
            The current working directory. Default is the current directory.
        """
        self.cwd = Path(cwd)
        self.top_level_history_path = self.cwd / history_name
        self.work_history_path = self.cwd / work_history_name
        self.orig_inputs_dir = (self.cwd / DEFAULT_SUPP
                                / ORIGINAL_INPUTS_DIR_NAME)

        self.job_name = job_name
        # attach a stream handler to logger if not already present
        if not any(isinstance(h, logging.StreamHandler)
                   for h in logger.handlers):
            logger.addHandler(logging.StreamHandler())
        logger.setLevel(logging.INFO)
        logger.propagate = True

        # Make top level history folder if not there yet
        try:
            self.top_level_history_path.mkdir(exist_ok=True)
        except OSError:
            logger.error('Error creating history folder.')
            raise

        # history.info handler - creates file if not yet there
        self.history_info = HistoryInfoFile(self.cwd / HISTORY_INFO_NAME,
                                            create_new=True)

        self.update_from_cwd()

    def update_from_cwd(self):
        """Updates timestamp, tensor number and log lines, etc. from cwd.
        
        This method is called in __init__, but can also be called manually if
        a new run happens during the lifetime of the bookkeeper.
        """
        # Figure out the number of the tensor (it's the most recent one)
        # and the highest run number currently stored for each tensor in
        # history_path
        self.tensor_number = getMaxTensorIndex(home=self.cwd, zip_only=True)
        self.max_job_for_tensor = _find_max_run_per_tensor(
            self.top_level_history_path)

        # Infer timestamp from log file, if possible
        self.timestamp, self.last_log_lines = _read_most_recent_log(self.cwd)

        # get history dir to deal with
        self.history_dir = (self.top_level_history_path /
                            self._get_new_history_directory_name())

    @property
    def base_history_dir_name(self):
        """The name of the history directory based on timestam, number and name.

        This is not necessarily the name of the history directory! Use
        self.history_dir.name instead."""
        suffix = (self.timestamp +
                  ('' if self.job_name is None else f'_{self.job_name}'))
        job_number = self.max_job_for_tensor[self.tensor_number]
        dir_name = lambda job_number: (
            f't{self.tensor_number:03d}.r{job_number:03d}_{suffix}')
        # if there is already a folder with the same name and correct
        # timestamp/suffix, we take that, otherwise, we increase the job number
        if not (self.top_level_history_path / dir_name(job_number)).is_dir():
            job_number += 1
        return dir_name(job_number)


    def _get_new_history_directory_name(self):
        """Return the name of a history directory for a given run.

        Folder has form 'tXXX.rYYY_<suffix>'. <suffix> can vary. It may
        be:
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
        if (self.top_level_history_path / dir_name).is_dir():
            bookkeeper_timestamp = time.strftime('%y%m%d-%H%M%S',
                                                 time.localtime())
            dir_name = f'{dir_name}_moved-{bookkeeper_timestamp}'
        return dir_name


    @property
    def cwd_ori_files(self):
        return [self.cwd / f'{file}_ori' for file in STATE_FILES]

    @property
    def cwd_logs(self):
        calc_logs, other_logs = [], []
        for file in self.cwd.glob('*.log'):
            if not file.is_file():
                continue
            container = (calc_logs if file.name.startswith(_CALC_LOG_PREFIXES)
                        else other_logs)
            container.append(file)
        return calc_logs, other_logs

    @property
    def all_cwd_logs(self):
        return self.cwd_logs[0] + self.cwd_logs[1]


    @property
    def files_needs_archiving(self):
        """Check if there are any files that need archiving."""
        # check for OUT and SUPP
        files_to_archive = [self.cwd / DEFAULT_OUT, self.cwd / DEFAULT_SUPP]
        # any calc logs
        files_to_archive.extend(self.cwd_logs[0])
        # check workhistory
        if self.work_history_path.is_dir():
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
        """Runs the bookkeeper in the given mode.

        Parameters
        ----------
        mode : str or BookkeeperMode
            Which bookkeeper mode to use. See help(BookkeeperMode).

        Returns
        -------
        exit_code : int
            0 -> success
            1 -> No files to be copied

        Raises
        ------
        OSError
            If creation of the history folder or any of the subfolders
            where results are to be stored fails.
        """
        mode = BookkeeperMode(mode)

        # ARCHIVE, CLEAR, DISCARD do archiving
        if mode is BookkeeperMode.ARCHIVE:
            logger.info("Running bookkeeper in ARCHIVE mode.")
            return self._run_archive_mode()
        elif mode is BookkeeperMode.CLEAR:
            logger.info("Running bookkeeper in CLEAR mode.")
            return self._run_clear_mode()
        elif mode is BookkeeperMode.DISCARD:
            logger.info("Running bookkeeper in DISCARD mode.")
            return self._run_discard_mode()
        elif mode is BookkeeperMode.DISCARD_FULL:
            logger.info("Running bookkeeper in DISCARD_FULL mode.")
            return self._run_discard_full_mode()
        else:
            raise ValueError(f'Unknown mode {mode}')

    def _run_archive_mode(self):
        if self.history_with_same_base_name_exists:
            logger.info(
                f'History directory for run {self.base_history_dir_name} '
                'exists. Exiting without doing anything.')
            return 1
        if not self.files_needs_archiving:
            logger.info('No files to be moved to history. Exiting without doing'
                           ' anything.')
            return 1
        self._make_and_copy_to_history(use_ori=False)

        # move old files to _ori and replace from OUT
        _update_state_files_from_out(self.cwd)

        # workhistory and history.info
        self._deal_with_workhistory_and_history_info(discard=False)

        return 0

    def _run_clear_mode(self):
        if (not self.history_with_same_base_name_exists
            and self.files_needs_archiving):
            logger.info(f'History folder {self.history_dir} does not yet exist.'
                        ' Running archive mode first.')
            self._make_and_copy_to_history(use_ori=True)

            # workhistory and history.info
            self._deal_with_workhistory_and_history_info(discard=False)

        # remove OUT, SUPP, logs and _ori files
        self.remove_log_files()
        self.remove_out_and_supp()
        self.remove_ori_files()

        return 0

    def _run_discard_mode(self):
        if (not self.history_with_same_base_name_exists
            and self.files_needs_archiving):
            logger.info(f'History folder {self.history_dir} does not yet exist.'
                        ' Running archive mode first.')
            self._make_and_copy_to_history(use_ori=False)

            # workhistory and history.info
            self._deal_with_workhistory_and_history_info(discard=True)

        self._discard_common()
        try:
            self.history_info.discard_last_entry()
        except ValueError:
            logger.warning('Error: Failed to mark last entry as discarded in '
                         f'{HISTORY_INFO_NAME}.')

    def _run_discard_full_mode(self):
        # check for notes in history.info
        if self.history_info.last_entry_has_notes:
            logger.warning(f'The last entry in {HISTORY_INFO_NAME} has user notes. '
                           'If you really want to purge the last run, remove '
                           'the notes first.')
            return 1

        # the directory we want to remove is not self.history_dir (since that
        # would be _moved-<timestamp>), but the one with the same base name!
        dir_to_remove = (self.top_level_history_path
                         / self.base_history_dir_name)
        # remove history entry
        if dir_to_remove.is_dir():
            try:
                shutil.rmtree(dir_to_remove)
            except OSError:
                logger.error(f'Error: Failed to delete {dir_to_remove}.')
                return 1
            self._discard_common()
        else:
            logger.error(f'FULL_DISCARD mode failed: could not identify '
                         'directory to remove. Please proceed manually.')
        # remove history entry from history.info
        try:
            self.history_info.remove_last_entry()
        except ValueError:
            logger.warning('Error: Failed to remove last entry from '
                         f'{HISTORY_INFO_NAME}.')


    def _discard_common(self):
        """Removes files that get discarded for both DISCARD and DISCARD_FULL"""
        # replace input files from _ori
        self.replace_state_files_from_ori()

        if self.tensor_number not in self.max_job_for_tensor:
            self.remove_tensors_and_deltas()

        # remove OUT, SUPP and logs
        self.remove_log_files()
        self.remove_out_and_supp()

    def _make_and_copy_to_history(self, use_ori=False):
                # Copy everything to history
        _create_new_history_dir(self.history_dir)
        self.copy_out_and_supp()
        self.copy_input_files_from_original_inputs_and_cwd(use_ori)
        self.copy_log_files_to_history()

    def _deal_with_workhistory_and_history_info(self, discard=False):
        tensor_nums = _move_and_cleanup_workhistory(self.work_history_path,
                                                    self.top_level_history_path,
                                                    self.timestamp,
                                                    self.max_job_for_tensor,
                                                    discard)
        tensor_nums.add(self.tensor_number)
        
        run_info, r_ref, r_super = self._infer_run_info_from_log()

        self.history_info.append_entry(HistoryInfoEntry(
            tensor_nums=list(tensor_nums),
            job_nums=[self.max_job_for_tensor[tensor] + 1
                      for tensor in tensor_nums],
            timestamp=self.timestamp,
            discarded=discard,
            folder_name=self.history_dir.name,
            notes=_read_and_clear_notes_file(self.cwd),
            # optionals
            job_name=self.job_name,
            run_info=run_info,
            r_ref=r_ref,
            r_super=r_super,
            ))

    def _infer_run_info_from_log(self):
        run_info, r_ref, r_super = None, None, None
        for line in self.last_log_lines:
            if line.startswith('Executed segments: '):
                run_info = line.split('Executed segments: ')[1].strip()
                break
        for line in self.last_log_lines:
            if 'refcalc' in line:
                r_ref = float(line.split(':', maxsplit=1)[1].strip())
                break
        for line in self.last_log_lines:
            if 'superpos' in line:
                r_super = float(line.split(':', maxsplit=1)[1].strip())
                break
        return run_info, r_ref, r_super


    def remove_ori_files(self):
        _move_or_discard_files(self.cwd_ori_files, self.cwd, True)

    def remove_out_and_supp(self):
        _move_or_discard_files([self.cwd / DEFAULT_OUT, self.cwd / DEFAULT_SUPP],
                               self.cwd, True)

    def remove_log_files(self):
        _move_or_discard_files(self.all_cwd_logs, self.cwd, True)

    def remove_tensors_and_deltas(self):
        tensor_file = self.cwd / 'Tensors' / f'Tensors_{self.tensor_number:03d}.zip'
        delta_file = self.cwd / 'Deltas' / f'Deltas_{self.tensor_number:03d}.zip'
        _move_or_discard_files((tensor_file, delta_file), self.cwd, True)

    def replace_state_files_from_ori(self):
        for file in STATE_FILES:
            ori_file = self.cwd / f'{file}_ori'
            if ori_file.is_file():
                try:
                    shutil.move(ori_file, self.cwd / file)
                except OSError:
                    logger.error(f'Error: failed to move {ori_file} to {file}.')
                    raise

    def copy_input_files_from_original_inputs_and_cwd(self, use_ori=False):
        """Copy files from original_inputs_path to target_path."""
        for file in ALL_INPUT_FILES:
            original_file = self.orig_inputs_dir / file
            cwd_file = self.cwd / file
            if file in STATE_FILES and use_ori:
                cwd_file = self.cwd / f'{file}_ori'
            if original_file.is_file() and cwd_file.is_file():
                # copy original, but warn if cwd is newer
                _copy_one_file_to_history(original_file, self.history_dir)
                original_timestamp = original_file.stat().st_mtime
                cwd_timestamp = cwd_file.stat().st_mtime
                if original_timestamp < cwd_timestamp:
                    logger.warning(
                        f'File {file} from {ORIGINAL_INPUTS_DIR_NAME} was '
                        'copied to history, but the file in the input '
                        'directory is newer.')
            elif original_file.is_file():
                # just copy original
                _copy_one_file_to_history(original_file, self.history_dir)
            elif cwd_file.is_file():
                # copy cwd and warn
                _copy_one_file_to_history(
                    self.cwd / file, self.history_dir)
                logger.warning(f'File {file} not found in '
                            f'{ORIGINAL_INPUTS_DIR_NAME}. Using file from root '
                            'directory instead and renaming to '
                            f'{cwd_file.name}_from_root.')
                try:
                    os.rename(self.history_dir / file,
                            self.history_dir / f'{cwd_file.name}_from_root')
                except OSError:
                    logger.error(f'Failed to rename {file} to '
                                f'{cwd_file.name}_from_root.')

    def copy_out_and_supp(self):
        """Copy OUT and SUPP directories to history."""
        for name in (DEFAULT_SUPP, DEFAULT_OUT):
            dir = self.cwd / name
            if not dir.is_dir():
                logger.warning(f'Could not find {name} directory in '
                               f'{self.cwd}. It will not be copied to history.')
                continue
            try:
                shutil.copytree(dir, self.history_dir / name)
            except OSError:
                logger.error(f'Failed to copy {name} directory to history.')

    def copy_log_files_to_history(self):
        """Copy log files to history."""
        calc_logs, other_logs = self.cwd_logs
        log_files = calc_logs + other_logs
        for file in log_files:
            _copy_one_file_to_history(file, self.history_dir)


def store_input_files_to_history(root_path, history_path):
    """Find input files in `root_path` and copy them to `history_path`.

    Parameters
    ----------
    root_path : pathlike
        Root directory from which to take files. Should be cwd.
    history_path : pathlike
        Path to the history directory in which the files should
        be stored.
    """
    root_path, history_path = Path(root_path), Path(history_path)
    original_inputs_path = root_path / 'SUPP' / ORIGINAL_INPUTS_DIR_NAME
    if original_inputs_path.is_dir():
        input_origin_path = original_inputs_path
    else:
        input_origin_path = root_path
        print(f'Could not find directory {ORIGINAL_INPUTS_DIR_NAME!r} with '
              'unaltered input files. Files will instead be copied from the '
              'root directory.')

    # Only files, no directories
    files_to_copy = (file for file in input_origin_path.iterdir()
                     if file.is_file() and file.name in ALL_INPUT_FILES)
    for file in files_to_copy:
        try:
            shutil.copy2(file, history_path/file.name)
        except OSError as exc:
            print(f'Failed to copy file {file} to history: {exc}')

class HistoryInfoFile:
    """Deals with the history.info file in a history directory."""

    def __init__(self, file_path, create_new=False):
        self.path = file_path
        if not self.path.is_file():
            if not create_new:
                raise FileNotFoundError("history.info file not found at "
                                        f"{self.path}.")
            # create new file
            self.path.touch()

    @property
    def raw_contents(self):
        with open(self.path, 'r', encoding='utf-8') as f:
            raw_contents = f.read()
        return raw_contents

    @property
    def last_entry(self):
        if not self.raw_contents:
            return None
        return self._parse_entry(
            self.raw_contents.split(HISTORY_INFO_SEPARATOR.strip())[-1])

    @property
    def last_entry_has_notes(self):
        if self.last_entry is None:
            return False
        return bool(self.last_entry.notes)

    @property
    def last_entry_was_discarded(self):
        if self.last_entry is None:
            return False
        return self.last_entry.discarded

    def append_entry(self, entry):
        skip_separator = (
            self.raw_contents.strip().endswith(HISTORY_INFO_SEPARATOR.strip())
            or not self.raw_contents.strip())
        with open(self.path, 'a', encoding='utf-8') as f:
            if not skip_separator:
                f.write(HISTORY_INFO_SEPARATOR)
            f.write(str(entry))

    def discard_last_entry(self):
        """Mark the last entry in the history.info file as discarded."""
        if self.last_entry is None:
            raise ValueError("No entries to discard.")
        last_entry = self.last_entry
        if last_entry.discarded:
            logger.warning('Last entry is already discarded.')
        last_entry.discarded = True
        self.remove_last_entry()
        self.append_entry(last_entry)

    def remove_last_entry(self):
        """Discard the last entry form the history.info file."""
        if self.last_entry is None:
            raise ValueError("No entries to remove.")
        if HISTORY_INFO_SEPARATOR.strip() in self.raw_contents:
            content_without_last = self.raw_contents.rsplit(
                HISTORY_INFO_SEPARATOR.strip(), 1)[0]
        else:
            # only one entry
            content_without_last = ''
        if content_without_last.endswith('\n'):
            content_without_last = content_without_last[:-len('\n')]
        # clear file and write back entries
        self.path.write_text(content_without_last)

    def _parse_entry(self, entry_str):
        # remove leading and trailing whitespace
        entry_str = entry_str.strip()
        # check for 'DISCARDED' at the end
        if entry_str.endswith('DISCARDED'):
            discarded = True
            entry_str = entry_str[:-len('DISCARDED')]
        else:
            discarded = False
        # split at 'Notes: ' because that is the only one that can be multiline
        general_info, notes = entry_str.split('Notes:', 1)
        notes = notes.replace('Notes: ', '').strip()
        # parse general_info
        general_info = iter(general_info.split('\n'))
        tensors = [int(num) for num in
                   next(general_info).replace('# TENSORS ', '').split()[1:]]
        jobs = [int(num) for num in
                next(general_info).replace('# JOB ID ', '').split()[1:]]
        the_rest = next(general_info)
        # optionals
        job_name, run_info, r_ref, r_super = None, None, None, None
        if the_rest.startswith('# JOB NAME '):
            job_name = the_rest.replace('# JOB NAME ', '').strip()
            the_rest = next(general_info)
        if the_rest.startswith('# RUN '):
            run_info = the_rest.replace('# RUN ', '').strip()
            the_rest = next(general_info)
        if the_rest.startswith('# TIME'):
            time = the_rest.replace('# TIME ', '').strip()
            the_rest = next(general_info)
        if the_rest.startswith('# R REF'):
            r_ref = float(the_rest.replace('# R REF ', ''))
            the_rest = next(general_info)
        if the_rest.startswith('# R SUPER'):
            r_super = float(the_rest.replace('# R SUPER ', ''))
            the_rest = next(general_info)
        folder = the_rest.replace('# FOLDER ', '').strip()

        # this should be all, if there is more, we have a problem
        try:
            # there may be empty lines at the end
            while not next(general_info):
                pass
        except StopIteration:
            pass
        else:
            raise ValueError("Error parsing history.info file.")

        return HistoryInfoEntry(tensors, jobs, time, folder, notes, discarded,
                                job_name, run_info, r_ref, r_super)

@dataclass
class HistoryInfoEntry:
    tensor_nums: List[int]
    job_nums: List[int]
    timestamp: str
    folder_name: str
    notes: str
    discarded: bool
    job_name: Optional[str] = None
    run_info: Optional[str] = None
    r_ref: Optional[float] = None
    r_super: Optional[float] = None


    def __str__(self):
        tensor_str = 'None' if self.tensor_nums == [0] else str(self.tensor_nums)[1:-1]
        job_str = str(self.job_nums)[1:-1]
        # translate timestamp if necessary
        time_str = (_translate_timestamp(self.timestamp)
                    if '-' in self.timestamp else self.timestamp)
        folder_str = self.folder_name
        return (
            '\n'
            + '# TENSORS '.ljust(_HISTORY_INFO_SPACING) + tensor_str + '\n'
            + '# JOB ID '.ljust(_HISTORY_INFO_SPACING) + job_str + '\n'
            + ('# JOB NAME '.ljust(_HISTORY_INFO_SPACING) + self.job_name + '\n'
                if self.job_name else '')
            + ('# RUN '.ljust(_HISTORY_INFO_SPACING) + self.run_info + '\n'
                if self.run_info else '')
            + '# TIME '.ljust(_HISTORY_INFO_SPACING) + time_str + '\n'
            + ('# R REF '.ljust(_HISTORY_INFO_SPACING) +f'{self.r_ref:.4f}\n'
                if self.r_ref is not None else '')
            + ('# R SUPER '.ljust(_HISTORY_INFO_SPACING) + f'{self.r_super:.4f}\n'
                if self.r_super is not None else '')
            + '# FOLDER '.ljust(_HISTORY_INFO_SPACING) + folder_str + '\n'
            + 'Notes: ' + self.notes.strip()
            + ('\nDISCARDED' if self.discarded else '')
            )


def _create_new_history_dir(new_history_path):
    try:
        new_history_path.mkdir()
    except OSError:
        logger.error('Error: Could not create target directory '
                     f'{new_history_path}\n Stopping...')
        raise


def _discard_workhistory_previous(work_history_path):
    """Remove 'previous'-labelled directories in `work_history_path`."""
    work_hist_prev = _get_workhistory_directories(work_history_path,
                                                  contains=PREVIOUS_LABEL)
    for directory in work_hist_prev:
        try:
            shutil.rmtree(directory)
        except OSError:
            print(f'Failed to delete {directory} directory '
                  f'from {work_history_path}')


def _find_max_run_per_tensor(history_path):
    """Return maximum run numbers for all directories in `history_path`.

    Parameters
    ----------
    history_path : Path
        The path to the history folder in which
        tensor-run directories are looked up.

    Returns
    -------
    max_nums : defaultdict(int)
        The maximum run number currently stored in `history_path`
        for each of the tensors in `history_path`.
    """
    max_nums = defaultdict(int)  # max. job number per tensor number
    for directory in history_path.iterdir():
        match = HIST_FOLDER_RE.match(directory.name)
        if not directory.is_dir() or not match:
            continue
        tensor_num = int(match['tensor_num'])
        job_num = int(match['job_num'])
        max_nums[tensor_num] = max(max_nums[tensor_num], job_num)
    return max_nums


def _get_current_workhistory_directories(work_history_path, contains=''):
    """Return a generator of subfolders in work_history_path.

    As compared with _get_workhistory_directories, only those
    directories not marked as 'previous' are returned.

    Parameters
    ----------
    work_history_path : Path
        Path to the workhistory folder where
        subdirectories are looked up.
    contains : str, optional
        Select only those subdirectories whose name contains this
        string. Default is an empty string, corresponding to no
        filtering other than the one described in Returns.

    Returns
    -------
    subfolders : generator
        When iterated over, it yields paths to the immediate
        subdirectories of `work_history_path` whose name matches
        the HIST_FOLDER_RE regular expression, whose name
        include `contains`, but does not contain 'previous'.
    """
    directories = _get_workhistory_directories(work_history_path,
                                               contains=contains)
    return (d for d in directories if PREVIOUS_LABEL not in d.name)


def _get_workhistory_directories(work_history_path, contains=''):
    """Return a generator of subfolders in `work_history_path`.

    Parameters
    ----------
    work_history_path : Path
        Path to the workhistory folder where
        subdirectories are looked up.
    contains : str, optional
        Select only those subdirectories whose name contains this
        string. Default is an empty string, corresponding to no
        filtering other than the one described in Returns.

    Returns
    -------
    subfolders : generator
        When iterated over, it yields paths to the immediate
        subdirectories of work_history_path whose name matches
        the HIST_FOLDER_RE regular expression, and whose name
        include `contains`.
    """
    globbed = (work_history_path.glob(f'*{contains}*') if contains
               else work_history_path.iterdir())
    return (d for d in globbed if d.is_dir() and HIST_FOLDER_RE.match(d.name))


def _move_and_cleanup_workhistory(work_history_path,
                                  history_path,
                                  timestamp,
                                  max_job_for_tensor,
                                  discard):
    """Move files from `work_history_path` to `history_path`, then clean up.

    If `work_history_path` is empty, it is always deleted. Otherwise
    only if `discard` is True.

    Parameters
    ----------
    work_history_path : Path
        Path to the workhistory directory, moved in root by
        sections.cleanup, as it is included in manifest in
        move_oldruns.
    history_path : Path
        The target history folder where files should be copied.
    timestamp : str
        The time-stamp selecting the subfolders to be moved. Only
        those folders whose name contains timestamp are copied over.
    max_job_for_tensor : defaultdict(int)
        The maximum job number known for each tensor number. Used
        to generate new folder names (incrementing by one unit)
        in `history_path`.
    discard : bool
        Whether files should only be discarded without copying them.

    Returns
    -------
    tensor_nums : set of int
        Indices of tensors found in `work_history_path` that have
        been moved to `history_path` as new history entries.
    """
    tensor_nums = set()
    if not work_history_path.is_dir():
        return tensor_nums

    # Always remove any 'previous'-labelled folders
    _discard_workhistory_previous(work_history_path)
    if not discard:
        tensor_nums = _move_workhistory_folders(work_history_path,
                                                history_path,
                                                timestamp,
                                                max_job_for_tensor)
    is_empty = not any(work_history_path.iterdir())
    if is_empty or discard:
        try:
            shutil.rmtree(work_history_path)
        except OSError as exc:
            err_ = f'Failed to {{}} {work_history_path.name} folder: {exc}'
            print(err_.format('discard' if discard else 'delete empty'))
    return tensor_nums


def _move_or_discard_files(file_paths, target_folder, discard):
    """Move file_paths to target_folder or delete them."""
    for file in file_paths:
        _move_or_discard_one_file(file, target_folder, discard)


def _move_or_discard_one_file(file, target_folder, discard):
    """Move file to target_folder or delete it."""
    if not file.exists():
        return
    if discard and file.is_file():
        try:
            file.unlink()
        except OSError:
            print(f'Failed to discard file {file.name}.')
        return
    if discard:  # Should be a directory
        try:
            shutil.rmtree(file)
        except OSError:
            print(f'Failed to discard directory {file.name}.')
        return
    # Move it
    try:
        shutil.move(file, target_folder / file.name)
    except OSError:
        logger.error(f'Failed to move {file.name}.')



def _copy_one_file_to_history(file_path, history_path):
    """Copy file_path to history_path."""
    try:
        shutil.copy2(file_path, history_path / file_path.name)
    except OSError:
        logger.error(f'Failed to copy {file_path} to history.')



def _move_workhistory_folders(work_history_path, history_path,
                              timestamp, max_job_for_tensor):
    """Move timestamp-labelled folders in work_history_path to history_path.

    Parameters
    ----------
    work_history_path : Path
        Base path to the workhistory folder from
        which stuff should be moved.
    history_path : Path
        Destination path for the subfolders of
        `work_history_path` that will be moved.
    timestamp : str
        Selects which sub-folders of `work_history_path` should be
        moved. Only those whose name contains `timestamp` are moved.
    max_job_for_tensor : defaultdict(int)
        The maximum job number known for each tensor number. Used
        to generate new folder names (incrementing by one unit)
        in `history_path` for each of the directories moved from
        `work_history_path`.

    Returns
    -------
    tensor_nums : set
        The tensor numbers of the folders that have been moved.
    """
    tensor_nums = set()
    directories = _get_current_workhistory_directories(work_history_path,
                                                       contains=timestamp)
    for directory in directories:
        match = HIST_FOLDER_RE.match(directory.name)
        if not match:  # Should always match
            continue
        tensor_num = int(match['tensor_num'])
        try:
            search_num = int(directory.name[6:9])                               # TODO: is this right? 6:9 should be job_num
        except (ValueError, IndexError):
            continue

        tensor_nums.add(tensor_num)
        job_num = max_job_for_tensor[tensor_num] + 1
        newname = (
            f't{tensor_num:03d}.r{job_num:03d}.{search_num:03d}'
            + directory.name[9:]
            )
        try:
            shutil.move(directory, history_path / newname)
        except OSError:
            print(f'Error: Failed to move {directory}.')
    return tensor_nums


def _read_and_clear_notes_file(cwd):
    """Return notes read from file. Clear the file contents."""
    notes_path = next(cwd.glob('notes*'), None)
    if notes_path is None:
        return ''
    try:
        notes = notes_path.read_text(encoding='utf-8')
    except OSError:
        print(f'Error: Failed to read {notes_path.name} file.')
        return ''
    try:  # pylint: disable=too-many-try-statements
        with notes_path.open('w', encoding='utf-8'):
            pass
    except OSError:
        logger.error(f'Failed to clear the {notes_path.name} '
                    'file after reading.')
    return notes


def _read_most_recent_log(cwd):
    """Return timestamp and lines from the most-recent log file in cwd."""
    last_log_lines = []
    try:
        most_recent_log = max(
            (log_file
             for prefix in _CALC_LOG_PREFIXES
             for log_file in cwd.glob(f'{prefix}*.log')
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
    except OSError:
        pass
    return old_timestamp, last_log_lines


def _update_state_files_from_out(cwd):
    """Moved states files from OUT to root if available and rename old to _ori."""
    out_path = cwd / 'OUT'
    if not out_path.is_dir():
        return
    for file in STATE_FILES:
        out_file = out_path / f'{file}_OUT'
        if out_file.is_file():
            try:
                shutil.move(cwd / file, cwd / f'{file}_ori')
                shutil.move(out_file, cwd / file)
            except OSError:
                print(f'Error: failed to copy {out_file} as new {file}.')



# This could probably be done by datetime.strptime
def _translate_timestamp(time_stamp):
    """Return a 'DD.MM.YY hh:mm:ss' timestamp from a YYMMDD-hhmmss one."""
    time_stamp = time_stamp.replace('moved-', '')
    if len(time_stamp) != 13:
        print(time_stamp)
        raise ValueError('Error translating timestamp: Invalid length '
                         f'{len(time_stamp)}. Expected 13 characters.')
    year, month, day = time_stamp[:2], time_stamp[2:4], time_stamp[4:6]
    hour, minutes, secs = time_stamp[7:9], time_stamp[9:11], time_stamp[11:13]
    return f'{day}.{month}.{year} {hour}:{minutes}:{secs}'



class BookkeeperCLI(ViPErLEEDCLI, cli_name='bookkeeper'):
    """The main command-line interface for the bookkeeper utility."""

    def add_parser_arguments(self, parser):
        """Add bookkeeper arguments to parser."""
        super().add_parser_arguments(parser)
        what_next = parser.add_mutually_exclusive_group()
        what_next.add_argument(
            '-a', '--archive',
            help=('Store last run in history. Overwrite PARAMETERS, POSCAR &'
                  'VIBROCC from OUT. Runs after viperleed.calc by default.'),
            action='store_true',
            )
        what_next.add_argument(
            '-c', '--clear',
            help=('Clear the input directory and add last run to history if not'
                  ' already there. Runs before viperleed.calc by default.'),
            action='store_true',
            )
        what_next.add_argument(
            '-d', '--discard',
            help=('Discard all results from the last run, and restore the '
                  'previous inputs. The discarded run is kept in history.'),
            action='store_true',
            )
        what_next.add_argument(
            '-df', '--discard-full',
            help=('Discard all results from the last run as if it never '
                  'happened. The discarded run is removed from history.'),
            action='store_true',
            )
        parser.add_argument(
            '-j', '--job-name',
            help=('define a string to be appended to the name of the history '
                  f'folder that is created, and is logged in {HISTORY_INFO_NAME}'),
            type=str
            )
        parser.add_argument(
            '--history-name',
            help=('define the name of the history folder that is '
                  f'created/used. Default is {DEFAULT_HISTORY!r}'),
            type=str,
            default=DEFAULT_HISTORY
            )
        parser.add_argument(
            '--work-history-name',
            help=('define the name of the workhistory folder that is '
                  f'created/used. Default is {DEFAULT_WORK_HISTORY!r}'),
            type=str,
            default=DEFAULT_WORK_HISTORY
            )


    def __call__(self, args=None):
        """Call the bookkeeper with command-line args."""
        parsed_args = self.parse_cli_args(args)

        bookkeeper = Bookkeeper(job_name=parsed_args.job_name,
                                history_name=parsed_args.history_name,
                                work_history_name=parsed_args.work_history_name,
                                cwd=Path.cwd().resolve())
        # Select mode
        if parsed_args.clear:
            mode = BookkeeperMode.CLEAR
        elif parsed_args.discard:
            mode = BookkeeperMode.DISCARD
        elif parsed_args.discard_full:
            mode = BookkeeperMode.DISCARD_FULL
        else:
            mode = BookkeeperMode.ARCHIVE
        return bookkeeper.run(mode)


if __name__ == '__main__':
    BookkeeperCLI.run_as_script()
