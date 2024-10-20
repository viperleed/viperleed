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
from pathlib import Path
import shutil

from viperleed.calc import DEFAULT_HISTORY
from viperleed.calc import ORIGINAL_INPUTS_DIR_NAME
from viperleed.calc.lib.leedbase import getMaxTensorIndex
from viperleed.calc.lib.log_utils import logging_silent
from viperleed.calc.lib.time_utils import DateTimeFormat
from viperleed.calc.sections.calc_section import ALL_INPUT_FILES
from viperleed.calc.sections.cleanup import DEFAULT_OUT
from viperleed.calc.sections.cleanup import DEFAULT_SUPP

from . import log
from .constants import HISTORY_FOLDER_RE
from .constants import STATE_FILES
from .history.constants import HISTORY_INFO_NAME
from .history.entry.entry import HistoryInfoEntry
from .history.errors import CantDiscardEntryError
from .history.errors import CantRemoveEntryError
from .history.errors import NoHistoryEntryError
from .history.info import HistoryInfoFile
from .history.meta import BookkeeperMetaFile
from .log import LOGGER
from .mode import BookkeeperMode
from .root_explorer import RootExplorer
from .utils import discard_files
from .utils import make_property
from .utils import needs_update_for_attr


# Input files that may be generated at runtime - don't warn if missing
RUNTIME_GENERATED_INPUT_FILES = ('IVBEAMS', 'PHASESHIFTS')


# Suffix for files moved from root rather than original_inputs
_FROM_ROOT = '_from_root'


class _FileNotOlderError(Exception):
    """Exception used internally for file age checks."""


class BookkeeperExitCode(IntEnum):
    """Exit code of the bookkeeper."""
    SUCCESS = 0
    NOTHING_TO_DO = -1
    FAIL = 1


class Bookkeeper:
    """Bookkeeper to archive or discard the most recent viperleed calc run."""

    def __init__(self, cwd=Path.cwd()):
        """Initialize the bookkeeper using `cwd` as the root folder."""
        self._folder_names = {
            'history_dir_base_name': None,   # Set in update_from_cwd
            'history_dir': None,             # Set in update_from_cwd
            }
        self._history_info = None   # Set in update_from_cwd
        self._mode = None           # Set in run
        self._root = RootExplorer(path=cwd, bookkeeper=self)
        self._state_info = {
            'logger_prepared': False,
            # The next ones are set in update_from_cwd
            'max_job_for_tensor': None,   # defaultdict[int]
            'tensor_number': None,        # int
            'timestamp': None,            # str
            'history_with_same_base_name_exists': None,  # bool
            }

    # Simple dynamic @properties
    cwd = make_property('_root.path')
    files_need_archiving = make_property('_root.needs_archiving')
    history_dir_base_name = make_property(
        # Notice that this is not necessarily the name of the
        # history subfolder that is created. Use .history_dir
        # for that.
        '_folder_names[history_dir_base_name]',
        needs_update=True,
        )
    history_info = make_property('_history_info', needs_update=True)
    history_with_same_base_name_exists = make_property(
        '_state_info[history_with_same_base_name_exists]',
        needs_update=True,
        )
    max_job_for_tensor = make_property('_state_info[max_job_for_tensor]',
                                       needs_update=True)
    tensor_number = make_property('_state_info[tensor_number]',
                                  needs_update=True)
    timestamp = make_property('_state_info[timestamp]', needs_update=True)
    _workhistory = make_property('_root.workhistory')

    @property
    def archiving_required(self):
        """Check if archiving is required."""
        return (self.files_need_archiving
                and not self.history_with_same_base_name_exists)

    @property
    @needs_update_for_attr('_folder_names[history_dir]')
    def history_dir(self):
        """Return the path to the history subfolder that we deal with.

        This is the directory to which files will be moved upon run().

        Returns
        -------
        Path
        """
        return self.top_level_history_path / self._folder_names['history_dir']

    @property
    def orig_inputs_dir(self):
        """Return the path to the folder containing untouched input files."""
        return self.cwd / DEFAULT_SUPP / ORIGINAL_INPUTS_DIR_NAME

    @property
    def top_level_history_path(self):
        """Return the path to the 'history' folder."""
        return self.cwd / DEFAULT_HISTORY

    def run(self, mode):
        """Run the bookkeeper in the given mode.

        Parameters
        ----------
        mode : str or BookkeeperMode
            Which bookkeeper mode to use. See help(BookkeeperMode).

        Returns
        -------
        exit_code : BookkeeperExitCode

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
        except ValueError as exc:
            raise ValueError(f'Unknown mode {mode}') from exc

        try:
            runner = getattr(self, f'_run_{mode.name.lower()}_mode')
        except AttributeError as exc:
            raise NotImplementedError from exc

        LOGGER.info(f'Running bookkeeper in {mode.name} mode.')
        self._mode = mode
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

            # Collect information from the root directory:
            # log files, SUPP, OUT, and workhistory.
            self._root.collect_info()
            self._state_info['timestamp'] = (
                self._root.calc_timestamp  # If a log file was found
                or f'moved-{DateTimeFormat.FILE_SUFFIX.now()}'
                )

            # Infer which history subfolder we should be handling
            self._find_base_name_for_history_subfolder()
            self._find_new_history_directory_name()

            # history.info handler (may create history.info file)
            self._history_info = HistoryInfoFile(self.cwd, create_new=True)
            self.history_info.read()

    def _add_history_info_entry(self, tensor_nums):
        """Add a new entry to the history.info file.

        Parameters
        ----------
        tensor_nums : iterable
            Progressive identification numbers of the Tensor files
            for which new history subfolders were created during
            this bookkeeper run.

        Returns
        -------
        None.
        """
        sorted_tensors = sorted(tensor_nums)
        with logging_silent():
            # It's OK to silence the logger here, as we will surely
            # replace the timestamp format in append_entry. In fact,
            # self.timestamp is just the one derived from the log file
            # and needs to be updated to a writable one. When we will
            # replace the format, logging messages for other faulty
            # stuff will pop up.
            new_info_entry = HistoryInfoEntry(
                tensor_nums=sorted_tensors,
                job_nums=[self.max_job_for_tensor[tensor] + 1
                          for tensor in sorted_tensors],
                timestamp=self.timestamp,
                folder_name=self.history_dir.name,
                notes=self._root.read_and_clear_notes_file(),
                discarded_=self._mode is BookkeeperMode.DISCARD,
                # Optional ones
                **self._root.infer_run_info()
                )
        self.history_info.append_entry(new_info_entry, fix_time_format=True)

    def _archive_to_history_and_add_info_entry(self):
        """Store all relevant files in root to history, add a new entry.

        This method does not check whether any archiving is indeed
        necessary. Make sure to check if self.archiving_required
        beforehand. After this method executes, the root folder
        still contains the output produced by viperleed.calc,
        except for those workhistory folders that were archived
        to history. The workhistory folder itself is cleaned-up
        of "previous" runs, and may not exist anymore if it was
        empty a the end of the archiving. A new entry is added
        to the history.info file.

        Returns
        -------
        None.
        """
        # Create the new 'primary' history directory...
        metadata = self._make_and_copy_to_history()
        # ...move workhistory folders...
        tensor_nums = self._workhistory.move_and_cleanup(metadata)
        tensor_nums.add(self.tensor_number)
        # ...and add a history.info entry
        self._add_history_info_entry(tensor_nums)

    def _clean_state(self):
        """Rebuild (almost) the same state as after __init__.

        The only surviving attributes are the ones given
        at __init__ and the state that concerns logging.

        Returns
        -------
        None.
        """
        logger_prepared = self._state_info['logger_prepared']

        # Note on the disable: while it is indeed not very elegant to
        # call __init__ again, it is the simplest way to (i) make it
        # clear in __init__ which attributes we have, and (ii) avoid
        # code repetition here. The alternative would be to set all
        # attributes and keys to None in __init__, and call another
        # method both in __init__ and here. This seems a lot of
        # complication just to reset (almost) everything to None.
        # pylint: disable-next=unnecessary-dunder-call
        self.__init__(cwd=self.cwd)
        self._state_info['logger_prepared'] = logger_prepared

    def _collect_max_run_per_tensor(self):
        """Find the maximum run numbers for all directories in history."""
        history_path = self.top_level_history_path
        max_jobs = defaultdict(int)  # Max. job nr. per tensor number
        for directory in history_path.iterdir():
            match = HISTORY_FOLDER_RE.match(directory.name)
            if not directory.is_dir() or not match:
                continue
            tensor_num = int(match['tensor_num'])
            job_num = int(match['job_num'])
            max_jobs[tensor_num] = max(max_jobs[tensor_num], job_num)
        self._state_info['max_job_for_tensor'] = max_jobs

    @needs_update_for_attr('_mode', updater='run')
    def _copy_input_files_from_original_inputs_or_cwd(self):
        """Copy input files to the history subfolder.

        The files are preferentially collected from the original_inputs
        subfolder of SUPP. Only if they are not found there, they are
        taken from the root directory (unless they are auto-generated
        ones). If the root directory contains a newer version of the
        file than in original_inputs, warn the user (but still copy
        the one from original_inputs).

        Returns
        -------
        None.
        """
        use_ori = self._mode.uses_ori_files_as_fallback
        for file in ALL_INPUT_FILES:
            original_file = self.orig_inputs_dir / file
            cwd_file = self.cwd / file
            if file in STATE_FILES and use_ori:
                cwd_file = self.cwd / f'{file}_ori'

            copy_file, with_name = None, file
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
            elif cwd_file.is_file():  # Copy cwd and warn
                copy_file, with_name = cwd_file, f'{cwd_file.name}{_FROM_ROOT}'
                LOGGER.warning(
                    f'File {file} not found in {ORIGINAL_INPUTS_DIR_NAME}. '
                    'Using file from root directory instead and renaming to '
                    f'{with_name}.'
                    )
            if copy_file:
                self._copy_one_file_to_history(copy_file, with_name=with_name)

    def _copy_log_files(self):
        """Copy all the log files in root to history."""
        for file in self._root.logs.files:
            self._copy_one_file_to_history(file)

    @needs_update_for_attr('_folder_names[history_dir]')
    def _copy_one_file_to_history(self, file_path, with_name=None):
        """Copy `file_path` to history, optionally renaming."""
        dest_name = with_name or file_path.name
        try:
            shutil.copy2(file_path, self.history_dir / dest_name)
        except OSError:
            LOGGER.error(f'Failed to copy {file_path} to history.')

    @needs_update_for_attr('_folder_names[history_dir]')
    def _copy_out_and_supp(self):
        """Copy OUT and SUPP directories to history."""
        for name in (DEFAULT_SUPP, DEFAULT_OUT):
            directory = self.cwd / name
            if not directory.is_dir():
                LOGGER.warning(f'Could not find {name} directory in '
                               f'{self.cwd}. It will not be copied '
                               'to history.')
                continue
            try:
                shutil.copytree(directory, self.history_dir / name)
            except OSError:
                LOGGER.error(f'Failed to copy {name} directory to history.')

    def _do_prerun_archiving(self):
        """If needed, archive files from root to history.

        This is a simple wrapper around the
        _archive_to_history_and_add_info_entry
        method that does nothing if no archiving
        is needed, logs an info message otherwise.

        Returns
        -------
        did_archive : bool
            Whether anything was archived at all.
        """
        if not self.archiving_required:
            return False
        LOGGER.info(f'History folder {self.history_dir} does not '
                    'exist yet. Running archive mode first.')
        self._archive_to_history_and_add_info_entry()
        return True

    def _find_base_name_for_history_subfolder(self):
        """Store internally the potential name of the history subfolder.

        This is the base name for the folder that will contain the
        final output and SUPP files. Folders coming from the work-
        history folder have their own naming. Notice also that this
        is not necessarily the name of the history directory that is
        actually created. Use self.history_dir.name instead. See also
        help(self._find_new_history_directory_name).

        Returns
        -------
        None.
        """
        suffix = self.timestamp
        job_number = self.max_job_for_tensor[self.tensor_number]
        dir_name_fmt = f't{self.tensor_number:03d}.r{{job:03d}}_{suffix}'
        # If there is already a folder with the same name and correct
        # timestamp/suffix, we take that, otherwise, we increase the
        # job number: it's a new run.
        base_name = dir_name_fmt.format(job=job_number)
        if not (self.top_level_history_path / base_name).is_dir():
            base_name = dir_name_fmt.format(job=job_number + 1)
        self._folder_names['history_dir_base_name'] = base_name

        # Now that we have a base name, we can also check if there
        # are already history folders with the same base name.
        self._state_info['history_with_same_base_name_exists'] = any(
            self._get_conflicting_history_subfolders()
            )

    def _find_new_history_directory_name(self):
        """Store internally the name of the history directory to be used.

        The path to the directory is then available via self.history_dir.

        Folder name has form 'tXXX.rYYY_<suffix>'. <suffix> can vary.
        It may be:
        - if there was a log file, and the folder is not
          present in history:
              <log_timestamp>
        - if there was a log file, but there is already a folder
          in history with the same timestamp (typically that means
          that bookkeeper was called multiple times on the same run):
              <log_timestamp>_moved-<bookie_timestamp>
        - if there was no log file:
              moved-<bookie_timestamp>
        - in the unlikely event that the folder in the previous point
          already exists (this would mean that the bookkeeper is called
          twice within one second):
              moved-<bookie_timestamp>_moved-<bookie_timestamp2>
          where <bookie_timestamp2> may be slightly later than
          <bookie_timestamp>, but is most likely the same.

        Returns
        -------
        None.
        """
        dir_name = self.history_dir_base_name
        if (self.top_level_history_path / dir_name).is_dir():
            bookkeeper_timestamp = DateTimeFormat.FILE_SUFFIX.now()
            dir_name = f'{dir_name}_moved-{bookkeeper_timestamp}'
        self._folder_names['history_dir'] = dir_name

    def _get_conflicting_history_subfolders(self):
        """Return an iterator of history subfolders with the same base name."""
        folders = self.top_level_history_path.glob(
            f'{self.history_dir_base_name}*'
            )
        return (f for f in folders if f.is_dir())

    def _make_and_copy_to_history(self):
        """Create new history subfolder and copy all files there.

        Returns
        -------
        meta : BookkeeperMetaFile
            The handler to the metadata file created with the
            new history subfolder.

        Raises
        ------
        OSError
            If creating the new history subfolder fails.
        """
        try:
            self.history_dir.mkdir()
        except OSError:
            LOGGER.error('Error: Could not create target directory '
                         f'{self.history_dir}\n Stopping...')
            raise
        self._copy_out_and_supp()
        self._copy_input_files_from_original_inputs_or_cwd()
        self._copy_log_files()

        # Now that all files are in place, add the metadata file
        meta = BookkeeperMetaFile(self.history_dir)
        meta.compute_hash()
        meta.write()
        return meta

    def _make_history_and_prepare_logger(self):
        """Make history folder and add handlers to the bookkeeper logger."""
        if self._state_info['logger_prepared']:
            return

        # Attach a stream handler to logger if not already present
        log.ensure_has_stream_handler()

        try:  # Make top level history folder if not there yet
            self.top_level_history_path.mkdir(exist_ok=True)
        except OSError:
            LOGGER.error('Error creating history folder.')
            raise

        # Attach file handler for history/bookkeeper.log
        log.add_bookkeeper_logfile(self.top_level_history_path)
        LOGGER.info(  # Log only once per instance
            '\n### Bookeeper running at '
            f'{DateTimeFormat.LOG_CONTENTS.now()} ###'
            )
        self._state_info['logger_prepared'] = True

    def _remove_tensors_and_deltas(self):
        """Delete the most recent tensor and delta files."""
        tensor_file = f'Tensors/Tensors_{self.tensor_number:03d}.zip'
        delta_file = f'Deltas/Deltas_{self.tensor_number:03d}.zip'
        discard_files(self.cwd / tensor_file, self.cwd / delta_file)

    def _run_archive_mode(self):
        """Execute bookkeeper in ARCHIVE mode."""
        if self.history_with_same_base_name_exists:
            LOGGER.info(
                f'History directory for run {self.history_dir_base_name} '
                'exists. Exiting without doing anything.'
                )
            return BookkeeperExitCode.NOTHING_TO_DO
        if not self.files_need_archiving:
            LOGGER.info('No files to be moved to history. Exiting '
                        'without doing anything.')
            return BookkeeperExitCode.NOTHING_TO_DO

        # Copy all relevant files to history, add an entry...
        self._archive_to_history_and_add_info_entry()
        # ...and move old inputs to _ori, replacing them from OUT
        self._root.update_state_files_from_out()
        return BookkeeperExitCode.SUCCESS

    def _run_clear_mode(self):
        """Execute bookkeeper in CLEAR mode."""
        self._do_prerun_archiving()
        self._root.clear_for_next_calc_run()
        return BookkeeperExitCode.SUCCESS

    def _run_discard_full_mode(self):
        """Execute bookkeeper in DISCARD_FULL mode."""
        try:
            self.history_info.may_remove_last_entry()
        except NoHistoryEntryError as exc:
            LOGGER.warning('Error: Failed to remove last entry '
                           f'from {HISTORY_INFO_NAME}: {exc}')
            return BookkeeperExitCode.NOTHING_TO_DO
        except CantRemoveEntryError as exc:
            LOGGER.error(str(exc))
            return BookkeeperExitCode.FAIL

        # The directory we want to remove is not self.history_dir (since that   # TODO: don't we also want to purge all the intermediate ones from workhistory?
        # would be _moved-<timestamp>), but the one with the same base name!
        dir_to_remove = (self.top_level_history_path
                         / self.history_dir_base_name)
        if not dir_to_remove.is_dir():
            LOGGER.error('FULL_DISCARD mode failed: could not identify '
                         'directory to remove. Please proceed manually.')
            return BookkeeperExitCode.FAIL

        # Clean up workhistory in root
        self._workhistory.discard_workhistory_root()

        # Remove history folder
        try:
            shutil.rmtree(dir_to_remove)
        except OSError:
            LOGGER.error(f'Error: Failed to delete {dir_to_remove}.')
            return BookkeeperExitCode.FAIL
        self._root.revert_to_previous_calc_run()

        # Tensors and Deltas
        if self.tensor_number not in self.max_job_for_tensor:
            # Tensor is new.
            self._remove_tensors_and_deltas()

        # And the history entry from history.info
        self.history_info.remove_last_entry()
        return BookkeeperExitCode.SUCCESS

    def _run_discard_mode(self):
        """Execute bookkeeper in DISCARD mode."""
        did_archive = self._do_prerun_archiving()
        self._root.revert_to_previous_calc_run()
        if did_archive:
            # We had to archive the contents of the root directory.
            # This means that we have already marked the entry as
            # discarded and we don't have to bother.
            return BookkeeperExitCode.SUCCESS
        try:
            self.history_info.discard_last_entry()
        except (NoHistoryEntryError, CantDiscardEntryError) as exc:
            emit_log = (LOGGER.error if isinstance(exc, NoHistoryEntryError)
                        else LOGGER.warning)
            emit_log('Failed to mark last entry as discarded '
                     f'in {HISTORY_INFO_NAME}: {exc}')
            return BookkeeperExitCode.FAIL
        return BookkeeperExitCode.SUCCESS


def _check_newer(older, newer):
    """Raise if file `older` is newer than file `newer`."""
    newer_timestamp = newer.stat().st_mtime
    older_timestamp = older.stat().st_mtime
    if newer_timestamp < older_timestamp:
        raise _FileNotOlderError
