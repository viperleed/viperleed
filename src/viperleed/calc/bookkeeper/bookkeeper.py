"""ViPErLEED bookkeeper module of package calc."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2020-01-30'
__license__ = 'GPLv3+'

from contextlib import nullcontext
from enum import IntEnum
from pathlib import Path

from viperleed.calc.constants import DEFAULT_OUT
from viperleed.calc.constants import DEFAULT_SUPP
from viperleed.calc.constants import ORIGINAL_INPUTS_DIR_NAME
from viperleed.calc.lib.log_utils import logging_silent
from viperleed.calc.lib.time_utils import DateTimeFormat
from viperleed.calc.sections.calc_section import ALL_INPUT_FILES

from . import log
from .errors import _FileNotOlderError
from .history.constants import HISTORY_INFO_NAME
from .history.entry.entry import HistoryInfoEntry
from .history.errors import CantDiscardEntryError
from .history.errors import CantRemoveEntryError
from .history.errors import MetadataMismatchError
from .history.errors import NoHistoryEntryError
from .history.meta import BookkeeperMetaFile
from .log import LOGGER
from .mode import BookkeeperMode
from .root_explorer import RootExplorer
from .utils import make_property


# Input files that may be generated at runtime - don't warn if missing
RUNTIME_GENERATED_INPUT_FILES = ('IVBEAMS', 'PHASESHIFTS')


# Suffix for files moved from root rather than original_inputs
_FROM_ROOT = '_from_root'


class BookkeeperExitCode(IntEnum):
    """Exit code of the bookkeeper."""
    SUCCESS = 0
    NOTHING_TO_DO = -1
    FAIL = 1


class Bookkeeper:
    """Bookkeeper to archive or discard the most recent viperleed calc run."""

    def __init__(self, cwd=Path.cwd()):
        """Initialize the bookkeeper using `cwd` as the root folder."""
        self._root = RootExplorer(path=cwd, bookkeeper=self)
        self._mode = None                        # Set in run
        self._requires_user_confirmation = None  # Set in run
        self._state_info = {
            'logger_prepared': False,
            # The next ones are set in update_from_cwd
            'tensor_number': None,  # int
            'timestamp': None,      # str
            }

    # Simple dynamic @properties
    cwd = make_property('_root.path')
    files_need_archiving = make_property('_root.needs_archiving')
    history = make_property('_root.history')
    max_job_for_tensor = make_property('history.max_run_per_tensor')
    timestamp = make_property('_state_info[timestamp]', needs_update=True)
    _workhistory = make_property('_root.workhistory')

    @property
    def archiving_required(self):
        """Check if archiving is required."""
        return self.files_need_archiving and not self.history.new_folder.exists

    @property
    def orig_inputs_dir(self):
        """Return the path to the folder containing untouched input files."""
        return self.cwd / DEFAULT_SUPP / ORIGINAL_INPUTS_DIR_NAME

    @property
    def tensor_number(self):
        """Return the index of the tensor currently handled."""
        tensor_number = self._state_info['tensor_number']
        return (tensor_number if tensor_number is not None
                else self._root.tensors.most_recent)

    def run(self, mode, requires_user_confirmation=True):
        """Run the bookkeeper in the given mode.

        Parameters
        ----------
        mode : str or BookkeeperMode
            Which bookkeeper mode to use. See help(BookkeeperMode).
        requires_user_confirmation : bool, optional
            Whether user confirmation is necessary before proceeding
            with destructive actions. Only used in DISCARD_FULL mode.
            Default is True.

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

        LOGGER.info(f'\nRunning bookkeeper in {mode.name} mode in {self.cwd}.')
        self._mode = mode
        self._requires_user_confirmation = requires_user_confirmation
        # Do not bother logging messages when the user asked to
        # --fix: some warnings will be fixed anyway. Others will
        # re-appear at the next run.
        self.update_from_cwd(silent=mode is BookkeeperMode.FIX)
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

            # Collect information from the root directory:
            # log files, SUPP, OUT, workhistory, and Tensors, as well
            # as information from the history folder (especially the
            # highest run number currently stored for each tensor).
            self._root.collect_info()
            self._state_info['timestamp'] = (
                self._root.calc_timestamp  # If a log file was found
                or f'moved-{DateTimeFormat.FILE_SUFFIX.now()}'
                )

            # Infer which history subfolder we should be handling
            self.history.find_new_history_directory(self.tensor_number,
                                                    suffix=self.timestamp)

            # history.info handler (may create history.info file)
            self.history.prepare_info_file()
            self.history.info.read()

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
        job_per_tensor = self.max_job_for_tensor
        with logging_silent():
            # It's OK to silence the logger here, as we will surely
            # replace the timestamp format in append_entry. In fact,
            # self.timestamp is just the one derived from the log file
            # and needs to be updated to a writable one. When we will
            # replace the format, logging messages for other faulty
            # stuff will pop up.
            new_info_entry = HistoryInfoEntry(
                tensor_nums=sorted_tensors,
                job_nums=[job_per_tensor[t] for t in sorted_tensors],
                timestamp=self.timestamp,
                folder_name=self.history.new_folder.name,
                notes=self._root.read_and_clear_notes_file(),
                discarded_=self._mode is BookkeeperMode.DISCARD,
                # Optional ones
                **self._root.infer_run_info()
                )
        self.history.info.append_entry(new_info_entry, fix_time_format=True)

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
        tensor_nums = self._workhistory.move_current_and_cleanup(metadata)
        tensor_nums.add(self.tensor_number)
        # ...and add a history.info entry
        self._add_history_info_entry(tensor_nums)
        LOGGER.info('Done archiving the current directory '
                    f'to {self.history.path.name}')

    def _check_may_discard_full(self):
        """Log and raise if it is not possible to DISCARD_FULL."""
        if self.archiving_required:
            LOGGER.error(f'Cannot {self._mode.long_flag} when the root '
                         'directory has not been archived yet. Run '
                         f'bookkeeper {BookkeeperMode.ARCHIVE.long_flag} '
                         'before.')
            raise CantRemoveEntryError
        try:
            self.history.info.may_remove_last_entry()
        except NoHistoryEntryError as exc:
            LOGGER.warning('Error: Failed to remove last entry '
                           f'from {HISTORY_INFO_NAME}: {exc}')
            raise
        except CantRemoveEntryError as exc:
            LOGGER.error(str(exc))
            raise

        # Make sure that we're going to remove consistent
        # stuff from history and history.info
        try:
            self.history.check_last_folder_consistent()
        except FileNotFoundError:
            LOGGER.error(f'{self._mode.name} mode failed: could not identify '
                         'directory to remove. Please proceed manually.')
            raise
        except MetadataMismatchError as exc:
            LOGGER.error(f'Error: {exc} Please proceed manually.')
            raise
        except CantRemoveEntryError as exc:
            LOGGER.error(
                f'Error: the most recent {self.history.path.name} folder is '
                f'inconsistent with the last entry in {HISTORY_INFO_NAME}. '
                f'{exc} Please proceed manually.'
                )
            raise

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
        history_folder = self.history.new_folder
        for file in ALL_INPUT_FILES:
            original_file = self.orig_inputs_dir / file
            cwd_file = self.cwd / file
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
                # by calc to original_inputs, and they should be
                # already handled by "elif original_file.is_file()".
                continue
            elif cwd_file.is_file():  # Copy cwd and warn
                copy_file, with_name = cwd_file, f'{cwd_file.name}{_FROM_ROOT}'
                LOGGER.warning(
                    f'File {file} not found in {ORIGINAL_INPUTS_DIR_NAME}. '
                    'Using file from root directory instead and renaming to '
                    f'{with_name}.'
                    )
            if copy_file:
                history_folder.copy_file_or_directory(copy_file,
                                                      with_name=with_name)

    def _copy_log_files(self):
        """Copy all the log files in root to history."""
        for file in self._root.logs.files:
            self.history.new_folder.copy_file_or_directory(file)

    def _copy_out_and_supp(self):
        """Copy OUT and SUPP directories to history."""
        history_folder = self.history.new_folder
        for name in (DEFAULT_SUPP, DEFAULT_OUT):
            directory = self.cwd / name
            if not directory.is_dir():
                LOGGER.warning(f'Could not find {name} directory in '
                               f'{self.cwd}. It will not be copied '
                               'to history.')
                continue
            history_folder.copy_file_or_directory(directory)

    def _do_prerun_archiving_and_mark_edited(self):
        """If needed, archive files from root to history.

        If any archiving is needed, input files in the current
        working directory that have been modified by the user
        since calc started and before bookkeeper ran are marked
        as '_edited'.

        Returns
        -------
        did_archive : bool
            Whether anything was archived at all.
        """
        in_archive_mode = self._mode is BookkeeperMode.ARCHIVE
        if in_archive_mode and self.history.new_folder.exists:
            LOGGER.info(
                f'History directory for run {self.history.new_folder.name} '
                'exists. Exiting without doing anything.'
                )
            return False
        if in_archive_mode and not self.files_need_archiving:
            LOGGER.info(f'No files to be moved to {self.history.path.name}. '
                        'Exiting without doing anything.')
            return False
        if not self.archiving_required:
            return False
        if not in_archive_mode:
            LOGGER.info(f'History folder {self.history.new_folder.name} '
                        'does not exist yet. Running archive mode first.')
        self._archive_to_history_and_add_info_entry()

        # Notice that it is OK to decide now whether files should
        # be marked as edited rather than before archiving. In
        # fact, during archiving we recognize the right files to
        # pull (in _copy_input_files_from_original_inputs_or_cwd):
        # always copy from original_inputs if possible, otherwise
        # copy with a _FROM_ROOT suffix.
        self._root.mark_edited_files()
        return True

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
            self.history.new_folder.path.mkdir()
        except OSError:
            LOGGER.error('Error: Could not create target directory '
                         f'{self.history.new_folder.name}\n Stopping...')
            raise
        LOGGER.info(f'Created history folder {self.history.new_folder.name} '
                    'for storing results of the most-recent viperleed.calc '
                    'execution')
        self._copy_out_and_supp()
        self._copy_input_files_from_original_inputs_or_cwd()
        self._copy_log_files()

        # Now that all files are in place, add the metadata file
        meta = BookkeeperMetaFile(self.history.new_folder.path)
        meta.compute_hash()
        meta.write()

        # Finally, register the new folder in history
        self.history.register_folder(self.history.new_folder.path)
        return meta

    def _make_history_and_prepare_logger(self):
        """Make history folder and add handlers to the bookkeeper logger."""
        if self._state_info['logger_prepared']:
            return

        # Attach a stream handler to logger if not already present
        log.ensure_has_stream_handler()

        try:  # Make top level history folder if not there yet
            self.history.path.mkdir(exist_ok=True)
        except OSError:
            LOGGER.error('Error creating history folder.')
            raise

        # Attach file handler for history/bookkeeper.log
        log.add_bookkeeper_logfile(self.history.path)
        LOGGER.info(  # Log only once per instance
            '\n### Bookeeper running at '
            f'{DateTimeFormat.LOG_CONTENTS.now()} ###'
            )
        self._state_info['logger_prepared'] = True

    def _print_discard_info(self):
        """Emit logging messages with files/folders/entry to be deleted."""
        LOGGER.info(f'About to delete files/folders from {self.cwd}')
        # First, files and folders that are simply deleted
        paths_to_discard = (
            *self._workhistory.list_paths_to_discard(),
            *self.history.list_paths_to_discard(),
            *self._root.list_paths_to_discard(),
            )
        if paths_to_discard:
            LOGGER.info('The following files/folders will '
                        'be permanently deleted:')
        for path in paths_to_discard:
            LOGGER.info(f'    {path.relative_to(self.cwd)}')

        # Then files that are replaced by others
        files_to_replace = self._root.list_files_to_replace()
        if files_to_replace:
            LOGGER.info('The following files will be permanently replaced:')
        for deleted, renamed in files_to_replace:
            LOGGER.info(f'    replace {deleted.relative_to(self.cwd)} '
                        f'with {renamed.relative_to(self.cwd)}')

        # Finally, the history.info entry
        LOGGER.info('The following entry will be deleted '
                    f'from {self.history.info.path.name}:'
                    f'{self.history.info.last_entry}'.rstrip())
        LOGGER.info('Calculation data will be lost irreversibly.')

    def _run_archive_mode(self):
        """Execute bookkeeper in ARCHIVE mode."""
        # Copy all the relevant files to history, add
        # a history.info entry, mark edited input files...
        did_archive = self._do_prerun_archiving_and_mark_edited()
        if not did_archive:
            return BookkeeperExitCode.NOTHING_TO_DO

        # ...mark non-edited inputs as _ori, and prepare inputs for
        # the next execution of calc (from OUT or original_inputs)
        try:
            self._root.prepare_for_next_calc_run()
        except OSError:
            return BookkeeperExitCode.FAIL
        return BookkeeperExitCode.SUCCESS

    def _run_clear_mode(self):
        """Execute bookkeeper in CLEAR mode."""
        self._do_prerun_archiving_and_mark_edited()
        self._root.clear_for_next_calc_run()
        return BookkeeperExitCode.SUCCESS

    def _run_discard_full_mode(self):
        """Execute bookkeeper in DISCARD_FULL mode."""
        try:
            self._check_may_discard_full()
        except NoHistoryEntryError:
            return BookkeeperExitCode.NOTHING_TO_DO
        except (CantRemoveEntryError,
                FileNotFoundError,
                MetadataMismatchError):
            return BookkeeperExitCode.FAIL

        self._print_discard_info()
        if not self._user_confirmed():
            return BookkeeperExitCode.NOTHING_TO_DO

        # Delete the history folders, stuff in workhistory,
        # and output files/folders in the root directory
        try:
            self.history.discard_most_recent_run()
        except OSError:
            return BookkeeperExitCode.FAIL
        self._workhistory.discard_workhistory_root()
        self._root.revert_to_previous_calc_run()

        # Tensors and Deltas, if created during the last run
        self._root.remove_tensors_and_deltas()

        # And the history entry from history.info
        self.history.info.remove_last_entry()
        return BookkeeperExitCode.SUCCESS

    def _run_discard_mode(self):
        """Execute bookkeeper in DISCARD mode."""
        did_archive = self._do_prerun_archiving_and_mark_edited()
        self._root.revert_to_previous_calc_run()
        if did_archive:
            # We had to archive the contents of the root directory.
            # This means that we have already marked the entry as
            # discarded and we don't have to bother.
            return BookkeeperExitCode.SUCCESS
        try:
            self.history.info.discard_last_entry()
        except (NoHistoryEntryError, CantDiscardEntryError) as exc:
            no_entry = isinstance(exc, NoHistoryEntryError)
            emit_log = LOGGER.error if no_entry else LOGGER.warning
            emit_log('Failed to mark last entry as discarded '
                     f'in {HISTORY_INFO_NAME}: {exc}')
            return (BookkeeperExitCode.NOTHING_TO_DO if no_entry
                    else BookkeeperExitCode.FAIL)
        return BookkeeperExitCode.SUCCESS

    def _run_fix_mode(self):
        """Fix format inconsistencies found in history and history.info."""
        self.history.fix()
        return BookkeeperExitCode.SUCCESS

    def _user_confirmed(self):
        """Return whether the user wants to proceed with discarding."""
        if not self._requires_user_confirmation:
            return True
        while True:
            reply = input('Are you sure you want to proceed (y/N)?')
            reply = reply.lower()
            if not reply or reply.startswith('n'):
                return False
            if reply.startswith('y'):
                return True


def _check_newer(older, newer):
    """Raise if file `older` is newer than file `newer`."""
    newer_timestamp = newer.stat().st_mtime
    older_timestamp = older.stat().st_mtime
    if newer_timestamp < older_timestamp:
        raise _FileNotOlderError
