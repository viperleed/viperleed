"""Module bookkeeper of viperleed.calc.bookkeeper."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2020-01-30'
__license__ = 'GPLv3+'

from contextlib import nullcontext
import logging
from pathlib import Path

from viperleed import __version__
from viperleed.calc.bookkeeper import log
from viperleed.calc.bookkeeper.domain_finder import DomainFinder
from viperleed.calc.bookkeeper.domain_finder import MainPathNotFoundError
from viperleed.calc.bookkeeper.errors import _FileNotOlderError
from viperleed.calc.bookkeeper.errors import BookkeeperUnexpectedError
from viperleed.calc.bookkeeper.errors import NotAnInteractiveShellError
from viperleed.calc.bookkeeper.exit_code import BookkeeperExitCode
from viperleed.calc.bookkeeper.history.constants import HISTORY_INFO_NAME
from viperleed.calc.bookkeeper.history.entry.entry import HistoryInfoEntry
from viperleed.calc.bookkeeper.history.errors import CantDiscardEntryError
from viperleed.calc.bookkeeper.history.errors import CantRemoveEntryError
from viperleed.calc.bookkeeper.history.errors import MetadataError
from viperleed.calc.bookkeeper.history.errors import MetadataMismatchError
from viperleed.calc.bookkeeper.history.errors import NoHistoryEntryError
from viperleed.calc.bookkeeper.history.meta import BookkeeperMetaFile
from viperleed.calc.bookkeeper.mode import BookkeeperMode
from viperleed.calc.bookkeeper.root_explorer import DomainRootExplorer
from viperleed.calc.bookkeeper.root_explorer import RootExplorer
from viperleed.calc.bookkeeper.utils import ask_user_confirmation
from viperleed.calc.bookkeeper.utils import file_contents_identical
from viperleed.calc.bookkeeper.utils import make_property
from viperleed.calc.constants import DEFAULT_OUT
from viperleed.calc.constants import DEFAULT_SUPP
from viperleed.calc.lib.log_utils import logging_silent
from viperleed.calc.lib.string_utils import harvard_commas
from viperleed.calc.lib.time_utils import DateTimeFormat
from viperleed.calc.sections.calc_section import ALL_INPUT_FILES
from viperleed.calc.sections.cleanup import MOVED_LABEL


LOGGER = logging.getLogger(__name__)


# Input files that may be generated at runtime - don't warn if missing
RUNTIME_GENERATED_INPUT_FILES = ('IVBEAMS', 'PHASESHIFTS')


# Suffix for files moved from root rather than original_inputs
_FROM_ROOT = '_from_root'


# Log message emitted when we detect that we're running
# in the tree created by a pre-0.13.0 calc version.
_MIN_CALC_WARN = '0.13.0'
_WARN_OLD_CALC = f'''\
Bookkeeper v{__version__} is running in a folder created by an older version of
viperleed (v%s). Please double-check that the correct files have been
processed, as bookkeeper is not entirely backward compatible. See the
documentation at https://viperleed.org/content/calc/sections/bookkeeper.html
for details on the changes introduced in bookkeeper since v{_MIN_CALC_WARN}.'''


# TODO: this module is pretty long with more than 1000 lines. It would
# perhaps be better to split off all the methods that concern archiving
# the contents of current directory into a dedicated class (e.g.,
# RootArchiver). This way, Bookkeeper would only concern itself
# with what to do in which mode. Portions of this class that would go
# there:
# _make_and_copy_to_history
# _archive_*
class Bookkeeper:
    """Bookkeeper to archive or discard the most recent viperleed calc run."""

    def __init__(self, cwd=None):
        """Initialize the bookkeeper using `cwd` as the root folder."""
        self._root = RootExplorer(path=cwd or Path.cwd(),
                                  bookkeeper=self)
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
        return self.files_need_archiving and not self._has_archived_folder

    @property
    def orig_inputs_dir(self):
        """Return the path to the folder containing untouched input files."""
        return self._root.orig_inputs_dir

    @property
    def tensor_number(self):
        """Return the index of the tensor currently handled."""
        tensor_number = self._state_info['tensor_number']
        return (tensor_number if tensor_number is not None
                else self._root.tensors.most_recent)

    def run(self, mode, requires_user_confirmation=True, domains=None):
        """Run the bookkeeper in the given mode.

        Parameters
        ----------
        mode : str or BookkeeperMode
            Which bookkeeper mode to use. See help(BookkeeperMode).
        requires_user_confirmation : bool, optional
            Whether user confirmation is necessary before proceeding.
            If False, the behavior is as if the user replied "yes".
            Used in DISCARD_FULL mode before performing destructive
            actions and in domain runs in case any inconsistencies
            are found in the information stored in history. Default
            is True.
        domains : Sequence or None, optional
            Paths to domain folders in which bookkeeper is executed
            in addition to self.cwd. The same `mode` is used for all.
            If not given or None, domain subfolders are searched in
            self.cwd, based on the metadata stored in the last history
            folder in self.cwd. Default is None.

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
        BookkeeperUnexpectedError
            If, when running in ARCHIVE mode, no unlabeled version
            is found for any of the _ori-suffixed files. This should
            never happen. See also Issue #353.
        """
        try:
            mode = BookkeeperMode(mode)
        except ValueError as exc:
            raise ValueError(f'Unknown mode {mode}') from exc

        kwargs = {
            'mode': mode,
            'requires_user_confirmation': requires_user_confirmation,
            }
        if domains is None:
            try:
                domains, path_to_main = self._find_domains(**kwargs)
            except NotAnInteractiveShellError:
                LOGGER.info('')
                return BookkeeperExitCode.FAIL
            except (MetadataError, MainPathNotFoundError):
                LOGGER.error('Please proceed manually.')
                LOGGER.info('')
                return BookkeeperExitCode.FAIL
            if path_to_main:
                # This is a subdomain. Make sure we pull
                # information out of the main log file.
                main_root = RootExplorer(path_to_main, self)
                main_root.collect_info(silent=True)
                self._root = DomainRootExplorer(self.cwd, self, main_root)

        if domains:
            exit_code = self._run_in_root_and_subdomains(domains, **kwargs)
        else:
            exit_code, _ = self._run_one_domain(**kwargs)
        self._check_exit_code(exit_code, domains)
        return exit_code

    def update_from_cwd(self, silent=False):
        """Update timestamp, tensor number, log lines, etc. from root.

        This method is always called automatically before every run
        of the bookkeeper.

        Parameters
        ----------
        silent : bool, optional
            Whether logging messages should be silenced or not.
            The default is False. One basic log message is always
            emitted, irrespective of the value of `silent`, if
            this instance accesses the log file for the first time.

        Returns
        -------
        None.
        """
        self._make_history_and_prepare_logger()
        context = logging_silent if silent else nullcontext
        with context():
            # Collect information from the root directory:
            # log files, SUPP, OUT, workhistory, and Tensors, as well
            # as information from the history folder (especially the
            # highest run number currently stored for each tensor).
            self._root.collect_info()
            self._state_info['timestamp'] = (
                self._root.calc_timestamp  # If a log file was found
                or f'{MOVED_LABEL}{DateTimeFormat.FILE_SUFFIX.now()}'
                )

            # Infer which history subfolder we should be handling
            self.history.find_new_history_directory(self.tensor_number,
                                                    suffix=self.timestamp)

            # history.info handler (may create history.info file)
            self.history.prepare_info_file()
            self.history.info.read()

    @property
    def _has_archived_folder(self):
        """Return whether a history folder already exists for this calc run."""
        # NB: It's enough to check the log timestamp, as it should be
        # a unique identifier. It is incorrect to include the tensor
        # number, as users may delete Tensors (they're large), but that
        # would be misinterpreted as requiring archiving (see #476).
        # We don't need to explicitly check also a job number, as those
        # are always "fresh" since we FIRST create the contents of the
        # main folder, and only LATER move the workhistory ones (in
        # _archive_to_history_and_add_info_entry).
        archived_name = rf'.*_{self.timestamp}'
        return self.history.has_subfolder(archived_name)

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
                **self._root.infer_run_info(),
                )
        self.history.info.append_entry(new_info_entry, fix_time_format=True)

    def _archive_input_files_from_original_inputs_or_cwd(self):
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
        orig_inputs_dir = self.orig_inputs_dir
        for file in ALL_INPUT_FILES:
            original_file = orig_inputs_dir / file
            cwd_file = self.cwd / file
            copy_file, with_name = None, file
            if original_file.is_file() and cwd_file.is_file():
                # Copy original, but warn if cwd is newer
                copy_file = original_file
                try:
                    _check_newer(older=cwd_file, newer=original_file)
                except _FileNotOlderError:
                    LOGGER.warning(
                        f'File {file} from {orig_inputs_dir.name} was '
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
                    f'File {file} not found in {orig_inputs_dir.name}. '
                    'Using file from root directory instead and renaming to '
                    f'{with_name}.'
                    )
            if copy_file:
                history_folder.copy_file_or_directory(copy_file,
                                                      with_name=with_name)

    def _archive_log_files(self):
        """Copy all the log files in root to history."""
        for file in self._root.logs.files:
            self.history.new_folder.copy_file_or_directory(file)

    def _archive_out_and_supp(self):
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

    def _archive_to_history_and_add_info_entry(self):
        """Store all relevant files in root to history, add a new entry.

        This method does not check whether any archiving is indeed
        necessary. Make sure to check if self.archiving_required
        beforehand. After this method executes, the root folder
        still contains the output produced by viperleed.calc,
        except for those workhistory folders that were archived
        to history. The workhistory folder itself is cleaned-up
        of "previous" runs, and may not exist anymore if it was
        empty at the end of the archiving. A new entry is added
        to the history.info file.

        Returns
        -------
        None.
        """
        # Create the new 'primary' history directory...
        main_history_subfolder = self._make_and_copy_to_history()
        # ...move workhistory folders...
        tensor_nums = self._workhistory.move_current_and_cleanup(
            main_history_subfolder,
            )
        tensor_nums.add(self.tensor_number)
        # ...and add a history.info entry
        self._add_history_info_entry(tensor_nums)
        LOGGER.info('Done archiving the current directory '
                    f'to {self.history.path.name}.')

    def _check_exit_code(self, exit_code, domains):
        """Raise if exit_code signals bugs in bookkeeper."""
        if exit_code is BookkeeperExitCode.MISSING_UNLABELED_FILES:
            err_msg = (
                'Some input files that have been renamed to *_ori do not have '
                'a corresponding unlabeled version. This is most likely a bug '
                'in bookkeeper. Please report this to the ViPErLEED team by '
                'opening an issue under '
                'https://github.com/viperleed/viperleed/issues attaching the '
                'contents of the most recent history folder'
                )
            if domains:
                err_msg += ' of both the main as well as all domain subfolders'
            err_msg += '.'
            LOGGER.error(err_msg)
            raise BookkeeperUnexpectedError(err_msg)

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

    def _check_may_discard_full_domains(self, domains):
        """Ensure bookkeeper can run in DISCARD_FULL in all domains."""
        msg = ('Checking whether it is possible to '
               f'{BookkeeperMode.DISCARD_FULL.name} '
               'the most-recent results %s.')
        self.update_from_cwd(silent=True)  # Only first log message
        LOGGER.info(msg, 'in the root directory')
        self._check_may_discard_full()     # Raises if not possible

        # If we manage to reach here, we should check that
        # it is also possible to execute in all domains.
        for domain in domains:
            # Notice that here, differently from what we do in
            # _run_in_root_and_subdomains, we can use self._root
            # as we don't run anything in the root directory.
            domain_bookie = DomainBookkeeper(self._root, cwd=domain)
            with logging_silent():  # No logs from the domains
                domain_bookie.update_from_cwd()
                log.remove_bookkeeper_logfile(domain_bookie.history.path)
            LOGGER.info(msg, f'of domain {domain.name}')
            domain_bookie.check_may_discard_full()

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
        if in_archive_mode and self._has_archived_folder:
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
        # pull (in _archive_input_files_from_original_inputs_or_cwd):
        # always copy from original_inputs if possible, otherwise
        # copy with a _FROM_ROOT suffix.
        self._root.mark_edited_files()
        return True

    def _find_domains(self, mode, requires_user_confirmation):
        """Find registered domain subfolders of the current directory.

        Only the last history entry in the current directory (as well
        as the last one in each domain subfolder) is checked.

        Parameters
        ----------
        mode : BookkeeperMode
            The mode bookkeeper is running in. Only used for
            error reporting.
        requires_user_confirmation : bool
            Whether users will be asked to confirm in case any
            inconsistencies are found in the metadata stored
            in history.

        Returns
        -------
        domain_paths : tuple
            Absolute paths to the subfolders of the current
            directory that are registered as domains.
        main_root_path : str or None
            Absolute path to the main root directory of the
            multi-domain calculation. Non-None only if the
            current directory is a registered domain subfolder
            without a viperleed.calc log file.
        """
        self.update_from_cwd(silent=True)  # Only the first log message
        finder = DomainFinder(self, requires_user_confirmation)
        finder.collect_info()
        try:
            domains = finder.find_domains(mode)
        except MetadataError as exc:
            LOGGER.error(str(exc))
            raise
        if finder.is_subdomain and self._root.logs.most_recent is None:
            # We're running in a domain subfolder without a log file
            main_root, _ = finder.domain_info
        else:
            main_root = None
        return tuple(self.cwd/d for d in domains), main_root

    def _make_and_copy_to_history(self):
        """Create new history subfolder and copy all files there.

        Returns
        -------
        main_history_folder : HistoryFolder
            The new folder added to history.

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
                    'for storing results of the most recent viperleed.calc '
                    'execution.')
        self._archive_out_and_supp()
        self._archive_input_files_from_original_inputs_or_cwd()
        self._archive_log_files()

        # Now that all files are in place, add the metadata file
        meta = BookkeeperMetaFile(self.history.new_folder.path)
        meta.compute_hash()
        meta.write()

        # Finally, register the new folder in history
        return self.history.register_folder(self.history.new_folder.path)

    def _make_history_and_prepare_logger(self):
        """Make history folder and add handlers to the bookkeeper logger."""
        if self._state_info['logger_prepared']:
            return

        # Attach a stream handler to the logger if not already present
        log.ensure_has_stream_handler()

        try:  # Make the top-level 'history' folder if not there yet
            self.history.path.mkdir(exist_ok=True)
        except OSError:
            LOGGER.error(f'Error creating {self.history.path.name} folder.')
            raise

        # Attach file handler for history/bookkeeper.log
        log.add_bookkeeper_logfile(self.history.path)
        LOGGER.info(  # Log only once per instance
            '### Bookkeeper running at '
            f'{DateTimeFormat.LOG_CONTENTS.now()} ###'
            )
        self._state_info['logger_prepared'] = True

    def _print_discard_info(self):
        """Emit logging messages with files/folders/entry to be deleted."""
        LOGGER.info(f'About to delete files/folders from {self.cwd}.')
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
                    f'from {self.history.info.path.name}: '
                    f'{self.history.info.last_entry}'.rstrip())
        LOGGER.info('Calculation data will be lost irreversibly.')

    def _run_archive_mode(self):
        """Execute bookkeeper in ARCHIVE mode."""
        # Copy all the relevant files to history, add
        # a history.info entry, mark edited input files...
        did_archive = self._do_prerun_archiving_and_mark_edited()
        if not did_archive:
            self._root.complain_about_edited_files()
            return BookkeeperExitCode.NOTHING_TO_DO, None

        last_folder = self.history.last_folder
        # ...mark non-edited inputs as _ori, and prepare inputs for
        # the next execution of calc (from OUT or original_inputs)
        try:
            self._root.prepare_for_next_calc_run()
        except OSError:
            return BookkeeperExitCode.FAIL, last_folder

        # Finally, make sure we always have unlabeled files for all
        # _ori suffixed ones.
        misses_unlabeled = self._root.ensure_has_unlabled_inputs()
        exit_code = (BookkeeperExitCode.SUCCESS if not misses_unlabeled
                     else BookkeeperExitCode.MISSING_UNLABELED_FILES)
        return exit_code, last_folder

    def _run_clear_mode(self):
        """Execute bookkeeper in CLEAR mode."""
        did_archive = self._do_prerun_archiving_and_mark_edited()
        did_clear_root = self._root.clear_for_next_calc_run()
        did_anything = did_archive or did_clear_root
        if did_clear_root:
            LOGGER.info('Successfully prepared the current directory '
                        'for the next run of viperleed.calc.')
        exit_code = (BookkeeperExitCode.SUCCESS if did_anything
                     else BookkeeperExitCode.NOTHING_TO_DO)
        last_folder = self.history.last_folder if did_archive else None
        return exit_code, last_folder

    # TODO: how to handle the TENSOR_INDEX case? Currently we attempt
    # removal of the MAXIMUM tensor number, not the LAST one that ran.
    # The latter would be available from .history.last_folder.tensor_num
    # if we would sort the history folders not by name but by creation
    # datetime (https://stackoverflow.com/questions/168409).
    def _run_discard_full_mode(self):
        """Execute bookkeeper in DISCARD_FULL mode."""
        try:
            self._check_may_discard_full()
        except NoHistoryEntryError:
            return BookkeeperExitCode.NOTHING_TO_DO, None
        except (CantRemoveEntryError,
                FileNotFoundError,
                MetadataMismatchError):
            return BookkeeperExitCode.FAIL, None

        self._print_discard_info()
        try:
            user_confirmed = self._user_confirmed()
        except NotAnInteractiveShellError:
            return BookkeeperExitCode.FAIL, None

        if not user_confirmed:
            return BookkeeperExitCode.NOTHING_TO_DO, None

        # Delete the history folders, stuff in workhistory,
        # and output files/folders in the root directory
        try:
            history_deleted = self.history.discard_most_recent_run()
        except OSError:
            return BookkeeperExitCode.FAIL, None
        self._workhistory.discard_workhistory_root()
        self._root.revert_to_previous_calc_run()

        # Tensors and deltas, if created during the last run
        self._root.remove_tensors_and_deltas(history_deleted)

        # And the history entry from history.info
        self.history.info.remove_last_entry()
        LOGGER.info('Successfully deleted the results of '
                    'the last viperleed.calc execution.')
        return BookkeeperExitCode.SUCCESS, None

    def _run_discard_mode(self):
        """Execute bookkeeper in DISCARD mode."""
        did_archive = self._do_prerun_archiving_and_mark_edited()
        self._root.revert_to_previous_calc_run()
        if not did_archive:
            # If we had to archive the contents of the root directory,
            # we'd have already marked the entry as discarded and we
            # wouldn't have to bother. Here we do it explicitly.
            try:
                self.history.info.discard_last_entry()
            except (NoHistoryEntryError, CantDiscardEntryError) as exc:
                no_entry = isinstance(exc, NoHistoryEntryError)
                emit_log = LOGGER.warning if no_entry else LOGGER.error
                emit_log('Failed to mark as discarded the last '
                         f'entry in {HISTORY_INFO_NAME}: {exc}')
                exit_code = (BookkeeperExitCode.NOTHING_TO_DO if no_entry
                             else BookkeeperExitCode.FAIL)
                return exit_code, None
        LOGGER.info('Successfully marked as discarded the '
                    f'last entry in {HISTORY_INFO_NAME}.')
        last_folder = self.history.last_folder if did_archive else None
        return BookkeeperExitCode.SUCCESS, last_folder

    def _run_fix_mode(self):
        """Fix format inconsistencies found in history and history.info."""
        did_fix = self.history.fix()
        exit_code = (BookkeeperExitCode.SUCCESS if did_fix
                     else BookkeeperExitCode.NOTHING_TO_DO)
        return exit_code, None

    def _run_in_root_and_subdomains(self, domains, mode, **kwargs):
        """Execute bookkeeper in its cwd as well as all subdomains."""
        if mode is BookkeeperMode.DISCARD_FULL:
            try:
                self._check_may_discard_full_domains(domains)
            except (CantRemoveEntryError,
                    FileNotFoundError,
                    MetadataMismatchError,
                    NoHistoryEntryError):
                LOGGER.error('Cannot safely run bookkeeper in mode '
                             f'{mode.name} on all domains. Please '
                             'proceed manually.')
                LOGGER.info('')
                return BookkeeperExitCode.FAIL

        # Store a RootExplorer instance BEFORE running on the main
        # domain, as _run_one_domain will collect information from
        # there and later on create a new instance (in _clean_state).
        # However, we need up-to-date information for running in the
        # subdomains, BEFORE any deleting/archiving is done, as the
        # log file may be deleted too. Also, we can't use self._root,
        # as clear_for_next_calc_run does internally re-collect info.
        main_root = RootExplorer(self.cwd, self)
        main_root.collect_info(silent=True)

        # Store information in the root explorers about the presence
        # of domains. This has impact on which input files we may store
        # in history as well as which ones will be kept as non-suffixed
        # in the root directory. We should not rely on DomainFinder,
        # whose find_potential_domains is bugged. See also #344.
        self._root.has_domains = True
        main_root.has_domains = True
        main_exit_code, main_folder = self._run_one_domain(mode, **kwargs)
        exit_codes = [
            main_exit_code,
            *self._run_subdomains(domains,
                                  main_root,
                                  main_folder,
                                  mode,
                                  **kwargs),
            ]
        exit_code = BookkeeperExitCode.from_codes(exit_codes)
        if exit_code is BookkeeperExitCode.FAIL and domains:
            LOGGER.warning('Bookkeeper failed and may not have processed '
                           'some domain directories. Make sure to invoke '
                           f'\'bookkeeper {mode.long_flag}\' in all domain '
                           'subfolders or try again at this path.')
            LOGGER.info('')
        return exit_code

    def _run_one_domain(self, mode, requires_user_confirmation=True):
        """Execute Bookkeeper in `mode` for a single domain.

        Parameters
        ----------
        mode : BookkeeperMode
            Which bookkeeper mode to use.
        requires_user_confirmation : bool, optional
            Whether user confirmation is necessary before proceeding
            with destructive actions. Only used in DISCARD_FULL mode.
            Default is True.

        Returns
        -------
        exit_code : BookkeeperExitCode
            The result of executing Bookkeeper in `mode`.
        last_folder : HistoryFolder or None
            The last folder that was archived, if any. None otherwise.

        Raises
        ------
        NotImplementedError
            If `mode` is a valid bookkeeper mode, but
            there is no method named _run_<mode>_mode.
        OSError
            If creation of the history folder or any of the subfolders
            where results are to be stored fails.
        """
        try:
            runner = getattr(self, f'_run_{mode.name.lower()}_mode')
        except AttributeError as exc:
            raise NotImplementedError from exc

        self._mode = mode
        self._requires_user_confirmation = requires_user_confirmation
        # Do not bother logging messages when the user asked to
        # --fix: some warnings will be fixed anyway. Others will
        # re-appear at the next run.
        self.update_from_cwd(silent=mode is BookkeeperMode.FIX)
        LOGGER.info(f'Running bookkeeper in {mode.name} mode in {self.cwd}.')
        self._warn_about_old_calc()
        try:
            exit_code, last_folder = runner()
        finally:
            # Clean up all the state attributes so we
            # don't risk giving the wrong information.
            self._clean_state()
        found_nothing = (
           exit_code is BookkeeperExitCode.NOTHING_TO_DO
           and mode not in (BookkeeperMode.ARCHIVE,)  # logs already
           )
        if found_nothing:
            LOGGER.info('Found nothing to do. Exiting...')
        LOGGER.info('')
        return exit_code, last_folder

    def _run_subdomains(self, domains, main_root, main_folder, mode, **kwargs):
        """Execute bookkeeper in `domains`.

        Parameters
        ----------
        domains : Sequence of Path
            Absolute path to each of the subdomain folders to be
            processed.
        main_root : RootExplorer
            Handler of the **main** root folder of the multi-domain
            calculation (i.e., where calc was originally invoked).
        main_folder : HistoryFolder or None
            The subfolder of `main_root`.path/'history' where the
            results of the main calc execution were stored (i.e.,
            not one of those coming from workhistory). Used for
            storing metadata information about subdomains.
        mode : BookkeeperMode
            The mode in which to run in all `domains`.
        **kwargs : object, optional
            Other keyword arguments, passed unaltered to each of
            the .run calls that execute in `domains`.

        Yields
        ------
        domain_exit_code : BookkeeperExitCode
            One exit code for each of the `domains` in which
            bookkeeper was executed in `mode`.

        Notes
        -----
        `main_folder` is only updated with domain information after all
            the `domains` have been processed. This means that it is
            necessary to **exhaust this generator function** before
            `main_folder` contains domain-related metadata information.
            Conversely, each of the `domains` subfolders contain
            domain markings as soon as they have been processed.
        """
        assert domains
        domain_rel_paths = []
        for path in domains:
            try:
                path = path.relative_to(self.cwd).as_posix()
            except ValueError:
                pass
            domain_rel_paths.append(path)
        LOGGER.info('Running bookkeeper in domain folders %s:',
                    harvard_commas(*domain_rel_paths))
        domain_folders = []    # Archived HistoryFolder for each domain
        for path in domains:
            domain_bookie = DomainBookkeeper(main_root, cwd=path)
            dom_exit, folder = domain_bookie.run_in_subdomain(main_folder,
                                                              mode,
                                                              **kwargs)
            domain_folders.append(folder)
            yield dom_exit
        if main_folder:
            main_folder.mark_as_domains_main(domain_rel_paths, domain_folders)
            main_folder.metadata.write()
        LOGGER.info('Done processing domain folders %s.',
                    harvard_commas(*domain_rel_paths))
        LOGGER.info('')

    def _user_confirmed(self):
        """Return whether the user wants to proceed with discarding."""
        if not self._requires_user_confirmation:
            return True
        return ask_user_confirmation(self._mode)

    def _warn_about_old_calc(self):
        """Emit warnings when this tree was created by an early calc."""
        if self._mode is BookkeeperMode.FIX:
            # FIX already emits log messages if it finds funky stuff.
            return
        version = self._root.logs.version
        # See if the last history folder was created with an old calc
        # version. All modes (except --fix) may be affected, as the
        # handling of files has changed since v0.13.0.
        if not version and self.history.last_folder:
            version = self.history.last_folder.logs.version
        if version and version < _MIN_CALC_WARN:
            LOGGER.warning(_WARN_OLD_CALC, version)


class DomainBookkeeper(Bookkeeper):
    """A Bookkeeper for handling a domain subfolder.

    The primary reason for having an explicit class for this is that
    domain subfolders do not contain a log file. Hence information
    derived from the log file should be taken from the root folder
    of the main calculation.
    """

    def __init__(self, main_root, cwd=None):
        """Initialize instance.

        Parameters
        ----------
        main_root : RootExplorer
            Information from the root folder of the main calculation
            directory. It must already be up to date by calling
            main_root.collect_info before any editing occurs in this
            root directory. This instance will never collect_info
            again on `main_root`.
        cwd : Pathlike, optional
            Path to the domain subfolder to be handled. If not given
            or None, take the current directory. Default is None.

        Returns
        -------
        None.
        """
        super().__init__(cwd)
        # IMPORTANT for developers: NEVER call collect_info on this
        # _main_root RootExplorer in this class, as the root folder
        # may have been modified by the main bookkeeper by that time.
        # Similarly, DO NOT re-create this attribute from its .path.
        # This means that properties such as _main_root.needs_archiving
        # should NEVER be considered up to date.
        self._main_root = main_root
        self._root = DomainRootExplorer(self._root.path, self, main_root)

    def check_may_discard_full(self):
        """Log and raise unless DISCARD_FULL is possible."""
        # This is just a wrapper around the private method in order
        # to make it public for DomainBookkeeper (as opposed to
        # Bookkeeper).
        self._check_may_discard_full()

    def run_in_subdomain(self, main_folder, *args, **kwargs):
        """Execute Bookkeeper in this subdomain.

        Parameters
        ----------
        main_folder : HistoryFolder
            The main folder created in the history of the root
            of the multi-domain calculation. Used for storing
            metadata information in history folders that may
            be created in this subdomain.
        *args : object
            Positional arguments, passed on unaltered to the
            "run" call for this subdomain.
        **kwargs : object, optional
            Optional keyword arguments, passed on unaltered
            to the "run" call for this subdomain.

        Returns
        -------
        exit_code : BookkeeperExitCode
            The exit code resulting from running bookkeeper
            in this subdomain.
        archived_folder : HistoryFolder or None
            The "main" folder that was added to the history of
            this subdomain as a result of running bookkeeper.

        Raises
        ------
        NotImplementedError
            If `mode` is a valid bookkeeper mode, but
            there is no method named _run_<mode>_mode.
        OSError
            If creation of the history folder or any of the subfolders
            where results are to be stored fails.
        """
        try:
            exit_code, archived_folder = self._run_one_domain(*args, **kwargs)
        finally:
            log.remove_bookkeeper_logfile(self.history.path)
        if archived_folder:
            archived_folder.mark_as_domain(self._main_root.path, main_folder)
            archived_folder.metadata.write()
        return exit_code, archived_folder

    # The next method needs to be overridden because the signature
    # of __init__ is changed compared to the one of Bookkeeper.
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
        self.__init__(self._main_root, cwd=self.cwd)
        self._state_info['logger_prepared'] = logger_prepared


def _check_newer(older, newer):
    """Raise if file `older` is newer than file `newer`."""
    newer_timestamp = newer.stat().st_mtime
    older_timestamp = older.stat().st_mtime
    if newer_timestamp >= older_timestamp:
        return
    if file_contents_identical(older, newer):
        return
    raise _FileNotOlderError
