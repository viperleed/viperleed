"""Module root_explorer of viperleed.calc.bookkeeper.

Defines objects useful for identifying, reading, and modifying files
and folders present in the directory in which Bookkeeper runs.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2024-10-14'
__license__ = 'GPLv3+'

from collections import defaultdict
from contextlib import nullcontext
from datetime import datetime
from functools import partial
import logging
from operator import attrgetter
from pathlib import Path
import shutil

from viperleed.calc.bookkeeper.constants import EDITED_SUFFIX
from viperleed.calc.bookkeeper.constants import ORI_SUFFIX
from viperleed.calc.bookkeeper.constants import STATE_FILES
from viperleed.calc.bookkeeper.domain_finder import DomainFinder
from viperleed.calc.bookkeeper.errors import FileOperationFailedError
from viperleed.calc.bookkeeper.history.explorer import HistoryExplorer
from viperleed.calc.bookkeeper.history.workhistory import WorkhistoryHandler
from viperleed.calc.bookkeeper.log import LogFiles
from viperleed.calc.bookkeeper.utils import discard_files
from viperleed.calc.bookkeeper.utils import file_contents_identical
from viperleed.calc.bookkeeper.utils import make_property
from viperleed.calc.bookkeeper.utils import needs_update_for_attr
from viperleed.calc.constants import DEFAULT_DELTAS
from viperleed.calc.constants import DEFAULT_OUT
from viperleed.calc.constants import DEFAULT_SUPP
from viperleed.calc.constants import DEFAULT_TENSORS
from viperleed.calc.constants import ORIGINAL_INPUTS_DIR_NAME
from viperleed.calc.constants import SKIP_IN_DOMAIN_MAIN
from viperleed.calc.lib.leedbase import getMaxTensorIndex
from viperleed.calc.lib.log_utils import logging_silent
from viperleed.calc.lib.string_utils import harvard_commas
from viperleed.calc.lib.time_utils import DateTimeFormat


LOGGER = logging.getLogger(__name__)


# TODO: read_and_clear_notes_file  -- support multiple files
class RootExplorer:
    """A class for handling files and folders in the bookkeeper's CWD."""

    _needs_collect = partial(needs_update_for_attr, updater='collect_info')

    def __init__(self, path, bookkeeper):
        """Initialize this explorer at `path` for a `bookkeeper`."""
        self._path = Path(path).resolve()
        self.has_domains = False       # Is this the root of a DOMAINS?
        self._logs = None              # LogFiles, set in collect_info
        self._files_to_archive = None  # See _collect_files_to_archive
        self.tensors = TensorAndDeltaInfo(self.path)
        self.history = HistoryExplorer(self.path)
        self.workhistory = WorkhistoryHandler(root=self,
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
        bool
            Whether any file operation took place.
        """
        return self._clean_up_root(deal_with_ori_files=self._remove_ori_files)

    def collect_info(self, silent=False):
        """Collect information from the root directory."""
        log_context = logging_silent if silent else nullcontext
        with log_context():
            self._logs = LogFiles(self.path)
            self._logs.collect()
            self.tensors.collect()
            self._collect_files_to_archive()
            self.history.collect_subfolders()
            # TODO: consider whether we can partly undo the changes
            # that fixed #344 when we fix #215. Then we would likely
            # be reading the PARAMETERS file and know whether we're
            # running in a DOMAINS tree without the need of having
            # Bookkeeper let us know by setting .has_domains.
            # self._find_potential_domain_subfolders()

    def complain_about_edited_files(self):
        """Log warnings if any _edited file is found in root."""
        edited_files = tuple(f'{self._relative_path(f)}'
                             for f in self.path.glob(f'*{EDITED_SUFFIX}*'))
        if not edited_files:
            return
        edited = ', '.join(edited_files)
        inputs = ', '.join(f.split(EDITED_SUFFIX)[0] for f in edited_files)
        LOGGER.warning(f'Found user-edited files {edited}. Make sure to port '
                       'any desired changes to the corresponding input files '
                       f'(i.e., {inputs}) or to delete the *{EDITED_SUFFIX} '
                       'files.')

    def ensure_has_unlabled_inputs(self):
        """Make sure there is always unlabeled versions of any _ori file.

        The unlabeled versions are created as copies of the _ori-suffixed
        ones (unless there is also an _edited version).

        Returns
        -------
        any_missing : bool
            Whether any unlabeled file was added.
        """
        any_missing = False
        for ori_file in self.path.glob(f'*{ORI_SUFFIX}'):
            unlabeled = ori_file.name.replace(ORI_SUFFIX, '')
            unlabeled_file = ori_file.with_name(unlabeled)
            edited_file = ori_file.with_name(f'{unlabeled}{EDITED_SUFFIX}')
            if edited_file.is_file() or unlabeled_file.is_file():
                continue
            any_missing = True
            try:
                shutil.copy2(ori_file, unlabeled_file)
            except OSError as exc:
                LOGGER.error(f'Could not create file '
                             f'{self._relative_path(unlabeled_file)} '
                             f'from its *{ORI_SUFFIX} version.')
        return any_missing

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
        ori_files = {self.path/file: self.path/f'{file}{ORI_SUFFIX}'
                     for file in STATE_FILES}
        ori_files = {file: ori_file for file, ori_file in ori_files.items()
                     if ori_file.exists()}
        return tuple(ori_files.items())

    def list_paths_to_discard(self):
        """Return a tuple of files and folders that will be discarded."""
        to_discard = (
            self.path / DEFAULT_OUT,
            self.path / DEFAULT_SUPP,
            *self.logs.files,
            *self.tensors.list_paths_to_discard(self.history),
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
        for filename in STATE_FILES:
            cwd_file = self.path / filename
            try:
                mod_time = datetime.fromtimestamp(cwd_file.stat().st_mtime)
            except FileNotFoundError:
                continue
            if mod_time < calc_started:
                continue
            # The input file was modified. See whether its contents
            # differ from the one we have copied in original_inputs.
            # If they differ (or there is no file in original_inputs),
            # label this file as having been user-edited.
            if file_contents_identical(cwd_file, ori_path / filename):
                continue
            edited = filename + EDITED_SUFFIX
            LOGGER.warning(f'File {self._relative_path(cwd_file)} was modified'
                           f' after calc started. Renaming it to {edited}.')
            try:
                cwd_file.replace(self.path / edited)
            except OSError:
                LOGGER.error(
                    f'Failed to rename {self._relative_path(cwd_file)} to '
                    f'{edited}. It may be overwritten by the corresponding '
                    f'file in {DEFAULT_OUT}.'
                    )

    def prepare_for_next_calc_run(self):
        """Rename inputs to _ori, pull new ones from OUT or original_inputs."""
        errors = []
        try:
            self._mark_state_files_as_ori()
        except FileOperationFailedError as exc:  # Re-raise it later
            errors.append(exc)

        # At this point, all STATE_FILES should have been renamed
        # to either _ori (above) or _edited (via an external call
        # to mark_edited_files).
        try:  # Pull in the inputs for the next calc run.
            self._copy_state_files_from_out_or_original_inputs()
        except FileOperationFailedError as exc:
            errors.append(exc)

        self.complain_about_edited_files()
        if errors:
            raise OSError('\n'.join(str(e) for e in errors))

    def read_and_clear_notes_file(self):
        """Return notes read from file. Clear the file contents."""
        notes_path = next(self.path.glob('notes*'), None)
        if notes_path is None:
            return ''
        try:
            notes = notes_path.read_text(encoding='utf-8')
        except OSError:
            LOGGER.error(
                f'Failed to read {self._relative_path(notes_path)} file.',
                exc_info=True,
                )
            return ''
        try:
            notes_path.write_text('', encoding='utf-8')
        except OSError:
            LOGGER.error(
                f'Failed to clear file {self._relative_path(notes_path)} '
                'after reading.',
                exc_info=True,
                )
        return notes

    def remove_tensors_and_deltas(self, history_folders):
        """Delete Tensors and Deltas archives related to `history_folders`."""
        self.tensors.discard(self.history, history_folders)

    def revert_to_previous_calc_run(self):
        """Remove files generated by the previous calc run in self.path.

        Notice that files in the workhistory directory, as
        well as the workhistory directory are not touched.

        Returns
        -------
        bool
            Whether any file operation took place.

        Raises
        ------
        OSError
            If replacing input files with their '_ori'-suffixed
            version fails.
        """
        return self._clean_up_root(
            deal_with_ori_files=self._replace_state_files_from_ori
            )

    def _clean_up_root(self, deal_with_ori_files):
        """Clean the root directory from calc results.

        Notice that workhistory files and the workhistory
        directory itself are left untouched.

        Parameters
        ----------
        deal_with_ori_files : callable
            A callable that handles _ori-suffixed files.

        Returns
        -------
        bool
            Whether any file operation took place.
        """
        deleted_logs = self.logs.discard()
        deleted_out_supp = self._remove_out_and_supp()
        dealt_with_ori = deal_with_ori_files()
        self.complain_about_edited_files()
        # Let's not complain again about funny stuff. We have already
        # done so when collect_info was called the first time. Notice
        # that this must have happened, as this method relies on the
        # logs attribute to be up to date.
        self.collect_info(silent=True)
        return any((deleted_logs, deleted_out_supp, dealt_with_ori))

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

    def _copy_state_files_from(self, source, *name_fmts, only_files=None):
        """Copy input files from `source` to the root directory.

        Parameters
        ----------
        source : Path
            The path where state files should be copied from.
            It must be a subfolder of self.path.
        *name_fmts : str
            Strings that will be .format()ted with the name of
            each state file to pick which file should be copied
            from `source`. The order of `name_fmts` corresponds
            to the priority in which files are searched. If no
            format is given, files are assumed to be named
            identically to their destination name. `source`
            is glob()-bed with each of the `name_fmts`.
        only_files : iterable or None
            Only consider these file names when copying. If
            None or not given, all the "state" input files are
            considered. The default is None.

        Raises
        ------
        FileOperationFailedError
            If any copying fails.
        """
        failed = {}
        name_fmts = name_fmts or ('{}',)
        if only_files is None:
            only_files = STATE_FILES
        dont_complain = set(SKIP_IN_DOMAIN_MAIN if self.has_domains else ())
        for file in only_files:
            cwd_file = self.path / file
            patterns = [fmt.format(file) for fmt in name_fmts]
            new_inputs = (f for p in patterns for f in source.glob(p))
            try:
                new_input = next(f for f in new_inputs if f.is_file())
            except StopIteration:
                if file in dont_complain:
                    continue
                failed[file] = [source/p for p in patterns]
                continue
            try:
                shutil.copy2(new_input, cwd_file)
            except OSError as exc:
                failed[file] = exc
        if failed:
            raise FileOperationFailedError(failed)

    def _copy_state_files_from_out_or_original_inputs(self):
        """Copy to root OUT or original_inputs files as new calc inputs."""
        # Try OUT first
        failed_out = {}
        try:
            # Prefer those without an _OUT suffix, but fall back onto
            # _OUT-suffixed ones to ensure backward compatibility. NB:
            # The old-style _OUT-suffixed files used to also have a
            # timestamp following '_OUT_'. Be explicit about the
            # presence of a digit not fetch POSCAR_OUT_mincell.
            self._copy_state_files_from(self.path/DEFAULT_OUT,
                                        '{}',
                                        '{}_OUT_[0-9]*')
        except FileOperationFailedError as exc:
            failed_out.update(exc.failures)

        # Then original_inputs, for those that were not found in OUT
        failed_ori = {}
        try:
            self._copy_state_files_from(self.orig_inputs_dir,
                                        only_files=failed_out)
        except FileOperationFailedError as exc:
            failed_ori.update(exc.failures)

        # Consider failed those that failed in both attempts.
        # Report exceptions with higher priority than files we tried.
        failed = {}
        for file in failed_ori:  # pylint: disable=C0206
            reasons = failed_out[file], failed_ori[file]
            try:
                failed[file] = next(r for r in reasons
                                    if isinstance(r, Exception))
            except StopIteration:  # Both are lists of files
                failed[file] = sum(reasons, [])
        if failed:
            self._log_failures_when_copying_input_files(failed)
            raise FileOperationFailedError(failed)

    def _find_potential_domain_subfolders(self):
        """Find all subfolders of self.path that look like subdomains."""
        finder = DomainFinder(self.workhistory.bookkeeper)
        self.has_domains = any(finder.find_potential_domains())

    def _log_failures_when_copying_input_files(self, failed):
        """Emit logging errors about failures to pull input files.

        Parameters
        ----------
        failed : dict
            Keys are names of input files that we failed
            to copy to the root folder. Values are either
            exceptions that occurred while copying, or list
            of paths of the files that we attempted to copy.

        Returns
        -------
        None.
        """
        LOGGER.error('Failed to collect input files for the next calc run:')
        for file, info in failed.items():
            if isinstance(info, Exception):
                reason = f'raised {type(info).__name__} - {info}'
            else:  # List of paths to files that were tried
                files = harvard_commas(
                    *(f'{self._relative_path(f)}' for f in info),
                    sep='or',
                    )
                reason = f'No {files} found'
            LOGGER.error(f'    {file}: {reason}')

    def _mark_state_files_as_ori(self):
        """Suffix input files in root as _ori. Raise on failure."""
        failed = {}
        for file in STATE_FILES:
            file_ori = f'{file}{ORI_SUFFIX}'
            cwd_file = self.path / file
            try:
                cwd_file.replace(self.path / file_ori)
            except FileNotFoundError:  # May have been _edited
                continue
            except OSError as exc:
                LOGGER.error(
                    f'Failed to rename {self._relative_path(cwd_file)} '
                    f'to {file_ori}.'
                    )
                failed[file] = exc
        if failed:
            raise FileOperationFailedError(failed)

    def _relative_path(self, path):
        """Return a version of `path` relative to this root folder."""
        return path.relative_to(self.path).as_posix()

    def _remove_ori_files(self):
        """Delete '_ori'-suffixed files from root."""
        ori_files = (self.path / f'{file}{ORI_SUFFIX}' for file in STATE_FILES)
        return discard_files(*ori_files)

    def _remove_out_and_supp(self):
        """Delete the SUPP and OUT directories from root."""
        return discard_files(self.path / DEFAULT_OUT, self.path / DEFAULT_SUPP)

    def _replace_state_files_from_ori(self):
        """Replace input files with their '_ori'-suffixed version."""
        to_replace = self.list_files_to_replace()
        for state_file, ori_file in to_replace:
            try:
                ori_file.replace(state_file)
            except OSError:
                LOGGER.error(
                    f'Failed to rename {self._relative_path(ori_file)} '
                    f'to {state_file.name}.'
                    )
                raise
        return any(to_replace)


class DomainRootExplorer(RootExplorer):
    """A RootExplorer for a domain subfolder."""

    def __init__(self, path, bookkeeper, main):
        """Initialize instance.

        Parameters
        ----------
        path : Pathlike
            Path to the domain subfolder handled.
        bookkeeper : Bookkeeper
            A Bookkeeper instance. May be the one running in the
            main folder, or the one running in this domain subfolder.
        main : RootExplorer
            The handler of the root folder of the "main" calculation,
            i.e., the one in which calc was started during a DOMAINS
            calculation.

        Returns
        -------
        None.
        """
        super().__init__(path, bookkeeper)
        self._main = main

    @property
    def calc_timestamp(self):
        """Return the timestamp of the main calc log file."""
        # There usually is no calc log file in the domain subfolders,
        # unless calc was explicitly invoked in this subfolder.
        return super().calc_timestamp or self._main.calc_timestamp

    def infer_run_info(self):
        """Return a dictionary of information read from the newest calc log."""
        return super().infer_run_info() or self._main.infer_run_info()

    def _relative_path(self, path):
        """Return a relative version of `path`."""
        try:
            return path.relative_to(self._main.path).as_posix()
        except ValueError:
            pass
        return super()._relative_path(path)


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

    def discard(self, history, history_folders):
        """Delete Tensor and Delta zip files related to `history_folders`.

        Parameters
        ----------
        history : HistoryExplorer
            The handler of the history folder. Tensor/Delta files are
            removed only if they have: (i) at most one history run left
            for existing `history_folders`, or (ii) no run at all for
            non-existing (i.e., typically already deleted) ones.
        history_folders : Sequence of HistoryFolder
            The history folders whose Tensor/Delta files may be
            removed. Only the files with Tensor indices included
            in `history_folders` are considered.

        Returns
        -------
        None.
        """
        to_discard = self._list_tensors_deltas_to_discard(history,
                                                          history_folders)
        discard_files(*to_discard)

    def list_paths_to_discard(self, history):
        """Return a tuple of paths to Tensor/Delta files to be discarded."""
        # The tensors that may be removed are those of
        # the last run archived in the history folder
        last_run = history.last_folder_and_siblings
        return self._list_tensors_deltas_to_discard(history, last_run)

    def _list_tensors_deltas_to_discard(self, history, history_folders):
        """Return paths to Tensor/Delta files to be discarded.

        Parameters
        ----------
        history : HistoryExplorer
            The handler of the history folder.
        history_folders : Sequence of HistoryFolder
            The history folders whose Tensor/Delta files may be
            removed. Only the files with Tensor indices included
            in `history_folders` are considered.

        Returns
        -------
        paths_to_discard : tuple
            Absolute path to Tensor/Delta archives to be deleted.
        """
        tensors_to_folders = defaultdict(list)
        for folder in history_folders:
            tensors_to_folders[folder.tensor_num].append(folder)

        # However, we should only discard the "new" ones, i.e., those
        # created during the last calc execution. These will have only
        # one run left (if the history folders have not been deleted
        # yet) or zero (if they have been deleted already).
        max_run = history.max_run_per_tensor
        tensor_nums = [
            t
            for t, folders in tensors_to_folders.items()
            if max_run[t] <= (1 if any(f.exists for f in folders) else 0)
            ]

        tensor_fmt = f'{DEFAULT_TENSORS}/{DEFAULT_TENSORS}_{{ind:03d}}.zip'
        delta_fmt = f'{DEFAULT_DELTAS}/{DEFAULT_DELTAS}_{{ind:03d}}.zip'
        to_discard = [tensor_fmt.format(ind=t) for t in tensor_nums]
        to_discard.extend(delta_fmt.format(ind=t) for t in tensor_nums)

        # Make them into paths, keeps those that exist, sort by name
        to_discard_paths = (self._root / f for f in to_discard)
        to_discard_paths = (f for f in to_discard_paths if f.is_file())
        return tuple(sorted(to_discard_paths, key=attrgetter('name')))
