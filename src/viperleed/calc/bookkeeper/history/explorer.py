"""Module explorer of viperleed.calc.bookkeeper.history.

Defines functionality useful for navigating the history folder.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2024-10-14'
__license__ = 'GPLv3+'

from collections import Counter
from collections import defaultdict
import logging
from operator import attrgetter
from pathlib import Path
import re
import shutil

from viperleed.calc.bookkeeper.history.folder import FolderFixAction
from viperleed.calc.bookkeeper.history.folder import HistoryFolder
from viperleed.calc.bookkeeper.history.folder import IncompleteHistoryFolder
from viperleed.calc.bookkeeper.history.info import HistoryInfoFile
from viperleed.calc.bookkeeper.utils import make_property
from viperleed.calc.bookkeeper.utils import needs_update_for_attr
from viperleed.calc.constants import DEFAULT_HISTORY
from viperleed.calc.sections.cleanup import MOVED_LABEL


LOGGER = logging.getLogger(__name__)


class HistoryExplorer:
    """A file-system explorer for the 'history' folder."""

    _updater = 'collect_subfolders'

    def __init__(self, root):
        """Initialize a history explorer at a `root` path."""
        self._path = root / DEFAULT_HISTORY
        # While the history.info file is technically in the root
        # directory and not in the history directory, it makes
        # sense to keep it in the same object that handles the
        # history directory since they are so much interrelated.
        self._info = None
        self._subfolders = []  # Always sorted by name
        self._maps = {
            'main_hash_to_folders': None,  # {parent_hash: {folders}}
            'hash_to_parent': None,   # {folder.hash_: parent_folder}
            'jobs_for_tensor': None,
            }
        self._new_calc_run_folder = None    # IncompleteHistoryFolder

    path = make_property('_path')        # The 'history' folder
    info = make_property('_info',        # The history.info file
                         needs_update=True,
                         updater='prepare_info_file')
    root = make_property('path.parent')  # The one containing 'history'
    new_folder = make_property('_new_calc_run_folder',
                               needs_update=True,
                               updater='find_new_history_directory')

    @property
    @needs_update_for_attr('_maps[hash_to_parent]', updater=_updater)
    def last_folder(self):
        """Return the HistoryFolder of the most-recent "main" calc run."""
        try:
            last_appended = self._subfolders[-1]
        except IndexError:
            return None
        # pylint: disable-next=unsubscriptable-object  # Decorated
        return self._maps['hash_to_parent'][last_appended.hash_]

    @property
    def last_folder_and_siblings(self):
        """Return all the HistoryFolder(s) created during the last run."""
        last_folder = self.last_folder
        if not last_folder or not last_folder.exists:
            return tuple()
        # pylint: disable-next=unsubscriptable-object     # It's a dict
        return tuple(self._maps['main_hash_to_folders'][last_folder.hash_])

    @property
    @needs_update_for_attr('_maps[jobs_for_tensor]', updater=_updater)
    def max_run_per_tensor(self):
        """Return a {tensor_num: max_job_num} dictionary."""
        map_ = self._maps['jobs_for_tensor']
        max_jobs = {t_num: max(job_nums)
                    for t_num, job_nums in map_.items()
                    if job_nums}
        return defaultdict(int, max_jobs)

    def check_last_folder_consistent(self):
        """Raise if the last folder in history has some inconsistency.

        Raises
        ------
        FileNotFoundError
            If there is no last folder.
        MetadataMismatchError
            If the hash computed from the contents of the last folder
            and the one stored as metadata do not match.
        CantRemoveEntryError
            If some inconsistency is found between the last folder and
            the last entry in the history.info file.
        """
        last_folder = self.last_folder
        if not last_folder or not last_folder.exists:
            raise FileNotFoundError
        last_folder.check_metadata()
        last_folder.check_consistent_with_entry(self.info.last_entry)

    def collect_subfolders(self):
        """Navigate the main history folder and pull subfolder information."""
        self._subfolders.clear()
        self._maps['jobs_for_tensor'] = defaultdict(set)
        self._maps['main_hash_to_folders'] = defaultdict(list)
        self._maps['hash_to_parent'] = {}
        if not self.path.is_dir():
            LOGGER.warning(f'No {DEFAULT_HISTORY!r} folder in {self.root}.')
            return
        # pathlib.Path.iterdir has no guarantees on the sorting order
        # with which the directories are returned. We use alphabetical
        # sorting by name. This way, the most recent folders are read
        # last.
        contents = sorted(self.path.iterdir(), key=attrgetter('name'))
        for directory in contents:
            try:
                self._append_existing_folder(directory, insert_sorted=False)
            except ValueError:  # Not a directory or not a valid name
                pass
        self._update_maps()

    def discard_most_recent_run(self):
        """Delete all subfolders created during the last calc run."""
        subfolders = self.list_paths_to_discard()
        for folder_path in subfolders:
            try:
                shutil.rmtree(folder_path)
            except OSError:
                folder_name = folder_path.relative_to(self.root)
                LOGGER.error(f'Error: Failed to delete {folder_name}.')
                raise

    def find_new_history_directory(self, tensor_number, suffix):
        """Store information to create a history folder for a new calc run."""
        dir_name = self._find_name_for_new_history_subfolder(tensor_number,
                                                             suffix)
        new_folder_path = self.path / dir_name
        self._new_calc_run_folder = IncompleteHistoryFolder(new_folder_path)

    def fix(self):
        """Fix folders in history and the history.info file.

        Returns
        -------
        bool
            Whether any fixing was actually performed.
        """
        fixed_info = self._fix_info_file()
        fixed_folders = self._fix_subfolders()
        fixed_anything = fixed_info or fixed_folders
        if fixed_anything:
            LOGGER.info(f'Successfully fixed {self.info.path.name} '
                        f'file and {self.path.name} folder.')
        return fixed_anything

    def has_subfolder(self, name_rgx):
        """Return whether a a subfolder matching `name_rgx` is present."""
        return any(re.fullmatch(name_rgx, f.name) for f in self._subfolders)

    def list_paths_to_discard(self):
        """Return a tuple of paths to folders that will be discarded."""
        paths = (f.path for f in self.last_folder_and_siblings)
        return tuple(sorted(paths, key=attrgetter('name')))

    def prepare_info_file(self):
        """Prepare a history.info file in the root directory."""
        self._info = HistoryInfoFile(self.root, create_new=True)

    @needs_update_for_attr('_maps[main_hash_to_folders]', updater=_updater)
    def register_folder(self, path_to_folder):
        """Register a new folder as part of the history tree."""
        appended = self._append_existing_folder(path_to_folder,
                                                insert_sorted=True)
        self._update_maps()
        return appended

    @needs_update_for_attr('_maps[hash_to_parent]', updater=_updater)
    def subfolder_from_hash(self, hash_):
        """Return the subfolder with a given `hash_`.

        Parameters
        ----------
        hash_ : str
            Hash value of the folder to be returned.

        Returns
        -------
        main_folder : HistoryFolder or None
            The "main" folder added by bookkeeper in history
            together with the one whose hash equals `hash_`.
            None if no subfolder with such hash exists.
        """
        try:
            # pylint: disable-next=unsubscriptable-object   # Inference
            return self._maps['hash_to_parent'][hash_]
        except (KeyError, TypeError):
            return None

    def _append_existing_folder(self, path_to_folder, insert_sorted=True):
        """Register a new folder, without updating the parent mapping.

        Parameters
        ----------
        path_to_folder : Path
            Path to the folder to be added.
        insert_sorted : bool, optional
            Whether `path_to_folder` should be added in a sorted
            manner, based on its .name. Passing False speeds up
            insertion. IMPORTANT NOTE FOR DEVELOPERS: it is critical
            not to pass False here unless you are sure that you call
            _append_existing_folder already in a sorted manner.
            Otherwise, the .last_folder property will NOT be the
            most-recent "main" calc run. Default is True.

        Returns
        -------
        appended_folder : HistoryFolder
            The history folder that was added from `path_to_folder`.

        Raises
        ------
        ValueError
            If `path_to_folder` is not a subfolder of history.
        ValueError
            If `path_to_folder` is not an existing directory.
        """
        if path_to_folder.parent != self.path:
            raise ValueError(f'Not a subfolder of {self.path}.')
        folder = HistoryFolder(path_to_folder)
        parent_hash = folder.parent or folder.hash_
        self._subfolders.append(folder)
        if insert_sorted:
            self._subfolders.sort(key=attrgetter('path.name'))

        # pylint: disable=unsubscriptable-object  # Called in decorated
        self._maps['main_hash_to_folders'][parent_hash].append(folder)
        self._maps['jobs_for_tensor'][folder.tensor_num].add(folder.job_num)
        return folder

    def _backup_info_file(self):
        """Create a numbered duplicate of history.info as a backup."""
        info = self.info.path
        if not info.is_file():
            return None
        i = 1
        while True:
            backup_info = Path(f'{info}.bak{i}')
            if not backup_info.is_file():
                info.replace(backup_info)
                return backup_info
            i += 1

    @needs_update_for_attr('_maps[jobs_for_tensor]', updater=_updater)
    def _find_name_for_new_history_subfolder(self, tensor_number, suffix):
        """Return the potential name of a new history subfolder.

        This is the name of .new_folder, i.e., the one that will
        collect the final output and SUPP files. Folders coming from
        the workhistory folder have their own naming.

        Returns
        -------
        str
        """
        jobs = self._maps['jobs_for_tensor']
        # pylint: disable-next=unsubscriptable-object   # Inference
        job_number = max(jobs[tensor_number], default=0) + 1
        return f't{tensor_number:03d}.r{job_number:03d}_{suffix}'

    def _fix_info_file(self):
        """Fix format problems in the history.info file."""
        # Create some backup of history.info for now, as we don't want
        # to risk completely messing with the user's file. This may be
        # removed in the future once we're sure we do the right thing.
        backup = self._backup_info_file()
        self.info.fix()
        if not backup:  # There was no history.info to back up
            return False
        if backup.read_text(encoding='utf-8') == self.info.raw_contents:
            # We haven't modified anything. No need to keep the backup.
            backup.unlink()
            return False
        LOGGER.info(
            f'The original {self.info.path.name} file has been renamed to '
            f'{backup.name}. You can safely delete it if no unexpected '
            'changes were introduced. Please open an Issue under '
            'https://github.com/viperleed/viperleed/issues if you '
            'experience any unexpected change.'
            )
        return True

    def _fix_subfolders(self):
        """Fix issues found in all history subfolders."""
        # Fix all folders: those that don't need a fix return empty
        folder_fix = {folder: folder.fix() for folder in self._subfolders}
        fixed_anything = any(actions for actions in folder_fix.values())

        # Now all folders should have metadata. Those to which
        # metadata was added must also be correctly marked into
        # parent/child relationships. Do so now.
        added_metadata = [f for f, actions in folder_fix.items()
                          if FolderFixAction.ADD_METADATA in actions]
        for folder in self._mark_subfolders_as_siblings(added_metadata):
            folder_fix[folder].remove(FolderFixAction.ADD_METADATA)

        folder_fix = {f: a for f, a in folder_fix.items() if a}
        folder_fix_actions = {a for actions in folder_fix.values()
                              for a in actions}
        for action in sorted(folder_fix_actions, key=attrgetter('name')):
            if action.value:
                LOGGER.info(action.value)
        return fixed_anything

    def _mark_subfolders_as_siblings(self, folders):
        """Infer relations between folders based on their timestamp."""
        # Pick first all the "main" folders: they have names with
        # format tTTT.rRRR_* and don't have an additional .SSS
        main_folder_re = re.compile(r't\d{3,}\.r\d{3,}_')

        # We can safely assign relationships only for the folders that
        # clearly come from calc, i.e., those that have no 'moved' in
        # their timestamp, as these are marked as such by bookkeeper.
        main_folders = [f for f in folders
                        if main_folder_re.match(f.name)
                        and MOVED_LABEL not in f.timestamp]

        # Also, we can only process main folders as long as there
        # is a unique timestamp: having two folders with the same
        # timestamp makes it impossible to discern which other
        # folders were created with these.
        n_timestamps = Counter(f.timestamp for f in main_folders)
        main_folder_timestamps = {
            f.timestamp: f
            for f in main_folders
            # pylint: disable-next=magic-value-comparison  # Clear
            if n_timestamps[f.timestamp] < 2
            }

        # Now we can mark siblings of these main folders
        marked_folders = set(main_folder_timestamps.values())
        for folder in folders:
            if folder in marked_folders:
                continue
            try:
                parent = main_folder_timestamps[folder.timestamp]
            except KeyError:
                continue
            folder.metadata.collect_from(parent.metadata)
            folder.metadata.write()
            marked_folders.add(folder)
            # Update also the maps
            # pylint: disable-next=unsubscriptable-object   # Inference
            self._maps['main_hash_to_folders'][parent.hash_].append(folder)
            # pylint: disable-next=unsupported-delete-operation
            del self._maps['main_hash_to_folders'][folder.hash_]
            # pylint: disable-next=unsupported-assignment-operation
            self._maps['hash_to_parent'][folder.hash_] = parent
        return marked_folders

    def _update_maps(self):
        """Update inner mappings for fast access."""
        self._maps['hash_to_parent'] = to_parent = {}
        for folder_group in self._maps['main_hash_to_folders'].values():
            parent = next(f for f in folder_group if not f.parent)
            hashes = (f.hash_ for f in folder_group)
            to_parent.update(zip(hashes, (parent,)*len(folder_group)))
