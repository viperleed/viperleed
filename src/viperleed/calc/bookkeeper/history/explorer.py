"""Module explorer of viperleed.calc.bookkeeper.history.

Defines functionality useful for navigating the history folder.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-10-14'
__license__ = 'GPLv3+'

from collections import defaultdict
from operator import attrgetter

from viperleed.calc import DEFAULT_HISTORY

from ..utils import make_property
from ..utils import needs_update_for_attr
from .folder import HistoryFolder
from .folder import IncompleteHistoryFolder
from .info import HistoryInfoFile


class HistoryExplorer:
    """A file-system explorer for the 'history' folder."""

    _updater = 'collect_subfolders'

    def __init__(self, root):
        """Initialize a history explorer at a root path."""
        self._path = root / DEFAULT_HISTORY
        # While the history.info file is technically in the root
        # directory and not in the history directory, it makes
        # sense to keep it in the same object that handles the
        # history directory since they are so much interrelated.
        self._info = None
        self._subfolders = []
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
    @needs_update_for_attr('_maps[jobs_for_tensor]', updater=_updater)
    def max_run_per_tensor(self):
        """Return a {tensor_num: max_job_num} dictionary."""
        map_ = self._maps['jobs_for_tensor']
        max_jobs = {t_num: max(job_nums)
                    for t_num, job_nums in map_.items()
                    if job_nums}
        return defaultdict(int, max_jobs)

    def collect_subfolders(self):
        """Navigate the main history folder and pull subfolder information."""
        self._subfolders.clear()
        self._maps['jobs_for_tensor'] = defaultdict(set)
        self._maps['main_hash_to_folders'] = defaultdict(list)
        self._maps['hash_to_parent'] = {}
        # pathlib.Path.iterdir has no guarantees on the sorting order
        # with which the directories are returned. We use alphabetical
        # sorting by name. This way, the most recent folders are read
        # last.
        contents = sorted(self.path.iterdir(), key=attrgetter('name'))
        for directory in contents:
            try:
                self._append_existing_folder(directory, insert_sorted=False)
            except ValueError: # Not a directory or not a valid name
                pass
        self._update_maps()

    def find_new_history_directory(self, tensor_number, suffix):
        """Store information to create a history folder for a new calc run."""
        dir_name = self._find_name_for_new_history_subfolder(tensor_number,
                                                             suffix)
        new_folder_path = self.path / dir_name
        self._new_calc_run_folder = IncompleteHistoryFolder(new_folder_path)

    def prepare_info_file(self):
        """Prepare a history.info file in the root directory."""
        self._info = HistoryInfoFile(self.root, create_new=True)

    @needs_update_for_attr('_maps[main_hash_to_folders]', updater=_updater)
    def register_folder(self, path_to_folder):
        """Register a new folder as part of the history tree."""
        self._append_existing_folder(path_to_folder, insert_sorted=True)
        self._update_maps()

    def _append_existing_folder(self, path_to_folder, insert_sorted=True):
        """Register a new folder, without updating the parent mapping."""
        if path_to_folder.parent != self.path:
            raise ValueError(f'Not a subfolder of {self.path}')
        folder = HistoryFolder(path_to_folder)
        parent_hash = folder.parent or folder.hash_
        self._subfolders.append(folder)
        if insert_sorted:
            self._subfolders.sort(key=attrgetter('path.name'))

        # pylint: disable=unsubscriptable-object  # Called in decorated
        self._maps['main_hash_to_folders'][parent_hash].append(folder)
        self._maps['jobs_for_tensor'][folder.tensor_num].add(folder.job_num)

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
        try:
            # pylint: disable-next=unsubscriptable-object   # Inference
            job_number = max(jobs[tensor_number])
        except ValueError:  # No jobs for this tensor yet
            job_number = 0
        dir_name_fmt = f't{tensor_number:03d}.r{{job:03d}}_{suffix}'
        # If there is already a folder with the same name and correct
        # suffix, we take that, otherwise, we increase the job number:
        # it's a new run.
        base_name = dir_name_fmt.format(job=job_number)
        if not (self.path / base_name).is_dir():
            base_name = dir_name_fmt.format(job=job_number + 1)
        return base_name

    def _update_maps(self):
        """Update inner mappings for fast access."""
        self._maps['hash_to_parent'] = to_parent = {}
        for folder_group in self._maps['main_hash_to_folders'].values():
            parent = next(f for f in folder_group if not f.parent)
            hashes = (f.hash_ for f in folder_group)
            to_parent.update(zip(hashes, (parent,)*len(folder_group)))
