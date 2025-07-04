"""Module workhistory of viperleed.calc.bookkeeper.history.

Defines the WorkhistoryHandler class that takes care of cleaning up
folders found in the workhistory folder, i.e., those that calc has
moved out of work as a result of repeated runs (e.g., looped searches
that required the production of multiple delta-amplitude files, or
a RUN = 1-3 1).

This module comes from the refactor of the bookkeeper module by
@michele-riva in October 2024.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2020-01-30'
__license__ = 'GPLv3+'

import logging
from operator import attrgetter
import shutil

from viperleed.calc.bookkeeper.constants import HISTORY_FOLDER_RE
from viperleed.calc.bookkeeper.history.meta import BookkeeperMetaFile
from viperleed.calc.bookkeeper.utils import make_property
from viperleed.calc.constants import DEFAULT_WORK_HISTORY
from viperleed.calc.sections.cleanup import PREVIOUS_LABEL


LOGGER = logging.getLogger(__name__)


class WorkhistoryHandler:
    """A class that takes care of the workhistory folder."""

    def __init__(self, root, bookkeeper):
        """Initialize instance.

        Parameters
        ----------
        root : RootExplorer
            The handler of the root folder that contains the
            'workhistory' directory on which this instance
            will work.
        bookkeeper : Bookkeeper
            The Bookkeeper instance associated
            with this WorkhistoryHandler.

        Returns
        -------
        None.
        """
        self._root = root
        self._path = root.path / DEFAULT_WORK_HISTORY
        self.bookkeeper = bookkeeper

    path = make_property('_path')
    root = make_property('_root.path')

    @property
    def history(self):
        """Return the path to the 'history' directory."""
        return self._root.history.path

    @property
    def timestamp(self):
        """Return the time-stamp of the current run."""
        return self.bookkeeper.timestamp

    def discard_workhistory_root(self):
        """Remove the top-level workhistory folder and all its contents."""
        if not self.path.is_dir():
            return
        try:
            shutil.rmtree(self.path)
        except OSError:
            is_empty = not any(self.path.iterdir())
            err_ = f'Failed to {{}} {self.path.name} folder'
            err_ = err_.format('delete empty' if is_empty else 'discard')
            LOGGER.error(err_, exc_info=True)

    def find_current_directories(self, contains=''):
        """Return a generator of subfolders in the current workhistory folder.

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
            matches the HISTORY_FOLDER_RE regular expression, contains
            `contains`, but does not contain 'previous'.
        """
        directories = self._find_directories(contains=contains)
        return (d for d in directories if PREVIOUS_LABEL not in d.name)

    def list_paths_to_discard(self):
        """Return a tuple of paths to folders that will be discarded."""
        if not self.path.is_dir():
            return tuple()
        subfolders = sorted(self.path.iterdir(), key=attrgetter('name'))
        return (self.path, *subfolders)

    def move_current_and_cleanup(self, main_history_subfolder):
        """Move files from the current workhistory folder, then clean up.

        Any subfolder of workhistory that is labeled as "previous"
        (i.e., corresponding to an earlier run of viperleed.calc) is
        always removed. If the current workhistory folder is empty
        (after removal of the "previous" folders and moving the new
        ones), it is deleted.

        Parameters
        ----------
        main_history_subfolder : HistoryFolder
            The primary history subfolder created by bookkeeper. Used
            to determine the appropriate name for the the workhistory
            folders moved to history as well as marking them as being
            created 'together' with the `main_history_subfolder`.

        Returns
        -------
        tensor_nums : set of int
            Indices of tensors found in the current workhistory folder
            that have been moved to history as new history entries.
        """
        if not self.path.is_dir():
            return set()

        # Always remove any 'previous'-labeled folders
        self._discard_previous()
        tensor_nums = self._move_folders_to_history(main_history_subfolder)
        is_empty = not any(self.path.iterdir())
        if is_empty:
            self.discard_workhistory_root()
        return tensor_nums

    def _discard_previous(self):
        """Remove 'previous'-labeled directories in workhistory."""
        previous = self._find_directories(contains=PREVIOUS_LABEL)
        for directory in previous:
            try:
                shutil.rmtree(directory)
            except OSError:
                LOGGER.error(f'Failed to delete {directory.name} directory '
                             f'from {self.path.name}', exc_info=True)

    def _find_directories(self, contains=''):
        """Return a generator of subfolders in the current workhistory folder.

        Parameters
        ----------
        contains : str, optional
            Select only those subfolders whose name contains this
            string. Default is an empty string, corresponding to no
            filtering other than the one described in Returns.

        Returns
        -------
        subfolders : generator
            When iterated over, it yields paths to the immediate
            subfolders of the current workhistory folder whose
            name matches the HISTORY_FOLDER_RE regular expression,
            and whose name includes `contains`.
        """
        globbed = (self.path.glob(f'*{contains}*') if contains
                   else self.path.iterdir())
        return (d for d in globbed
                if d.is_dir() and HISTORY_FOLDER_RE.match(d.name))

    def _move_folders_to_history(self, main_history_subfolder):
        """Move relevant folders from the current workhistory to history.

        Parameters
        ----------
        main_history_subfolder : HistoryFolder
            The primary history subfolder created by bookkeeper. Used
            to determine the appropriate name for the the workhistory
            folders moved to history as well as marking them as being
            created 'together' with the `main_history_subfolder`.

        Returns
        -------
        tensor_nums : set
            The tensor numbers of the relevant workhistory folders.

        Raises
        ------
        FileExistsError
            If moving one of the workhistory directories fails because
            there already is a history directory with the same name.
        """
        tensor_nums = set()
        # NB: it is important to create the tensor-to-job map before
        # we start adding more history folders: all of those added
        # at the same time must have the same job_num.
        max_job_for_tensor = self._root.history.max_run_per_tensor
        directories = self.find_current_directories(contains=self.timestamp)
        # Workhistory directories have the following naming convention
        # [see cleanup.move_oldruns(rpars, prerun=False)]:
        #    tTTT.rSSS[_<short_labels_of_sections>]_<log_timestamp>
        # Here we 'shift' the SSS part further down the line, to make
        # folders of the form
        #    tTTT.rRRR.SSS[_<short_labels_of_sections>]_<log_timestamp>
        for directory in directories:
            match_ = HISTORY_FOLDER_RE.match(directory.name)
            tensor_num = int(match_['tensor_num'])
            search_num = int(match_['job_num'])  # Misuse the job_num
            job_num = max_job_for_tensor[tensor_num]
            if tensor_num != main_history_subfolder.tensor_num:
                # Workhistory is always processed after the primary
                # history folder, so if tensor_num is the same as for
                # the main folder, also job_num should be, i.e., we
                # should not increase job_num by one.
                job_num += 1
            newname = (
                f't{tensor_num:03d}.r{job_num:03d}.{search_num:03d}'
                + match_['rest']
                )
            target = self.history / newname
            try:
                # NB: using fs_util.move would not help to discern
                # the FileExistsError case below. If target exists,
                # directory would be moved inside target instead of
                # replacing it.
                directory.replace(target)
            except OSError as exc:
                err_msg = ('Error: Failed to move '
                           f'{directory.relative_to(self.root)} '
                           f'to {target.relative_to(self.root)}')
                if target.is_dir():
                    err_msg += ': Target path already exists. Stopping...'
                    LOGGER.error(err_msg)
                    raise FileExistsError(str(exc)) from exc
                LOGGER.error(err_msg, exc_info=True)
                continue
            tensor_nums.add(tensor_num)

            # Add metadata file to the folder we have just moved.
            meta = BookkeeperMetaFile(target)
            meta.compute_hash()
            meta.collect_from(main_history_subfolder.metadata)
            meta.write()

            # And register the folder in history. Notice that we need
            # the metadata file to already be written to successfully
            # register a folder.
            self._root.history.register_folder(target)
        return tensor_nums
