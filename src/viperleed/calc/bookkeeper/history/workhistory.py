"""Module workhistory.py of viperleed.calc.bookkeeper.history.

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
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2020-01-30'
__license__ = 'GPLv3+'

import shutil

from viperleed.calc.sections.cleanup import PREVIOUS_LABEL

from ..constants import HISTORY_FOLDER_RE
from ..log import LOGGER


class WorkhistoryHandler:
    """A class that takes care of the workhistory folder."""

    def __init__(self, work_history_path, bookkeeper):
        """Initialize an instance to handle `work_history_path`.

        Parameters
        ----------
        work_history_path : Path
            The path to the 'workhistory' folder
            that this instance will work on.
        bookkeeper : Bookkeeper
            The Bookkeeper instance associated
            with this WorkhistoryHandler.

        Returns
        -------
        None.
        """
        self.path = work_history_path
        self.bookkeeper = bookkeeper

    @property
    def history(self):
        """Return the path to the 'history' directory."""
        return self.bookkeeper.top_level_history_path

    @property
    def timestamp(self):
        """Return the time-stamp of the current run."""
        return self.bookkeeper.timestamp

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

    def move_and_cleanup(self, discard):
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
        tensor_nums = set()
        if not self.path.is_dir():
            return tensor_nums

        # Always remove any 'previous'-labelled folders
        self._discard_previous()
        if not discard:
            tensor_nums = self._move_folders_to_history()
        is_empty = not any(self.path.iterdir())
        if is_empty or discard:
            try:
                shutil.rmtree(self.path)
            except OSError:
                err_ = f'Failed to {{}} {self.path.name} folder'
                err_ = err_.format('discard' if discard else 'delete empty')
                LOGGER.error(err_, exc_info=True)
        return tensor_nums

    def _discard_previous(self):
        """Remove 'previous'-labelled directories in work history."""
        previous = self._find_directories(contains=PREVIOUS_LABEL)
        for directory in previous:
            try:
                shutil.rmtree(directory)
            except OSError:
                LOGGER.error(f'Failed to delete {directory} directory '
                             f'from {self.path}', exc_info=True)

    def _find_directories(self, contains=''):
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
            name matches the HISTORY_FOLDER_RE regular expression,
            and whose name include `contains`.
        """
        globbed = (self.path.glob(f'*{contains}*') if contains
                   else self.path.iterdir())
        return (d for d in globbed
                if d.is_dir() and HISTORY_FOLDER_RE.match(d.name))

    def _move_folders_to_history(self):
        """Move relevant folders from the current work history to history.

        Returns
        -------
        tensor_nums : set
            The tensor numbers of the folders that have been moved.

        Raises
        ------
        FileExistsError
            If moving one of the work-history directories fails because
            there already is a history directory with the same name.
        """
        tensor_nums = set()
        max_job_for_tensor = self.bookkeeper.max_job_for_tensor
        directories = self.find_current_directories(contains=self.timestamp)
        # Work-history directories have the following naming convention
        # [see cleanup.move_oldruns(rp, prerun=False)]:
        #    tTTT.rSSS[_<short_labels_of_sections>]_<log_timestamp>
        # Here we 'shift' the SSS part further down the line, to make
        # folders of the form
        #    tTTT.rRRR.SSS[_<short_labels_of_sections>]_<log_timestamp>
        for directory in directories:
            match = HISTORY_FOLDER_RE.match(directory.name)
            tensor_num = int(match['tensor_num'])
            search_num = int(match['job_num'])  # Misuse the job_num
            job_num = max_job_for_tensor[tensor_num] + 1
            newname = (
                f't{tensor_num:03d}.r{job_num:03d}.{search_num:03d}'
                + match['rest']
                )
            target = self.history / newname
            try:
                directory.replace(target)
            except FileExistsError:
                LOGGER.error(f'Error: Failed to move {directory} to {target}: '
                             'Target path already exists. Stopping...')
                raise
            except OSError:
                LOGGER.error(f'Error: Failed to move {directory}.',
                             exc_info=True)
                continue
            tensor_nums.add(tensor_num)
        return tensor_nums
