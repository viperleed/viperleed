"""Module domain_finder of viperleed.calc.bookkeeper."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-04-10'
__license__ = 'GPLv3+'

import logging
from pathlib import Path

from viperleed.calc.bookkeeper.constants import CALC_LOG_PREFIXES
from viperleed.calc.bookkeeper.history.constants import HISTORY_INFO_NAME
from viperleed.calc.bookkeeper.history.explorer import HistoryExplorer
from viperleed.calc.bookkeeper.history.errors import MetadataMismatchError
from viperleed.calc.bookkeeper.utils import ask_user_confirmation
from viperleed.calc.bookkeeper.utils import make_property
from viperleed.calc.bookkeeper.utils import needs_update_for_attr
from viperleed.calc.constants import DEFAULT_DOMAIN_FOLDER_PREFIX
from viperleed.calc.constants import DEFAULT_HISTORY
from viperleed.calc.constants import DEFAULT_OUT
from viperleed.calc.constants import DEFAULT_SUPP
from viperleed.calc.constants import DEFAULT_WORK_HISTORY
from viperleed.calc.lib.log_utils import logging_silent
from viperleed.calc.lib.string_utils import harvard_commas


LOGGER = logging.getLogger(__name__)

# The following are keys in the 'domains' section of the
# metadata files that signal whether a folder belongs to
# the history of the main folder in a domains calculation
# (_MAIN_DOMAIN_KEY) or to that of a subdomain thereof
# (_SUBDOMAIN_KEY). If neither key exists, the folder
# belongs to a "normal" single-domain calc run.
# See history.meta.BookkeeperMetaFile.domains for the
# definition of the keys.
_MAIN_DOMAIN_KEY = 'domains'
_SUBDOMAIN_KEY = 'main'


# Names of folders that are never considered to be the root of
# domain subfolders when searching potential ones.
_NOT_A_DOMAIN_SUBFOLDER = {
    DEFAULT_OUT,
    DEFAULT_SUPP,
    DEFAULT_HISTORY,
    DEFAULT_WORK_HISTORY,
    }


class DomainFinderError(Exception):
    """Base class for errors in DomainFinder."""


class FolderNotFoundError(DomainFinderError):
    """Some expected folder was not found."""

    def __init__(self, path):
        """Initialize instance from a path."""
        self.path = path
        super().__init__(f'Could not find folder at {path}.')


class DomainPathNotFoundError(FolderNotFoundError):
    """The path of a domain subfolder was not found from the main."""


class MainPathNotFoundError(FolderNotFoundError):
    """The main path of the calculation was not found, or is mismatched."""


class DomainFinder:
    """A class for finding domain subfolders of a root directory."""

    def __init__(self, bookkeeper):
        """Initialize an instance from a `bookkeeper`."""
        self._bookkeeper = bookkeeper
        self._domain_info = None

    path = make_property('_bookkeeper.cwd')

    @property
    @needs_update_for_attr('_domain_info', updater='collect_info')
    def is_subdomain(self):
        """Return whether self.path is a subdomain of a main one."""
        return _SUBDOMAIN_KEY in self._domain_info

    @property
    def domain_info(self):
        """Return domain information found in self.path."""
        key = _SUBDOMAIN_KEY if self.is_subdomain else _MAIN_DOMAIN_KEY
        return self._domain_info.get(key)

    def collect_info(self):
        """Collect domain information from self.path."""
        self._bookkeeper.update_from_cwd(silent=True)
        last_folder = self._bookkeeper.history.last_folder
        self._domain_info = (
            {} if last_folder is None
            else last_folder.metadata.domains
            )

    def find_potential_domains(self):
        """Return paths to subfolders of self.path that may be domains."""
        subfolders = (d for d in self.path.iterdir()
                      if d.is_dir()
                      and d.name not in _NOT_A_DOMAIN_SUBFOLDER)
        return tuple(d for d in subfolders
                     if self._is_potential_domain_subfolder(d))

    @needs_update_for_attr('_domain_info', updater='collect_info')
    def find_domains(self, mode):
        """Return domains found in self.path.

        Domain information is collected from the metadata files
        found in the most-recent subfolder of the 'history' directory
        of self.path.

        Parameters
        ----------
        mode : BookkeeperMode
            The mode in which self.bookkeeper is currently executed.
            Only used for error reporting.

        Returns
        -------
        domain_paths : tuple
            Paths to all subdomains of self.path.

        Raises
        ------
        MetadataError
            If any metadata file contains corrupted domain information.
        MainPathNotFoundError
            If we're running in the main directory of a domains
            calculation, but this path is inconsistent with the
            one stored in any of the domain subfolders. This
            exception is not raised if users explicitly confirm
            their will to proceed anyway.
        """
        domain_info = self.domain_info
        if not domain_info:
            return tuple()
        if self.is_subdomain:
            main_path, _ = domain_info
            self._warn_running_in_subdomain(mode, main_path)
            return tuple()
        return self._find_domains_from_main()

    @staticmethod
    def _ask_user_confirmation_on_main_path_mismatch(already_confirmed):
        """Main path is mismatched. Ask confirmation to proceed."""
        return already_confirmed or ask_user_confirmation()

    def _check_subdomain_info_consistent(self, domain_path, domain_hash):
        """Raise if domain information is not up to date."""
        domain_path = Path(domain_path)
        try:
            folder = self._get_history_folder(domain_path, domain_hash)
        except FileNotFoundError:
            LOGGER.warning(f'Could not find expected domain '
                           f'folder {domain_path}.')
            raise DomainPathNotFoundError(domain_path) from None
        if not folder:  # No history folder with that hash
            err_msg = ('Could not find history folder with hash '
                       f'{domain_hash} in domain {domain_path}.')
            LOGGER.warning(err_msg)
            raise MetadataMismatchError(err_msg)
        main_path, _ = folder.metadata.domains[_SUBDOMAIN_KEY]
        if Path(main_path) != self.path:
            LOGGER.error(f'Domain {domain_path}: Inconsistent path to the '
                         f'main folder of the domain calculation. Expected '
                         f'{main_path}, found {self.path} instead.')
            raise MainPathNotFoundError(main_path)

    def _find_domains_from_main(self):
        """Return paths to domain subfolders of the current directory."""
        domain_paths = []
        user_confirmed = False
        for domain in self.domain_info:
            try:
                self._check_subdomain_info_consistent(*domain)
            except DomainPathNotFoundError:
                # TODO: once --fix can take care of renaming of the
                # domain subfolders, we can suggest running it here.
                continue
            except MainPathNotFoundError:
                # We have found the domain subfolder (albeit in the
                # wrong path).
                # TODO: once --fix can take care of renaming of the
                # main folder, we can suggest running it here and
                # raise again. For now make sure users are sure.
                user_confirmed = (
                    self._ask_user_confirmation_on_main_path_mismatch(
                        user_confirmed,
                        )
                    )
                if not user_confirmed:
                    raise
            except MetadataMismatchError:
                # Users must have edited some files in the history
                # folder. We do not allow updating the hashes in the
                # metadata files, not even via --fix. We have found
                # the domain subfolder, though, so we can still
                # return it.
                pass
            domain_paths.append(Path(domain.path))
        return tuple(domain_paths)

    def _find_domains_from_subdomain(self):
        """Return relative domain paths when self.path is a subdomain."""
        main_info = self.domain_info
        try:
            main_folder = self._get_history_folder(*main_info)
        except FileNotFoundError:  # Probably main path renamed
            raise MainPathNotFoundError(main_info.path) from None
        if not main_folder:
            return tuple()
        domain_info = main_folder.metadata.domains.get(_MAIN_DOMAIN_KEY, ())
        return tuple(Path(d.path) for d in domain_info)

    @staticmethod
    def _get_history_folder(root_path, folder_hash):
        """Return the `root_path`/'history' folder with `folder_hash`."""
        root_path = Path(root_path)
        if not root_path.is_dir():
            raise FileNotFoundError
        history = HistoryExplorer(root_path)
        with logging_silent():
            history.collect_subfolders()
        return history.subfolder_from_hash(folder_hash)

    def _is_potential_domain_subfolder(self, path):
        """Return whether `path` points to a (likely) domain subfolder."""
        if path.name.startswith(DEFAULT_DOMAIN_FOLDER_PREFIX):
            return True
        contents = (
            # A folder is likely a domain subfolder if it contains:
            DEFAULT_OUT,            # an OUT subfolder
            DEFAULT_SUPP,           # or a SUPP subfolder
            DEFAULT_WORK_HISTORY,   # or a workhistory subfolder
            # or it has a calc log file (this means calc
            # was run explicitly in the subfolder)
            *(f'{prefix}*.log' for prefix in CALC_LOG_PREFIXES),
            # or it has been processed before by bookkeeper
            DEFAULT_HISTORY,
            HISTORY_INFO_NAME,
            )
        return any(f for pattern in contents for f in path.glob(pattern))

    def _warn_running_in_subdomain(self, mode, main_path):
        """Emit warnings when calling find_domains in a subdomain."""
        log_msg = (
            f'Running bookkeeper in mode {mode.name} in a domain subfolder. '
            )
        try:
            domains = self._find_domains_from_subdomain()
        except MainPathNotFoundError:
            # TODO: add here a suggestion to run in --fix mode in
            # the parent path once we can fix renaming of main
            log_msg += ('The path to the corresponding main directory could '
                        f'not be identified. It should have been {main_path}, '
                        f'but the parent directory is now {self.path.parent}.')
            LOGGER.error(log_msg)
            return

        siblings = (d for d in domains if d.name != self.path.name)
        LOGGER.warning(f'{log_msg} Make sure to also execute \'bookkeeper '
                       f'{mode.long_flag}\' in the main path ({main_path}) '
                       'or in other domain subfolders (%s).',
                       harvard_commas(*siblings))
