"""Module errors of viperleed.calc.bookkeeper.

Defines the base BookkeeperError exception and its subclasses.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2024-11-29'
__license__ = 'GPLv3+'


class BookkeeperError(Exception):
    """Base class of all bookkeeper-related errors."""


class FileOperationFailedError(BookkeeperError):
    """Something went wrong when moving/copying files."""

    def __init__(self, files_and_info):
        """Initialize instance with a dict of failure information."""
        self.failures = files_and_info
        msgs = (f'{file}: {reason}' for file, reason in self.failures.items())
        super().__init__('\n'.join(msgs))


class NotAnInteractiveShellError(BookkeeperError):
    """Tried to prompt the user in a non-interactive shell."""


class _FileNotOlderError(BookkeeperError):
    """Exception used internally for file-age checks."""
