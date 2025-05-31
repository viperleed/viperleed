"""Module errors of viperleed.gui.dialogs."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    'Florian Dörr (@FlorianDoerr)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2024-10-15'
__license__ = 'GPLv3+'


class DialogError(Exception):
    """Base exception for all dialog errors."""


class DialogDismissedError(DialogError):
    """Raised if the dialog was discarded without a conclusive result."""
