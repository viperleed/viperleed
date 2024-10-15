"""Module errors of viperleed.guilib.dialogs.

========================================
   ViPErLEED Graphical User Interface
========================================

Created: 2024-10-15
Author: Michele Riva
Author: Florian Doerr

Defines exceptions used in the other ViPErLEED modules.
"""


class DialogError(Exception):
    """Base exception for all dialog errors."""


class DialogDismissedError(DialogError):
    """Raised if the dialog was discarded without a conclusive result."""
