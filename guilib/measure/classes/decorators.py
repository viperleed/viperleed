"""Module decorators of viperleed.guilib.measure.classes

========================================
   ViPErLEED Graphical User Interface
========================================

Created: 2024-10-16
Author: Michele Riva
Author: Florian Doerr

Defines decorators used in the other ViPErLEED modules.
"""

from functools import wraps

from viperleed.guilib.measure.classes.abc import QObjectSettingsErrors
from viperleed.guilib.measure.classes.settings import DefaultSettingsError
from viperleed.guilib.measure.hardwarebase import emit_error


def emit_default_faulty(func):
    """Emit an error_occurred when a _defaults settings file has problems."""
    @wraps(func)
    def _wrapper(self, *args, **kwargs):
        try:
            return func(self, *args, **kwargs)
        except DefaultSettingsError as exc:
            emit_error(self, QObjectSettingsErrors.DEFAULT_SETTINGS_CORRUPTED,
                       exc)
            raise
    return _wrapper
