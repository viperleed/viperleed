"""Module decorators of viperleed.gui.measure.classes

Defines decorators used in the other ViPErLEED modules.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    'Florian Dörr (@FlorianDoerr)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2024-10-16'
__license__ = 'GPLv3+'

from functools import wraps

from viperleed.gui.measure.classes.abc import QObjectSettingsErrors
from viperleed.gui.measure.classes.settings import DefaultSettingsError
from viperleed.gui.measure.hardwarebase import emit_error


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
