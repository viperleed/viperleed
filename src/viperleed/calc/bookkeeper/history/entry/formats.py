"""Module entry of viperleed.calc.bookkeeper.history.entry.

Defines classes that handle the contents of a single 'block'
of the history.info file.
"""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2020-01-30'
__license__ = 'GPLv3+'

from enum import Enum

from viperleed.calc.lib.time_utils import DateTimeFormat
from viperleed.calc.lib.time_utils import now_


class TimestampFormat(Enum):
    """A collection of known formats for timestamps."""

    # Underscore names are used only for conversion purposes
    # while parsing strings, and never for formatting.
    GERMAN = '%d.%m.%y %H:%M:%S'              # < v1.0.0
    ISO = DateTimeFormat.ISO.value            # >= v1.0.0
    _CALC = DateTimeFormat.FILE_SUFFIX.value  # Format of log names
    DEFAULT = ISO  # Just an alias of one of the above!

    @property
    def writable(self):
        """Return whether self is suitable for writing to history.info."""
        return not self.name.startswith('_')

    def now(self):
        """Return the local time formatted according to this format."""
        return now_(self.value, use_gmt=False)
