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


class TimestampFormat(Enum):
    """A collection of known formats for timestamps."""

    # Underscore names are used only for conversion purposes
    # while parsing strings, and never for formatting.
    GERMAN = '%d.%m.%y %H:%M:%S'  # Format used for < v1.0.0
    ISO = '%Y-%m-%d %H:%M:%S'     # Format used for >= v1.0.0
    _CALC = '%y%m%d-%H%M%S'       # Format used by calc in log files
    DEFAULT = ISO

    @property
    def writable(self):
        """Return whether this format can be used for writing."""
        return not self.name.startswith('_')
