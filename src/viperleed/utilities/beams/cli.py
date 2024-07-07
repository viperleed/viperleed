"""Module cli of viperleed.utilities.beams.

Defines the main command-line interface for utilities that convert
various experimental beam formats into the EXPBEAMS.csv files used
by ViPErLEED.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-03-27'
__license__ = 'GPLv3+'


from viperleed.cli_base import ViPErLEEDCLIWithAutoChildren


_HELP = 'convert experimental beams file to EXPBEAMS.csv format'

class BeamsUtilCLI(ViPErLEEDCLIWithAutoChildren,
                   cli_name='beams', help_=_HELP):
    """The main command-line interface for the beams-conversion utility."""

    def __init__(self):
        """Initialize CLI with case-insensitive aliases."""
        super().__init__()
        self.add_child_aliases('from-SATLEED', *_make_aliases('from-SATLEED'))
        self.add_child_aliases('from-TensErLEED',
                               *_make_aliases('from-TensErLEED'),
                               'from-ErLEED',
                               *_make_aliases('from-ErLEED'))


def _make_aliases(name):
    """Return aliases for name."""
    underscore = name.replace('-', '_')
    return name.lower(), underscore, underscore.lower()
