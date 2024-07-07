"""Module cli of viperleed.utilities.

Defines the command-line interface for the ViPErLEED utilities.
"""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-08-02'
__license__ = 'GPLv3+'


from viperleed.cli_base import ViPErLEEDCLIWithAutoChildren


class ViPErLEEDUtilitiesCLI(ViPErLEEDCLIWithAutoChildren,
                            cli_name='utilities',
                            help_='run ViPErLEED utilities'):
    """The main command-line interface for the ViPErLEED utilities."""

    def __init__(self):
        """Initialize CLI with some aliases."""
        super().__init__()
        self.add_child_aliases('rearrange-phaseshifts',
                               'rearrange-phaseshift',
                               'rearrange_phaseshifts',
                               'rearrange_phaseshift')
