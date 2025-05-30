"""Module cli of ViPErLEED POSCAR utilities.

Defines the main-entry-point function for the POSCAR utilities
and other command-line related functionality.

The functionality in this module used to be in __main__.py.
"""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-08-04'
__license__ = 'GPLv3+'


from viperleed.cli_base import ViPErLEEDCLIWithAutoChildren


class PoscarUtilsCLI(ViPErLEEDCLIWithAutoChildren,
                     cli_name='poscar',
                     help_='utilities for POSCAR files'):
    """Command-line interface for all POSCAR utilities."""

    def add_parser_arguments(self, parser):
        """Add POSCAR arguments to `parser`."""
        self.add_verbose_option(parser)
        super().add_parser_arguments(parser)
