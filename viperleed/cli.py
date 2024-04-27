"""Module cli of viperleed.

Defines the main command-line interface for running viperleed
and its subpackages.
"""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-08-01'
__license__ = 'GPLv3+'

from viperleed.cli_base import ViPErLEEDCLI


class ViPErLEEDMain(ViPErLEEDCLI, cli_name='viperleed'):
    """The main CLI interface of viperleed."""

    def __init__(self, *args, **kwargs):
        """Initialize instance by registering sub-utilities."""
        super().__init__(*args, **kwargs)
        children = (
            # Notice that, because of a current fuckery with the way
            # the imports are in guilib, viperleed.gui should be the
            # first one, so that it is not imported later. Otherwise,
            # "python viperleed gui" thinks it is in command-line
            # mode, while "python viperleed.gui" works fine.
            'viperleed.gui',                                                    # TODO: gui arguments
            'viperleed.calc.bookkeeper',
            'viperleed.calc.cli',
            'viperleed.utilities.cli',
            'viperleed.utilities.poscar.cli',
            )
        self.register_valid_children(children)
        self.add_child_aliases('utilities', 'utils', 'util')
