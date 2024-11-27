"""Module cli of viperleed.calc.bookkeeper.

Defines the BookkeeperCLI class, the main command-line interface
for executing the bookkeeper.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2020-01-30'
__license__ = 'GPLv3+'

from argparse import Action
from pathlib import Path

from viperleed.calc.constants import DEFAULT_OUT
from viperleed.calc.constants import DEFAULT_HISTORY
from viperleed.cli_base import ViPErLEEDCLI

from .bookkeeper import Bookkeeper
from .history.constants import HISTORY_INFO_NAME
from .mode import BookkeeperMode


# All the pylint disables below are needed as Action subclasses
# must have this interface. The problem rather lies in argparse.
# pylint: disable-next=too-few-public-methods
class StoreBookkeeperMode(Action):
    """An action to select the correct BookkeeperMode.

    It behaves as 'store_true', making the arguments to which
    it is associated not require any other information.
    """

    # pylint: disable-next=too-many-arguments
    def __init__(self,
                 option_strings,
                 dest,
                 default=False,
                 required=False,
                 help=None):  # pylint: disable=redefined-builtin
        """Initialize instance to behave the same as 'store_true'."""
        super().__init__(
            option_strings=option_strings,
            dest=dest,
            nargs=0,
            const=True,
            default=default,
            type=None,
            choices=None,
            required=required,
            help=help,
            metavar=None,
            )

    def __call__(self, parser, args, values, option_string=None):
        """Set args.mode to the right BookkeeperMode."""
        try:
            mode = BookkeeperMode(self.dest)
        except ValueError:
            parser.error(f'Unknown bookkeeper mode {self.dest!r}')
        setattr(args, 'mode', mode)
        setattr(args, self.dest, self.const)


class BookkeeperCLI(ViPErLEEDCLI, cli_name='bookkeeper'):
    """The main command-line interface for the bookkeeper utility."""

    def add_parser_arguments(self, parser):
        """Add bookkeeper arguments to parser."""
        super().add_parser_arguments(parser)
        what_next = parser.add_mutually_exclusive_group()
        what_next.add_argument(
            *BookkeeperMode.ARCHIVE.flags,
            help=('Store last run in history. Overwrite PARAMETERS, POSCAR &'
                  f'VIBROCC from {DEFAULT_OUT}. Runs after viperleed.calc by '
                  'default.'),
            action=StoreBookkeeperMode,
            )
        what_next.add_argument(
            *BookkeeperMode.CLEAR.flags,
            help=('Clear the input directory and add last run '
                  'to history if not already there. Runs before '
                  'viperleed.calc by default.'),
            action=StoreBookkeeperMode,
            )
        what_next.add_argument(
            *BookkeeperMode.DISCARD.flags,
            help=('Discard all results from the last run, and restore the '
                  'previous inputs. The discarded run is kept in history.'),
            action=StoreBookkeeperMode,
            )
        what_next.add_argument(
            *BookkeeperMode.DISCARD_FULL.flags,
            help=('Discard all results from the last run as if it never '
                  'happened. The discarded run is removed from history.'),
            action=StoreBookkeeperMode,
            )
        what_next.add_argument(
            *BookkeeperMode.FIX.flags,
            help=(f'Automatically fix problems found in the {DEFAULT_HISTORY} '
                  f'directory and in the {HISTORY_INFO_NAME} file. This mode '
                  'can be used, for example, to upgrade the state of your '
                  f'{DEFAULT_HISTORY} to the format of the most recent '
                  'version of the bookkeeper.'),
            action=StoreBookkeeperMode,
            )
        parser.add_argument(
            '-y',
            help=('Do not ask for confirmation when running in '
                  f'{BookkeeperMode.DISCARD_FULL.long_flag} mode.'),
            action='store_true',
            dest='skip_confirmation',
            )

    def __call__(self, args=None):
        """Call the bookkeeper with command-line args."""
        parsed_args = self.parse_cli_args(args)
        bookkeeper = Bookkeeper(cwd=Path.cwd().resolve())
        mode = getattr(parsed_args, 'mode', BookkeeperMode.ARCHIVE)
        kwargs = {
            'requires_user_confirmation': not parsed_args.skip_confirmation,
            }
        return bookkeeper.run(mode, **kwargs)
