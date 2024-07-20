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

from viperleed.calc import DEFAULT_HISTORY
from viperleed.calc import DEFAULT_WORK_HISTORY
from viperleed.cli_base import ViPErLEEDCLI

from .bookkeeper import Bookkeeper
from .constants import HISTORY_INFO_NAME
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
        except ValueError:                                                      # TODO: untested
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
            '-a', '--archive',
            help=('Store last run in history. Overwrite PARAMETERS, POSCAR &'
                  'VIBROCC from OUT. Runs after viperleed.calc by default.'),
            action=StoreBookkeeperMode,
            )
        what_next.add_argument(
            '-c', '--clear',
            help=('Clear the input directory and add last run '
                  'to history if not already there. Runs before '
                  'viperleed.calc by default.'),
            action=StoreBookkeeperMode,
            )
        what_next.add_argument(
            '-d', '--discard',
            help=('Discard all results from the last run, and restore the '
                  'previous inputs. The discarded run is kept in history.'),
            action=StoreBookkeeperMode,
            )
        what_next.add_argument(
            '-df', '--discard-full',
            help=('Discard all results from the last run as if it never '
                  'happened. The discarded run is removed from history.'),
            action=StoreBookkeeperMode,
            )
        parser.add_argument(
            '-j', '--job-name',
            help=('define a string to be appended to the name '
                  'of the history folder that is created, and '
                  f'is logged in {HISTORY_INFO_NAME}'),
            type=str
            )
        parser.add_argument(
            '--history-name',
            help=('define the name of the history folder that is '
                  f'created/used. Default is {DEFAULT_HISTORY!r}'),
            type=str,
            default=DEFAULT_HISTORY
            )
        parser.add_argument(
            '--work-history-name',
            help=('define the name of the workhistory folder that is '
                  f'created/used. Default is {DEFAULT_WORK_HISTORY!r}'),
            type=str,
            default=DEFAULT_WORK_HISTORY
            )

    def __call__(self, args=None):                                              # TODO: untested
        """Call the bookkeeper with command-line args."""
        parsed_args = self.parse_cli_args(args)
        bookkeeper = Bookkeeper(
            job_name=parsed_args.job_name,
            history_name=parsed_args.history_name,
            work_history_name=parsed_args.work_history_name,
            cwd=Path.cwd().resolve()
            )
        mode = getattr(parsed_args, 'mode', BookkeeperMode.ARCHIVE)
        return bookkeeper.run(mode)
