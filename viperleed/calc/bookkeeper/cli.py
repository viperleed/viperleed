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

from pathlib import Path

from viperleed.calc import DEFAULT_HISTORY
from viperleed.calc import DEFAULT_WORK_HISTORY
from viperleed.cli_base import ViPErLEEDCLI

from .bookkeeper import Bookkeeper
from .constants import HISTORY_INFO_NAME
from .mode import BookkeeperMode


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
            action='store_true',
            )
        what_next.add_argument(
            '-c', '--clear',
            help=('Clear the input directory and add last run to history if not'
                  ' already there. Runs before viperleed.calc by default.'),
            action='store_true',
            )
        what_next.add_argument(
            '-d', '--discard',
            help=('Discard all results from the last run, and restore the '
                  'previous inputs. The discarded run is kept in history.'),
            action='store_true',
            )
        what_next.add_argument(
            '-df', '--discard-full',
            help=('Discard all results from the last run as if it never '
                  'happened. The discarded run is removed from history.'),
            action='store_true',
            )
        parser.add_argument(
            '-j', '--job-name',
            help=('define a string to be appended to the name of the history '
                  f'folder that is created, and is logged in {HISTORY_INFO_NAME}'),
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


    def __call__(self, args=None):
        """Call the bookkeeper with command-line args."""
        parsed_args = self.parse_cli_args(args)

        bookkeeper = Bookkeeper(job_name=parsed_args.job_name,
                                history_name=parsed_args.history_name,
                                work_history_name=parsed_args.work_history_name,
                                cwd=Path.cwd().resolve())
        # Select mode
        if parsed_args.clear:
            mode = BookkeeperMode.CLEAR
        elif parsed_args.discard:
            mode = BookkeeperMode.DISCARD
        elif parsed_args.discard_full:
            mode = BookkeeperMode.DISCARD_FULL
        else:
            mode = BookkeeperMode.ARCHIVE
        return bookkeeper.run(mode)
