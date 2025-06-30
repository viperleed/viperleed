"""ViPErLEED utility: Sort slab by z."""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-08-03'
__license__ = 'GPLv3+'

from viperleed.utilities.poscar.base import _PoscarStreamCLI


class SortByZCLI(_PoscarStreamCLI, cli_name='sort_by_z'):
    """Sort atoms by ascending or descending z coordinate."""

    long_name = 'sort slab by z'

    def add_parser_arguments(self, parser):
        """Add optional --reversed argument."""
        super().add_parser_arguments(parser)
        parser.add_argument(
            '-r', '--reversed',
            help='sort from bottom to top. Default is from top to bottom.',
            action='store_true',
            )

    def process_slab(self, slab, args):
        """Return a slab z-sorted according to args."""
        slab.sort_by_z(bottom_to_top=args.reversed)
        return slab


if __name__ == '__main__':
    SortByZCLI.run_as_script()
