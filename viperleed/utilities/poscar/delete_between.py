"""ViPErLEED utility: Delete atoms between certain c fractions."""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Florian Kraushofer (@fkraushofer)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-08-03'
__license__ = 'GPLv3+'

from viperleed.cli_base import float_in_zero_one
from viperleed.utilities.poscar.base import _RemoveAtomsCLI


class DeleteBetweenCLI(_RemoveAtomsCLI, cli_name='delete_between'):
    """Remove atoms between a pair of c fractions."""

    long_name = 'delete atoms between'

    def add_remove_condition_arguments(self, parser):
        """Add c fraction limits to parser."""
        parser.add_argument(
            'c',
            help='delete all atoms at or between these c fractions',
            type=float_in_zero_one,
            nargs=2,
            )

    def process_slab(self, slab, args):
        """Remove atoms above the c fraction specified in the CLI."""
        modified_slab = super().process_slab(slab, args)
        logger = self.get_logger()
        min_c, max_c = args.c
        logger.debug(f'Deleted {slab.n_atoms - modified_slab.n_atoms} atoms '
                     f'in the range c = [{min_c:5.3}, {max_c:5.3}].')
        return modified_slab

    def select_surviving_atoms(self, slab, args):
        """Return atoms of slab with c position inside the args.c interval."""
        min_c, max_c = args.c
        return [atom for atom in slab if min_c <= atom.pos[2] <= max_c]

    def parse_cli_args(self, args):
        """Validate c fractions passed after parsing."""
        parsed_args = super().parse_cli_args(args)
        min_c, max_c = parsed_args.c
        if min_c > max_c:
            self.parser.error('c fractions must be sorted '
                              'from smaller to larger')
        return parsed_args


if __name__ == '__main__':
    DeleteBetweenCLI.run_as_script()
