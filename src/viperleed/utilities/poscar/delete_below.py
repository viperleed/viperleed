"""ViPErLEED utility: Delete atoms below a certain height."""

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


class DeleteBelowCLI(_RemoveAtomsCLI, cli_name='delete_below'):
    """Remove atoms below a certain c fraction."""

    long_name = 'delete atoms below'

    def add_remove_condition_arguments(self, parser):
        """Add maximum c fraction to parser."""
        parser.add_argument(
            'c',
            help='delete all atoms strictly below this c fraction',
            type=float_in_zero_one,
            )

    def process_slab(self, slab, args):
        """Remove atoms below the c fraction specified in the CLI."""
        modified_slab = super().process_slab(slab, args)
        logger = self.get_logger()
        logger.debug(f'Deleted {slab.n_atoms - modified_slab.n_atoms} '
                     f'atoms at c < {args.c}.')
        return modified_slab

    def select_surviving_atoms(self, slab, args):
        """Return atoms of slab with c position of at least args.c."""
        return [atom for atom in slab if atom.pos[2] >= args.c]


if __name__ == '__main__':
    DeleteBelowCLI.run_as_script()
