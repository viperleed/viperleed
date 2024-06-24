"""ViPErLEED utility: Rescale cell.

This utility takes a slab in POSCAR format and rescales the unit cell.
"""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-08-03'
__license__ = 'GPLv3+'

from viperleed.cli_base import length_choices
from viperleed.utilities.poscar.base import _PoscarStreamCLI


class RescaleCellCLI(_PoscarStreamCLI, cli_name='rescale_cell'):
    """Apply scaling factors to the unit cell of a POSCAR."""

    long_name = 'rescale unit cell'

    def add_parser_arguments(self, parser):
        """Add mandatory scaling factor(s) argument."""
        super().add_parser_arguments(parser)
        parser.add_argument(
            'scaling',
            help=('One or three scaling factors for the unit cell. If three '
                  'values are given, the scaling factors are applied to the '
                  'a, b, and c vectors, respectively. If only one value is '
                  'given, an isotropic scaling is applied.'),
            type=float,
            nargs='+',
            action=length_choices(1, 3),
            )

    def process_slab(self, slab, args):
        """Return a slab with a unit cell scaled according to args."""
        slab.apply_scaling(*args.scaling)
        return slab


if __name__ == '__main__':
    RescaleCellCLI.run_as_script()
