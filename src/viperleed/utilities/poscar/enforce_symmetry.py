"""ViPErLEED utility: Enforce Symmetry."""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-08-03'
__license__ = 'GPLv3+'

from argparse import ArgumentTypeError

from viperleed.calc import symmetry
from viperleed.gui.classes.planegroup import _KNOWN_GROUPS
from viperleed.utilities.poscar.base import _PoscarSymmetryCLI


def validate_planegroup(plane_group):
    """Raise ArgumentTypeError if plane_group is not acceptable."""
    if plane_group not in _KNOWN_GROUPS:
        raise ArgumentTypeError(f'Invalid planegroup: {plane_group}')
    return plane_group


class SymmetrizeSlabCLI(_PoscarSymmetryCLI, cli_name='enforce_symmetry'):
    """Symmetrize a slab according to an optionally specified plane group."""

    long_name = 'symmetrize slab'

    def add_parser_arguments(self, parser):
        """Add an optional --plane-group argument."""
        super().add_parser_arguments(parser)
        parser.add_argument(
            '-g', '--plane-group',
            help=('Plane group to enforce. Default: detected automatically '
                  'from the slab. Use this option to override the automatic '
                  'detection and manually lower the symmetry.'),
            type=validate_planegroup,
            )

    def process_slab(self, slab, args):
        """Find slab symmetry and symmetrize it according to `args`."""
        rpars = self.prepare_rpars(slab, args)
        found_symmetry = symmetry.findSymmetry(slab, rpars)
        plane_group = args.plane_group or found_symmetry
        symmetry.enforceSymmetry(slab, rpars, plane_group)
        return slab


if __name__ == '__main__':
    SymmetrizeSlabCLI.run_as_script()
