"""ViPErLEED utility: Merge POSCAR files."""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Florian Kraushofer (@fkraushofer)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-08-03'
__license__ = 'GPLv3+'

from copy import deepcopy

import numpy as np
from scipy.spatial.distance import cdist

from viperleed.calc.files import poscar

from viperleed.cli_base import minimum_length, positive_float
from viperleed.utilities.poscar.base import _PoscarStreamCLI


def merge_slabs(slabs, check_collisions=True, eps=1e-5):
    """Return a single slab containing all the atoms from slabs.

    Parameters
    ----------
    slabs : Sequence of Slab
        The slab objects to be merged.
    check_collisions : bool, optional
        Whether collisions between atoms in the slabs to be merged
        should be reported. Default is True.
    eps : float, optional
        Cartesian tolerance (in angstrom) on the distances between
        pairs of atoms to be considered colliding. Default is 1e-5.

    Returns
    -------
    merged_slab : Slab
        A single slab with all the atoms from `slabs`.

    Raises
    ------
    RuntimeError
        If check_collisions is True, and there are atoms closer
        than `eps` when merging `slabs`.

    Notes
    -----
    The atoms of the slabs are currently NOT CORRECTLY UPDATED:
        Their slab, layer, etc. remain those of `slabs`, and
        not the corresponding properties of `merged_slab`.
    """
    merged_slab = deepcopy(slabs[0])
    # Allow duplicate atom.num till we fix
    # them later with update_atom_numbers()
    merged_slab.atlist.strict = False
    for slab in slabs[1:]:
        if check_collisions:
            # check for collisions
            dists = cdist([at.cartpos for at in merged_slab],
                          [at.cartpos for at in slab])
            if np.any(dists < eps):
                raise RuntimeError('Distance between atoms in different '
                                   f'slabs is smaller than eps={eps}.')
        merged_slab.atlist.extend(slab)                                         # TODO: Potential issue: The slab (and layer, etc.) of each added atom is NOT updated. Probably better to use duplicates if this function is to be used anywhere else.

    merged_slab.update_atom_numbers()
    merged_slab.update_element_count()
    merged_slab.update_cartesian_from_fractional()
    merged_slab.atlist.strict = True
    return merged_slab


_EPS_CELL_DEFAULT = 1e-3
_EPS_COLLISION_DEFAULT = 0.1


class MergePoscarsCLI(_PoscarStreamCLI, cli_name='merge'):
    """Merge two or more POSCAR files into one."""

    long_name = 'merge POSCAR files'

    def add_parser_arguments(self, parser):
        """Add positional list of files and optional arguments."""
        # This version of the overridden add_parser_arguments is
        # a bit tedious: we do not want to have a single stdin as
        # an input stream, but still want to support:
        # (1) the base arguments of ViPErLEEDCLI: we do this via
        #     super, skipping _PoscarStreamCLI in the inheritance
        super(_PoscarStreamCLI, self).add_parser_arguments(parser)
        # (2) the verbose option for the logger.
        self.add_verbose_option(parser)
        parser.add_argument('files', nargs='+', action=minimum_length(2))
        parser.add_argument(
            '--eps-cell',
            help=('maximum element-wise allowed difference in lattice '
                  'parameters (in angstrom) between the unit cells of '
                  'the POSCAR files to be merged. Default is '
                  f'{_EPS_CELL_DEFAULT:g}'),
            type=positive_float,
            default=_EPS_CELL_DEFAULT,
            )
        parser.add_argument(
            '--eps-collision',
            help=('minimum distance (in angstrom) between atoms in '
                  'different slabs before a collision is detected. '
                  f'Default is {_EPS_COLLISION_DEFAULT}'),
            type=positive_float,
            default=_EPS_COLLISION_DEFAULT
            )

    # DISABLE: the positional argument is renamed to slabs on purpose
    # as this makes it clearer that this method expects multiple slabs,
    # in contrast with the one of the base class. This also will raise
    # if anyone should pass it along as a named argument.
    # pylint: disable-next=arguments-renamed
    def process_slab(self, slabs, args):
        """Merge slabs into a new one."""
        return merge_slabs(slabs,
                           check_collisions=True,
                           eps=args.eps_collision)

    def read_poscar(self, args):
        """Read POSCAR files, checking that they are consistent."""
        slabs = [poscar.read(file) for file in args.files]

        # Check that all pairs of unit cells match within the tolerance
        for i, slab in enumerate(slabs[:-1]):
            this_ucell = slab.ucell
            other_ucells = [other_slab.ucell for other_slab in slabs[i:]]
            if not np.allclose(other_ucells, this_ucell, atol=args.eps_cell):
                self.parser.error('Files do not have the same unit cell. '
                                  'Use --eps-cell to increase tolerance.')
        return slabs


if __name__ == '__main__':
    MergePoscarsCLI.run_as_script()
