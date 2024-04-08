"""ViPErLEED utility: Modify vacuum gap.

This utility takes a slab in POSCAR format and modifies the vacuum gap.
"""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-08-03'
__license__ = 'GPLv3+'

from copy import deepcopy
import logging

import numpy as np

from viperleed.calc.classes.slab.errors import (NotEnoughVacuumError,
                                                WrongVacuumPositionError)
from viperleed.utilities.poscar.base import _PoscarStreamCLI


logger = logging.getLogger(__name__)


def modify_vacuum(slab, vacuum_gap_size, absolute=False):
    """Modify the vacuum gap size of a slab.

    Parameters
    ----------
    slab : Slab
        Slab to prepare for relaxation,
    vacuum_gap_size : float
        Size of the vacuum gap in Angstroms. If absolute is True, this
        is the absolute size of the vacuum gap, else it is the amount
        of vacuum to add or remove from the slab.
    absolute : bool, optional
        If True, vacuum_gap_size is the absolute size of the vacuum gap.

    Returns
    -------
    Slab
        The slab with the modified vacuum gap.

    Raises
    ------
    RuntimeError
        If the resulting vacuum gap size is negative.
    """
    processed_slab = deepcopy(slab)
    processed_slab.check_a_b_in_plane()
    processed_slab.project_c_to_z()
    processed_slab.update_cartesian_from_fractional()

    slab_thickness = processed_slab.thickness
    current_gap_size = processed_slab.vacuum_gap

    vacuum_gap_size += current_gap_size if not absolute else 0

    # check there is enough space on top of the slab to reduce the vacuum gap
    if vacuum_gap_size < current_gap_size:
        top_atom_positions = max(at.cartpos[2] for at in slab)
        gap_on_top = slab.c_vector[2] - top_atom_positions
        if gap_on_top < (current_gap_size - vacuum_gap_size):
            raise RuntimeError('The resulting vacuum gap size would be larger '
                               'than the distance between the topmost atom and '
                               'the top of the cell.')

    logger.debug(f'Current vacuum gap size:\t{current_gap_size:9.3f}')
    logger.debug(f'New vacuum gap size:\t\t{vacuum_gap_size:9.3f}')

    if vacuum_gap_size < 0:
        raise RuntimeError('The resulting vacuum gap size would be negative.')

    new_c_vector = (
        processed_slab.c_vector
        / np.linalg.norm(processed_slab.c_vector)
        * (vacuum_gap_size + slab_thickness)
    )
    print(f'New c vector: {new_c_vector}')
    print(f'Old c vector: {processed_slab.c_vector}')

    processed_slab.c_vector[:] = new_c_vector
    processed_slab.collapse_cartesian_coordinates()
    processed_slab.update_fractional_from_cartesian()
    try:
        processed_slab.check_vacuum_gap()
    except NotEnoughVacuumError:
        raise RuntimeError('The resulting vacuum gap would be too small.')
    except WrongVacuumPositionError:
        raise RuntimeError('Cannot modify the vaccum gap as requested. Check '
                           'that there already is a vacuum gap in the POSCAR.')
    return processed_slab


class ModifyVacuumCLI(_PoscarStreamCLI, cli_name='modify_vacuum'):
    """Change the vacuum gap of a POSCAR."""

    long_name = 'modify vacuum gap'

    def add_parser_arguments(self, parser):
        """Add mandatory vacuum and optional --absolute arguments."""
        super().add_parser_arguments(parser)
        parser.add_argument(
            'vacuum',
            help=('Add (or remove, if negative) vacuum on top of '
                  'the slab. Value in angstrom. Default is zero.'),
            type=float,
            default=0.0,
            )
        parser.add_argument(
            '-a', '--absolute',
            help=('If set, the value given to vacuum is the absolute size '
                  'of the vacuum gap, else it is the amount of vacuum to add '
                  'or remove from the slab.'),
            action='store_true'
            )

    def parse_cli_args(self, args):
        """Check consistency of vacuum and absolute."""
        parsed_args = super().parse_cli_args(args)
        if parsed_args.absolute and parsed_args.vacuum < 0:
            self.parser.error('Absolute vacuum gap must be non-negative, '
                              f'but was {parsed_args.vacuum}.')
        return parsed_args

    def process_slab(self, slab, args):
        """Return a new slab with modified vacuum."""
        if args.absolute:
            logger.debug('Using absolute vacuum gap size.')
        try:
            return modify_vacuum(slab, args.vacuum, absolute=args.absolute)
        except RuntimeError as exc:
            self.parser.error(str(exc))
        return slab  # This is unreachable as parser.error sys-exits


if __name__ == '__main__':
    ModifyVacuumCLI.run_as_script()
