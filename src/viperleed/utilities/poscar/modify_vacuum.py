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
from dataclasses import dataclass
import logging

import numpy as np

from viperleed.calc.classes.slab.errors import NotEnoughVacuumError
from viperleed.calc.classes.slab.errors import VacuumError
from viperleed.calc.classes.slab.errors import WrongVacuumPositionError
from viperleed.utilities.poscar.base import _PoscarStreamCLI


logger = logging.getLogger(__name__)


@dataclass
class VacuumGapInfo:
    """Information about a slab's desired vacuum gap.

    Attributes
    ----------
    size : float
        Size or size change of the desired vacuum gap, in angstrom.
    absolute : bool, optional
        Whether `size` should be considered as the absolute value of
        the new vacuum gap rather than its change. Default is False.
    accept_small_gap : bool, optional
        Whether the script should tolerate a resulting vacuum
        gap smaller than the minimum value for a viperleed.calc
        run, as long as it is non-negative. Default is False.
    """

    size: float
    absolute: bool = False
    accept_small_gap: bool = False


def modify_vacuum(slab, vacuum_gap_info):
    """Modify the vacuum gap size of a slab.

    Parameters
    ----------
    slab : Slab
        Slab to prepare for relaxation,
    vacuum_gap_info : VacuumGapInfo
        Information about the desired vacuum gap.

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
    vacuum_gap_size = vacuum_gap_info.size

    vacuum_gap_size += current_gap_size if not vacuum_gap_info.absolute else 0

    if vacuum_gap_size < 0:
        raise NotEnoughVacuumError(
            'The resulting vacuum gap size would be negative.',
            None,
            )

    logger.debug(f'Current vacuum gap size:\t{current_gap_size:9.3f}')
    logger.debug(f'New vacuum gap size:\t\t{vacuum_gap_size:9.3f}')

    new_c_vector = (
        processed_slab.c_vector
        / np.linalg.norm(processed_slab.c_vector)
        * (vacuum_gap_size + slab_thickness)
        )
    processed_slab.c_vector[:] = new_c_vector
    processed_slab.collapse_cartesian_coordinates(update_origin=True)

    try:
        processed_slab.check_vacuum_gap()
    except NotEnoughVacuumError as err:
        if not vacuum_gap_info.accept_small_gap:
            raise NotEnoughVacuumError(
                'The resulting vacuum gap would be too small.',
                None,
                ) from err
    except WrongVacuumPositionError as err:
        if not vacuum_gap_info.accept_small_gap:
            raise WrongVacuumPositionError(
                'Cannot modify the vacuum gap as requested. Check '
                'that there already is a vacuum gap in the POSCAR.',
                None,
                ) from err
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
        parser.add_argument(
            '-f', '--force',
            help=('If set, the script will not check if the resulting '
                  'vacuum gap is valid as long as it is non-negative.'),
            action='store_true',
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
        vacuum_gap_info = VacuumGapInfo(
            size=args.vacuum,
            absolute=args.absolute,
            accept_small_gap=args.force,
            )
        if args.absolute:
            logger.debug('Using absolute vacuum gap size.')
        try:
            return modify_vacuum(slab, vacuum_gap_info)
        except VacuumError as exc:
            self.parser.error(str(exc))
        return slab  # This is unreachable as parser.error sys-exits


if __name__ == '__main__':
    ModifyVacuumCLI.run_as_script()
