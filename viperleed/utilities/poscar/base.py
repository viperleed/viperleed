"""Module base of viperleed.utilities.poscar.

Defines base CLI classes and functions useful to limit code
repetition in poscar utilities.
"""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-09-27'
__license__ = 'GPLv3+'

from abc import ABC, abstractmethod
from copy import deepcopy
import logging
import sys

from viperleed.calc.classes.rparams import Rparams, SymmetryEps
from viperleed.calc.classes.rparams.defaults import DEFAULTS as PARAM_DEFAULTS
from viperleed.calc.files import poscar
from viperleed.cli_base import ViPErLEEDCLI, positive_float


# ##########################   IMPORTANT  #############################
# All new base classes should be named with a leading underscore to   #
# prevent collection of this module as a valid poscar utility!        #
# #####################################################################


class _PoscarStreamCLI(ViPErLEEDCLI, ABC, cli_name=None):
    """A POSCAR utilities that reads/writes from/to the terminal."""

    long_name = ''

    def __call__(self, args=None):
        """Call this utility."""
        logger = self.get_logger()
        parsed_args = self.parse_cli_args(args)
        if parsed_args.verbose:
            logger.setLevel(logging.DEBUG)
        logger.info(f'ViPErLEED POSCAR utility: {self.long_name}\n')
        slab = self.read_poscar(parsed_args)
        processed_slab = self.process_slab(slab, parsed_args)
        self.write_to_stdout(processed_slab, parsed_args)
        return 0

    def add_parser_arguments(self, parser):
        """Add generic arguments to this CLI."""
        super().add_parser_arguments(parser)
        self.add_verbose_option(parser)

    @abstractmethod
    def process_slab(self, slab, args):
        """Return a processed slab from slab.

        This is an abstract method that must be
        extended in concrete subclasses.

        Parameters
        ----------
        slab : Slab
            The slab that was read from stdin and should be processed
            by this utility.
        args : argparse.Namespace
            The parsed command-line arguments.

        Returns
        -------
        processed_slab : Slab
            A new slab processed by this CLI utility. Will be
            written to stout exactly as it is returned.
        """

    # DISABLE: args is unused in this implementation, but it is
    # needed for the interface, as subclasses may need it to
    # decide what/how to read.
    # pylint: disable-next=unused-argument
    def read_poscar(self, args):
        """Return a slab from args.

        The return value of this method is passed on, together
        with args, to self.process_slab. The default implementation
        reads a slab from stdin.

        Parameters
        ----------
        args : argparse.Namespace
            The processed CLI arguments.

        Returns
        -------
        slab : Slab
            A slab read from the terminal.
        """
        return poscar.read(sys.stdin)

    # DISABLE: args is unused in this implementation, but it is
    # needed for the interface, as subclasses may need it to
    # decide what/how to read.
    # pylint: disable-next=unused-argument
    def write_to_stdout(self, processed_slab, args):
        """Write output to the terminal.

        This is the last method that is called by this utility
        before terminating. The default implementation writes
        processed_slab to the terminal.

        Parameters
        ----------
        processed_slab : Slab
            The slab processed by this utility.
        args : argparse.Namespace
            The processed command-line arguments.
        """
        self.write_processed_slab(processed_slab)

    def write_processed_slab(self, slab):
        """Write slab to the terminal."""
        log_level = self.get_logger().getEffectiveLevel()
        poscar.write(slab,
                     filename=sys.stdout,
                     comments='none',
                     silent=log_level<=logging.DEBUG)


EPS_DEFAULT = PARAM_DEFAULTS['SYMMETRY_EPS']

class _PoscarSymmetryCLI(_PoscarStreamCLI, ABC, cli_name=None):
    """Base CLI class for utilities dealing with slab symmetry."""

    def add_parser_arguments(self, parser):
        """Add SYMMETRY_EPS optional arguments."""
        parser.add_argument(
            '-e', '--symmetry-eps',
            help=('Epsilon for symmetry detection in angstrom. '
                  f'Default: {EPS_DEFAULT} A'),
            type=positive_float,
            default=EPS_DEFAULT,
            )
        parser.add_argument(
            '--symmetry-eps-z',
            help=('Epsilon for symmetry detection in z in angstrom. If '
                  'not provided, the value of --symmetry-eps is used.'),
            type=positive_float,
            )

    def prepare_rpars(self, slab, args):
        """Return an Rparams with symmetry information, and update slab."""
        rpars = Rparams()
        rpars.SYMMETRY_EPS = SymmetryEps(args.symmetry_eps,
                                         args.symmetry_eps_z)
        rpars.SYMMETRY_FIND_ORI = True
        slab.full_update(rpars)
        return rpars


class _RemoveAtomsCLI(_PoscarStreamCLI, ABC, cli_name=None):
    """A command-line interface that removes atoms from a slab."""

    def add_parser_arguments(self, parser):
        """Add the condition for atom removal."""
        super().add_parser_arguments(parser)
        self.add_remove_condition_arguments(parser)

    @abstractmethod
    def add_remove_condition_arguments(self, parser):
        """Add arguments that define the condition for atom removal."""

    def process_slab(self, slab, args):
        """Return a new slab with atoms from slab removed according to args."""
        modified_slab = deepcopy(slab)
        atoms = self.select_surviving_atoms(modified_slab, args)
        modified_slab.atlist.clear()
        modified_slab.extend(atoms)
        modified_slab.update_element_count()
        return modified_slab

    @abstractmethod
    def select_surviving_atoms(self, slab, args):
        """Return the atoms of slab that should survive the removal."""
