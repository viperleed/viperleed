"""Module base of viperleed.utilities.poscar.

Defines base CLI classes and functions useful to limit code
repetition in poscar utilities.
"""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-03-27'
__license__ = 'GPLv3+'

from abc import ABC
from abc import abstractmethod
from copy import deepcopy
import logging
import sys

from viperleed.calc.classes.rparams.defaults import DEFAULTS as PARAM_DEFAULTS
from viperleed.calc.classes.rparams.rparams import Rparams
from viperleed.calc.classes.rparams.special.symmetry_eps import SymmetryEps
from viperleed.calc.files import poscar
from viperleed.calc.lib.log_utils import debug_or_lower
from viperleed.cli_base import ViPErLEEDCLI
from viperleed.cli_base import StreamArgument
from viperleed.cli_base import positive_float


# ##########################   IMPORTANT  #############################
# All new base classes should be named with a leading underscore to   #
# prevent collection of this module as a valid poscar utility!        #
# #####################################################################


class _PoscarStreamCLI(ViPErLEEDCLI, ABC, cli_name=None):
    """A POSCAR utility that reads/writes from/to files or the terminal."""

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
        self.write_output(processed_slab, parsed_args)
        return 0

    def add_infile_argument(self, parser):
        """Add an optional --infile/-i argument to parser.

        This method is automatically called by the base implementation
        of add_parser_arguments. It can be used by subclasses that
        intend to override (and not extend) add_parser_arguments but
        still desire to provide an --infile/-i argument.

        **IMPORTANT**: Subclasses that override add_parser_arguments
        but do not override read_poscar **must** call this method
        in add_parser_arguments.

        Parameters
        ----------
        parser : argparse.ArgumentParser
            The parser to which the optional argument should be
            added. The added argument defaults to the standard-in
            stream, i.e., the terminal.

        Returns
        -------
        None.
        """
        help_ = ('Name of the POSCAR input file. Default: read text '
                 'from the standard-input stream (i.e., the terminal)')
        stream = StreamArgument('r')
        parser.add_argument('--infile', '-i',
                            type=stream,
                            help=help_,
                            default=stream(sys.stdin))

    def add_outfile_argument(self, parser):
        """Add an optional --outfile/-o argument to parser.

        This method is automatically called by the base implementation
        of add_parser_arguments. It can be used by subclasses that
        intend to override (and not extend) add_parser_arguments but
        still desire to provide an --outfile/-o argument.

        ''IMPORTANT**: Subclasses that override add_parser_arguments
        but do not override write_output **must** call this method
        in add_parser_arguments.

        Parameters
        ----------
        parser : argparse.ArgumentParser
            The parser to which the optional argument should be
            added. The added argument defaults to the standard-out
            stream, i.e., the terminal.

        Returns
        -------
        None.
        """
        help_ = ('Name of the POSCAR output file. Default: write text '
                 'to the standard-output stream (i.e., the terminal)')
        stream = StreamArgument('w')
        parser.add_argument('--outfile', '-o',
                            type=stream,
                            help=help_,
                            default=stream(sys.stdout))

    def add_parser_arguments(self, parser):
        """Add generic arguments to this CLI.

        The base implementation adds:
        - all the arguments of ancestor classes.
        - two optional arguments for the input and output files
          (--infile/-i and --outfile/-o, respectively) which
          default to standard-in and standard-out.
        - an optional argument to turn logging verbose.

        **IMPORTANT**: Subclasses may extend this method via super()
        if they want to inherit the arguments of this base class. If
        not, they **must** provide outfile (positional) or --outfile
        (optional) CLI arguments (typically with sys.stdout default),
        unless they also override write_output. If, additionally, they
        do not override read_poscar, they must provide infile/--infile
        CLI arguments (typically with sys.stdin default). Subclasses
        may use the add_infile_argument/add_outfile_argument methods
        to add optional --infile/--outfile arguments.

        Parameters
        ----------
        parser : argparse.ArgumentParser
            The parser to which arguments are added.

        Returns
        -------
        None.
        """
        super().add_parser_arguments(parser)
        self.add_infile_argument(parser)
        self.add_outfile_argument(parser)
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

    def read_poscar(self, args):
        """Return a slab from args.

        The return value of this method is passed on, together with
        args, to self.process_slab. The default implementation reads
        a slab from args.infile, which defaults to the terminal.

        Parameters
        ----------
        args : argparse.Namespace
            The processed CLI arguments.

        Returns
        -------
        slab : Slab
            A slab read using information in args.

        Raises
        ------
        AttributeError
            If subclasses have overridden (and not extended) the
            add_parser_arguments method, have not overridden this
            method and forgot to call add_infile_argument in their
            modified add_parser_arguments.
        """
        try:
            _ = args.infile
        except AttributeError:
            raise AttributeError('--infile argument is missing. If you have '
                                 'overridden add_parser_arguments without '
                                 'also overriding read_poscar, make sure to '
                                 'call add_infile_argument in your overridden '
                                 'add_parser_arguments') from None
        try:
            return self._read_poscar_from_infile(args)
        except (ValueError, poscar.POSCARError, FileNotFoundError) as exc:
            self.parser.error(f'Failed to read POSCAR. Stopping. Info: {exc}')
        return None  # Unreachable as .error does SystemExit

    def write_output(self, processed_slab, args):
        """Write output to an output file or the terminal.

        This is the last method that is called by this utility
        before terminating. The default implementation writes
        processed_slab to args.outfile, which defaults to the
        terminal.

        Parameters
        ----------
        processed_slab : Slab
            The slab processed by this utility.
        args : argparse.Namespace
            The processed command-line arguments.

        Raises
        ------
        AttributeError
            If subclasses have overridden (and not extended) the
            add_parser_arguments method, have not overridden this
            method and forgot to call add_outfile_argument in their
            modified add_parser_arguments.
        """
        try:
            _ = args.outfile
        except AttributeError:
            raise AttributeError('--outfile argument is missing. If you have '
                                 'overridden add_parser_arguments without '
                                 'also overriding write_output, make sure to '
                                 'call add_outfile_argument in your overridden'
                                 ' add_parser_arguments') from None
        with args.outfile as outfile:
            poscar.write(processed_slab,
                         filename=outfile,
                         comments='none',
                         silent=debug_or_lower(self.get_logger()))

    def _read_poscar_from_infile(self, args):
        """Return a POSCAR read from args.infile."""
        if args.infile.is_interactive:
            sys.stderr.write('Please input the contents of a POSCAR file:\n')
        with args.infile as infile:
            return poscar.read(infile)



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
        super().add_parser_arguments(parser)

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
        modified_slab.atlist.extend(atoms)
        modified_slab.update_element_count()
        return modified_slab

    @abstractmethod
    def select_surviving_atoms(self, slab, args):
        """Return the atoms of slab that should survive the removal."""
