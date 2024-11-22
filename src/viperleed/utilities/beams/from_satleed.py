"""ViPErLEED utility to convert experimental input for SATLEED to EXPBEAMS.csv.

Reads a file containing experimental I(V) curves as used by SATLEED
and converts it into the standard CSV format used by ViPErLEED.
"""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-03-20'
__license__ = 'GPLv3+'

import logging
from pathlib import Path

import fortranformat as ff
import numpy as np

from viperleed.calc.classes.beam import Beam
from viperleed.calc.files.beams import averageBeams
from viperleed.calc.files.beams import writeOUTBEAMS
from viperleed.cli_base import ViPErLEEDCLI
from viperleed.guilib.base import BeamIndex
from viperleed.utilities.beams import EXPBEAMS_DEFAULT

logger = logging.getLogger(__name__)


_HELP = ('convert a SATLEED-style experimental-'
         'beams input file to ViPErLEED format')
_SCALING_FACTOR_DEVIATION_FROM_ONE = 1e-3  # Warn if > than this


class FromSATLEEDCLI(ViPErLEEDCLI, cli_name='from-SATLEED', help_=_HELP):
    """Utility to convert SATLEED beams file to EXPBEAMS format."""

    def add_parser_arguments(self, parser):
        """Add command-line arguments to parser."""
        super().add_parser_arguments(parser)
        parser.add_argument(
            'input',
            help='input file containing the I(V) curves as used by SATLEED',
            type=str,
            )
        parser.add_argument(
            '-o', '--output',
            help=('output file to write the beams to. '
                  f'Default is {EXPBEAMS_DEFAULT}'),
            type=str,
            default=EXPBEAMS_DEFAULT,
            )
        parser.add_argument(
            '--average',
            help=('Average beams according to the averaging '
                  'scheme in the input file. Default is False.'),
            action='store_true',
            default=False,
            )

    def __call__(self, args=None):
        """Parse CLI arguments, then convert file."""
        args = self.parse_cli_args(args)
        # Read input file
        try:
            iv_beams_data = _read_satleed_file(Path(args.input),
                                               average=args.average)
        except ValueError as exc:
            self.parser.error(str(exc))
            raise  # Unreachable as .error does sys.exit
        except StopIteration:
            self.parser.error(f'Not enough lines in {args.input}.')
            raise  # Unreachable as .error does sys.exit
        # Check if out_file exists...
        out_file = args.output
        if Path(out_file).exists():
            err_ = f'Output file {out_file} already exists.'
            self.parser.error(err_)
            raise FileExistsError(err_)
        # ...and write to it
        writeOUTBEAMS(filename=out_file, beams=iv_beams_data)


def _average_beams_if_requested(beams, averaging_scheme, average):
    """Return a list of averaged beams.

    Parameters
    ----------
    beams : list of Beam
        The beams to be averaged.
    averaging_scheme : Sequence of int
        Grouping indices for averaging `beams`. Same length as
        `beams`. Elements are group 'ids'. Beams with the same
        group id are averaged together.
    average: bool
        Whether any averaging should be carried out.

    Returns
    -------
    averaged_beams : list of Beam
        The beams averaged according to `averaging_scheme`.

    Raises
    ------
    ValueError
        If no averaging was requested via `average`, but the
        `averaging_scheme` would suggest averaging is needed.
    """
    averaging_groups = np.unique(averaging_scheme)
    if not average and averaging_groups.size != len(beams):
        raise ValueError(
            'Averaging scheme contains duplicate values. '
            'Use the --average option to average the beams if '
            'this is intended.'
            )
    if not average:
        return beams

    averaged_beams = []
    for group_id in averaging_groups:
        group_beams = [
            [beams[i]]
            for i in np.argwhere(averaging_scheme == group_id).flatten()        # TODO: probably nicer to do this with a boolean mask after converting beams into an array. Something along the lines of group_beams = beams[averaging_scheme == group_id]
            ]
        try:
            averaged_beams.append(averageBeams(group_beams, weights=None))
        except ValueError as exc:
            logger.warning(f'Failed to average group {group_id}: {exc}')
    return averaged_beams


def _read_averaging_scheme(line, n_beams):
    """Return a list of indices for beams to be averaged."""
    averaging_scheme = [int(x) for x in line.split()]
    if len(averaging_scheme) != n_beams:
        raise ValueError('Averaging scheme (line 3) has length '
                         f'{len(averaging_scheme)}, expected {n_beams}.')
    return averaging_scheme


def _read_beams(lines, n_beams, energy_increment, reader):
    """Return a list of `n_beams` Beam objects read from lines."""
    beam_data = {}
    for _ in range(n_beams):
        # read beam hk
        beam_id = next(lines).strip()
        # strip brackets left and right an convert to BeamIndex
        beam_id = beam_id[1:-1] if beam_id.startswith('(') else beam_id
        beam_id = BeamIndex(beam_id)

        # number of energies and scaling factor
        n_energies, scaling_factor = next(lines).split()
        n_energies = int(n_energies)
        scaling_factor = float(scaling_factor)

        if abs(scaling_factor - 1.0) > _SCALING_FACTOR_DEVIATION_FROM_ONE:
            logger.warning(
                f'Scaling factor for beam {beam_id} is {scaling_factor}.'
                )

        block_lines = (next(lines) for _ in range(n_energies))
        block_values = [reader.read(line) for line in block_lines]
        energies = np.array([value[0] for value in block_values])
        intensities = np.array([value[1] for value in block_values])

        # check that energies are equally spaced
        # this also checks that the energies are sorted
        if not np.allclose(np.diff(energies), energy_increment):
            raise ValueError(
                f'Energies for beam {beam_id} are not spaced '
                'according to the energy increment specified on '
                'line 5. Interpolation is not supported.'
                )

        # store into the dict and process the next beam
        beam_data[beam_id] = (energies, intensities)

    processed_beams = []
    for beam_id, (energies, intensities) in beam_data.items():
        hk_beam = Beam(beam_id)
        hk_beam.intens = dict(zip(energies, intensities))
        processed_beams.append(hk_beam)
    return processed_beams


def _read_satleed_file(file_name, average=False):
    """Read `file_name` and return the contents as a list of Beam objects.

    Parameters
    ----------
    file_name : Path
        Path to the input file.
    average : bool, optional
        If True, average beams according to the averaging
        scheme in `file_name`. Default is False.

    Returns
    -------
    list of Beam
        List of Beam objects containing the data from the input file.

    Raises
    ------
    ValueError
        If the file is not formatted correctly, or if
        averaging would be needed but `average` was False.
    StopIteration
        If there are not enough lines in `file_name`.
    """
    with file_name.open('r', encoding='utf-8') as file:
        lines = iter(file.readlines())

    # First line is name of the system
    _ = next(lines).strip()
    # number of beams
    n_beams = int(next(lines))

    # Averaging scheme
    averaging_scheme = _read_averaging_scheme(next(lines), n_beams)

    # Formatting
    data_format = next(lines).strip()
    reader = ff.FortranRecordReader(data_format)

    # Energy increment
    energy_increment = float(next(lines).strip())

    # Beam energies and intensities for each beam index
    processed_beams = _read_beams(lines, n_beams, energy_increment, reader)

    # Check that we have read all lines
    for surplus_line in lines:
        if surplus_line.strip():
            raise ValueError('Found surplus data after reading all beams.')

    # Average if requested
    return _average_beams_if_requested(processed_beams,
                                       averaging_scheme,
                                       average)


if __name__ == '__main__':
    FromSATLEEDCLI.run_as_script()
