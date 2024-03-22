"""ViPErLEED utility to convert experimental input for SATLEED to EXPBEAMS.csv.

Reads a file containing experimental I(V) curves as used by SATLEED
and converts it into the standard CSV format used by ViPErLEED.
"""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-03-20'
__license__ = 'GPLv3+'

import argparse
import logging
from pathlib import Path

import fortranformat as ff
import numpy as np

from viperleed.calc.classes.beam import Beam
from viperleed.calc.files.beams import writeOUTBEAMS, averageBeams
from viperleed.guilib.base import BeamIndex

logger = logging.getLogger(__name__)


def add_cli_parser_arguments(parser):
    parser.add_argument(
        "input",
        help="Input file containing the I(V) curves as used by SATLEED",
        type=str,
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Output file to write the beams to",
        type=str,
        default="EXPBEAMS.csv",
    )
    parser.add_argument(
        "--average",
        help="Average beams according to the averaging scheme in the input file",
        action="store_true",
        default=False,
    )


def read_file(file, average=False):
    """Read the file and return the contents as a list of Beam objects.

    Parameters
    ----------
    file : Pathlike
        Path to the input file.
    average : bool
        If True, average beams according to the averaging scheme in the input
        file.

    Raises
    ------
    ValueError
        If the file is not formatted correctly.

    Returns
    -------
    list of Beam
        List of Beam objects containing the data from the input file.
    """

    with open(file, "r") as f:
        lines = iter(f.readlines())

    # first line is name of the system
    name = next(lines).strip()
    # number of beams
    n_beams = int(next(lines))

    # averaging scheme
    averaging_scheme = next(lines)
    averaging_scheme = [int(x) for x in averaging_scheme.split()]

    if len(averaging_scheme) != n_beams:
        _err = (
            f"Averaging scheme (line 3) has length {len(averaging_scheme)},"
            f" expected {n_beams}."
        )
        raise ValueError(_err)

    # formatting
    data_format = next(lines).strip()
    reader = ff.FortranRecordReader(data_format)

    # energy increment
    energy_increment = float(next(lines).strip())

    beam_data = {}

    for _ in range(n_beams):
        # read beam hk
        beam_id = next(lines).strip()
        # strip brackets left and right an convert to BeamIndex
        beam_id = beam_id[1:-1] if beam_id.startswith("(") else beam_id
        beam_id = BeamIndex(beam_id)

        # number of energies and scaling factor
        n_energies, scaling_factor = next(lines).split()
        n_energies = int(n_energies)
        scaling_factor = float(scaling_factor)

        if abs(scaling_factor - 1.0) > 1e-3:
            logger.warning(
                "Scaling factor for beam {beam_index_str} is " "{scaling_factor}."
            )

        block_lines = [next(lines) for _ in range(n_energies)]
        block_values = [reader.read(line) for line in block_lines]
        energies = np.array([value[0] for value in block_values])
        intensities = np.array([value[1] for value in block_values])

        # check that energies are equally spaced
        # this also checks that the energies are sorted
        if not np.allclose(np.diff(energies), energy_increment):
            raise ValueError(
                f"Energies for beam {beam_id} are not spaced "
                "according to the energy increment specified on "
                " line 5. Interpolation is not supported."
            )

        # store into the dict and process the next beam
        beam_data[beam_id] = (energies, intensities)

    # check that we have read all lines; i.e there are no more lines with data
    try:
        while True:
            surplus_line = next(lines)
            if surplus_line.strip():
                raise ValueError("Found surplus data after reading all beams.")
    except StopIteration:
        pass

    processed_beams = []
    for beam_id, (energies, intensities) in beam_data.items():
        hk_beam = Beam(beam_id)
        hk_beam.intens = {en: intens for en, intens in zip(energies, intensities)}
        processed_beams.append(hk_beam)

    # average if requested
    if average and np.unique(averaging_scheme).size != n_beams:
        averaged_beams = []
        for group_id in np.unique(averaging_scheme):
            group_beams = [
                [processed_beams[i]]
                for i in np.argwhere(averaging_scheme == group_id).flatten()
            ]
            try:
                averaged_beams.append(averageBeams(group_beams, weights=None))
            except ValueError as e:
                logger.warning(f"Failed to average group {group_id}: {e}")
        processed_beams = averaged_beams
    elif np.unique(averaging_scheme).size != n_beams:
        raise ValueError(
            "Averaging scheme contains duplicate values. "
            "Use the --average option to average the beams if "
            "this is intended."
        )

    return processed_beams


def main(args=None):
    if args is None:
        parser = argparse.ArgumentParser()
        add_cli_parser_arguments(parser)
        args = parser.parse_args()

    # Read the file
    iv_beams_data = read_file(Path(args.input), average=args.average)

    out_file = args.output if args.output else "EXPBEAMS.csv"

    # check if out_file exists
    if Path(out_file).exists():
        raise FileExistsError(f"Output file {out_file} already exists.")

    # write the output file
    writeOUTBEAMS(filename=out_file, beams=iv_beams_data)


if __name__ == "__main__":
    main()
