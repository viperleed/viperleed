"""ViPErLEED utility to convert TensErLEED-style beams to EXPBEAMS."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2019-12-16'
__license__ = 'GPLv3+'

from viperleed.calc.files.beams import readAUXEXPBEAMS, writeOUTBEAMS
from viperleed.cli_base import ViPErLEEDCLI
from viperleed.utilities.beams import EXPBEAMS_DEFAULT


_DEFAULT_INPUT = 'AUXEXPBEAMS'
_HEADER = '''
This utility reads an AUXEXPBEAMS file (i.e., beams formatted as TensErLEED
experimental input) and writes the contents in the csv formatting applied
in viperleed.calc for THEOBEAMS.csv and EXPBEAMS.csv files.
'''
_HELP = ('convert a TensErLEED-style experimental-'
         'beams input file to ViPErLEED format')


class FromTensErLEEDBeamsCLI(ViPErLEEDCLI,
                             cli_name='from-TensErLEED',
                             help_=_HELP):
    """Utility to convert TensErLEED beams file to EXPBEAMS format."""

    def __call__(self, args=None):
        """Execute utility."""
        # print some info
        print(_HEADER)

        # read the AUXEXPBEAMS file
        filename = ''
        while not filename:
            filename = input('Enter name of TensErLEED-beams file '
                             f'(default: {_DEFAULT_INPUT!r}): ')
            filename = filename or _DEFAULT_INPUT
            try:
                beams = readAUXEXPBEAMS(filename, interactive=True)
            except FileNotFoundError:
                print(f'File {filename} not found.')
                filename = ''

        if not beams:
            print('Error reading AUXEXPBEAMS file.')
            return 1
        # print some info
        print(f'Found file with {len(beams)} beams.\n')

        # relabel beams
        label_width = max(beam.getLabel()[1] for beam in beams)                 # TODO: repeated in many places
        for beam in beams:
            beam.label, _ = beam.getLabel(lwidth=label_width)

        # get output file name
        filename = input('Enter output file name '
                         f'(default: {EXPBEAMS_DEFAULT!r}): ')
        filename = filename or EXPBEAMS_DEFAULT

        # write new file
        try:
            writeOUTBEAMS(beams, filename)
        except OSError:
            print(f'Error writing new file {filename}')
            return 1
        print(f'Wrote output as {filename}')
        return 0


if __name__ == '__main__':
    FromTensErLEEDBeamsCLI.run_as_script()
