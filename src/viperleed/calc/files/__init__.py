"""Module files of viperleed.calc.

Contains functionality for reading and writing files. This
includes both input/output files for the whole viperleed.calc
as well as input/output files for the TensErLEED FORTRAN code.

Modules
-------
beamgen
    Functions to generate the BEAMLIST file, a TensErLEED input.
beams
    Functions for reading and writing various beams files.
displacements
    Functions for reading and interpreting the DISPLACEMENTS file.
iodeltas
    Reading and writing files for the delta-amplitudes calculation.
ioerrorcalc
    Reading and writing files for the error calculation.
iofdopt
    Functions for writing output from full-dynamic optimization.
iorefcalc
    Reading and writing files for the reference calculation.
iorfactor
    Reading and writing files for the R-factor calculation.
iosearch
    Reading, processing and writing files for the structure search.
iosuperpos
    Writing files relevant to the superpos calculation.
ivplot
    Functions for producing plots of I(V) curves.
patterninfo
    Functions for producing input files usable in the pattern
    simulator plug-in of the ViPErLEED GUI.
poscar
    Functions for reading and writing POSCAR files.
searchpdf
    Functions for writing SearchProgress/Report.pdf files.
vibrocc
    Functions for reading and writing the VIBROCC file.

Packages
--------
parameters
    Reading, writing, editing, and interpreting a PARAMETERS file.
"""
# Currently undocumented:
# delta_intensities
#     Too early, too messy.
# new_search
#     Too early, too messy.

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2020-08-19'
__license__ = 'GPLv3+'
