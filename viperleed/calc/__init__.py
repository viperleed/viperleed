"""Package calc of viperleed.

This package contains all the main functionality of the TensErLEED
manager part of ViPErLEED (tleedm), i.e., the part of ViPErLEED
dedicated to the calculation of theoretical I(V) curves as well
as the optimization of structural models such that their calculated
I(V) curves best fit experimental ones.

Modules
-------
base:
    Non-LEED-related classes and functions used throughout tleedmlib.
beamgen:
    Handle input/output to the FORTRAN program that generates the list
    of beams used internally by TensErLEED
checksums:
    Functions used to prevent user tinkering with FORTRAN source files.
    The actual checksums are contained in _checksums.dat and should
    only be edited by expert contributors.
leedbase:
    LEED- and TensErLEED-related base functions used in tleedmlib
periodic_table:
    Collection of basi atom data
psgen:
    Handle input/output to the FORTRAN program(s) that generate atomic
    phase-shifts and information concerning the energy dependence of
    the inner potential of the crystal
symmetry:
    Function for detecting symmetry of 2D- and 3D-periodic slabs

Packages
--------
classes:
    Definition of LEED-related classes used throughout viperleed.calc.
files:
    Input/output handling for TensErLEED sections
sections:
    Functionality to run the various, logically different parts of
    TensErLEED
wrapped:
    Python extensions written in FORTRAN
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__created__ = '2020-08-03'

import logging

LOGGER = logging.getLogger(__name__)
LOG_PREFIX = 'viperleed-calc'
