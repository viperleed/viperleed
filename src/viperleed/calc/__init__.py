"""Package calc of viperleed.

This package contains all the main functionality of the TensErLEED
manager part of ViPErLEED, i.e., the part of ViPErLEED dedicated to
the calculation of theoretical I(V) curves as well as the optimization
of structural models such that their calculated I(V) curves best
fit experimental ones.

This package can be invoked as
>>> /path/to/python3 -m viperleed.calc <options>

Modules
-------
bookkeeper (script)
    Functionality for keeping track of
    calculation history, as in a logbook.
from_ase:
    API for executing viperleed.calc from an ase.Atoms object.
psgen
    Handle input/output to the FORTRAN program(s) that generate
    atomic phase-shifts and information concerning the energy
    dependence of the inner potential of the crystal.
run
    Main functionality for running calc from a set of input files.
symmetry
    Function for detecting symmetry of 2D- and 3D-periodic slabs.

Packages
--------
classes
    Definition of LEED-related classes used throughout viperleed.calc.
files
    Input/output handling for both TensErLEED sections and calc.
lib
    Functions and classes used in various parts of calc.
sections
    Functionality to run the various, logically different parts
    of TensErLEED.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2020-08-03'
__license__ = 'GPLv3+'

import logging
import multiprocessing as mp


LOGGER = logging.getLogger(__name__)

if not mp.parent_process():  # Only in the main process
    try:
        mp.set_start_method("spawn", force=False)                               # TODO: we could add a command-line option that allows the use of force=True
    except RuntimeError:
        # don't do anything for now, will be checked later
        pass
