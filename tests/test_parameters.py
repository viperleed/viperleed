"""Test module parameters of viperleed.tests.

Created on 2023-06-09

@author: Alexander M. Imre

"""

import pytest
import shutil, tempfile
import sys
import os
from pathlib import Path
from zipfile import ZipFile
from copy import deepcopy
import numpy as np

vpr_path = str(Path(__file__).parent.parent.parent)
if os.path.abspath(vpr_path) not in sys.path:
    sys.path.append(os.path.abspath(vpr_path))

import viperleed.tleedmlib.files.parameters as parameters
from viperleed.tleedmlib.files.parameters import (readPARAMETERS,
                                                  interpretPARAMETERS,
                                                  modifyPARAMETERS)
from viperleed.tleedmlib.files.poscar import readPOSCAR
from viperleed.tleedmlib.classes.rparams import Rparams
from viperleed.tleedmlib.files.parameter_errors import ParameterError

@pytest.fixture()
def ag100_parameters_example():
    # read Ag(100) POSCAR and PARAMETERS files
    slab = readPOSCAR('tests/fixtures/Ag(100)/initialization/POSCAR')
    rpars = readPARAMETERS('tests/fixtures/Ag(100)/initialization/PARAMETERS')
    # interpret PARAMETERS file
    interpretPARAMETERS(rpars, slab)
    return rpars, slab

def test_read_parameters_for_ag100():
    # just check that readPARAMETERS does not crash; not interpreted yet
    filename = 'tests/fixtures/parameters/PARAMETERS_1'
    rpars = readPARAMETERS(filename)
    assert rpars


def test_interpret_parameters_for_ag100(ag100_parameters_example):
    rpars, slab = ag100_parameters_example
    # check that the parameters are interpreted correctly based on a few examples
    assert rpars.V0_IMAG == pytest.approx(5.0)
    assert rpars.THEO_ENERGIES == pytest.approx([50, 350, 3])

# unit tests for parameter interpretation

# _interpret_intpol_deg()
def test_interpret_intpol_deg():
    rpars = Rparams()
    param = 'INTPOL_DEG'

    for val in rpars.get_limits(param):
        parameters._interpret_intpol_deg(rpars, param, val)
        assert rpars.INTPOL_DEG == int(val)
    incompatible_values = ['1', 'text']
    for val in incompatible_values:
        with pytest.raises(ParameterError):
            parameters._interpret_intpol_deg(rpars, param, val)