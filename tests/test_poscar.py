"""Test module file/poscar.

Created on 2023-07-12

@author: Alexander M. Imre

"""

import pytest
from pathlib import Path
import os, sys
from copy import deepcopy
import numpy as np

vpr_path = str(Path(__file__).parent.parent.parent)
if os.path.abspath(vpr_path) not in sys.path:
    sys.path.append(os.path.abspath(vpr_path))

import viperleed.tleedmlib.files.parameters as parameters
from viperleed.tleedmlib.files.parameters import (readPARAMETERS,
                                                  ParameterInterpreter,
                                                  Assignment, NumericBounds)
from viperleed.tleedmlib.files.poscar import readPOSCAR, writePOSCAR
from viperleed.tleedmlib.classes.rparams import Rparams
from viperleed.tests.helpers import (ag100_parameters_example,
                                     _FIXTURES_PATH,
                                     _POSCARs_PATH,
                                     _EXAMPLE_POSCARs,
                                     example_poscars,
                                     slab_and_expectations,
)


class Test_readPOSCAR:
    def test_readPOSCAR_slab_exists(self, example_poscars):
        slab = example_poscars
        assert slab

    def test_read_POSCAR_slab_has_atoms(self, example_poscars):
        slab = example_poscars
        assert len(slab.atlist) > 0

    def test_readPOSCAR_slab_n_atom_correct(self, slab_and_expectations):
        slab, expected_n_atoms, *_ = slab_and_expectations
        assert len(slab.atlist) == expected_n_atoms

class Test_writePOSCAR:
    # TODO: Test that the written POSCAR is the same as the read one.
    def test_writePOSCAR_creates_files(self, example_poscars, tmp_path):
        slab = example_poscars
        writePOSCAR(slab, tmp_path / 'POSCAR')
        assert (tmp_path / 'POSCAR').exists()