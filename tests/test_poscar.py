"""Test module poscar.

Created on 2023-07-12

@author: Alexander M. Imre
"""

from pathlib import Path
import os
import sys
from unittest.mock import mock_open, patch

import pytest

vpr_path = str(Path(__file__).parent.parent.parent)
if os.path.abspath(vpr_path) not in sys.path:
    sys.path.append(os.path.abspath(vpr_path))

from viperleed.tleedmlib.files.parameters import (readPARAMETERS,
                                                  ParameterInterpreter,
                                                  Assignment, NumericBounds)
from viperleed.tleedmlib.files.poscar import (readPOSCAR,
                                              writePOSCAR,
                                              ensure_away_from_c_edges,
                                              POSCARReader,
                                              POSCARSyntaxError,
                                              InvalidUnitCellError)
from viperleed.tleedmlib.classes.rparams import Rparams
from viperleed.tests.helpers import (ag100_parameters_example,
                                     example_poscars,
                                     slab_and_expectations)


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

    def test_writePOSCAR_does_not_modify_slab(self, example_poscars, tmp_path):
        slab = example_poscars
        n_atoms_original = len(slab.atlist)
        writePOSCAR(slab, tmp_path / 'POSCAR')
        assert slab == example_poscars
        assert len(slab.atlist) == n_atoms_original

    def test_writePOSCAR_writes_correct_number_of_atoms(self, slab_and_expectations, tmp_path):
        slab, expected_n_atoms, *_ = slab_and_expectations
        writePOSCAR(slab, tmp_path / 'POSCAR')
        # read the written POSCAR and check that the number of atoms is correct
        written_slab = readPOSCAR(tmp_path / 'POSCAR')
        assert len(written_slab.atlist) == expected_n_atoms
