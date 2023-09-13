"""Test module poscar.

Created on 2023-07-12

@author: Alexander M. Imre
"""

from pathlib import Path
import os
import sys

import pytest

from viperleed.calc.classes.rparams import Rparams
from viperleed.calc.files import poscar
from viperleed.calc.lib import symmetry


class TestReadPOSCAR:
    def test_read_slab_exists(self, example_poscars):
        slab = example_poscars
        assert slab

    def test_read_slab_has_atoms(self, example_poscars):
        slab = example_poscars
        assert len(slab.atlist) > 0

    def test_read_slab_n_atom_correct(self, slab_and_expectations):
        slab, expected_n_atoms, *_ = slab_and_expectations
        assert len(slab.atlist) == expected_n_atoms


class TestWritePOSCAR:
    # TODO: Test that the written POSCAR is the same as the read one.
    def test_write_creates_files(self, example_poscars, tmp_path):
        slab = example_poscars
        poscar.write(slab, tmp_path / 'POSCAR')
        assert (tmp_path / 'POSCAR').exists()

    def test_write_does_not_modify_slab(self, example_poscars, tmp_path):
        slab = example_poscars
        n_atoms_original = len(slab.atlist)
        poscar.write(slab, tmp_path / 'POSCAR')
        assert slab == example_poscars                                          # TODO: will always succeed because tests "is" if there's no __eq__
        assert len(slab.atlist) == n_atoms_original

    def test_write_writes_correct_number_of_atoms(self, slab_and_expectations, tmp_path):
        slab, expected_n_atoms, *_ = slab_and_expectations
        poscar.write(slab, tmp_path / 'POSCAR')
        # read the written POSCAR and check that the number of atoms is correct
        written_slab = poscar.read(tmp_path / 'POSCAR')
        assert len(written_slab.atlist) == expected_n_atoms

    def test_write_slab_with_symmetry(self, slab_pg_rp, tmp_path):
        slab, *_ = slab_pg_rp
        poscar.write(slab, tmp_path / 'POSCAR')
        assert (tmp_path / 'POSCAR').exists()

    def test_write_with_symmetry_n_atoms(self, slab_pg_rp, tmp_path):
        slab, *_ = slab_pg_rp
        poscar.write(slab, tmp_path / 'POSCAR')
        written_slab = poscar.read(tmp_path / 'POSCAR')
        assert len(written_slab.atlist) == len(slab.atlist)

    def test_write_with_symmetry_pg(self, slab_pg_rp, tmp_path):
        slab, pg, rp = slab_pg_rp
        poscar.write(slab, tmp_path / 'POSCAR')
        written_slab = poscar.read(tmp_path / 'POSCAR')
        written_pg = symmetry.findSymmetry(written_slab, rp)
        assert written_pg == pg
