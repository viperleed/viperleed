"""Tests for the viperleed poscar delete above/below/between utilities."""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-04-03'
__license__ = 'GPLv3+'

import argparse

import numpy as np

import pytest
import pytest_cases

from viperleed.calc.bookkeeper import BookkeeperCLI
from viperleed.calc.cli import ViPErLEEDCalcCLI


from viperleed.utilities.poscar.delete_above import DeleteAboveCLI
from viperleed.utilities.poscar.delete_below import DeleteBelowCLI
#from viperleed.utilities.poscar.delete_between import DeleteBetweenCLI

from ... import poscar_slabs

TEST_CUT_FRACTIONS = [0.1234, 0.5, 0.75]

class TestDeleteAbove:
    """Tests for the delete_above utility."""

    @pytest.mark.parametrize('c', TEST_CUT_FRACTIONS)
    @pytest_cases.parametrize_with_cases('test_slab', cases=poscar_slabs)
    def test_delete_above_cli(self, test_slab, c):
        """Test the DeleteAboveCLI class."""
        parser = DeleteAboveCLI().parser
        slab, *_ = test_slab
        original_ucell = slab.ucell.copy()
        args = parser.parse_args([str(c)])
        n_atoms_below_c = sum(1 for atom in slab if atom.pos[2] <= c)
        # right number of atoms
        assert DeleteAboveCLI().process_slab(slab, args).n_atoms == n_atoms_below_c
        # unit cell is unchanged
        assert np.allclose(slab.ucell, original_ucell)

class TestDeleteBelow:
    """Tests for the delete_below utility."""

    @pytest.mark.parametrize('c', TEST_CUT_FRACTIONS)
    @pytest_cases.parametrize_with_cases('test_slab', cases=poscar_slabs)
    def test_delete_below_cli(self, test_slab, c):
        """Test the DeleteBelowCLI class."""
        parser = DeleteBelowCLI().parser
        slab, *_ = test_slab
        original_ucell = slab.ucell.copy()
        args = parser.parse_args([str(c)])
        n_atoms_above_c = sum(1 for atom in slab if atom.pos[2] >= c)
        # right number of atoms
        assert DeleteBelowCLI().process_slab(slab, args).n_atoms == n_atoms_above_c
        # unit cell is unchanged
        assert np.allclose(slab.ucell, original_ucell)
