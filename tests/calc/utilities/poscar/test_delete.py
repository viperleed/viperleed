"""Tests for the viperleed poscar delete above/below/between utilities."""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-04-03'
__license__ = 'GPLv3+'

import numpy as np

import pytest
import pytest_cases

from viperleed.utilities.poscar.delete_above import DeleteAboveCLI
from viperleed.utilities.poscar.delete_below import DeleteBelowCLI
from viperleed.utilities.poscar.delete_between import DeleteBetweenCLI

from ... import poscar_slabs

VALID_CUT_FRACTIONS = [0.1234, 0.5, 0.75]
INVALID_CUT_FRACTIONS = [-0.1, 1.1, 0.0, 'abc']

class TestDeleteAbove:
    """Tests for the delete_above utility."""

    @pytest.mark.parametrize('c', VALID_CUT_FRACTIONS)
    @pytest_cases.parametrize_with_cases('test_slab', cases=poscar_slabs,
                                         has_tag=poscar_slabs.Tag.NO_INFO)
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

    @pytest.mark.parametrize('c', INVALID_CUT_FRACTIONS)
    def test_invalid_c(self, ag100, c):
        """Test the DeleteAboveCLI class."""
        parser = DeleteAboveCLI().parser
        slab, *_ = ag100
        with pytest.raises(SystemExit):
            args = parser.parse_args([str(c)])


class TestDeleteBelow:
    """Tests for the delete_below utility."""

    @pytest.mark.parametrize('c', VALID_CUT_FRACTIONS)
    @pytest_cases.parametrize_with_cases('test_slab', cases=poscar_slabs,
                                         has_tag=poscar_slabs.Tag.NO_INFO)
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

    @pytest.mark.parametrize('c', INVALID_CUT_FRACTIONS)
    def test_invalid_c(self, ag100, c):
        """Test the DeleteAboveCLI class."""
        parser = DeleteBelowCLI().parser
        slab, *_ = ag100
        with pytest.raises(SystemExit):
            args = parser.parse_args([str(c)])


class TestDeleteBetween:
    """Tests for the delete_between utility."""

    @pytest.mark.parametrize('c1, c2', [(0.2, 0.8), (0.4, 0.6)])
    @pytest_cases.parametrize_with_cases('test_slab', cases=poscar_slabs,
                                         has_tag=poscar_slabs.Tag.NO_INFO)
    def test_delete_between_cli(self, test_slab, c1, c2):
        """Test the DeleteBetweenCLI class."""
        parser = DeleteBetweenCLI().parser
        slab, *_ = test_slab
        original_ucell = slab.ucell.copy()
        args = parser.parse_args([str(c1), str(c2)])
        n_atoms_between_c1_c2 = sum(1 for atom in slab if atom.pos[2] < c1
                                    or atom.pos[2] > c2)
        # right number of atoms
        assert DeleteBetweenCLI().process_slab(slab, args).n_atoms == n_atoms_between_c1_c2
        # unit cell is unchanged
        assert np.allclose(slab.ucell, original_ucell)
