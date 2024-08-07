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
from ...tags import CaseTag as Tag

VALID_CUT_FRACTIONS = [0.1234, 0.5, 0.75]
INVALID_CUT_FRACTIONS = [-0.1, 1.1, 0.0, 'abc']

infoless = pytest_cases.parametrize_with_cases('test_slab',
                                               cases=poscar_slabs,
                                               has_tag=Tag.NO_INFO)

class TestDeleteAbove:
    """Tests for the delete_above utility."""

    @pytest.mark.parametrize('c', VALID_CUT_FRACTIONS)
    @infoless
    def test_delete_above_cli(self, test_slab, c):
        """Test the DeleteAboveCLI class."""
        parser = DeleteAboveCLI().parser
        slab, *_ = test_slab
        original_ucell = slab.ucell.copy()
        args = parser.parse_args([str(c)])
        n_atoms_below_c = sum(1 for atom in slab if atom.pos[2] <= c)
        processed = DeleteAboveCLI().process_slab(slab, args)
        # right number of atoms
        assert processed.n_atoms == n_atoms_below_c
        # unit cell is unchanged
        assert np.allclose(slab.ucell, original_ucell)
        assert np.allclose(processed.ucell, original_ucell)

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
    @infoless
    def test_delete_below_cli(self, test_slab, c):
        """Test the DeleteBelowCLI class."""
        parser = DeleteBelowCLI().parser
        slab, *_ = test_slab
        original_ucell = slab.ucell.copy()
        args = parser.parse_args([str(c)])
        n_atoms_above_c = sum(1 for atom in slab if atom.pos[2] >= c)
        processed = DeleteBelowCLI().process_slab(slab, args)
        # right number of atoms
        assert processed.n_atoms == n_atoms_above_c
        # unit cell is unchanged
        assert np.allclose(slab.ucell, original_ucell)
        assert np.allclose(processed.ucell, original_ucell)

    @pytest.mark.parametrize('c', INVALID_CUT_FRACTIONS)
    def test_invalid_c(self, ag100, c):
        """Test the DeleteBelowCLI class."""
        parser = DeleteBelowCLI().parser
        slab, *_ = ag100
        with pytest.raises(SystemExit):
            args = parser.parse_args([str(c)])


class TestDeleteBetween:
    """Tests for the delete_between utility."""

    @pytest.mark.parametrize('c_min, c_max', [(0.2, 0.8), (0.4, 0.6)])
    @infoless
    def test_delete_between_cli(self, test_slab, c_min, c_max):
        """Test the DeleteBetweenCLI class."""
        parser = DeleteBetweenCLI().parser
        slab, *_ = test_slab
        original_ucell = slab.ucell.copy()
        args = parser.parse_args([str(c_min), str(c_max)])
        n_atoms_in_c_range = sum(1 for atom in slab
                                 if atom.pos[2] < c_min
                                 or atom.pos[2] > c_max)
        processed = DeleteBetweenCLI().process_slab(slab, args)
        # right number of atoms
        assert processed.n_atoms == n_atoms_in_c_range
        # unit cell is unchanged
        assert np.allclose(slab.ucell, original_ucell)
        assert np.allclose(processed.ucell, original_ucell)
