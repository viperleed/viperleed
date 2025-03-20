"""Tests for the viperleed poscar delete above/below/between utilities."""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-04-03'
__license__ = 'GPLv3+'

import numpy as np
import pytest
from pytest_cases import parametrize
from pytest_cases import parametrize_with_cases

from viperleed.utilities.poscar.delete_above import DeleteAboveCLI
from viperleed.utilities.poscar.delete_below import DeleteBelowCLI
from viperleed.utilities.poscar.delete_between import DeleteBetweenCLI

from ... import poscar_slabs
from ...tags import CaseTag as Tag

VALID_CUT_FRACTIONS = [0.1234, 0.5, 0.75]
INVALID_CUT_FRACTIONS = [-0.1, 1.1, 0.0, 'abc']

infoless = parametrize_with_cases('test_slab',
                                  cases=poscar_slabs,
                                  has_tag=Tag.NO_INFO)

class TestDeleteAbove:
    """Tests for the delete_above utility."""

    @parametrize(c_max=VALID_CUT_FRACTIONS)
    @infoless
    def test_delete_above_cli(self, test_slab, c_max):
        """Test the DeleteAboveCLI class."""
        slab, *_ = test_slab
        original_ucell = slab.ucell.copy()
        n_atoms_below_c = sum(1 for atom in slab if atom.pos[2] <= c_max)
        cli = DeleteAboveCLI()
        args = cli.parse_cli_args([str(c_max)])
        processed = cli.process_slab(slab, args)
        # right number of atoms
        assert processed.n_atoms == n_atoms_below_c
        # unit cell is unchanged
        assert np.allclose(slab.ucell, original_ucell)
        assert np.allclose(processed.ucell, original_ucell)

    @parametrize(invalid_c=INVALID_CUT_FRACTIONS)
    def test_invalid_c(self, invalid_c):
        """Test the DeleteAboveCLI class."""
        cli = DeleteAboveCLI()
        with pytest.raises(SystemExit):
            cli.parse_cli_args([str(invalid_c)])


class TestDeleteBelow:
    """Tests for the delete_below utility."""

    @parametrize(c_min=VALID_CUT_FRACTIONS)
    @infoless
    def test_delete_below_cli(self, test_slab, c_min):
        """Test the DeleteBelowCLI class."""
        slab, *_ = test_slab
        original_ucell = slab.ucell.copy()
        n_atoms_above_c = sum(1 for atom in slab if atom.pos[2] >= c_min)
        cli = DeleteBelowCLI()
        args = cli.parse_cli_args([str(c_min)])
        processed = cli.process_slab(slab, args)
        # right number of atoms
        assert processed.n_atoms == n_atoms_above_c
        # unit cell is unchanged
        assert np.allclose(slab.ucell, original_ucell)
        assert np.allclose(processed.ucell, original_ucell)

    @parametrize(invalid_c=INVALID_CUT_FRACTIONS)
    def test_invalid_c(self, invalid_c):
        """Test the DeleteBelowCLI class."""
        cli = DeleteBelowCLI()
        with pytest.raises(SystemExit):
            cli.parse_cli_args([str(invalid_c)])


class TestDeleteBetween:
    """Tests for the delete_between utility."""

    @parametrize('c_min, c_max', [(0.2, 0.8), (0.4, 0.6)])
    @infoless
    def test_delete_between_cli(self, test_slab, c_min, c_max):
        """Test the DeleteBetweenCLI class."""
        slab, *_ = test_slab
        original_ucell = slab.ucell.copy()
        n_atoms_in_c_range = sum(1 for atom in slab
                                 if not c_min <= atom.pos[2] <= c_max)
        cli = DeleteBetweenCLI()
        args = cli.parse_cli_args([str(c_min), str(c_max)])
        processed = cli.process_slab(slab, args)
        # right number of atoms
        assert processed.n_atoms == n_atoms_in_c_range
        # unit cell is unchanged
        assert np.allclose(slab.ucell, original_ucell)
        assert np.allclose(processed.ucell, original_ucell)

    @parametrize(invalid_c_min=INVALID_CUT_FRACTIONS)
    @parametrize(invalid_c_max=INVALID_CUT_FRACTIONS)
    def test_invalid_c(self, invalid_c_min, invalid_c_max):
        """Test complaints when c_min/c_max are out of range."""
        cli = DeleteBelowCLI()
        with pytest.raises(SystemExit):
            cli.parse_cli_args([str(invalid_c_min), str(invalid_c_max)])
