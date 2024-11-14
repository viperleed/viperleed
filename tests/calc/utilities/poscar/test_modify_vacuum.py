"""Tests for the viperleed poscar utility modify_vacuum."""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-04-05'
__license__ = 'GPLv3+'

import numpy as np
import pytest
from pytest_cases import parametrize
from pytest_cases import parametrize_with_cases

from viperleed.calc.classes.slab.errors import NotEnoughVacuumError
from viperleed.calc.classes.slab.surface_slab import _MIN_VACUUM
from viperleed.utilities.poscar import modify_vacuum
from viperleed.utilities.poscar.modify_vacuum import ModifyVacuumCLI
from viperleed.utilities.poscar.modify_vacuum import VacuumGapInfo

from ... import poscar_slabs
from ...tags import CaseTag as Tag

infoless = parametrize_with_cases('test_slab',
                                  cases=poscar_slabs,
                                  has_tag=Tag.NO_INFO)


class TestModifyVacuum:
    """Tests for the modify_vacuum utility."""

    @parametrize('vacuum_gap_size', [1.0, 3.14, 15.0, -1.0])
    @infoless
    def test_gap_relative(self, test_slab, vacuum_gap_size):
        """Test the ModifyVacuumCLI class with relative vacuum gap."""
        parser = ModifyVacuumCLI().parser
        slab, *_ = test_slab
        args = parser.parse_args([str(vacuum_gap_size)])
        original_gap = slab.vacuum_gap
        if vacuum_gap_size + original_gap < _MIN_VACUUM:
            with pytest.raises(SystemExit):
                modified_slab = ModifyVacuumCLI().process_slab(slab, args)
        else:
            modified_slab = ModifyVacuumCLI().process_slab(slab, args)
            # check if vacuum gap size is modified correctly
            assert (modified_slab.vacuum_gap
                    == pytest.approx(vacuum_gap_size + original_gap))

    @parametrize('vacuum_gap_size', [15.0, 100.0])
    @infoless
    def test_gap_absolute(self, test_slab, vacuum_gap_size):
        """Test the ModifyVacuumCLI class with absolute vacuum gap size."""
        parser = ModifyVacuumCLI().parser
        slab, *_ = test_slab
        args = parser.parse_args([str(vacuum_gap_size), '-a'])
        modified_slab = ModifyVacuumCLI().process_slab(slab, args)
        # check if vacuum gap size is modified correctly
        assert modified_slab.vacuum_gap == pytest.approx(vacuum_gap_size)
        assert modified_slab.thickness == pytest.approx(slab.thickness)

    @parametrize('vacuum_gap_size', [1.0])
    @infoless
    def test_gap_too_small(self, test_slab, vacuum_gap_size):
        """Test the ModifyVacuumCLI class with absolute vacuum gap size."""
        parser = ModifyVacuumCLI().parser
        slab, *_ = test_slab
        args = parser.parse_args([str(vacuum_gap_size), '-a'])
        with pytest.raises(SystemExit):
            modified_slab = ModifyVacuumCLI().process_slab(slab, args)

    def test_parse_cli_args_absolute_negative_gap(self):
        """Test that parse_cli_args raises for a negative absolute gap."""
        parser = ModifyVacuumCLI()
        with pytest.raises(SystemExit):
            parser.parse_cli_args(['-1.0', '--absolute'])

    @parametrize('vacuum_gap_size', [-1.0, -5.0])
    @infoless
    def test_negative_gap_raises(self, test_slab, vacuum_gap_size):
        """Check complaints when the vacuum gap size is negative."""
        slab, *_ = test_slab
        vacuum_gap_info = VacuumGapInfo(size=vacuum_gap_size, absolute=True)
        with pytest.raises(NotEnoughVacuumError):
            modified_slab = modify_vacuum.modify_vacuum(slab, vacuum_gap_info)

    @infoless
    def test_not_enough_vacuum_error(self, test_slab):
        """Test that NotEnoughVacuumError is raised correctly."""
        slab, *_ = test_slab
        gap = VacuumGapInfo(size=0.5, absolute=True, accept_small_gap=False)
        with pytest.raises(RuntimeError):
            modified_slab = modify_vacuum.modify_vacuum(slab, gap)

    @infoless
    def test_zero_gap_raises(self, test_slab):
        """Check complaints when the vacuum gap size is negative."""
        slab, *_ = test_slab
        vacuum_gap_info = VacuumGapInfo(size=0., absolute=True)
        if slab.thickness <= 1e-8:  # monolayer; would lead to zero volume cell
            return
        with pytest.raises(RuntimeError):
            modified_slab = modify_vacuum.modify_vacuum(slab, vacuum_gap_info)
