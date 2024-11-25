"""Tests for the viperleed poscar utility modify_vacuum."""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-04-05'
__license__ = 'GPLv3+'

import pytest
from pytest_cases import parametrize
from pytest_cases import parametrize_with_cases

from viperleed.calc.classes.slab.errors import VacuumError
from viperleed.calc.classes.slab.surface_slab import _MIN_VACUUM
from viperleed.utilities.poscar.modify_vacuum import ModifyVacuumCLI
from viperleed.utilities.poscar.modify_vacuum import VacuumGapInfo
from viperleed.utilities.poscar.modify_vacuum import modify_vacuum

from ... import poscar_slabs
from ...tags import CaseTag as Tag

infoless = parametrize_with_cases('test_slab',
                                  cases=poscar_slabs,
                                  has_tag=Tag.NO_INFO)
SINGLE_LAYER = 1e-8  # Slab thickness of a one-layer slab


class TestModifyVacuum:
    """Tests for the modify_vacuum utility."""

    @infoless
    def test_accept_small_gap(self, test_slab):
        """Test that NotEnoughVacuumError is raised correctly."""
        slab, *_ = test_slab
        if slab.thickness <= SINGLE_LAYER:
            pytest.skip('Single layer; would lead to zero volume cell')
        gap = VacuumGapInfo(size=0.0, absolute=True, accept_small_gap=True)
        # Gap size recognition likely does not work for slabs with
        # gaps < 5AA so we shouldn't test for the exact gap size
        modify_vacuum(slab, gap)

    @parametrize(vacuum_gap_size=[1.0, 3.14, 15.0, -1.0])
    @infoless
    def test_gap_relative(self, test_slab, vacuum_gap_size):
        """Test the ModifyVacuumCLI class with relative vacuum gap."""
        cli = ModifyVacuumCLI()
        slab, *_ = test_slab
        original_gap = slab.vacuum_gap
        desired_gap = vacuum_gap_size + original_gap

        args = cli.parse_cli_args([str(vacuum_gap_size)])
        if desired_gap < _MIN_VACUUM:
            with pytest.raises(SystemExit):
                modified_slab = cli.process_slab(slab, args)
        else:
            modified_slab = cli.process_slab(slab, args)
            # Check that vacuum gap size is modified correctly
            assert modified_slab.vacuum_gap == pytest.approx(desired_gap)

    @parametrize(vacuum_gap_size=[15.0, 100.0])
    @infoless
    def test_gap_absolute(self, test_slab, vacuum_gap_size):
        """Test the ModifyVacuumCLI class with absolute vacuum gap size."""
        cli = ModifyVacuumCLI()
        slab, *_ = test_slab
        args = cli.parse_cli_args([str(vacuum_gap_size), '-a'])
        modified_slab = cli.process_slab(slab, args)
        # Check that vacuum gap size is modified correctly
        assert modified_slab.vacuum_gap == pytest.approx(vacuum_gap_size)
        assert modified_slab.thickness == pytest.approx(slab.thickness)

    @parametrize(vacuum_gap_size=[1.0])
    @infoless
    def test_gap_too_small(self, test_slab, vacuum_gap_size):
        """Test the ModifyVacuumCLI class with absolute vacuum gap size."""
        cli = ModifyVacuumCLI()
        slab, *_ = test_slab
        args = cli.parse_cli_args([str(vacuum_gap_size), '-a'])
        with pytest.raises(SystemExit):
            cli.process_slab(slab, args)

    @parametrize(vacuum_gap_size=[-1.0, -5.0])
    @infoless
    def test_negative_gap_raises(self, test_slab, vacuum_gap_size):
        """Check complaints when the vacuum gap size is negative."""
        slab, *_ = test_slab
        vacuum_gap_info = VacuumGapInfo(size=vacuum_gap_size, absolute=True)
        with pytest.raises(VacuumError):
            modify_vacuum(slab, vacuum_gap_info)

    @infoless
    def test_not_enough_vacuum_error(self, test_slab):
        """Test that NotEnoughVacuumError is raised correctly."""
        slab, *_ = test_slab
        gap = VacuumGapInfo(size=0.5, absolute=True, accept_small_gap=False)
        with pytest.raises(VacuumError):
            modify_vacuum(slab, gap)

    @parametrize(vacuum_gap_size=[1.0, 3.14, 15.0, -1.0])
    @infoless
    def test_shifted_slab_down(self, test_slab, vacuum_gap_size):
        """Test that we get the same result if we shift the slab down."""
        cli = ModifyVacuumCLI()
        slab, *_ = test_slab
        min_z = min(at.pos[2] for at in slab)
        for at in slab:
            at.pos[2] -= min_z
        args = cli.parse_cli_args([str(vacuum_gap_size)])
        original_gap = slab.vacuum_gap
        desired_gap = vacuum_gap_size + original_gap
        if desired_gap < _MIN_VACUUM:
            with pytest.raises(SystemExit):
                modified_slab = cli.process_slab(slab, args)
        else:
            modified_slab = cli.process_slab(slab, args)
            # Check that vacuum gap size is modified correctly
            assert modified_slab.vacuum_gap == pytest.approx(desired_gap)

    @infoless
    def test_zero_gap_raises(self, test_slab):
        """Check complaints when the vacuum gap size is negative."""
        slab, *_ = test_slab
        vacuum_gap_info = VacuumGapInfo(size=0., absolute=True)
        if slab.thickness <= SINGLE_LAYER:
            pytest.skip('Single layer; would lead to zero-volume cell')
        with pytest.raises(VacuumError):
            modify_vacuum(slab, vacuum_gap_info)


class TestParseCLIArgs:
    """Tests for the ModifyVacuumCLI.parse_cli_args method."""

    def test_absolute_negative_gap_raises(self):
        """Test that parse_cli_args raises for a negative absolute gap."""
        cli = ModifyVacuumCLI()
        with pytest.raises(SystemExit):
            cli.parse_cli_args(['-1.0', '--absolute'])

    @parametrize(vacuum_gap_size=[1.0, 3.14, 15.0])
    def test_gap_absolute(self, vacuum_gap_size):
        """Test the ModifyVacuumCLI class with relative vacuum gap."""
        cli = ModifyVacuumCLI()
        parsed_args = cli.parse_cli_args([str(vacuum_gap_size), '-a'])
        assert parsed_args.vacuum == vacuum_gap_size
        assert parsed_args.absolute is True

    @parametrize(vacuum_gap_size=[1.0, 3.14, 15.0, -1.0])
    def test_gap_relative(self, vacuum_gap_size):
        """Test the ModifyVacuumCLI class with relative vacuum gap."""
        cli = ModifyVacuumCLI()
        parsed_args = cli.parse_cli_args([str(vacuum_gap_size)])
        assert parsed_args.vacuum == vacuum_gap_size
        assert parsed_args.absolute is False
