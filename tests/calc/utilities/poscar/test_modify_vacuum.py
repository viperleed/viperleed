"""Tests for the viperleed poscar utility modify_vacuum."""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-04-05'
__license__ = 'GPLv3+'

import numpy as np
import pytest
import pytest_cases

from viperleed.utilities.poscar.modify_vacuum import ModifyVacuumCLI

from ... import poscar_slabs
from ...tags import CaseTag as Tag

infoless = pytest_cases.parametrize_with_cases('test_slab',
                                               cases=poscar_slabs,
                                               has_tag=Tag.NO_INFO)


class TestModifyVacuum:
    """Tests for the modify_vacuum utility."""

    @pytest.mark.parametrize('vacuum_gap_size', [1.0, 3.14, 15.0, -1.0])
    @infoless
    def test_modify_vacuum_relative(self, test_slab, vacuum_gap_size):
        """Test the ModifyVacuumCLI class with relative vacuum gap."""
        parser = ModifyVacuumCLI().parser
        slab, *_ = test_slab
        args = parser.parse_args([str(vacuum_gap_size)])
        original_gap = slab.vacuum_gap
        if vacuum_gap_size + original_gap < 5.0:
            with pytest.raises(SystemExit):
                modified_slab = ModifyVacuumCLI().process_slab(slab, args)
        else:
            modified_slab = ModifyVacuumCLI().process_slab(slab, args)
            # check if vacuum gap size is modified correctly
            assert (modified_slab.vacuum_gap
                    == pytest.approx(vacuum_gap_size + original_gap))

    @pytest.mark.parametrize('vacuum_gap_size', [15.0, 100.0])
    @infoless
    def test_modify_vacuum_gap_absolute(self, test_slab, vacuum_gap_size):
        """Test the ModifyVacuumCLI class with absolute vacuum gap size."""
        parser = ModifyVacuumCLI().parser
        slab, *_ = test_slab
        args = parser.parse_args([str(vacuum_gap_size), '-a'])
        modified_slab = ModifyVacuumCLI().process_slab(slab, args)
        # check if vacuum gap size is modified correctly
        assert modified_slab.vacuum_gap == pytest.approx(vacuum_gap_size)
        assert modified_slab.thickness == pytest.approx(slab.thickness)

    @pytest.mark.parametrize('vacuum_gap_size', [1.0])
    @infoless
    def test_modify_vacuum_gap_too_small(self, test_slab, vacuum_gap_size):
        """Test the ModifyVacuumCLI class with absolute vacuum gap size."""
        parser = ModifyVacuumCLI().parser
        slab, *_ = test_slab
        args = parser.parse_args([str(vacuum_gap_size), '-a'])
        with pytest.raises(SystemExit):
            modified_slab = ModifyVacuumCLI().process_slab(slab, args)
