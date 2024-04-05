"""Tests for the viperleed poscar utility modify_vacuum."""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-04-05'
__license__ = 'GPLv3+'

import numpy as np

import pytest_cases

from viperleed.utilities.poscar.modify_vacuum import ModifyVacuumCLI

from ... import poscar_slabs
import pytest
import numpy as np
from viperleed.utilities.poscar.modify_vacuum import ModifyVacuumCLI

class TestModifyVacuum:
    """Tests for the modify_vacuum utility."""

    @pytest.mark.parametrize('vacuum_gap_size', [1.0, 3.14, 10.0, -1.0])
    @pytest_cases.parametrize_with_cases('test_slab', cases=poscar_slabs,
                                         has_tag=poscar_slabs.Tag.NO_INFO)
    def test_modify_vacuum_relative(self, test_slab, vacuum_gap_size):
        """Test the ModifyVacuumCLI class."""
        parser = ModifyVacuumCLI().parser
        slab, *_ = test_slab
        args = parser.parse_args([str(vacuum_gap_size)])
        original_gap = slab.vacuum_gap
        modified_slab = ModifyVacuumCLI().process_slab(slab, args)
        # check if vacuum gap size is modified correctly
        assert (modified_slab.vacuum_gap
                == pytest.approx(vacuum_gap_size + original_gap))

    @pytest.mark.parametrize('vacuum_gap_size', [1.0, 3.14, 10.0,])
    @pytest_cases.parametrize_with_cases('test_slab', cases=poscar_slabs,
                                         has_tag=poscar_slabs.Tag.NO_INFO)
    def test_modify_vacuum_gap_absolute(self, test_slab, vacuum_gap_size):
        """Test the ModifyVacuumCLI class with invalid vacuum gap size."""
        parser = ModifyVacuumCLI().parser
        slab, *_ = test_slab
        args = parser.parse_args([str(vacuum_gap_size), '-a'])
        modified_slab = ModifyVacuumCLI().process_slab(slab, args)
        # check if vacuum gap size is modified correctly
        assert modified_slab.vacuum_gap == pytest.approx(vacuum_gap_size)
