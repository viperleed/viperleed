"""Tests for the viperleed poscar find symmetry utility."""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-11-20'
__license__ = 'GPLv3+'

import pytest_cases

from viperleed.utilities.poscar.find_symmetry import FindSymmetryCLI

from ... import poscar_slabs

with_info = pytest_cases.parametrize_with_cases(
    'test_slab',
    cases=poscar_slabs.CasePOSCARSlabs.case_poscar,
)


@with_info
def test_find_symmetry_parser(test_slab):
    slab, _, info = test_slab
    expected_plane_group = info.symmetry.hermann
    parser = FindSymmetryCLI().parser
    args = parser.parse_args([])

    # Check if the correct symmetry group is found
    found_plane_group = (
        FindSymmetryCLI().process_slab(slab, args).foundplanegroup)
    assert found_plane_group == expected_plane_group

