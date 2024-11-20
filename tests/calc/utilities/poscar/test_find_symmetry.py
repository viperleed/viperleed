"""Tests for the viperleed poscar find symmetry utility."""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-11-20'
__license__ = 'GPLv3+'

import pytest
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
    find_symmetry_cli = FindSymmetryCLI()
    parser = find_symmetry_cli.parser
    args = parser.parse_args([])

    # Check if the correct symmetry group is found
    found_plane_group = (
        find_symmetry_cli.process_slab(slab, args).foundplanegroup)
    assert found_plane_group == expected_plane_group

@with_info
def test_find_symmetry_to_outfile(test_slab, tmp_path):
    slab, _, info = test_slab
    out_path = tmp_path / f'out_{info.poscar.name}.txt'
    expected_plane_group = info.symmetry.hermann
    cli = FindSymmetryCLI()
    parser = cli.parser
    args = parser.parse_args(['-o', str(out_path)])

    # process the slab and write the output to the outfile
    processed_slab = cli.process_slab(slab, args)
    cli.write_output(processed_slab, args)
    args.outfile.close()

    # Check if the correct symmetry group is written to the outfile
    with open(out_path, 'r') as f:
        read_plane_group = f.read().strip()
    assert read_plane_group == expected_plane_group
