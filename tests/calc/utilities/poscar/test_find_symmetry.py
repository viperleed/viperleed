"""Tests for the viperleed poscar find symmetry utility."""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-11-20'
__license__ = 'GPLv3+'

from io import StringIO

from pytest_cases import parametrize_with_cases

from viperleed.utilities.poscar.find_symmetry import FindSymmetryCLI

from ....helpers import POSCAR_PATH
from ... import poscar_slabs

with_info = parametrize_with_cases(
    'test_case',
    cases=poscar_slabs.CasePOSCARSlabs.case_poscar,
    )


class TestFindSymmetry:
    """Tests for the find_symmetry utility."""

    @with_info
    def test_plane_group_to_outfile(self, test_case, tmp_path):
        """Check identification of plane group with an explicit outfile."""
        *_, info = test_case
        poscar_path = POSCAR_PATH/info.poscar.name
        out_path = tmp_path / f'out_{info.poscar.name}.txt'
        expected_plane_group = info.symmetry.hermann
        find_symmetry_cli = FindSymmetryCLI()
        find_symmetry_cli([
            '-o', str(out_path),
            '-i', str(poscar_path),
            ])

        # Check that the correct symmetry group is written to file
        read_plane_group = out_path.read_text(encoding='utf-8').strip()
        assert read_plane_group == expected_plane_group

    @with_info
    def test_plane_group_to_stdout(self, test_case, capsys):
        """Check that the correct plane group is printed to STDOUT."""
        *_, info = test_case
        poscar_path = POSCAR_PATH/info.poscar.name
        expected_plane_group = info.symmetry.hermann
        find_symmetry_cli = FindSymmetryCLI()
        find_symmetry_cli([
            '-i', str(poscar_path),
            ])
        captured = capsys.readouterr()
        assert captured.out.rstrip() == expected_plane_group
