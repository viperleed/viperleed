"""Tests for module viperleed.tleedmlib.files.poscar.

Created on 2023-07-12

@author: Alexander M. Imre (@amimre)
@author: Michele Riva (@michele-riva)
"""

from pathlib import Path
import sys

from pytest_cases import parametrize_with_cases, fixture

VPR_PATH = str(Path(__file__).resolve().parents[3])
if VPR_PATH not in sys.path:
    sys.path.append(VPR_PATH)

# pylint: disable=wrong-import-position
# Will be fixed in installable version
from viperleed.tleedmlib import symmetry
from viperleed.tleedmlib.files import poscar

from ..poscar_slabs import CasePOSCARSlabs
from ..helpers import exclude_tags, CaseTag
# pylint: enable=wrong-import-position


_WITH_INFO = {'cases': CasePOSCARSlabs,
              'filter': exclude_tags(CaseTag.NO_INFO)}


class TestPOSCARRead:
    """Collection of tests for reading a POSACR file."""

    @parametrize_with_cases('args', cases=CasePOSCARSlabs)
    def test_read_slab_has_atoms(self, args):
        """Ensure that some atoms were read."""
        slab, *_ = args
        assert slab.atlist

    @parametrize_with_cases('args', **_WITH_INFO)
    def test_nr_atom_correct(self, args):
        """Ensure that the correct number of atoms was read."""
        slab, *_, info = args
        assert len(slab.atlist) == info.poscar.n_atoms


# TODO: Test that the written POSCAR is the same as the read one.
class TestPOSCARWrite:
    """Collection of tests for writing a POSACR file."""

    @fixture(name='write_and_read')
    def factory_write_and_read(self, tmp_poscar):
        """Return the slab read after having written it to file."""
        def write_then_read(slab, **kwargs):
            poscar.write(slab, tmp_poscar, **kwargs)
            return poscar.read(tmp_poscar)
        return write_then_read

    @staticmethod
    def check_not_empty(file_path):
        """Ensure that file_path points to a non-empty file."""
        assert file_path.exists()
        with file_path.open('r', encoding='utf-8') as written:
            assert written.readline()

    @parametrize_with_cases('args', cases=CasePOSCARSlabs)
    def test_file_written(self, args, tmp_poscar):
        """Ensure that a POSCAR file exists after writing."""
        slab, *_ = args
        poscar.write(slab, tmp_poscar)
        self.check_not_empty(tmp_poscar)

    @parametrize_with_cases('args', cases=CasePOSCARSlabs)
    def test_ori_slab_nr_atoms_unchanged(self, args, tmp_poscar):
        """Ensure that writing does not change the number of atoms."""
        slab, *_ = args
        n_atoms_original = len(slab.atlist)
        poscar.write(slab, tmp_poscar)
        assert len(slab.atlist) == n_atoms_original

    @parametrize_with_cases('args', **_WITH_INFO)
    def test_nr_atoms_expected(self, args, write_and_read):
        """Check that a written slab has the expected n_atoms."""
        slab, *_, info = args
        written_slab = write_and_read(slab)
        assert len(written_slab.atlist) == info.poscar.n_atoms

    def test_write_slab_with_symmetry(self, poscar_with_group, tmp_poscar):
        """Ensure symmetry information are written without errors."""
        slab, *_ = poscar_with_group
        poscar.write(slab, tmp_poscar, comments='all')
        self.check_not_empty(tmp_poscar)
        with tmp_poscar.open('r', encoding='utf-8') as written:
            sym_lines = (line for line in written if 'group' in line.lower())
            sym_line = next(sym_lines, '')
            assert str(slab.planegroup) in sym_line

    def test_nr_atoms_conserved(self, poscar_with_group, write_and_read):
        """Ensure that writing a slab conserves the number of atoms."""
        slab, *_ = poscar_with_group
        written_slab = write_and_read(slab, comments='all')
        assert len(written_slab.atlist) == len(slab.atlist)

    def test_group_conserved(self, poscar_with_group, write_and_read):
        """Ensure that writing a slab conserves its detected plane group."""
        slab, rpars, *_ = poscar_with_group
        written_slab = write_and_read(slab)
        symmetry.findSymmetry(written_slab, rpars)
        assert written_slab.foundplanegroup == slab.foundplanegroup
