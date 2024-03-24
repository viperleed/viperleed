"""Tests for module viperleed.calc.files.poscar."""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-07-12'
__license__ = 'GPLv3+'

import io

import pytest
from pytest_cases import parametrize_with_cases, fixture

from viperleed.calc import symmetry
from viperleed.calc.files import poscar
from viperleed.calc.lib.base import pairwise

from ...helpers import exclude_tags
from ..poscar_slabs import CasePOSCARSlabs
from ..tags import CaseTag
from . import case_poscar_files


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
        assert slab.n_atoms == info.poscar.n_atoms

    @parametrize_with_cases('args', **_WITH_INFO)
    def test_nr_atom_by_elem_correct(self, args):
        """Ensure that the correct number of atoms was read."""
        slab, *_, info = args
        if not info.poscar.n_atoms_by_elem:
            pytest.skip('No element-resolved atom counts available')
        assert slab.n_per_elem == info.poscar.n_atoms_by_elem


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
        n_atoms_original = slab.n_atoms
        poscar.write(slab, tmp_poscar)
        assert slab.n_atoms == n_atoms_original

    @parametrize_with_cases('args', **_WITH_INFO)
    def test_nr_atoms_expected(self, args, write_and_read):
        """Check that a written slab has the expected n_atoms."""
        slab, *_, info = args
        written_slab = write_and_read(slab)
        assert written_slab.n_atoms == info.poscar.n_atoms

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
        assert written_slab.n_atoms == slab.n_atoms

    def test_group_conserved(self, poscar_with_group, write_and_read):
        """Ensure that writing a slab conserves its detected plane group."""
        slab, rpars, *_ = poscar_with_group
        written_slab = write_and_read(slab)
        symmetry.findSymmetry(written_slab, rpars)
        assert written_slab.foundplanegroup == slab.foundplanegroup

    def test_reorder(self, poscar_with_group, write_and_read):
        """Ensure that a sorted slab is written in the same order."""
        slab, *_ = poscar_with_group
        slab.sort_by_z()
        written_slab = write_and_read(slab, comments='all')
        assert all(a1.pos[2] <= a2.pos[2] or a1.el != a2.el
                   for a1, a2 in pairwise(written_slab))
        assert all(a1.el == a2.el for a1, a2 in zip(slab, written_slab))


@fixture(name='poscar_stream')
@parametrize_with_cases('args', cases=case_poscar_files)
def fixture_poscar_stream(args):
    """Return a POSCARStreamReader object."""
    file, contents = args
    stream = io.StringIO(contents)
    return file, poscar.POSCARStreamReader(stream)


class TestPOSCARStreamReader:
    """Collection of tests for reading a POSCAR file from a stream."""

    def test_read_slab_has_atoms(self, poscar_stream):
        """Ensure that some atoms were read."""
        _, reader = poscar_stream
        slab = reader.read()
        assert slab.atlist

    def test_slab_matches_read_from_file(self, poscar_stream):
        """Check consistency when reading the same file as a stream or file."""
        file, reader = poscar_stream
        stream_slab = reader.read()
        file_slab = poscar.read(file)
        assert stream_slab.n_atoms == file_slab.n_atoms
        assert stream_slab.n_per_elem == file_slab.n_per_elem


@fixture(name='output_stream')
def fixture_output_stream():
    """Return a StringIO object to write to."""
    return io.StringIO()


class TestPOSCARStreamWriter:
    """Collection of tests for writing a POSCAR file to a stream."""

    @parametrize_with_cases('args', cases=CasePOSCARSlabs)
    def test_stream_written(self, args, output_stream):
        """Ensure that a POSCAR file exists after writing."""
        slab, *_ = args
        poscar_steam_writer = poscar.POSCARStreamWriter(output_stream)
        poscar_steam_writer.write(slab)
        written_contents = output_stream.getvalue()
        assert written_contents

    @parametrize_with_cases('args', cases=CasePOSCARSlabs)
    def test_written_stream_is_valid(self, args, output_stream):
        """Ensure that the written stream is a valid POSCAR."""
        slab, *_ = args
        poscar_steam_writer = poscar.POSCARStreamWriter(output_stream)
        poscar_steam_writer.write(slab)
        written_contents = output_stream.getvalue()
        written_slab = poscar.read(io.StringIO(written_contents))
        assert written_slab.atlist
        assert written_slab.n_atoms == slab.n_atoms
        assert written_slab.n_per_elem == slab.n_per_elem


class TestPOSCARVaspWrite:
    """Collection of tests for writing a VASP-style POSCAR file."""

    @parametrize_with_cases('args', cases=CasePOSCARSlabs)
    def test_selective_dynamics_tag_is_added(self, args, output_stream):
        """Ensure that the Selective Dynamics line is added to the POSCAR."""
        slab, *_ = args
        vasp_writer = poscar.VASPPOSCARWriter(stream=output_stream,
                                              relax_info={'c_only': True})
        vasp_writer.write(slab)
        content = output_stream.getvalue()
        assert 'selective dynamics' in content.lower()

    @parametrize_with_cases('args', cases=CasePOSCARSlabs)
    def test_flag_c_only(self, args, output_stream):
        """Check correct constraints for displacing only along z.""" 
        slab, *_ = args
        vasp_writer = poscar.VASPPOSCARWriter(output_stream,
                                              relax_info={'c_only': True})
        vasp_writer.write(slab)
        content = output_stream.getvalue()
        assert '   T   T   T' not in content
        assert '   F   F   T' in content

    @parametrize_with_cases('args', cases=CasePOSCARSlabs)
    def test_flag_c_only_off(self, args, output_stream):
        """Check correct constraints for displacing in all 3 directions.""" 
        slab, *_ = args
        vasp_writer = poscar.VASPPOSCARWriter(output_stream,
                                              relax_info={'c_only': False})
        vasp_writer.write(slab)
        content = output_stream.getvalue()
        assert '   T   T   T' in content
        assert '   F   F   T' not in content

    @parametrize_with_cases('args', cases=CasePOSCARSlabs)
    def test_flag_above_c(self, args, output_stream):
        """Check correct constraints for displacing only above a given c.""" 
        slab, *_ = args
        above_c = 0.3
        vasp_writer = poscar.VASPPOSCARWriter(output_stream,
                                              relax_info={'c_only': True,
                                                          'above_c': above_c})
        vasp_writer.write(slab)
        content = output_stream.getvalue()
        content_lines = content.split('\n')
        for line in content_lines[9:]:
            if not line:
                continue  # Skip empty lines
            z_coord = float(line.split()[2])
            if z_coord > above_c:
                assert '   F   F   T' in line
            else:
                assert '   F   F   F' in line

    @parametrize_with_cases('args', cases=CasePOSCARSlabs)
    def test_flag_c_between_0_and_1(self, args, output_stream):
        """Check complaints for invalid c value.""" 
        slab, *_ = args
        with pytest.raises(ValueError):
            vasp_writer = poscar.VASPPOSCARWriter(output_stream,
                                                  relax_info={'c_only': True,
                                                              'above_c': 1.2})
            vasp_writer.write(slab)
