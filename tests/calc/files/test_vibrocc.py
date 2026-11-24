"""Tests for module viperleed.calc.files.vibrocc."""

__authors__ = (
    'Alexandra Mia Imre (@alexmiame)',
    'Michele Riva (@michele-riva)',
    'Florian Kraushofer (@fkraushofer)',
)
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2023-07-28'
__license__ = 'GPLv3+'

import logging

import pytest

from viperleed.calc.files.vibrocc import readVIBROCC

from ..poscar_slabs import CasePOSCARSlabs


# TODO: would need more extensive testing

class TestReadVIBROCC:
    """Collection of tests for reading a VIBROCC file."""

    def test_read_offset_occ(self, displaced_atom):
        """Read correctly an occupation offset."""
        atom = displaced_atom
        assert atom.offset_occ[atom.el] == pytest.approx(-0.1)

    def test_interpret_offset_allowed(self, displaced_atom):
        """Interpret correctly an allowed occupation offset."""
        slab = displaced_atom.slab
        for atom in slab:
            atom.dispInitialized = True
            atom.mergeDisp(atom.el)
        atom = displaced_atom
        expected = pytest.approx([0.4, 0.5, 0.6, 0.7, 0.8, 0.9])
        assert atom.disp_occ[atom.el] == expected

    def test_interpret_offset_not_allowed(self, displaced_atom):                 # TODO: This is kind of duplicated in test_atom
        """Interpret correctly a non-allowed occupation offset."""
        atom = displaced_atom
        atom.offset_occ[atom.el] = +0.2
        atom.dispInitialized = True
        atom.mergeDisp(atom.el)
        expected = pytest.approx([0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
        assert atom.disp_occ[atom.el] == expected


class TestGenerateVIBROCC:
    """Tests for generating vibration amplitudes from temperatures."""

    def test_empty_file(self, ag100, tmp_path):
        """Check amplitudes are generated for an empty VIBROCC."""
        slab, rpars, *_ = ag100
        rpars.T_DEBYE = 330
        rpars.T_EXPERIMENT = 100
        vibrocc = tmp_path / 'VIBROCC'
        vibrocc.touch()
        generated = readVIBROCC(rpars, slab, vibrocc)
        assert generated
        assert not rpars.halt
        for site in slab.sitelist:
            assert site.vibamp[site.el] > 0

    def test_generate_mixed_site(self, tmp_path):
        """Check amplitudes are generated for a mixed site only."""
        slab, rpars, _ = CasePOSCARSlabs().case_poscar_ptrh_100_1x3_o()
        slab.full_update(rpars)  # Chemical elements of ELEMENT_MIX
        rpars.T_DEBYE = 300
        rpars.T_EXPERIMENT = 100
        vibrocc = tmp_path / 'VIBROCC'
        vibrocc.write_text('= Vibrational Amplitudes\nO_def = 0.12\n',
                           encoding='utf-8')
        generated = readVIBROCC(rpars, slab, vibrocc)
        assert generated
        assert not rpars.halt
        o_site, mixed_site = (next(s for s in slab.sitelist if s.el == el)
                              for el in ('O', 'Me'))
        assert o_site.vibamp['O'] == pytest.approx(0.12)
        assert mixed_site.vibamp['Pt'] > 0
        assert mixed_site.vibamp['Rh'] > 0

    def test_no_file_no_temperatures(self, ag100, tmp_path, caplog):
        """Check the error for a missing VIBROCC suggests T_DEBYE."""
        slab, rpars, *_ = ag100
        caplog.set_level(logging.ERROR)
        with pytest.raises(FileNotFoundError):
            readVIBROCC(rpars, slab, tmp_path / 'VIBROCC')
        assert 'T_DEBYE and T_EXPERIMENT' in caplog.text
