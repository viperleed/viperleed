"""Tests for module viperleed.calc.files.vibrocc."""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-07-28'
__license__ = 'GPLv3+'

from pytest import approx


# TODO: would need more extensive testing

class TestReadVIBROCC:
    """Collection of tests for reading a VIBROCC file."""

    def test_read_offset_occ(self, displaced_atom):
        """Read correctly an occupation offset."""
        atom = displaced_atom
        assert atom.offset_occ[atom.el] == approx(-0.1)

    def test_interpret_offset_allowed(self, displaced_atom):
        """Interpret correctly an allowed occupation offset."""
        slab = displaced_atom.slab
        for atom in slab:
            atom.dispInitialized = True
            atom.mergeDisp(atom.el)
        atom = displaced_atom
        assert atom.disp_occ[atom.el] == approx([0.4, 0.5, 0.6, 0.7, 0.8, 0.9])

    def test_interpret_offset_not_allowed(self, displaced_atom):                 # TODO: This is kind of duplicated in test_atom
        """Interpret correctly a non-allowed occupation offset."""
        atom = displaced_atom
        atom.offset_occ[atom.el] = +0.2
        atom.dispInitialized = True
        atom.mergeDisp(atom.el)
        assert atom.disp_occ[atom.el] == approx([0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
