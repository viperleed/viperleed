"""Tests for readDISPLACEMENTS_block of viperleed.calc.files.displacements."""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2023-07-28'
__license__ = 'GPLv3+'

import numpy as np
from pytest import approx


class TestReadDISPLACEMENTSBlock:
    """Test successful interpretation of a DISPLACEMENTS block."""

    def test_geo(self, displaced_atom):
        """Test successful reading of a geometric DISPLACEMENTS block."""
        displ_read = np.array(displaced_atom.disp_geo['all'])
        assert displ_read[:, 2] == approx([0.2, 0.1, 0.0, -0.1, -0.2])

    def test_geo_inverted_ranges(self, ag100_with_displacements_and_offsets):
        """Test reading geometric DISPLACEMENTS with start > stop."""
        slab, *_ = ag100_with_displacements_and_offsets
        atom = slab.atlist[1]
        displ_read = np.array(atom.disp_geo['all'])
        assert displ_read[:, 2] == approx([-0.2, -0.1, 0.0, 0.1, 0.2])

    def test_occ(self, displaced_atom):
        """Test successful reading of an occupation DISPLACEMENTS block."""
        displ_read = displaced_atom.disp_occ[displaced_atom.el]
        assert displ_read == approx([0.5, 0.6, 0.7, 0.8, 0.9, 1.0])

    def test_vib(self, displaced_atom):
        """Test successful reading of a vibration DISPLACEMENTS block."""
        displ_read = displaced_atom.disp_vib['all']
        assert displ_read == approx([-0.1, 0.0, 0.1])
