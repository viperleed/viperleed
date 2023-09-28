"""Tests for viperleed.tleedmlib.files.displacements.

Created on 2023-07-28

@author: Alexander M. Imre (@amimre)
@author: Michele Riva (@michele-riva)
"""

import numpy as np
from pytest import approx


class TestReadDISPLACEMENTS:
    """Test successful reading of a DISPLACEMENTS block."""

    def test_read_geo(self, displaced_atom):
        """Test successful reading of a geometric DISPLACEMENTS block."""
        displ_read = np.array(displaced_atom.disp_geo['all'])
        assert displ_read[:, 2] == approx([0.2, 0.1, 0.0, -0.1, -0.2])

    def test_read_vib(self, displaced_atom):
        """Test successful reading of a vibration DISPLACEMENTS block."""
        displ_read = displaced_atom.disp_vib['all']
        assert displ_read == approx([-0.1, 0.0, 0.1])

    def test_read_occ(self, displaced_atom):
        """Test successful reading of an occupation DISPLACEMENTS block."""
        displ_read = displaced_atom.disp_occ[displaced_atom.el]
        assert displ_read == approx([0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
