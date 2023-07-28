"""classes/test_atom.py

Created on 2023-07-28

@author: Alexander M. Imre
"""

import pytest

import numpy as np


class Test_Atom_mergeDisp:
    def test_atom_mergeDisp_allowed(self, atom_with_disp_and_offset):
        """Test method mergeDisp of Atom.
        Offsets are allowed and should be combined with stored
        displacements.
        """
        atom = atom_with_disp_and_offset
        el = atom.el
        atom.offset_geo[el] = +0.1
        atom.offset_vib[el] = +0.1
        atom.offset_occ[el] = -0.1
        atom.mergeDisp(el)
        assert np.allclose(atom.disp_geo[el], [-0.1, 0.1, 0.3])
        assert np.allclose(atom.disp_vib[el], [0.0, 0.1, 0.2])
        assert np.allclose(atom.disp_occ[el], [0.6, 0.7, 0.8, 0.9])

    def test_atom_mergeDisp_not_allowed(self, atom_with_disp_and_offset):
        """Test method mergeDisp of Atom.
        Offsets are not allowed and should not be combined with stored
        displacements. Instead the final displacements should be unchanged.
        """
        atom = atom_with_disp_and_offset
        el = atom.el
        atom.offset_vib[el] = -0.1
        atom.offset_occ[el] = +0.2
        atom.mergeDisp(el)
        assert np.allclose(atom.disp_vib[el], [-0.1, 0.0, 0.1])
        assert np.allclose(atom.disp_occ[el], [0.7, 0.8, 0.9, 1.0])
