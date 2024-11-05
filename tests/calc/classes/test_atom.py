"""Tests for module viperleed.calc.classes.atom."""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-07-28'
__license__ = 'GPLv3+'

import numpy as np
import pytest
from pytest import approx, mark


# TODO: would need more extensive testing

class TestAtomMergeDisplacementOffset:
    """Test combination of displacements and offsets."""

    def test_add_displ_offset_runtime_error(self, manually_displaced_atom):
        """Test that a non-initialized displacement raises RuntimeError."""
        with pytest.raises(RuntimeError):
            manually_displaced_atom.mergeDisp(None)

    @mark.xfail(
        raises=KeyError,
        reason='Incorrect implementation in Atom for accessing site.vibamp'
        )
    def test_offset_merge_allowed(self, manually_displaced_atom, subtests):
        """Test successful combination of VIBROCC offset and displacements."""
        atom = manually_displaced_atom
        element = atom.el
        atom.offset_geo[element] = np.array([0.1, 0, 0])
        atom.offset_vib[element] = +0.1
        atom.offset_occ[element] = -0.1
        atom.dispInitialized = True
        atom.mergeDisp(element)
        with subtests.test('geo'):
            assert np.allclose(atom.disp_geo[element],
                               [[-0.1, 0, 0], [0.1, 0, 0], [0.3, 0, 0]])
        with subtests.test('vib'):
            assert atom.disp_vib[element] == approx([0.0, 0.1, 0.2])
        with subtests.test('occ'):
            assert atom.disp_occ[element] == approx([0.6, 0.7, 0.8, 0.9])

    @mark.xfail(
        raises=KeyError,
        reason='Incorrect implementation in Atom for accessing site.vibamp'
        )
    def test_displacement_out_of_range(self, manually_displaced_atom,
                                       subtests):
        """Check skipping of offsets that cause out-of-range displacements."""
        atom = manually_displaced_atom
        element = atom.el
        atom.offset_vib[element] = -0.1
        atom.offset_occ[element] = +0.2
        atom.dispInitialized = True
        atom.mergeDisp(element)
        with subtests.test('vib'):
            pytest.xfail(reason=(
                'The current implementation is wrong (does not check site) '
                'fixed by @fkraushofer in installable branch. See '
                'pull/94/commits/4a4a5a5236aa4a2cd03de4e4ec2a448e64cc63d1. '
                'However, the fix is also wrong, as it does not consider '
                'that site.vibamp may not be present. It seems much better '
                'to check for consistency only when we know that we have all '
                'information.'
                ))
            assert atom.disp_vib[element] == approx([-0.2, -0.1, 0.0])
        with subtests.test('occ'):
            assert atom.disp_occ[element] == approx([0.7, 0.8, 0.9, 1.0])


def test_displacement_ranges(manually_displaced_atom):
    atom = manually_displaced_atom
    disp_ranges = atom.disp_ranges

    # without any displacement, the range should be zero
    assert (disp_ranges.geo['all'][0] - disp_ranges.geo['all'][1]).sum() == approx(0)
    assert disp_ranges.vib['all'][0] - disp_ranges.vib['all'][1] == approx(0)
    assert disp_ranges.occ['Ag'][0] - disp_ranges.occ['Ag'][1] == approx(0)

    # set some displacements
    atom.disp_geo['Ag'] = np.array([0.1, 0, 0])
    atom.offset_vib['Ag'] = -0.2
    atom.offset_occ['Ag'] = +0.3
    atom.dispInitialized = True
    # merge the displacements which should update the ranges
    atom.mergeDisp('all')

    # check the new ranges
    assert (
        disp_ranges.geo['Ag'][1] - disp_ranges.geo['Ag'][0]
    ).sum() == approx(0.1)
    assert disp_ranges.vib['Ag'][1] - disp_ranges.vib['Ag'][0] == approx(
        0.2
    )
    assert disp_ranges.occ['Ag'][1] - disp_ranges.occ['Ag'][0] == approx(
        0.3
    )
