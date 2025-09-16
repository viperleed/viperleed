"""Tests for module viperleed.calc.classes.atom."""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-07-28'
__license__ = 'GPLv3+'

import numpy as np
import pytest
from pytest import approx, mark
from pytest_cases import parametrize


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


class TestDisplacementRanges:
    """Tests for the DisplacementRanges class."""

    @staticmethod
    def diff(range_):
        """Return max - min from `range_`."""
        return range_[1] - range_[0]

    def test_no_displacements(self, manually_displaced_atom):
        """Check that ranges are zero without any displacement."""
        atom = manually_displaced_atom
        ranges = atom.disp_ranges
        assert self.diff(ranges.geo['all']).sum() == approx(0)
        assert self.diff(ranges.vib['all']) == approx(0)
        assert self.diff(ranges.occ['Ag']) == approx(0)

    def test_with_displacements(self, manually_displaced_atom):
        """Check expected ranges with displacements set."""
        atom = manually_displaced_atom
        atom.disp_geo['Ag'] = np.array([0.1, 0, 0])
        atom.offset_vib['Ag'] = -0.2
        atom.offset_occ['Ag'] = +0.3
        atom.dispInitialized = True

        # Merge the displacements which should update the ranges
        atom.mergeDisp('all')

        ranges = atom.disp_ranges
        assert self.diff(ranges.geo['Ag']).sum() == approx(0.1)
        assert self.diff(ranges.vib['Ag']) == approx(0.2)
        assert self.diff(ranges.occ['Ag']) == approx(0.3)


class TestDistance:
    """Tests for the distance function."""

    distance_zero = {
        'trivial': ((0, 0, 0), False),
        '1D in a': ((1, 0, 0), False),
        '1D in c': ((0, 0, 1), True),
        '2D corner': ((1, 1, 0), False),
        '3D corner': ((1, 1, 1), True),
        }

    @parametrize('vec,include_c', distance_zero.values(), ids=distance_zero)
    def test_distance_zero(self, manually_displaced_atom, vec, include_c):
        """Check cases in which the atom matches a replica exactly."""
        atom = manually_displaced_atom
        ucell = atom.slab.ucell.T
        assert atom.distance(atom.cartpos + np.dot(vec, ucell),
                             include_c_replicas=include_c) == approx(0)

    def test_distance_non_zero(self, manually_displaced_atom):
        """Check cases in which the distance is not zero."""
        atom = manually_displaced_atom
        ucell = atom.slab.ucell.T
        abs_dist = np.linalg.norm(ucell[2])
        assert atom.distance(atom.cartpos + ucell[2],
                             include_c_replicas=False) == approx(abs_dist)
        cartpos = atom.cartpos + ucell[0]*0.75
        assert atom.distance(cartpos) <= 0.2501 * np.linalg.norm(ucell[0])

    def test_atom_as_cartpos(self, manually_displaced_atom):
        """Test that atoms are interpreted as their cartpos."""
        atom = manually_displaced_atom
        assert atom.distance(atom) == approx(0)
