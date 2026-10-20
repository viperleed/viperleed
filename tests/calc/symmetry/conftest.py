"""Test configuration for viperleed.tests.calc.symmetry."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-09-08'
__license__ = 'GPLv3+'

import numpy as np
import pytest

from pytest_cases import fixture, parametrize_with_cases

from viperleed.calc import symmetry

from ...helpers import duplicate_all
from ...helpers import flat_fixture
from .. import poscar_slabs
from . import simple_slabs


def get_cases(which='all'):
    """Return a sequence of test cases."""
    cases = []
    if which in 'all|simple':
        cases.append(simple_slabs)
    if which in 'all|poscar':
        cases.append(poscar_slabs)
    return cases


slab_p1 = flat_fixture(simple_slabs.CaseSimpleSlabs().case_p1)
slab_pmg = flat_fixture(simple_slabs.CaseSimpleSlabs().case_pmg)
slab_p6m = flat_fixture(simple_slabs.CaseSimpleSlabs().case_p6m)


@fixture(name='with_plane_group')
@parametrize_with_cases('args', cases=get_cases('all'))
def fixture_factory_with_plane_group(args, draws):
    """Find plane group, including a (random) rigid origin shift.

    The returned factory takes one of the `draws` and returns
    slab, Rparams, and info that belong to that draw only.
    """
    slab, param, info, *_ = args
    # Each draw gets its own pristine copy of the case, and its own
    # shift. Both are taken before anything modifies the objects.
    per_draw = {
        draw.index: (
            *duplicate_all(slab, param, info),
            draw.rng.uniform(-1, 1, size=2),
            )
        for draw in draws
        }

    def _make(draw, random_shifts=True):
        slab, param, info, shift = per_draw[draw.index]
        if random_shifts:
            slab.translate_atoms_2d(shift)
            slab.collapse_cartesian_coordinates()
        symmetry.findSymmetry(slab, param, forceFindOri=True)
        return slab, param, info
    return _make


@fixture(scope='session')
def displace_atoms():
    """Move atoms by displacements (if allowed) and update symmetry group."""
    _mode = 4

    def _is_displacement_forbidden(displacements, info):
        """Return whether a given displacement direction is forbidden."""
        if not info.atom_on_axis:
            return False
        displacements = np.asarray(displacements)
        return np.any(abs(displacements[:, :2]) > 1e-3)

    def _displace(slab, param, info, displacements, check_raises=True):
        for (displ_info, displ) in zip(info.displacements, displacements):
            displ = (displ,)
            atom = slab.atlist.get(displ_info.atom_nr)
            if check_raises and _is_displacement_forbidden(displ, displ_info):
                pytest.xfail('No complaint about invalid directions yet')
                with pytest.raises(Exception):
                    atom.assignDisp(_mode, displ)
            else:
                atom.assignDisp(_mode, displ)
        for atom in slab:
            try:
                disp = atom.disp_geo_offset[atom.el]
            except KeyError:
                disp = atom.disp_geo_offset['all']
            disp = disp[0]
            atom.cartpos += disp
        slab.update_fractional_from_cartesian()
        symmetry.findSymmetry(slab, param)
    return _displace


@fixture
def with_symmetry_constraints(with_plane_group):
    """Apply symmetry constraints, then return slab, Rparams, and info."""
    def _make(draw, random_shifts=True):
        slab, param, *rest = with_plane_group(draw, random_shifts)
        symmetry.enforceSymmetry(slab, param)
        return (slab, param, *rest)
    return _make
