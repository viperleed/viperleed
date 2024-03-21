"""Test configuration for viperleed.tests.symmetry.

Created on 2023-09-08

@author: Michele Riva (@michele-riva)
"""

from pathlib import Path
import sys

import numpy as np
import pytest

from pytest_cases import fixture, parametrize_with_cases

VPR_PATH = str(Path(__file__).resolve().parents[3])
if VPR_PATH not in sys.path:
    sys.path.append(VPR_PATH)

# pylint: disable=wrong-import-position
# Cannot do anything about it until we make viperleed installable
from viperleed.tleedmlib import symmetry

from ..helpers import flat_fixture
from . import simple_slabs
from .. import poscar_slabs
# pylint: enable=wrong-import-position


RANDOM = np.random.default_rng().uniform


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
def fixture_factory_with_plane_group(args):
    """Find plane group, including a (random) rigid origin shift."""
    slab, param, info, *_ = args
    shift = RANDOM(-1, 1, size=2)

    def _make(random_shifts=True):
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
    def _make(random_shifts=True):
        slab, param, *rest = with_plane_group(random_shifts)
        symmetry.enforceSymmetry(slab, param)
        return (slab, param, *rest)
    return _make
