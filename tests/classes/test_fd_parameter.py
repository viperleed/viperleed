"""Tests for viperleed.tleedmlib.classes.fd_parameter.

Created on 2024-01-05

@author: Alexander M. Imre (@amimre)
"""
from copy import deepcopy
import pytest

from viperleed.tleedmlib.classes.fd_optimizer.fd_parameter import FDParameter
from viperleed.tleedmlib.classes.fd_optimizer.fd_parameter import FD_PARAMETERS
from ..poscar_slabs import SLAB_Cu2O_111

def test_FDParameter():
    """Test that FDParameter is initialized correctly."""
    
    fd_parameter = FDParameter('v0i', 0.5, (0, 1), lambda r, s, v: setattr(r, "V0_IMAG", v))
    assert fd_parameter.name == 'v0i'
    assert fd_parameter.start_value == 0.5
    assert fd_parameter.bounds == (0, 1)
    assert fd_parameter.transform_func.__name__ == '<lambda>'


def test_FDParameter_v0i(make_poscar):
    """Test that FDParameter is initialized correctly for v0i."""
    fd_param_name = 'v0i'
    slab, rpars, *_ = make_poscar(SLAB_Cu2O_111)
    fd_parameter = _get_fd_parameter(fd_param_name, rpars)
    fd_parameter.apply(rpars, slab, 2.0)
    assert rpars.V0_IMAG == 2.0


def test_FDParameter_theta(make_poscar):
    """Test that FDParameter is initialized correctly for v0i."""
    fd_param_name = 'theta'
    slab, rpars, *_ = make_poscar(SLAB_Cu2O_111)
    fd_parameter = _get_fd_parameter(fd_param_name, rpars)
    fd_parameter.apply(rpars, slab, 5.0)
    assert rpars.THETA == 5.0


def test_FDParameter_phi(make_poscar):
    """Test that FDParameter is initialized correctly for v0i."""
    fd_param_name = 'phi'
    slab, rpars, *_ = make_poscar(SLAB_Cu2O_111)
    fd_parameter = _get_fd_parameter(fd_param_name, rpars)
    fd_parameter.apply(rpars, slab, 5.0)
    assert rpars.PHI == 5.0


@pytest.mark.parametrize('vector_elements',
                         (('a', ([0,0],)), ('b', ([1,1],)), ('c', ([2,2],)),
                          ('ab', ([0,0],[1,1])), ('bc', ([1,1],[2,2])),
                          ('abc', ([0,0],[1,1],[2,2]))))
def test_FDParameter_vector_scaling(vector_elements,make_poscar):
    """Test that FDParameter is initialized correctly for v0i."""
    fd_param_name, elements_to_check = vector_elements
    slab, rpars, *_ = make_poscar(SLAB_Cu2O_111)
    slab.bulkslab = slab.makeBulkSlab(rpars)
    orig_slab = deepcopy(slab)
    fd_parameter = _get_fd_parameter(fd_param_name, rpars)
    fd_parameter.apply(rpars, slab, 1.1)
    for element in elements_to_check:
        assert slab.ucell[element[0],element[1]] / orig_slab.ucell[element[0],element[1]] == pytest.approx(1.1)


def _get_fd_parameter(fd_param_name, rpars):
    """Get the FDParameter object for the given parameter name."""
    fd_parameter = FDParameter(
        fd_param_name,
        FD_PARAMETERS[fd_param_name]['x0'](rpars),
        FD_PARAMETERS[fd_param_name]['max_bounds'],
        FD_PARAMETERS[fd_param_name]['eval']
    )
    return fd_parameter
