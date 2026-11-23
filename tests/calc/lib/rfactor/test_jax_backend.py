"""Tests for the R-factor functions running with JAX."""

__authors__ = ('Florian Kraushofer (@fkraushofer)',)
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-10-08'
__license__ = 'GPLv3+'

from functools import partial

import numpy as np
import pytest

from viperleed.calc.lib import dynamic_numerical_lib as dnl
from viperleed.calc.lib.rfactor import pendry
from viperleed.calc.lib.rfactor import r_1
from viperleed.calc.lib.rfactor import r_2
from viperleed.calc.lib.rfactor import smooth
from viperleed.calc.lib.rfactor import zannazi_jona
from viperleed.calc.lib.rfactor.groups import group_rfactors

jax = pytest.importorskip('jax')
interpax = pytest.importorskip('interpax')

# Attributes of dnl as patched by viperleed-jax (tensor_calculator.py)
_JAX_BACKEND = {
    'xp': jax.numpy,
    'stop_gradient': jax.lax.stop_gradient,
    'bincount': jax.numpy.bincount,
    'vmap': jax.vmap,
    'segment_sum': jax.ops.segment_sum,
    'CubicSpline': partial(interpax.CubicSpline, check=False),
    'PPoly': interpax.PPoly,
    }
_NUM_GROUPS = 3
# interpax and scipy splines and the JAX and NumPy R-factor code paths
# agree to within a few ULPs (largest relative difference measured on
# the Fe2O3(012) data: 2.8e-15, for per-group R2).
_RTOL = 1e-12


# A fixture, as test_jit_matches_numpy must request it only after
# computing its NumPy reference, and x64 must be reset on teardown.
@pytest.fixture(name='jax_backend')
def fixture_jax_backend(monkeypatch):
    """Switch dnl to JAX as viperleed-jax does, with float64 on."""
    x64_enabled = jax.config.read('jax_enable_x64')
    jax.config.update('jax_enable_x64', val=True)
    for attr, value in _JAX_BACKEND.items():
        monkeypatch.setattr(dnl, attr, value)
    yield
    jax.config.update('jax_enable_x64', x64_enabled)


@pytest.fixture(name='theo_on_full_grid', scope='module')
def fixture_theo_on_full_grid(fe2o3_012_1x1_dataset):
    """Return Fe2O3(012) data restricted to where all THEO beams exist.

    A spline through intensities that contain NaN is NaN everywhere
    with interpax. Like the tensor-LEED intensities of viperleed-jax,
    the THEO beams are thus restricted to energies without NaN.
    """
    theo_grid = fe2o3_012_1x1_dataset['theo_grid']
    theo_beams = fe2o3_012_1x1_dataset['theo_beams']
    theo_rows = ~np.isnan(theo_beams).any(axis=1)
    theo_grid = theo_grid[theo_rows]
    energies = fe2o3_012_1x1_dataset['energies']
    in_range = energies >= theo_grid[0]
    return {
        'energy_step': fe2o3_012_1x1_dataset['energy_step'],
        'energies': energies[in_range],
        'exp_triplet': tuple(
            arr[in_range] for arr in fe2o3_012_1x1_dataset['exp_triplet']
            ),
        'theo_grid': theo_grid,
        'theo_beams': theo_beams[theo_rows],
        'v0i': fe2o3_012_1x1_dataset['v0i'],
        }


def _r_factor_from_theo_spline(rfactor_func,
                               data,
                               groups=None,
                               num_groups=None):
    """Return R as calculated by viperleed-jax's calc_r_factor."""
    theo_spline = dnl.CubicSpline(
        data['theo_grid'],
        data['theo_beams'],
        extrapolate=False,
        )
    return rfactor_func(
        data['v0i'],
        data['energy_step'],
        data['energies'],
        data_and_derivatives_1=data['exp_triplet'],
        data_spline_2=theo_spline,
        groups=groups,
        num_groups=num_groups,
        shift_2nd_spline=0.0,
        )


_jit_r_factor_from_theo_spline = jax.jit(
    _r_factor_from_theo_spline,
    static_argnames=('rfactor_func', 'groups', 'num_groups'),
    )


@pytest.mark.parametrize('per_group', [False, True], ids=('overall', 'group'))
@pytest.mark.parametrize('rfactor_func', [
    pendry.r_pendry,
    r_2.r_2,
    smooth.r_s,
    zannazi_jona.r_zjj,
    r_1.r_1,
    ])
def test_jit_matches_numpy(rfactor_func,
                           per_group,
                           theo_on_full_grid,
                           request):
    """Check that R under jax.jit with interpax matches NumPy."""
    group_kwargs = {}
    if per_group:
        n_beams = theo_on_full_grid['theo_beams'].shape[1]
        group_kwargs = {
            'groups': tuple(i % _NUM_GROUPS for i in range(n_beams)),
            'num_groups': _NUM_GROUPS,
            }
    expected = _r_factor_from_theo_spline(rfactor_func,
                                          theo_on_full_grid,
                                          **group_kwargs)
    assert np.all(np.isfinite(expected))
    # Switch only now: the reference above must use NumPy and SciPy
    request.getfixturevalue('jax_backend')
    jax_data = jax.tree_util.tree_map(jax.numpy.asarray, theo_on_full_grid)
    result = _jit_r_factor_from_theo_spline(rfactor_func,
                                            jax_data,
                                            **group_kwargs)
    assert isinstance(result, jax.Array)
    np.testing.assert_allclose(np.asarray(result), expected, rtol=_RTOL)


@pytest.mark.usefixtures('jax_backend')
@pytest.mark.parametrize(('group_kwargs', 'expected'), [
    ({'groups': 'beam'}, (0.5, 0.0, 1.0)),
    ({'groups': (0, 1, 0), 'num_groups': 2}, (1/3, 0.0, 1/3)),
    ], ids=('beam', 'group'))
def test_group_rfactors_gradient_zero_denominator(group_kwargs, expected):
    """Check a zero denominator gives zero, not NaN, gradients."""
    denominators = jax.numpy.asarray((2.0, 0.0, 1.0))

    def _nansum_rfactors(numerators):
        rfactors = group_rfactors(numerators,
                                  denominators,
                                  **group_kwargs)
        return jax.numpy.nansum(rfactors)

    gradient = jax.grad(_nansum_rfactors)(jax.numpy.asarray((1.0, 2.0, 0.5)))
    np.testing.assert_allclose(np.asarray(gradient), expected)
