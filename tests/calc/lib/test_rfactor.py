"""Test Module R-factor"""

import interpax
import numpy as np
import pytest
from jax import numpy as jnp
from pytest_cases import case, fixture, parametrize_with_cases
from viperleed.calc.files.beams import readOUTBEAMS
from viperleed.calc.files.iorfactor import beamlist_to_array

from viperleed_jax import rfactor
from viperleed_jax.interpolation import interpolate_ragged_array

_MOCK_ROUGH_E_AXIS = np.arange(15, 300, 3)
_MOCK_FINE_E_STEP = 0.5
_MOCK_FINE_E_AXIS = np.arange(15, 300, _MOCK_FINE_E_STEP)


_MOCK_SIN_INTENSITY = np.sin(_MOCK_ROUGH_E_AXIS / (3 * np.pi)) ** 2
_MOCK_COS_INTENSITY = np.cos(_MOCK_ROUGH_E_AXIS / (3 * np.pi)) ** 2
_MOCK_CONST_INTENSITY = np.ones_like(_MOCK_ROUGH_E_AXIS)

_FE2O3_0123_EXPECTED_VALUES = [
    (-1.0, 0.17540223),
    (-0.5, 0.1573342),
    (0.0, 0.16663656),
    (0.5, 0.19706298),
    (1.0, 0.24195227),
]


@fixture
def sin_spline():
    return interpax.CubicSpline(_MOCK_ROUGH_E_AXIS, _MOCK_SIN_INTENSITY)


@fixture
def cos_spline():
    return interpax.CubicSpline(_MOCK_ROUGH_E_AXIS, _MOCK_COS_INTENSITY)


@fixture
def const_spline():
    return interpax.CubicSpline(_MOCK_ROUGH_E_AXIS, _MOCK_CONST_INTENSITY)


# TODO: replace paths with info object when available
@fixture(scope='session')
def fe2o3_012_exp_spline_and_hk():
    exp_out_beams = readOUTBEAMS(
        '/Users/alexander/GitHub/on-the-fly-deltas/tests/test_data/Fe2O3_012/converged/EXPBEAMS.csv'
    )
    exp_energies, _, _, exp_beams = beamlist_to_array(exp_out_beams)
    exp_spline = interpolate_ragged_array(exp_energies, exp_beams)
    exp_hk = [b.hk for b in exp_out_beams]
    return exp_spline, exp_hk


@fixture(scope='session')
@pytest.mark.parametrize('v0r_shift,expected_R', _FE2O3_0123_EXPECTED_VALUES)
def fe2o3_012_theo_spline_exp_spline_and_R(
    v0r_shift, expected_R, fe2o3_012_exp_spline_and_hk
):
    exp_spline, exp_hk = fe2o3_012_exp_spline_and_hk
    theo_out_beams = readOUTBEAMS(
        '/Users/alexander/GitHub/on-the-fly-deltas/tests/test_data/Fe2O3_012/converged/OUT/THEOBEAMS.csv'
    )
    theo_hk = [b.hk for b in theo_out_beams]
    theo_out_dict = {b.hk: b for b in theo_out_beams}
    sorted_theo = [theo_out_dict[hk] for hk in exp_hk]
    theo_energies, _, _, theo_beams = beamlist_to_array(sorted_theo)
    theo_spline = interpolate_ragged_array(
        theo_energies + v0r_shift, theo_beams
    )

    return exp_spline, theo_spline, theo_energies, expected_R


class SplinesWithExpectedValues:
    @case(tags='pendry')
    def case_sin_correlated(self, sin_spline):
        theo_spline, exp_spline = sin_spline, sin_spline
        v0i, energy_step, energies = 3.0, 0.5, _MOCK_FINE_E_AXIS
        expected_R = 0.0
        return theo_spline, v0i, energy_step, energies, exp_spline, expected_R

    @case(tags='pendry')
    def case_cos_correlated(self, cos_spline):
        theo_spline, exp_spline = cos_spline, cos_spline
        v0i, energy_step, energies = 3.0, 0.5, _MOCK_FINE_E_AXIS
        expected_R = 0.0
        return theo_spline, v0i, energy_step, energies, exp_spline, expected_R

    @case(tags='pendry')
    def case_sin_anticorrelated(self, sin_spline, cos_spline):
        theo_spline, exp_spline = sin_spline, cos_spline
        v0i, energy_step, energies = 3.0, 0.5, _MOCK_FINE_E_AXIS
        expected_R = 1.9047921
        return theo_spline, v0i, energy_step, energies, exp_spline, expected_R

    @case(tags='pendry')
    def case_fe2o3_012_refcalc(self, fe2o3_012_theo_spline_exp_spline_and_R):
        exp_spline, theo_spline, energies, expected_R = (
            fe2o3_012_theo_spline_exp_spline_and_R
        )
        v0i, energy_step = 5.0, 0.5
        return theo_spline, v0i, energy_step, energies, exp_spline, expected_R


@parametrize_with_cases(
    'theo_spline, v0i, energy_step, energies, exp_spline, expected_R',
    cases=SplinesWithExpectedValues,
    has_tag='pendry',
)
def test_pendry_R_from_splines(
    theo_spline, v0i, energy_step, energies, exp_spline, expected_R
):
    """Test the R-factor calculation from spline objects."""
    calc_R = rfactor.pendry_R(
        theo_spline, v0i, energy_step, energies, exp_spline
    )
    assert calc_R == pytest.approx(expected_R, rel=1e-6)
