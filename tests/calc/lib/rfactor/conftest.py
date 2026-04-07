import numpy as np
from pytest_cases import fixture
from scipy.interpolate import CubicSpline

from viperleed.calc.files.beams import readOUTBEAMS
from viperleed.calc.files.iorfactor import beamlist_to_array
from viperleed.calc.lib.spline_interpolation import interpolate_ragged_array

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
    return CubicSpline(_MOCK_ROUGH_E_AXIS, _MOCK_SIN_INTENSITY)


@fixture
def cos_spline():
    return CubicSpline(_MOCK_ROUGH_E_AXIS, _MOCK_COS_INTENSITY)


@fixture
def const_spline():
    return CubicSpline(_MOCK_ROUGH_E_AXIS, _MOCK_CONST_INTENSITY)


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


class SplinesWithExpectedValues:
    def case_sin_correlated(self, sin_spline):
        theo_spline, exp_spline = sin_spline, sin_spline
        v0i, energy_step, energies = 3.0, 0.5, _MOCK_FINE_E_AXIS
        expected_R = {'pendry': 0.005350277}
        return theo_spline, v0i, energy_step, energies, exp_spline, expected_R

    def case_cos_correlated(self, cos_spline):
        theo_spline, exp_spline = cos_spline, cos_spline
        v0i, energy_step, energies = 3.0, 0.5, _MOCK_FINE_E_AXIS
        expected_R = {'pendry': 0.0}
        return theo_spline, v0i, energy_step, energies, exp_spline, expected_R

    def case_sin_anticorrelated(self, sin_spline, cos_spline):
        theo_spline, exp_spline = sin_spline, cos_spline
        v0i, energy_step, energies = 3.0, 0.5, _MOCK_FINE_E_AXIS
        expected_R = {'pendry': 1.9047921884401726}
        return theo_spline, v0i, energy_step, energies, exp_spline, expected_R
