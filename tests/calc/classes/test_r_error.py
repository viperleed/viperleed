"""Tests for viperleed.calc.classes.RError."""

__authors__ = (
    'Alexandra Mia Imre (@alexmiame)',
    'Michele Riva (@michele-riva)',
)
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2023-03-09'
__license__ = 'GPLv3+'

from types import SimpleNamespace

import numpy as np
import pytest

from viperleed.calc.classes import r_error
from viperleed.calc.classes.rparams.special.r_factor_type import RFactorType
from viperleed.calc.files import ioerrorcalc


class TestZeroCrossings:
    """Tests for methods detecting zeros."""

    _x_only = (  # x_values, nr_expected_crossings
        (np.array([-1, 1, -1]),       2),
        (np.array([0, 0, 0]),         1),
        ([-1, 0, -1],                 1),
        ([-1, 1, 10, 2, -50, 0.1, 0], 3),
        )

    @pytest.mark.parametrize('x_vals,n_crossings', _x_only)
    def test_nr_zero_crossings(self, x_vals, n_crossings):
        """Check that the number of detected crossings is correct."""
        assert r_error.get_n_zero_crossings(x_vals) == n_crossings

    _x_and_y = (  # (x_values, y_values), cross_x_position
        (([0, 1, 2],       [-1, 0, 1]),         1),
        (([0, 1, 2, 3, 4], [-1, -3, -2, 2, 5]), 2.5),
        (([-2, -1, 0, 1, 2], [-5, -4, -3, 1, 3]), 0.75),
        (([-2, -1, 0, 1, 2], [5, 4, 3, -1, -3]), 0.75),
        )

    @pytest.mark.parametrize('xy_vals,crossing_x', _x_and_y)
    def test_pos_zero_crossing(self, xy_vals, crossing_x):
        """Check that the x position of the detected crossing is correct."""
        assert r_error.get_zero_crossing(*xy_vals) == crossing_x


_V0I = 4.5
_E_RANGE = 500.0
_DISPLACEMENTS = [-0.2, -0.1, 0.0, 0.1, 0.2]
_RFACS = [0.5, 0.3, 0.2, 0.25, 0.45]
_VAR_R = np.sqrt(8 * _V0I / _E_RANGE) * min(_RFACS)
_BOUNDS = (  # Zero crossings of _RFACS - (R_min + var(R))
    0.1 * _VAR_R / (0.3 - 0.2),
    0.1 + 0.1 * (0.2 + _VAR_R - 0.25) / (0.45 - 0.25),
    )


def _make_geo_error(r_type):
    """Return an R_Error for geometric displacements of one atom."""
    atom = SimpleNamespace(num=1,
                           el='O',
                           disp_geo={'all': _DISPLACEMENTS})
    return r_error.R_Error(r_type,
                           [atom],
                           'geo',
                           _RFACS,
                           'z',
                           _DISPLACEMENTS,
                           v0i=_V0I,
                           energy_range=_E_RANGE)


@pytest.mark.parametrize('r_type', [RFactorType.PENDRY, RFactorType.SMOOTH])
def test_error_estimates_with_var_r(r_type):
    """Check var(R) and error bounds of the Pendry and smooth R."""
    error = _make_geo_error(r_type)
    assert error.var_r == pytest.approx(_VAR_R)
    assert error.get_error_estimates == pytest.approx(_BOUNDS)
    var_r_info = ioerrorcalc.extract_var_r([error])
    assert var_r_info == {'geo': error.var_r, 'vib': None, 'occ': None}


@pytest.mark.parametrize('r_type', [RFactorType.R2, RFactorType.ZJ])
def test_no_error_estimates_without_var_r(r_type):
    """Check that R2 and ZJ give no var(R) and no error bounds."""
    error = _make_geo_error(r_type)
    assert error.var_r is None
    assert error.get_error_estimates == (None, None)
    assert not any(ioerrorcalc.extract_var_r([error]).values())
    with pytest.raises(ValueError, match='Pendry or the smooth'):
        error.calc_var_r(_V0I, _E_RANGE)


_NON_FLOAT_ARGS = {  # v0i, energy_range
    'int v0i': (4, _E_RANGE),
    'int energy_range': (_V0I, 500),
    'str energy_range': (_V0I, '500'),
    }


@pytest.mark.parametrize(('v0i', 'energy_range'),
                         _NON_FLOAT_ARGS.values(),
                         ids=_NON_FLOAT_ARGS)
def test_var_r_rejects_non_float(v0i, energy_range):
    """Check that var(R) refuses a v0i or energy_range not float."""
    error = _make_geo_error(RFactorType.PENDRY)
    with pytest.raises(ValueError, match='must be float'):
        error.calc_var_r(v0i, energy_range)
