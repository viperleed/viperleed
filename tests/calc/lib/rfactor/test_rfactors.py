import numpy as np
import pytest
from pytest_cases import case, fixture, parametrize_with_cases

from .conftest import SplinesWithExpectedValues

from viperleed.calc.lib import rfactor


@parametrize_with_cases(
    'theo_spline, v0i, energy_step, energies, exp_spline, expected_R',
    cases=SplinesWithExpectedValues,
)
def test_pendry_R_from_splines(
    theo_spline, v0i, energy_step, energies, exp_spline, expected_R
):
    """Test the R-factor calculation from spline objects."""
    calc_R = rfactor.R_pendry(
        theo_spline, v0i, energy_step, energies, exp_spline
    )
    assert calc_R == pytest.approx(expected_R['pendry'], rel=1e-6)
