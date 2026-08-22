"""Unit tests for R-factor implementations on synthetic data."""

import numpy as np
import pytest

from viperleed.calc.lib import rfactor


RFACTOR_FUNCS = (
    rfactor.r_pendry,
    rfactor.r_1,
    rfactor.r_2,
    rfactor.r_s,
    rfactor.r_zjj,
)


@pytest.mark.parametrize("func", RFACTOR_FUNCS)
def test_rfactor_requires_input_data(func, synthetic_pair):
    with pytest.raises(TypeError):
        func(
            synthetic_pair["v0i"],
            synthetic_pair["energy_step"],
            synthetic_pair["energies"],
        )


@pytest.mark.parametrize("func", RFACTOR_FUNCS)
def test_rfactor_spline_and_precomputed_paths_agree(func, synthetic_pair):
    energies = synthetic_pair["energies"]
    spline_1 = synthetic_pair["data_spline_1"]
    spline_2 = synthetic_pair["data_spline_2"]

    deriv_1_1 = spline_1.derivative()
    deriv_1_2 = deriv_1_1.derivative()
    deriv_2_1 = spline_2.derivative()
    deriv_2_2 = deriv_2_1.derivative()

    data_and_derivatives_1 = (
        spline_1(energies),
        deriv_1_1(energies),
        deriv_1_2(energies),
    )
    data_and_derivatives_2 = (
        spline_2(energies),
        deriv_2_1(energies),
        deriv_2_2(energies),
    )

    result_from_splines = func(
        synthetic_pair["v0i"],
        synthetic_pair["energy_step"],
        energies,
        data_spline_1=spline_1,
        data_spline_2=spline_2,
    )
    result_from_arrays = func(
        synthetic_pair["v0i"],
        synthetic_pair["energy_step"],
        energies,
        data_and_derivatives_1=data_and_derivatives_1,
        data_and_derivatives_2=data_and_derivatives_2,
    )

    assert result_from_splines == pytest.approx(result_from_arrays, rel=1e-12, abs=1e-14)


@pytest.mark.parametrize("func", RFACTOR_FUNCS)
def test_rfactor_per_beam_shape(func, synthetic_pair):
    spline_1 = synthetic_pair["data_spline_1"]
    spline_2 = synthetic_pair["data_spline_2"]

    result = func(
        synthetic_pair["v0i"],
        synthetic_pair["energy_step"],
        synthetic_pair["energies"],
        data_spline_1=spline_1,
        data_spline_2=spline_2,
        groups="beam",
    )

    result = np.asarray(result)
    assert result.ndim == 0 or result.shape == ()  # scalar spline case


def test_pendry_matches_existing_reference_value(sin_spline, cos_spline):
    energies = np.arange(15.0, 300.0, 0.5)
    calc_r = rfactor.r_pendry(
        3.0,
        0.5,
        energies,
        data_spline_1=sin_spline,
        data_spline_2=cos_spline,
    )
    assert calc_r == pytest.approx(1.904792, rel=1e-6)
