"""Unit tests for R-factor implementations on synthetic data."""

import numpy as np
import pytest

from viperleed.calc.lib.rfactor import pendry
from viperleed.calc.lib.rfactor import r_1
from viperleed.calc.lib.rfactor import r_2
from viperleed.calc.lib.rfactor import smooth
from viperleed.calc.lib.rfactor import zannazi_jona


RFACTOR_FUNCS = (
    pendry.r_pendry,
    r_1.r_1,
    r_2.r_2,
    smooth.r_s,
    zannazi_jona.r_zjj,
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
    calc_r = pendry.r_pendry(
        3.0,
        0.5,
        energies,
        data_spline_1=sin_spline,
        data_spline_2=cos_spline,
    )
    assert calc_r == pytest.approx(1.904792, rel=1e-6)


NAN = np.nan


def test_pendry_from_y_weighs_overlap_edges_by_half():
    """Check trapezoid weights of Pendry's R at the overlap edges."""
    # Overlap is the three middle points; trapezoid weights 1/2, 1, 1/2.
    y_1 = np.array([NAN, 1.0, 2.0, 3.0, NAN])[:, None]
    y_2 = np.array([0.5, 0.5, 1.0, 1.0, 1.0])[:, None]
    numerator = 0.5 * 0.25 + 1.0 + 0.5 * 4.0
    denominator = 0.5 * 1.25 + 5.0 + 0.5 * 10.0
    result = pendry.r_pendry_from_y(y_1, y_2, 1.0)
    assert result == pytest.approx(numerator / denominator, rel=1e-12)


def test_zj_weighs_overlap_edges_by_half():
    """Check trapezoid weights of ZJ's R at the overlap edges."""
    # Overlap is the three middle points. With I_e = 2, I_t = 1 there,
    # the scale factor is 2. Weighted quotients are [0, 0, 1], so the
    # integral is 1/2; the largest |I'_e| is 4, and int I_e = 4.
    exp = (
        np.array([NAN, 2, 2, 2, NAN])[:, None],
        np.array([NAN, 1, 4, 2, NAN])[:, None],
        np.array([NAN, 1, 1, 1, NAN])[:, None],
    )
    theo = (
        np.array([NAN, 1, 1, 1, NAN])[:, None],
        np.array([NAN, 0.5, 2, 0, NAN])[:, None],
        np.array([NAN, 0.5, 1, 2, NAN])[:, None],
    )
    result = zannazi_jona.r_zjj(
        0.0, 1.0, np.arange(5.0),
        data_and_derivatives_1=exp,
        data_and_derivatives_2=theo,
    )
    assert result == pytest.approx(0.5 / (4 * 0.027), rel=1e-12)


@pytest.mark.parametrize("func", [r_1.r_1, r_2.r_2])
def test_r1_r2_vanish_for_scaled_copy_on_overlap(func):
    """Check R1 and R2 ignore theory outside the experimental range."""
    # The experiment covers only the middle of the theory range. On the
    # overlap, theory is exactly 3 times the experiment; outside it the
    # theory is arbitrary and must not affect scale factor or R.
    exp = np.array([NAN, 1.0, 2.0, 3.0, 1.5, NAN])[:, None]
    theo = np.array([9.0, 3.0, 6.0, 9.0, 4.5, 20.0])[:, None]
    result = func(
        0.0, 1.0, np.arange(6.0),
        data_and_derivatives_1=(exp, None, None),
        data_and_derivatives_2=(theo, None, None),
    )
    assert result == pytest.approx(0.0, abs=1e-12)
