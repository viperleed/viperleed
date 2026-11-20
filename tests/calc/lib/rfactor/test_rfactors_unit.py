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


@pytest.mark.parametrize('func', [r_1.r_1, r_2.r_2])
def test_r1_r2_vanish_for_scaled_copy_in_wide_experiment(func):
    """Check R1 and R2 ignore experiment outside the theory range."""
    # Mirror of the test above: the theory covers only the middle of
    # the experimental range. On the overlap, theory is exactly 3 times
    # the experiment; outside it the experiment is arbitrary and must
    # not affect scale factor or R.
    exp = np.array([9.0, 1.0, 2.0, 3.0, 1.5, 20.0])[:, None]
    theo = np.array([NAN, 3.0, 6.0, 9.0, 4.5, NAN])[:, None]
    result = func(
        0.0, 1.0, np.arange(6.0),
        data_and_derivatives_1=(exp, None, None),
        data_and_derivatives_2=(theo, None, None),
    )
    assert result == pytest.approx(0.0, abs=1e-12)


def test_r2_normalizes_by_experiment_on_overlap_only():
    """Check R2's denominator excludes experiment outside the theory."""
    # Overlap is the three middle points, with weights 1/2, 1, 1/2.
    # I_e = [1, 2, 3] and I_t = [2, 2, 2] have the same integral (4),
    # so the scale factor is 1. The numerator is 0.5 + 0 + 0.5 = 1 and
    # the denominator is 0.5 + 4 + 4.5 = 9.
    exp = np.array([9.0, 1.0, 2.0, 3.0, 20.0])[:, None]
    theo = np.array([NAN, 2.0, 2.0, 2.0, NAN])[:, None]
    result = r_2.r_2(
        0.0, 1.0, np.arange(5.0),
        data_and_derivatives_1=(exp, None, None),
        data_and_derivatives_2=(theo, None, None),
    )
    assert result == pytest.approx(1 / 9, rel=1e-12)


def test_pendry_from_y_ignores_nan_in_second_dataset():
    """Check that points where only y_2 is NaN leave the overlap."""
    # Overlap is the three middle points; trapezoid weights 1/2, 1, 1/2.
    y_1 = np.array([1.0, 2.0, 3.0, 4.0, 5.0])[:, None]
    y_2 = np.array([NAN, 0.5, 1.0, 1.0, NAN])[:, None]
    numerator = 0.5 * 1.5**2 + 2.0**2 + 0.5 * 3.0**2
    denominator = 0.5 * 4.25 + 10.0 + 0.5 * 17.0
    result = pendry.r_pendry_from_y(y_1, y_2, 1.0)
    assert result == pytest.approx(numerator / denominator, rel=1e-12)


def test_zj_ignores_experiment_outside_theory_range():
    """Check ZJ uses only the overlap when only the theory has NaN."""
    # Same data as the edge-weight test, but the experiment is finite
    # (and large) outside the overlap. Neither the scale factor, nor
    # epsilon (largest |I'_e| = 4), nor the integrals may change.
    exp = (
        np.array([9, 2, 2, 2, 9])[:, None],
        np.array([7, 1, 4, 2, 7])[:, None],
        np.array([7, 1, 1, 1, 7])[:, None],
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


def test_zj_epsilon_uses_largest_derivative_magnitude():
    """Check that ZJ's epsilon is the largest experimental |dI/dE|."""
    # Same data as the edge-weight test, but the experimental derivative
    # of largest magnitude (-4) is negative. epsilon is its modulus, so
    # the denominators are [5, 8, 6].
    exp = (
        np.array([NAN, 2, 2, 2, NAN])[:, None],
        np.array([NAN, 1, -4, 2, NAN])[:, None],
        np.array([NAN, 1, 1, 1, NAN])[:, None],
    )
    theo = (
        np.array([NAN, 1, 1, 1, NAN])[:, None],
        np.array([NAN, 0.5, -2, 0, NAN])[:, None],
        np.array([NAN, 0.5, 1, 2, NAN])[:, None],
    )
    result = zannazi_jona.r_zjj(
        0.0, 1.0, np.arange(5.0),
        data_and_derivatives_1=exp,
        data_and_derivatives_2=theo,
    )
    assert result == pytest.approx(0.5 / (4 * 0.027), rel=1e-12)


STEP = 0.5
# Three beams on a 6-point grid, overlapping on 6, 3 and 5 points.
OVERLAPS = ((0, 6), (1, 4), (0, 5))
EXP_I = ((1.0, 2.0, 3.0, 2.0, 1.0, 2.0),
         (9.9, 1.0, 3.0, 2.0, 9.9, 9.9),
         (2.0, 1.0, 1.5, 3.0, 2.5, 9.9))
EXP_D1 = ((1.0, 1.0, -1.0, -1.0, 0.5, 2.0),
          (9.9, 2.0, 0.0, -2.0, 9.9, 9.9),
          (0.5, 1.0, 1.5, -1.0, -3.0, 9.9))
EXP_D2 = ((0.5, -1.0, 1.0, 2.0, -2.0, 1.0),
          (9.9, 1.0, -1.0, 0.5, 9.9, 9.9),
          (1.0, 0.5, -1.0, 1.0, 2.0, 9.9))
THEO_I = ((1.5, 1.0, 2.0, 3.0, 1.0, 1.0),
          (9.9, 2.0, 2.0, 1.0, 9.9, 9.9),
          (1.0, 2.0, 1.0, 2.0, 3.0, 9.9))
THEO_D1 = ((0.5, 0.5, 1.0, 0.0, -1.0, 1.0),
           (9.9, 1.0, 0.5, -1.0, 9.9, 9.9),
           (1.0, 0.5, -0.5, 1.5, -1.0, 9.9))
THEO_D2 = ((1.0, 0.5, -1.0, 1.0, -0.5, 2.0),
           (9.9, 0.5, 0.5, -1.0, 9.9, 9.9),
           (-1.0, 1.0, 0.5, 0.5, 1.0, 9.9))
GROUPS = [0, 0, 1]


def _integral(values):
    """Return the trapezoidal integral of values on the test grid."""
    values = np.asarray(values)
    return STEP * np.sum(values[1:] + values[:-1]) / 2


def _beam_slices(beam):
    """Return intensities and derivatives of beam within its overlap."""
    first, last = OVERLAPS[beam]
    return [np.array(arr[beam][first:last], dtype=float)
            for arr in (EXP_I, EXP_D1, EXP_D2, THEO_I, THEO_D1, THEO_D2)]


def _r2_of_beam(beam):
    """Return the R2 of a single beam from first principles."""
    exp, _, _, theo, _, _ = _beam_slices(beam)
    scale = _integral(exp) / _integral(theo)
    return _integral((exp - scale*theo)**2) / _integral(exp**2)


def _zj_of_beam(beam):
    """Return the ZJ R factor of a single beam from first principles."""
    exp, exp_1, exp_2, theo, theo_1, theo_2 = _beam_slices(beam)
    scale = _integral(exp) / _integral(theo)
    quotient = (abs(scale*theo_2 - exp_2) * abs(scale*theo_1 - exp_1)
                / (abs(exp_1) + np.max(abs(exp_1))))
    return _integral(quotient) / _integral(exp) / 0.027


def _overlap_length(beam):
    """Return the energy length of the overlap of beam."""
    first, last = OVERLAPS[beam]
    return (last - first - 1) * STEP


def _weighted_mean(values, beams):
    """Return the mean of values weighted by overlap length of beams."""
    lengths = [_overlap_length(b) for b in beams]
    return sum(v*w for v, w in zip(values, lengths)) / sum(lengths)


def _nan_outside_overlap(rows):
    """Return rows as a (6, 3) array, NaN outside the overlaps."""
    data = np.full((6, 3), np.nan)
    for beam, (first, last) in enumerate(OVERLAPS):
        data[first:last, beam] = rows[beam][first:last]
    return data


def _call_with_ragged_beams(func, **kwargs):
    """Return func evaluated on the ragged test beams."""
    exp = tuple(_nan_outside_overlap(rows)
                for rows in (EXP_I, EXP_D1, EXP_D2))
    theo = tuple(_nan_outside_overlap(rows)
                 for rows in (THEO_I, THEO_D1, THEO_D2))
    return func(
        0.0, STEP, np.arange(6) * STEP,
        data_and_derivatives_1=exp,
        data_and_derivatives_2=theo,
        **kwargs,
    )


@pytest.mark.parametrize(
    ("func", "per_beam"),
    [(r_2.r_2, _r2_of_beam), (zannazi_jona.r_zjj, _zj_of_beam)],
    ids=["r2", "zj"],
)
class TestBeamAveraging:
    """Tests for averaging of R2 and ZJ over beams."""
    def test_per_beam_values(self, func, per_beam):
        """Check per-beam values."""
        result = _call_with_ragged_beams(func, groups="beam")
        expected = [per_beam(b) for b in range(3)]
        assert result == pytest.approx(expected, rel=1e-12)

    def test_overall_is_overlap_length_weighted_mean(self, func, per_beam):
        """Check the overall value is weighted by overlap length."""
        result = _call_with_ragged_beams(func)
        expected = _weighted_mean([per_beam(b) for b in range(3)],
                                  range(3))
        assert result == pytest.approx(expected, rel=1e-12)

    def test_groups_are_overlap_length_weighted_means(self, func, per_beam):
        """Check group values are means weighted by overlap length."""
        result = _call_with_ragged_beams(func, groups=GROUPS)
        expected = [
            _weighted_mean([per_beam(0), per_beam(1)], (0, 1)),
            per_beam(2),
        ]
        assert result == pytest.approx(expected, rel=1e-12)

    def test_beam_without_overlap_does_not_contribute(self, func, per_beam):
        """Check that a beam without overlap is left out."""
        exp = tuple(_nan_outside_overlap(rows)
                    for rows in (EXP_I, EXP_D1, EXP_D2))
        theo = tuple(_nan_outside_overlap(rows)
                     for rows in (THEO_I, THEO_D1, THEO_D2))
        for arr in (*exp, *theo):
            arr[:, 1] = np.nan
        result = func(
            0.0, STEP, np.arange(6) * STEP,
            data_and_derivatives_1=exp,
            data_and_derivatives_2=theo,
        )
        expected = _weighted_mean([per_beam(0), per_beam(2)], (0, 2))
        assert result == pytest.approx(expected, rel=1e-12)

    @pytest.mark.parametrize('groups', [None, 'beam'])
    def test_beam_with_zero_experimental_intensity_is_excluded(
        self, func, per_beam, groups
    ):
        """Check that a beam with undefined R (zero I_e) is left out."""
        exp, theo = (
            tuple(
                np.column_stack([_nan_outside_overlap(rows)[:, 0]] * 2)
                for rows in beams
            )
            for beams in ((EXP_I, EXP_D1, EXP_D2), (THEO_I, THEO_D1, THEO_D2))
        )
        exp[0][:, 1] = 0.0  # I_e = 0 on the whole overlap of beam 1
        result = func(
            0.0, STEP, np.arange(6) * STEP,
            data_and_derivatives_1=exp,
            data_and_derivatives_2=theo,
            groups=groups,
        )
        if groups is None:
            assert result == pytest.approx(per_beam(0), rel=1e-12)
        else:
            assert result[0] == pytest.approx(per_beam(0), rel=1e-12)
            assert np.isnan(result[1])
