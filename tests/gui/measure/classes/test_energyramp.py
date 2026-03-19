"""Tests for module energyramp of viperleed.gui.measure.classes."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    'Florian Dörr (@FlorianDoerr)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-03-17'
__license__ = 'GPLv3+'

from configparser import ConfigParser

import pytest
from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.gui.measure.classes.energyramp import ABRUPT
# from viperleed.gui.measure.classes.energyramp import DEFAULT_DELTA            # TODO: add tests for default values
# from viperleed.gui.measure.classes.energyramp import DEFAULT_END
# from viperleed.gui.measure.classes.energyramp import DEFAULT_START
from viperleed.gui.measure.classes.energyramp import LINEAR
from viperleed.gui.measure.classes.energyramp import ConstantEnergyRamp
from viperleed.gui.measure.classes.energyramp import EndlessLinearEnergyRamp
from viperleed.gui.measure.classes.energyramp import LinearEnergyRamp


@fixture(name='linear_ramp')
def fixture_linear_ramp():
    """Return a fresh LinearEnergyRamp instance."""
    return LinearEnergyRamp()


@fixture(name='endless_ramp')
def fixture_endless_ramp():
    """Return a fresh EndlessLinearEnergyRamp instance."""
    return EndlessLinearEnergyRamp()


@fixture(name='constant_ramp')
def fixture_constant_ramp():
    """Return a fresh ConstantEnergyRamp instance."""
    return ConstantEnergyRamp()


def _make_settings(section_dict):
    """Return a ConfigParser populated from a dictionary."""
    cfg = ConfigParser()
    cfg.read_dict(section_dict)
    return cfg


class TestReturnMatchingEnergyRamp:
    """Tests for EnergyRampABC.get_matching_energy_ramp."""

    def test_returns_linear_ramp_by_default(self):
        """Check that an empty settings object yields LinearEnergyRamp."""
        settings = _make_settings({})
        result = LinearEnergyRamp.get_matching_energy_ramp(settings)
        assert result is LinearEnergyRamp

    def test_returns_linear_ramp_when_constant_energy_and_endless_false(self):
        """Check that constant_energy/endless=false yield LinearEnergyRamp."""
        settings = _make_settings(
            {'energies': {'constant_energy': 'false', 'endless': 'false'}}
            )
        result = LinearEnergyRamp.get_matching_energy_ramp(settings)
        assert result is LinearEnergyRamp

    def test_returns_constant_ramp_when_constant_energy_true(self):
        """Check that constant_energy=true yields ConstantEnergyRamp."""
        settings = _make_settings({'energies': {'constant_energy': 'true'}})
        result = LinearEnergyRamp.get_matching_energy_ramp(settings)
        assert result is ConstantEnergyRamp

    def test_returns_endless_ramp_when_endless_true(self):
        """Check that endless=true yields EndlessLinearEnergyRamp."""
        settings = _make_settings({'energies': {'endless': 'true'}})
        result = LinearEnergyRamp.get_matching_energy_ramp(settings)
        assert result is EndlessLinearEnergyRamp

    def test_constant_energy_takes_priority_over_endless(self):
        """Check that constant_energy=true takes precedence over endless."""
        settings = _make_settings(
            {'energies': {'constant_energy': 'true', 'endless': 'true'}}
            )
        result = LinearEnergyRamp.get_matching_energy_ramp(settings)
        assert result is ConstantEnergyRamp

    _invalid = {
        'constant_energy not a bool': {'constant_energy': 'maybe'},
        'endless not a bool': {'endless': 'perhaps'},
    }

    @parametrize('section', _invalid.values(), ids=_invalid)
    def test_returns_linear_ramp_for_invalid_values(self, section):
        """Check that invalid boolean values fall back to LinearEnergyRamp."""
        settings = _make_settings({'energies': section})
        result = LinearEnergyRamp.get_matching_energy_ramp(settings)
        assert result is LinearEnergyRamp

# pylint: disable=protected-access
class TestLinearEnergyRampFinished:
    """Tests for LinearEnergyRamp.ramp_finished."""

    def test_not_finished_positive_delta_within_range(self, linear_ramp):
        """Check ramp not finished when next step is within range."""
        linear_ramp._delta_energy = 1.0
        linear_ramp._end_energy = 20.0
        linear_ramp.current_energy = 19.0
        assert not linear_ramp.ramp_finished()

    def test_finished_positive_delta_at_end(self, linear_ramp):
        """Check ramp finished when next step would exceed end energy."""
        linear_ramp._delta_energy = 1.0
        linear_ramp._end_energy = 20.0
        linear_ramp.current_energy = 20.0
        assert linear_ramp.ramp_finished()

    def test_finished_positive_delta_beyond_end(self, linear_ramp):
        """Check ramp finished when current energy is already past end."""
        linear_ramp._delta_energy = 1.0
        linear_ramp._end_energy = 20.0
        linear_ramp.current_energy = 25.0
        assert linear_ramp.ramp_finished()

    def test_not_finished_negative_delta_within_range(self, linear_ramp):
        """Check ramp not finished when next step is within range."""
        linear_ramp._delta_energy = -1.0
        linear_ramp._end_energy = 10.0
        linear_ramp.current_energy = 11.0
        assert not linear_ramp.ramp_finished()

    def test_finished_negative_delta_at_end(self, linear_ramp):
        """Check ramp finished when next step would go below end energy."""
        linear_ramp._delta_energy = -1.0
        linear_ramp._end_energy = 10.0
        linear_ramp.current_energy = 10.0
        assert linear_ramp.ramp_finished()

    def test_finished_negative_delta_beyond_end(self, linear_ramp):
        """Check ramp finished when current energy is already below end."""
        linear_ramp._delta_energy = -1.0
        linear_ramp._end_energy = 10.0
        linear_ramp.current_energy = 5.0
        assert linear_ramp.ramp_finished()

    def test_energy_steps_positive_delta(self, linear_ramp):
        """Check energy_steps computed correctly for positive delta."""
        linear_ramp._start_energy = 10.0
        linear_ramp._end_energy = 20.0
        linear_ramp._delta_energy = 2.0
        assert linear_ramp.energy_steps == 6  # 10, 12, 14, 16, 18, 20

    def test_energy_steps_negative_delta(self, linear_ramp):
        """Check energy_steps computed correctly for negative delta."""
        linear_ramp._start_energy = 20.0
        linear_ramp._end_energy = 10.0
        linear_ramp._delta_energy = -2.0
        assert linear_ramp.energy_steps == 6  # 20, 18, 16, 14, 12, 10

    def test_energy_steps_zero_delta(self, linear_ramp):
        """Check energy_steps is 0 for 0.0 delta."""
        linear_ramp._start_energy = 10.0
        linear_ramp._end_energy = 20.0
        linear_ramp._delta_energy = 0.0
        assert linear_ramp.energy_steps == 0


class TestEndlessLinearEnergyRampIncrementEnergy:
    """Tests for EndlessLinearEnergyRamp.increment_energy."""

    def test_wraps_to_start_on_positive_overflow(self, endless_ramp):
        """Check that energy wraps to start when next step exceeds end."""
        endless_ramp._start_energy = 10.0
        endless_ramp._end_energy = 20.0
        endless_ramp._delta_energy = 1.0
        endless_ramp.current_energy = 20.0
        endless_ramp.increment_energy()
        assert endless_ramp.current_energy == 10.0

    def test_wraps_to_start_on_negative_overflow(self, endless_ramp):
        """Check that energy wraps to start when next step goes below end."""
        endless_ramp._start_energy = 20.0
        endless_ramp._end_energy = 10.0
        endless_ramp._delta_energy = -1.0
        endless_ramp.current_energy = 10.0
        endless_ramp.increment_energy()
        assert endless_ramp.current_energy == 20.0

    def test_normal_increment_within_range_positive(self, endless_ramp):
        """Check that energy is incremented normally within range."""
        endless_ramp._start_energy = 10.0
        endless_ramp._end_energy = 20.0
        endless_ramp._delta_energy = 1.0
        endless_ramp.current_energy = 15.0
        endless_ramp.increment_energy()
        assert endless_ramp.current_energy == 16.0

    def test_normal_increment_within_range_negative(self, endless_ramp):
        """Check that energy is decremented normally within range."""
        endless_ramp._start_energy = 20.0
        endless_ramp._end_energy = 10.0
        endless_ramp._delta_energy = -1.0
        endless_ramp.current_energy = 15.0
        endless_ramp.increment_energy()
        assert endless_ramp.current_energy == 14.0

    def test_ramp_never_finished(self, endless_ramp):
        """Check that an endless ramp is never considered finished."""
        endless_ramp._delta_energy = 1.0
        endless_ramp._end_energy = 20.0
        endless_ramp.current_energy = 100.0
        assert not endless_ramp.ramp_finished()


class TestConstantEnergyRamp:
    """Tests for ConstantEnergyRamp."""

    def test_ramp_never_finished(self, constant_ramp):
        """Check that a constant-energy ramp is never considered finished."""
        assert not constant_ramp.ramp_finished()

    def test_increment_energy_is_no_op(self, constant_ramp):
        """Check that increment_energy does not change the energy."""
        constant_ramp.current_energy = 42.0
        constant_ramp.increment_energy()
        assert constant_ramp.current_energy == 42.0


class TestStepProfileParsing:
    """Tests for EnergyRampABC._set_step_profile and step_profile property."""

    def test_abrupt_profile_stored(self, linear_ramp):
        """Check that an ABRUPT profile is stored correctly."""
        linear_ramp._set_step_profile((ABRUPT,))
        assert linear_ramp._step_profile == (ABRUPT,)

    def test_abrupt_profile_yields_empty_tuple(self, linear_ramp):
        """Check that step_profile returns empty tuple for ABRUPT profile."""
        linear_ramp._set_step_profile((ABRUPT,))
        assert linear_ramp.step_profile == ()

    def test_linear_profile_valid(self, linear_ramp):
        """Check that a valid LINEAR profile is stored correctly."""
        linear_ramp._set_step_profile((LINEAR, '5', '100'))
        assert linear_ramp._step_profile == (LINEAR, '5', '100')

    def test_linear_profile_too_few_params(self, linear_ramp):
        """Check that too few params for LINEAR falls back to ABRUPT."""
        linear_ramp._set_step_profile((LINEAR, '5'))
        assert linear_ramp._step_profile == (ABRUPT,)

    def test_linear_profile_too_many_params(self, linear_ramp):
        """Check that too many params for LINEAR falls back to ABRUPT."""
        linear_ramp._set_step_profile((LINEAR, '5', '100', '200'))
        assert linear_ramp._step_profile == (ABRUPT,)

    def test_linear_profile_non_integer_params(self, linear_ramp):
        """Check that non-integer params for LINEAR falls back to ABRUPT."""
        linear_ramp._set_step_profile((LINEAR, 'abc', 'xyz'))
        assert linear_ramp._step_profile == (ABRUPT,)

    def test_linear_profile_zero_n_steps(self, linear_ramp):
        """Check that n_steps <= 0 for LINEAR falls back to ABRUPT."""
        linear_ramp._set_step_profile((LINEAR, '0', '100'))
        assert linear_ramp._step_profile == (ABRUPT,)

    def test_linear_profile_negative_time(self, linear_ramp):
        """Check that negative total_time for LINEAR falls back to ABRUPT."""
        linear_ramp._set_step_profile((LINEAR, '5', '-10'))
        assert linear_ramp._step_profile == (ABRUPT,)

    def test_custom_profile_valid(self, linear_ramp):
        """Check that a valid custom (numeric) profile is stored correctly."""
        linear_ramp._set_step_profile(('0.5', '50'))
        assert linear_ramp._step_profile == ('0.5', '50')

    def test_custom_profile_multiple_pairs(self, linear_ramp):
        """Check that a multi-pair custom profile is stored correctly."""
        linear_ramp._set_step_profile(('0.3', '30', '0.7', '70'))
        assert linear_ramp._step_profile == ('0.3', '30', '0.7', '70')

    def test_custom_profile_odd_length(self, linear_ramp):
        """Check that an odd-length custom profile falls back to ABRUPT."""
        linear_ramp._set_step_profile(('0.5', '50', '0.8'))
        assert linear_ramp._step_profile == (ABRUPT,)

    def test_custom_profile_negative_time(self, linear_ramp):
        """Check that negative time in custom profile falls back to ABRUPT."""
        linear_ramp._set_step_profile(('0.5', '-50'))
        assert linear_ramp._step_profile == (ABRUPT,)

    def test_unknown_profile_shape(self, linear_ramp):
        """Check that an unknown profile shape falls back to ABRUPT."""
        linear_ramp._set_step_profile(('UNKNOWN_SHAPE',))
        assert linear_ramp._step_profile == (ABRUPT,)

    def test_custom_step_profile_from_strings(self, linear_ramp):
        """Check step_profile property from a numeric custom profile."""
        linear_ramp._step_profile = ('0.5', '50')
        # delta = current - previous: set previous to 80, current to 100
        linear_ramp._previous_energy = 80.0
        linear_ramp._current_energy = 100.0
        # fraction=0.5 → (0.5 - 1) * 20 + 100 = 90.0, time=50
        result = linear_ramp.step_profile
        assert result == (pytest.approx(90.0), 50)

    def test_custom_step_profile_from_strings_no_delta(self, linear_ramp):
        """Check step_profile returns empty tuple when there is no delta."""
        linear_ramp._step_profile = ('0.5', '50')
        linear_ramp._previous_energy = 100.0
        linear_ramp._current_energy = 100.0
        assert linear_ramp.step_profile == ()
