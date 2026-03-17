"""Tests for module energyramp of viperleed.gui.measure.classes."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    'Florian Dörr (@FlorianDoerr)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-03-17'
__license__ = 'GPLv3+'

from pytest_cases import parametrize
import pytest

from viperleed.gui.measure.classes.energyramp import ABRUPT
from viperleed.gui.measure.classes.energyramp import ConstantEnergyRamp
from viperleed.gui.measure.classes.energyramp import EndlessLinearEnergyRamp
from viperleed.gui.measure.classes.energyramp import EnergyRampABC
from viperleed.gui.measure.classes.energyramp import LinearEnergyRamp
from viperleed.gui.measure.classes.settings import ViPErLEEDSettings


def _make_settings(**energies_options):
    """Return a ViPErLEEDSettings with the given energies section."""
    settings = ViPErLEEDSettings()
    settings.read_dict({'energies': energies_options})
    return settings


def _make_linear_ramp(start=10.0, end=100.0, delta=5.0):
    """Return a LinearEnergyRamp with the given start/end/delta."""
    ramp = LinearEnergyRamp()
    ramp._start_energy = start
    ramp._end_energy = end
    ramp._delta_energy = delta
    ramp.current_energy = start
    return ramp


class TestReturnMatchingEnergyRamp:
    """Tests for EnergyRampABC.return_matching_energy_ramp."""

    def test_default_returns_linear(self):
        """Check that an empty settings returns a LinearEnergyRamp."""
        settings = _make_settings()
        result = EnergyRampABC.return_matching_energy_ramp(settings)
        assert result is LinearEnergyRamp

    def test_constant_energy_flag(self):
        """Check that constant_energy=true returns ConstantEnergyRamp."""
        settings = _make_settings(constant_energy='true')
        result = EnergyRampABC.return_matching_energy_ramp(settings)
        assert result is ConstantEnergyRamp

    def test_endless_flag(self):
        """Check that endless=true returns EndlessLinearEnergyRamp."""
        settings = _make_settings(endless='true')
        result = EnergyRampABC.return_matching_energy_ramp(settings)
        assert result is EndlessLinearEnergyRamp

    def test_constant_energy_takes_precedence_over_endless(self):
        """Check that constant_energy=true overrides endless=true."""
        settings = _make_settings(constant_energy='true', endless='true')
        result = EnergyRampABC.return_matching_energy_ramp(settings)
        assert result is ConstantEnergyRamp

    @parametrize(flag=('constant_energy', 'endless'))
    def test_invalid_bool_flag_falls_back_to_linear(self, flag):
        """Check that an invalid boolean flag falls back to LinearEnergyRamp."""
        settings = _make_settings(**{flag: 'notabool'})
        result = EnergyRampABC.return_matching_energy_ramp(settings)
        assert result is LinearEnergyRamp

    def test_no_energies_section_returns_linear(self):
        """Check that missing energies section returns LinearEnergyRamp."""
        settings = ViPErLEEDSettings()
        result = EnergyRampABC.return_matching_energy_ramp(settings)
        assert result is LinearEnergyRamp


class TestLinearRampFinished:
    """Tests for LinearEnergyRamp.ramp_finished() with +/- deltas."""

    def test_positive_delta_not_finished_at_start(self):
        """Check ramp is not finished when at start energy."""
        ramp = _make_linear_ramp(start=10.0, end=100.0, delta=5.0)
        assert not ramp.ramp_finished()

    def test_positive_delta_not_finished_at_exact_end(self):
        """Check ramp is not finished when next step lands exactly on end."""
        ramp = _make_linear_ramp(start=10.0, end=100.0, delta=5.0)
        ramp.current_energy = 95.0  # 95 + 5 == 100, not > 100
        assert not ramp.ramp_finished()

    def test_positive_delta_finished_past_end(self):
        """Check ramp is finished when next step would exceed end."""
        ramp = _make_linear_ramp(start=10.0, end=100.0, delta=5.0)
        ramp.current_energy = 96.0  # 96 + 5 = 101 > 100
        assert ramp.ramp_finished()

    def test_negative_delta_not_finished_at_start(self):
        """Check ramp is not finished at start with negative delta."""
        ramp = _make_linear_ramp(start=100.0, end=10.0, delta=-5.0)
        assert not ramp.ramp_finished()

    def test_negative_delta_not_finished_at_exact_end(self):
        """Check ramp is not finished when next step lands exactly on end."""
        ramp = _make_linear_ramp(start=100.0, end=10.0, delta=-5.0)
        ramp.current_energy = 15.0  # 15 - 5 == 10, not < 10
        assert not ramp.ramp_finished()

    def test_negative_delta_finished_past_end(self):
        """Check ramp is finished when next step would go below end."""
        ramp = _make_linear_ramp(start=100.0, end=10.0, delta=-5.0)
        ramp.current_energy = 14.0  # 14 - 5 = 9 < 10
        assert ramp.ramp_finished()


class TestEndlessLinearEnergyRampWraparound:
    """Tests for EndlessLinearEnergyRamp wraparound behavior."""

    def _make_ramp(self, start=10.0, end=100.0, delta=5.0):
        """Return an EndlessLinearEnergyRamp with the given parameters."""
        ramp = EndlessLinearEnergyRamp()
        ramp._start_energy = start
        ramp._end_energy = end
        ramp._delta_energy = delta
        ramp.current_energy = start
        return ramp

    def test_ramp_never_finished(self):
        """Check that ramp_finished always returns False."""
        ramp = self._make_ramp()
        ramp.current_energy = 200.0  # way past end
        assert not ramp.ramp_finished()

    def test_positive_delta_normal_increment(self):
        """Check that energy increments normally when within range."""
        ramp = self._make_ramp(start=10.0, end=100.0, delta=5.0)
        ramp.current_energy = 50.0
        ramp.increment_energy()
        assert ramp.current_energy == pytest.approx(55.0)

    def test_positive_delta_wraparound_when_past_end(self):
        """Check wraparound when next step would exceed end energy."""
        ramp = self._make_ramp(start=10.0, end=100.0, delta=5.0)
        ramp.current_energy = 96.0  # 96 + 5 = 101 > 100 → wrap
        ramp.increment_energy()
        assert ramp.current_energy == pytest.approx(10.0)

    def test_positive_delta_no_wraparound_at_exact_end(self):
        """Check no wraparound when next step lands exactly on end."""
        ramp = self._make_ramp(start=10.0, end=100.0, delta=5.0)
        ramp.current_energy = 95.0  # 95 + 5 = 100, not > 100
        ramp.increment_energy()
        assert ramp.current_energy == pytest.approx(100.0)

    def test_negative_delta_normal_decrement(self):
        """Check that energy decrements normally when within range."""
        ramp = self._make_ramp(start=100.0, end=10.0, delta=-5.0)
        ramp.current_energy = 50.0
        ramp.increment_energy()
        assert ramp.current_energy == pytest.approx(45.0)

    def test_negative_delta_wraparound_when_past_end(self):
        """Check wraparound when next step would go below end energy."""
        ramp = self._make_ramp(start=100.0, end=10.0, delta=-5.0)
        ramp.current_energy = 14.0  # 14 - 5 = 9 < 10 → wrap
        ramp.increment_energy()
        assert ramp.current_energy == pytest.approx(100.0)

    def test_negative_delta_no_wraparound_at_exact_end(self):
        """Check no wraparound when next step lands exactly on end."""
        ramp = self._make_ramp(start=100.0, end=10.0, delta=-5.0)
        ramp.current_energy = 15.0  # 15 - 5 = 10, not < 10
        ramp.increment_energy()
        assert ramp.current_energy == pytest.approx(10.0)

    def test_energy_steps_is_zero(self):
        """Check that energy_steps is 0 for endless ramps."""
        ramp = self._make_ramp()
        assert ramp.energy_steps == 0


class TestConstantEnergyRamp:
    """Tests for ConstantEnergyRamp behavior."""

    def test_ramp_never_finished(self):
        """Check that a constant energy ramp is never finished."""
        ramp = ConstantEnergyRamp()
        assert not ramp.ramp_finished()

    def test_increment_does_not_change_energy(self):
        """Check that increment_energy is a no-op."""
        ramp = ConstantEnergyRamp()
        ramp.current_energy = 50.0
        ramp.increment_energy()
        assert ramp.current_energy == pytest.approx(50.0)

    def test_energy_steps_is_zero(self):
        """Check that energy_steps is 0 for constant ramps."""
        ramp = ConstantEnergyRamp()
        assert ramp.energy_steps == 0


class TestStepProfileParsing:
    """Tests for step_profile parsing and validation."""

    @pytest.fixture(autouse=True)
    def ramp_with_errors(self):
        """Create a LinearEnergyRamp and collect any emitted errors."""
        self.ramp = LinearEnergyRamp()
        self.ramp.set_ramp(None)
        self.errors = []
        self.ramp.error_occurred.connect(self.errors.append)

    def test_abrupt_profile(self):
        """Check that 'abrupt' sets the abrupt step profile."""
        self.ramp._set_step_profile(('abrupt',))
        assert self.ramp._step_profile == ('abrupt',)
        assert not self.errors

    def test_abrupt_profile_case_insensitive(self):
        """Check that 'ABRUPT' (uppercase) is also accepted."""
        self.ramp._set_step_profile(('ABRUPT',))
        assert self.ramp._step_profile == ('ABRUPT',)
        assert not self.errors

    def test_linear_profile_valid(self):
        """Check that a valid 'linear' profile is accepted."""
        self.ramp._set_step_profile(('linear', '5', '100'))
        assert self.ramp._step_profile == ('linear', '5', '100')
        assert not self.errors

    def test_linear_profile_too_few_params(self):
        """Check that a linear profile with too few params falls back."""
        self.ramp._set_step_profile(('linear', '5'))
        assert self.ramp._step_profile == (ABRUPT,)
        assert len(self.errors) == 1

    def test_linear_profile_too_many_params(self):
        """Check that a linear profile with too many params falls back."""
        self.ramp._set_step_profile(('linear', '5', '100', '200'))
        assert self.ramp._step_profile == (ABRUPT,)
        assert len(self.errors) == 1

    def test_linear_profile_non_integer_params(self):
        """Check that non-integer linear params fall back."""
        self.ramp._set_step_profile(('linear', 'x', 'y'))
        assert self.ramp._step_profile == (ABRUPT,)
        assert len(self.errors) == 1

    def test_linear_profile_zero_n_steps(self):
        """Check that n_steps=0 falls back for linear profile."""
        self.ramp._set_step_profile(('linear', '0', '100'))
        assert self.ramp._step_profile == (ABRUPT,)
        assert len(self.errors) == 1

    def test_linear_profile_negative_total_time(self):
        """Check that a negative total_time falls back for linear profile."""
        self.ramp._set_step_profile(('linear', '5', '-10'))
        assert self.ramp._step_profile == (ABRUPT,)
        assert len(self.errors) == 1

    def test_custom_profile_valid(self):
        """Check that a valid custom profile (fraction, time pairs) is set."""
        profile = ('0.3', '30', '0.7', '70')
        self.ramp._set_step_profile(profile)
        assert self.ramp._step_profile == profile
        assert not self.errors

    def test_custom_profile_single_pair(self):
        """Check that a single fraction/time pair is accepted."""
        profile = ('0.5', '50')
        self.ramp._set_step_profile(profile)
        assert self.ramp._step_profile == profile
        assert not self.errors

    def test_custom_profile_odd_number_of_entries(self):
        """Check that an odd number of custom entries falls back."""
        self.ramp._set_step_profile(('0.3', '30', '0.7'))
        assert self.ramp._step_profile == (ABRUPT,)
        assert len(self.errors) == 1

    def test_custom_profile_negative_time(self):
        """Check that a negative time in a custom profile falls back."""
        self.ramp._set_step_profile(('0.5', '-50'))
        assert self.ramp._step_profile == (ABRUPT,)
        assert len(self.errors) == 1

    def test_unknown_profile_shape(self):
        """Check that an unknown profile shape falls back to abrupt."""
        self.ramp._set_step_profile(('unknown_shape',))
        assert self.ramp._step_profile == (ABRUPT,)
        assert len(self.errors) == 1

    def test_step_profile_property_abrupt_returns_empty(self):
        """Check that the abrupt step_profile property returns an empty tuple."""
        self.ramp._step_profile = ('abrupt',)
        assert self.ramp.step_profile == ()

    def test_step_profile_property_custom_with_delta(self):
        """Check that a custom profile property returns scaled energies."""
        self.ramp._previous_energy = 10.0
        self.ramp._current_energy = 20.0
        self.ramp._step_profile = ('0.5', '50')
        result = self.ramp.step_profile
        # fraction=0.5 → (0.5 - 1) * delta + current = -0.5*10 + 20 = 15
        assert len(result) == 2
        assert result[0] == pytest.approx(15.0)
        assert result[1] == 50

    def test_step_profile_set_via_settings(self):
        """Check that step_profile is correctly parsed from settings."""
        settings = _make_settings(
            start_energy='10.0',
            end_energy='100.0',
            delta_energy='5.0',
            step_profile='(\'linear\', \'3\', \'90\')',
        )
        ramp = LinearEnergyRamp()
        ramp.set_ramp(settings)
        assert ramp._step_profile == ('linear', '3', '90')
