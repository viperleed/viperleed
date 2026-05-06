"""Tests for module energyramp of viperleed.gui.measure.classes."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    'Florian Dörr (@FlorianDoerr)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-03-17'
__license__ = 'GPLv3+'

import pytest
from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.gui.measure.classes.energyramp import ABRUPT
from viperleed.gui.measure.classes.energyramp import DEFAULT_DELTA
from viperleed.gui.measure.classes.energyramp import DEFAULT_END
from viperleed.gui.measure.classes.energyramp import DEFAULT_START
from viperleed.gui.measure.classes.energyramp import DELTA_E_NAME
from viperleed.gui.measure.classes.energyramp import LINEAR
from viperleed.gui.measure.classes.energyramp import MINIMUM_ENERGY
from viperleed.gui.measure.classes.energyramp import ConstantEnergyRamp
from viperleed.gui.measure.classes.energyramp import EnergyRampABC
from viperleed.gui.measure.classes.energyramp import SawtoothEnergyRamp
from viperleed.gui.measure.classes.energyramp import LinearEnergyRamp
from viperleed.gui.measure.classes.energyramp import get_ramp_from_settings
from viperleed.gui.measure.classes.settings import ViPErLEEDSettings


MODULE = 'viperleed.gui.measure.classes.energyramp'


class _ConcreteEnergyRamp(EnergyRampABC):
    """Concrete implementation used to test EnergyRampABC defaults."""

    @property
    def n_steps(self):
        """Return default energy steps from base class."""
        return super().n_steps

    @classmethod
    def get_settings_widgets(cls):
        """Return settings widgets from base class."""
        return super().get_settings_widgets()

    @classmethod
    def ramp_settings_ok(cls, energy_options):
        """Return settings validation result from base class."""
        return super().ramp_settings_ok(energy_options)

    def increment_energy(self):
        """Keep energy unchanged for tests."""

    def ramp_finished(self):
        """A concrete dummy ramp used for defaults tests is never done."""
        return False


@fixture(name='linear_ramp')
def fixture_linear_ramp():
    """Return a fresh LinearEnergyRamp instance."""
    return LinearEnergyRamp()


@fixture(name='sawtooth_ramp')
def fixture_sawtooth_ramp():
    """Return a fresh SawtoothEnergyRamp instance."""
    return SawtoothEnergyRamp()


@fixture(name='constant_ramp')
def fixture_constant_ramp():
    """Return a fresh ConstantEnergyRamp instance."""
    return ConstantEnergyRamp()


@fixture(name='concrete_ramp')
def fixture_concrete_ramp():
    """Return a concrete ramp exposing EnergyRampABC defaults."""
    return _ConcreteEnergyRamp()


def _make_settings(**kwargs):
    """Return a ViPErLEEDSettings."""
    settings_dict = {'energies': kwargs} if kwargs else {}
    settings = ViPErLEEDSettings()
    settings.read_dict(settings_dict)
    return settings


# pylint: disable-next=too-few-public-methods
class _FakeOption:
    """Minimal option stub for ramp_settings_ok tests."""

    def __init__(self, option_name, value):
        """Initialize name and value of fake option."""
        self.option_name = option_name
        self._value = value

    def get_(self):
        """Return the underlying option value."""
        return self._value


# pylint: disable-next=too-few-public-methods
class _FakeSpinBox:
    """Small spin-box stub used to avoid QWidget creation in tests."""

    def __init__(self, *_, **__):
        """Initialize fake spin box."""
        self.single_step = None

    # pylint: disable-next=invalid-name
    def setSingleStep(self, value):
        """Store single-step values for assertions."""
        self.single_step = value


# pylint: disable-next=too-few-public-methods
class _FakeSettingsDialogOption:
    """Small settings option stub used in widget-construction tests."""

    def __init__(self, option_name, handler_widget, *_args, **_kwargs):
        """Initialize fake settings dialog option."""
        self.option_name = option_name
        self.handler_widget = handler_widget


class TestReturnMatchingEnergyRamp:
    """Tests for get_ramp_from_settings."""

    def test_returns_linear_ramp_by_default(self):
        """Check that an empty settings object yields LinearEnergyRamp."""
        settings = _make_settings()
        result = get_ramp_from_settings(settings)
        assert result is LinearEnergyRamp

    def test_returns_linear_ramp_when_constant_energy_and_endless_false(self):
        """Check that constant_energy/endless=false yield LinearEnergyRamp."""
        settings = _make_settings(constant_energy='false', endless='false')
        result = get_ramp_from_settings(settings)
        assert result is LinearEnergyRamp

    def test_returns_constant_ramp_when_constant_energy_true(self):
        """Check that constant_energy=true yields ConstantEnergyRamp."""
        settings = _make_settings(constant_energy='true')
        result = get_ramp_from_settings(settings)
        assert result is ConstantEnergyRamp

    def test_returns_sawtooth_ramp_when_endless_true(self):
        """Check that endless=true yields SawtoothEnergyRamp."""
        settings = _make_settings(endless='true')
        result = get_ramp_from_settings(settings)
        assert result is SawtoothEnergyRamp

    def test_constant_energy_takes_priority_over_endless(self):
        """Check that constant_energy=true takes precedence over endless."""
        settings = _make_settings(constant_energy='true', endless='true')
        result = get_ramp_from_settings(settings)
        assert result is ConstantEnergyRamp

    _invalid = {
        'constant_energy not a bool': {'constant_energy': 'maybe'},
        'endless not a bool': {'endless': 'perhaps'},
    }

    @parametrize('section', _invalid.values(), ids=_invalid)
    def test_returns_linear_ramp_for_invalid_values(self, section):
        """Check that invalid boolean values fall back to LinearEnergyRamp."""
        settings = _make_settings(**section)
        result = get_ramp_from_settings(settings)
        assert result is LinearEnergyRamp


# pylint: disable=protected-access
# pylint: disable=magic-value-comparison
class TestEnergyRampABCDefaults:
    """Tests for default behavior implemented in EnergyRampABC."""

    def test_default_n_steps_returns_zero(self, concrete_ramp):
        """Check that the default abstract implementation yields 0 steps."""
        # pylint: disable-next=use-implicit-booleaness-not-comparison-to-zero
        assert concrete_ramp.n_steps == 0

    def test_default_ramp_settings_ok_returns_true(self):
        """Check that the default abstract implementation accepts settings."""
        assert _ConcreteEnergyRamp.ramp_settings_ok(()) == (True, '')

    def test_current_energy_setter_updates_previous(self, concrete_ramp):
        """Check that assigning current_energy preserves previous energy."""
        concrete_ramp.current_energy = 5.0
        concrete_ramp.current_energy = 7.0
        assert concrete_ramp.previous_energy == 5.0
        assert concrete_ramp.current_energy == 7.0

    def test_min_energy_property_reflects_internal_value(self, concrete_ramp):
        """Check min_energy property output."""
        concrete_ramp._min_energy = 4.2
        assert concrete_ramp.min_energy == 4.2


class TestLinearEnergyRampSetRamp:
    """Tests for LinearEnergyRamp.set_ramp and base set_ramp behavior."""

    def test_set_ramp_without_settings_uses_defaults(self, linear_ramp):
        """Check default ramp setup when no settings are provided."""
        linear_ramp._min_energy = 13.0
        linear_ramp._start_energy = 20.0
        linear_ramp._step_profile = (LINEAR, '2', '10')
        linear_ramp._delta_energy = 9.0
        linear_ramp._end_energy = 99.0

        linear_ramp.set_ramp(None)

        assert linear_ramp.min_energy == MINIMUM_ENERGY
        assert linear_ramp.start_energy == DEFAULT_START
        assert linear_ramp._step_profile == (ABRUPT,)
        assert linear_ramp._delta_energy == DEFAULT_DELTA
        assert linear_ramp._end_energy == DEFAULT_END

    def test_set_ramp_clamps_energies_to_minimum(self, linear_ramp):
        """Check that energies are never set below min_energy."""
        settings = _make_settings(
            min_energy='5.0',
            start_energy='2.0',
            delta_energy='0.5',
            end_energy='3.0',
            step_profile='abrupt',
        )

        linear_ramp.set_ramp(settings)

        assert linear_ramp.min_energy == 5.0
        assert linear_ramp.start_energy == 5.0
        assert linear_ramp._end_energy == 5.0
        assert linear_ramp._step_profile == (ABRUPT,)

    def test_set_ramp_converts_string_profile_to_tuple(self, linear_ramp):
        """Check that string profiles are converted to one-element tuples."""
        settings = _make_settings(
            min_energy='0.0',
            start_energy='1.0',
            delta_energy='1.0',
            end_energy='5.0',
            step_profile='abrupt',
        )

        linear_ramp.set_ramp(settings)

        assert linear_ramp._step_profile == (ABRUPT,)

    def test_set_ramp_handles_invalid_minimum_and_start_values(self,
                                                               linear_ramp):
        """Check fallback values when min/start energies are malformed."""
        settings = _make_settings(
            min_energy='bad min',
            start_energy='bad start',
            delta_energy='1.0',
            end_energy='5.0',
            step_profile='abrupt',
        )

        linear_ramp.set_ramp(settings)

        assert linear_ramp.min_energy == MINIMUM_ENERGY
        assert linear_ramp.start_energy == DEFAULT_START

    def test_set_ramp_handles_invalid_delta_and_end_values(self, linear_ramp):
        """Check fallback values when delta/end energies are malformed."""
        settings = _make_settings(
            min_energy='0.0',
            start_energy='1.0',
            delta_energy='bad delta',
            end_energy='bad end',
            step_profile='abrupt',
        )

        linear_ramp.set_ramp(settings)

        assert linear_ramp._delta_energy == DEFAULT_DELTA
        assert linear_ramp._end_energy == DEFAULT_END

    def test_set_ramp_handles_missing_delta_and_end(self, linear_ramp):
        """Check fallback values when delta/end options are missing."""
        settings = _make_settings(
            min_energy='0.0',
            start_energy='1.0',
            step_profile='abrupt',
        )

        linear_ramp.set_ramp(settings)

        assert linear_ramp._delta_energy == DEFAULT_DELTA
        assert linear_ramp._end_energy == DEFAULT_END

    def test_set_ramp_preserves_zero_delta(self, linear_ramp):
        """Check that a configured zero delta is kept by set_ramp."""
        settings = _make_settings(
            min_energy='0.0',
            start_energy='0.0',
            delta_energy='0.0',
            end_energy='5.0',
            step_profile='abrupt',
        )

        linear_ramp.set_ramp(settings)

        # pylint: disable-next=use-implicit-booleaness-not-comparison-to-zero
        assert linear_ramp._delta_energy == 0.0
        assert linear_ramp._end_energy == 5.0


class TestStepProfileFallbackAndGeneration:
    """Tests for fallback logic and linear interpolation profile generation."""

    def test_step_profile_falls_back_to_abrupt(self, linear_ramp):
        """Check ABRUPT fallback branch when parsing as fractions fails."""
        settings = _make_settings(step_profile='invalid')
        linear_ramp.set_ramp(settings)
        linear_ramp.current_energy = 3.0
        linear_ramp.current_energy = 6.0

        assert linear_ramp.step_profile == ()

    def test_step_profile_returns_linear_profile(self, linear_ramp):
        """Check linear-step generation when profile shape is LINEAR."""
        settings = _make_settings(step_profile=['linear', '5', '10'])
        linear_ramp.set_ramp(settings)
        linear_ramp.current_energy = 0.0
        linear_ramp.current_energy = 10.0

        assert linear_ramp.step_profile == pytest.approx(
            (1.0, 2, 3.0, 2, 5.0, 2, 7.0, 2, 9.0, 2)
            )

    def test_get_linear_step_returns_empty_when_time_is_zero(self,
                                                             linear_ramp):
        """Check that 0 ms time intervals do not produce steps."""
        linear_ramp.current_energy = 0.0
        linear_ramp.current_energy = 10.0

        # To calculate time intervals we divide the total time by the
        # number of steps. 5 // 20 then yields a time interval of 0 ms,
        # which should not produce any steps.
        assert linear_ramp._get_linear_step(20, 5) == ()
        assert linear_ramp._step_profile == (ABRUPT,)

    def test_get_linear_step_returns_empty_for_tiny_energy_delta(self,
                                                                 linear_ramp):
        """Check that tiny energy changes do not produce steps."""
        linear_ramp.current_energy = 1.0
        linear_ramp.current_energy = 1.0 + 1e-6

        assert linear_ramp._get_linear_step(2, 10) == ()
        assert linear_ramp._step_profile == (ABRUPT,)

    def test_set_step_profile_rejects_non_sequence_profile(self, linear_ramp):
        """Check that non-sequence profiles are rejected as invalid."""
        linear_ramp._set_step_profile(42)
        assert linear_ramp._step_profile == (ABRUPT,)


class TestRampSettingsWidgetsAndValidation:
    """Tests for get_settings_widgets and ramp_settings_ok methods."""

    # These are intentionally broken settings that we use to check that
    # ramp_settings_ok correctly indentifies invalid configurations.
    _invalid_linear_settings = {
        'zero delta': (10.0, 0.0, 20.0, f'{DELTA_E_NAME} cannot be zero.'),
        'positive delta descending range': (
            10.0,
            1.0,
            9.0,
            f'For positive {DELTA_E_NAME}',
            ),
        'negative delta ascending range': (
            10.0,
            -1.0,
            11.0,
            f'For negative {DELTA_E_NAME}',
            ),
        }

    def test_linear_settings_widgets_include_expected_options(self, mocker):
        """Check linear widget setup without creating real Qt widgets."""
        mocker.patch(f'{MODULE}.CoercingDoubleSpinBox', _FakeSpinBox)
        mocker.patch(f'{MODULE}.SettingsDialogOption',
                     _FakeSettingsDialogOption)
        widgets = LinearEnergyRamp.get_settings_widgets()

        assert [option.option_name for option in widgets] == [
            'start_energy',
            'delta_energy',
            'end_energy',
            ]
        assert widgets[1].handler_widget.single_step == 0.5

    def test_constant_settings_widgets_include_start_only(self, mocker):
        """Check constant widget setup without creating real Qt widgets."""
        mocker.patch(f'{MODULE}.CoercingDoubleSpinBox', _FakeSpinBox)
        mocker.patch(f'{MODULE}.SettingsDialogOption',
                     _FakeSettingsDialogOption)
        widgets = ConstantEnergyRamp.get_settings_widgets()

        assert [option.option_name for option in widgets] == ['start_energy']

    @parametrize('start, delta, end, reason_fragment',
                 _invalid_linear_settings.values(),
                 ids=_invalid_linear_settings)
    def test_linear_ramp_settings_validation_rejects_invalid_options(
            self, start, delta, end, reason_fragment):
        """Check invalid option combinations are rejected."""
        options = (
            _FakeOption('start_energy', start),
            _FakeOption('delta_energy', delta),
            _FakeOption('end_energy', end),
            )

        is_ok, reason = LinearEnergyRamp.ramp_settings_ok(options)

        assert not is_ok
        assert reason_fragment in reason

    def test_linear_ramp_settings_validation_accepts_valid_options(self):
        """Check a valid option combination is accepted."""
        options = (
            _FakeOption('start_energy', 10.0),
            _FakeOption('delta_energy', 0.5),
            _FakeOption('end_energy', 12.0),
            )

        assert LinearEnergyRamp.ramp_settings_ok(options) == (True, '')

    def test_constant_ramp_settings_validation_always_ok(self):
        """Check constant-energy validation always succeeds."""
        assert ConstantEnergyRamp.ramp_settings_ok(()) == (True, '')


# pylint: disable-next=too-few-public-methods
class TestLinearEnergyRampIncrement:
    """Tests for LinearEnergyRamp.increment_energy."""

    def test_increment_energy_updates_current_and_previous(self, linear_ramp):
        """Check that incrementing updates both current and previous values."""
        settings = _make_settings(delta_energy='1.5')
        linear_ramp.set_ramp(settings)
        linear_ramp.current_energy = 10.0

        linear_ramp.increment_energy()

        assert linear_ramp.previous_energy == 10.0
        assert linear_ramp.current_energy == pytest.approx(11.5)


class TestLinearEnergyRampFinished:
    """Tests for LinearEnergyRamp.ramp_finished."""

    def test_n_steps_positive_delta(self, linear_ramp):
        """Check n_steps computed correctly for positive delta."""
        settings = _make_settings(start_energy='10.0', end_energy='20.0',
                                  delta_energy='2.0')
        linear_ramp.set_ramp(settings)
        assert linear_ramp.n_steps == 6  # 10, 12, 14, 16, 18, 20

    def test_n_steps_negative_delta(self, linear_ramp):
        """Check n_steps computed correctly for negative delta."""
        settings = _make_settings(start_energy='20.0', end_energy='10.0',
                                  delta_energy='-2.0')
        linear_ramp.set_ramp(settings)
        assert linear_ramp.n_steps == 6  # 20, 18, 16, 14, 12, 10

    def test_n_steps_zero_delta(self, linear_ramp):
        """Check n_steps is 0 for 0.0 delta."""
        settings = _make_settings(start_energy='10.0', end_energy='20.0',
                                  delta_energy='0.0')
        linear_ramp.set_ramp(settings)

        # pylint: disable-next=use-implicit-booleaness-not-comparison-to-zero
        assert linear_ramp.n_steps == 0

    @parametrize('delta, end, current', ((1.0, 20.0, 20.0),
                                         (1.0, 20.0, 25.0),
                                         (-1.0, 10.0, 10.0),
                                         (-1.0, 10.0, 5.0)))
    def test_ramp_finished(self, linear_ramp, delta, end, current):
        """Check ramp correctly identifies when it is finished."""
        settings = _make_settings(end_energy=str(end), delta_energy=str(delta))
        linear_ramp.set_ramp(settings)
        linear_ramp.current_energy = current
        assert linear_ramp.ramp_finished()

    @parametrize('delta, end, current', ((1.0, 20.0, 19.0),
                                         (-1.0, 10.0, 11.0)))
    def test_ramp_not_finished(self, linear_ramp, delta, end, current):
        """Check ramp not finished when next step is within range."""
        settings = _make_settings(end_energy=str(end), delta_energy=str(delta))
        linear_ramp.set_ramp(settings)
        linear_ramp.current_energy = current
        assert not linear_ramp.ramp_finished()


# pylint: disable=too-many-arguments
# pylint: disable=too-many-positional-arguments
class TestSawtoothEnergyRampIncrementEnergy:
    """Tests for SawtoothEnergyRamp.increment_energy."""

    def test_n_steps_is_zero(self, sawtooth_ramp):
        """Check sawtooth-ramp n_steps value."""
        # pylint: disable-next=use-implicit-booleaness-not-comparison-to-zero
        assert sawtooth_ramp.n_steps == 0

    @parametrize('start, end, delta, current, expected',
                 (('10.0', '20.0', '1.0', 20.0, 10.0),
                  ('20.0', '10.0', '-1.0', 10.0, 20.0)))
    def test_wraps_to_start_on_overflow(self, sawtooth_ramp, start, end,
                                        delta, current, expected):
        """Check that energy wraps to start when limit is exceeded."""
        settings = _make_settings(start_energy=start, end_energy=end,
                                  delta_energy=delta)
        sawtooth_ramp.set_ramp(settings)
        sawtooth_ramp.current_energy = current
        sawtooth_ramp.increment_energy()
        assert sawtooth_ramp.current_energy == expected

    @parametrize('start, end, delta, current, expected',
                 (('10.0', '20.0', '1.0', 15.0, 16.0),
                  ('20.0', '10.0', '-1.0', 15.0, 14.0)))
    def test_normal_increment_within_range(self, sawtooth_ramp, start, end,
                                           delta, current, expected):
        """Check that energy is incremented normally within range."""
        settings = _make_settings(start_energy=start, end_energy=end,
                                  delta_energy=delta)
        sawtooth_ramp.set_ramp(settings)
        sawtooth_ramp.current_energy = current
        sawtooth_ramp.increment_energy()
        assert sawtooth_ramp.current_energy == expected

    def test_ramp_never_finished(self, sawtooth_ramp):
        """Check that a sawtooth ramp is never considered finished."""
        settings = _make_settings(end_energy='20.0', delta_energy='1.0')
        sawtooth_ramp.set_ramp(settings)
        sawtooth_ramp.current_energy = 100.0
        assert not sawtooth_ramp.ramp_finished()
# pylint: enable=too-many-arguments
# pylint: enable=too-many-positional-arguments


class TestConstantEnergyRamp:
    """Tests for ConstantEnergyRamp."""

    def test_n_steps_is_zero(self, constant_ramp):
        """Check constant-ramp n_steps value."""
        # pylint: disable-next=use-implicit-booleaness-not-comparison-to-zero
        assert constant_ramp.n_steps == 0

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

    def test_linear_profile_valid(self, linear_ramp):
        """Check that a valid LINEAR profile is stored correctly."""
        linear_ramp._set_step_profile((LINEAR, '5', '100'))
        assert linear_ramp._step_profile == (LINEAR, '5', '100')

    @parametrize('profile', ((LINEAR, '5'),
                             (LINEAR, '5', '100', '200'),
                             (LINEAR, 'abc', 'xyz'),
                             (LINEAR, '0', '100'),
                             (LINEAR, '5', '-10')))
    def test_linear_profile_invalid_falls_back_to_abrupt(self, linear_ramp,
                                                         profile):
        """Check that invalid LINEAR params fall back to ABRUPT."""
        linear_ramp._set_step_profile(profile)
        assert linear_ramp._step_profile == (ABRUPT,)

    def test_custom_profile_valid(self, linear_ramp):
        """Check that a valid custom (numeric) profile is stored correctly."""
        linear_ramp._set_step_profile(('0.5', '50'))
        assert linear_ramp._step_profile == ('0.5', '50')

    def test_custom_profile_multiple_pairs(self, linear_ramp):
        """Check that a multi-pair custom profile is stored correctly."""
        linear_ramp._set_step_profile(('0.3', '30', '0.7', '70'))
        assert linear_ramp._step_profile == ('0.3', '30', '0.7', '70')

    @parametrize('profile', (('0.5', '50', '0.8'),
                             ('0.5', '-50'),
                             ('UNKNOWN_SHAPE',)))
    def test_custom_profile_invalid_falls_back_to_abrupt(self, linear_ramp,
                                                         profile):
        """Check that invalid custom profiles fall back to ABRUPT."""
        linear_ramp._set_step_profile(profile)
        assert linear_ramp._step_profile == (ABRUPT,)

    def test_custom_step_profile_from_strings(self, linear_ramp):
        """Check step_profile property from a numeric custom profile."""
        settings = _make_settings(step_profile=['0.5', '50'])
        linear_ramp.set_ramp(settings)
        # delta = current - previous: set previous to 80, current to 100
        linear_ramp.current_energy = 80.0
        linear_ramp.current_energy = 100.0
        # fraction=0.5 → (0.5 - 1) * 20 + 100 = 90.0, time=50
        result = linear_ramp.step_profile
        assert result == (pytest.approx(90.0), 50)

    def test_custom_step_profile_from_strings_no_delta(self, linear_ramp):
        """Check step_profile returns empty tuple when there is no delta."""
        settings = _make_settings(step_profile=['0.5', '50'])
        linear_ramp.set_ramp(settings)
        linear_ramp.current_energy = 100.0
        linear_ramp.current_energy = 100.0
        assert linear_ramp.step_profile == ()
