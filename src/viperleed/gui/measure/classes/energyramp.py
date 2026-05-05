"""Module energyramp of viperleed.gui.measure.classes.

This module defines the EnergyRampABC class and various subclasses,
which decide which energy is next during a measurement.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    'Florian Dörr (@FlorianDoerr)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-02-27'
__license__ = 'GPLv3+'

from abc import abstractmethod
from configparser import NoSectionError, NoOptionError
import math

from viperleed.gui.measure.classes.abc import QMetaABC
from viperleed.gui.measure.classes.abc import QObjectSettingsErrors
from viperleed.gui.measure.classes.abc import QObjectWithError
from viperleed.gui.measure.classes.settings import NotASequenceError
from viperleed.gui.measure.dialogs.settingsdialog import SettingsDialogOption
from viperleed.gui.measure.dialogs.settingsdialog import SettingsTag
from viperleed.gui.measure.widgets.spinboxes import CoercingDoubleSpinBox


ABRUPT = 'abrupt'
LINEAR = 'linear'
DELTA_E_NAME = '\u0394E'
START_E_NAME = 'E start'
END_E_NAME = 'E end'
DEFAULT_DELTA = 0.5
DEFAULT_END = 0.0
DEFAULT_START = 0.0
MINIMUM_ENERGY = 0.0
MINIMUM_DELTA = 1e-4


def get_matching_energy_ramp(settings):
    """Determine and return the type of energy ramp."""
    constant_energy = False
    endless = False
    try:
        constant_energy = settings.getboolean('energies', 'constant_energy',
                                              fallback=False)
    except (TypeError, ValueError):
        pass
    try:
        endless = settings.getboolean('energies', 'endless', fallback=False)
    except (TypeError, ValueError):
        pass
    if constant_energy:
        return ConstantEnergyRamp
    if endless:
        return SawtoothEnergyRamp
    return LinearEnergyRamp


class EnergyRampABC(QObjectWithError, metaclass=QMetaABC):                      # TODO: Move profile settings over to controller settings.
    """Generic energy ramp class."""

    display_name = None
    info_text = None

    def __init__(self, *args, **kwargs):
        """Initialise generic energy ramp class."""
        super().__init__(*args, **kwargs)
        self._current_energy = 0.0
        self._min_energy = MINIMUM_ENERGY
        self._previous_energy = 0.0
        self._start_energy = DEFAULT_START
        self._step_profile = (ABRUPT,)

    @property
    @abstractmethod
    def n_steps(self):
        """Return the number of energy steps."""
        return 0

    @property
    def current_energy(self):
        """Return the current energy in electronvolts."""
        return self._current_energy

    @current_energy.setter
    def current_energy(self, new_energy):
        """Set a new value of the current_energy.

        Parameters
        ----------
        new_energy : float
            The new current energy
        """
        self._previous_energy = self.current_energy
        self._current_energy = new_energy

    @property
    def min_energy(self):
        """Return the minimum starting energy (in eV)."""
        return self._min_energy

    @property
    def previous_energy(self):
        """Return the previous energy in electronvolts."""
        return self._previous_energy

    @property
    def start_energy(self):
        """Return the first energy for the energy ramp.

        The returned value is limited below by a minimum energy
        (as found in 'energies/min_energy' if present, 0.0 eV otherwise).
        This is useful to avoid calibrating for the non-linearity of
        LEED electronics in the low-energy regime.

        Returns
        -------
        start_energy : float
            The first energy of the energy ramp.
        """
        return self._start_energy

    @property
    def step_profile(self):
        """Return a list of energies and times for setting the next energy.

        The returned value excludes the very last step, i.e.,
        self.current_energy and the settling time for it.
        A typical call to set_leed_energy would be
            measurement.set_leed_energy(*self._energy_ramp.step_profile,
                                        self.current_energy,
                                        last_settle_time,
                                        ...)

        Returns
        -------
        step_profile : tuple
            Sequence of energies and waiting intervals.
        """
        try:
            return self._step_profile_from_strings(self._step_profile)
        except (ValueError, TypeError):
            pass

        shape, *params = self._step_profile
        is_abrupt = shape.lower() == ABRUPT
        values = tuple() if is_abrupt else self._get_linear_step(*params)
        return values

    @abstractmethod
    def increment_energy(self):
        """Go to the next energy in the ramp."""

    @classmethod
    @abstractmethod
    def get_settings_widgets(cls):
        """Return the widgets necessary to set a ramp.

        Must be extended in subclasses to add the widgets necessary to
        set the settings for the ramp.

        Returns
        -------
        widgets : tuple
            A tuple of SettingsDialogOptions.
        """
        widget = CoercingDoubleSpinBox(decimals=1, soft_range=(0, 1000),
                                       suffix=' eV')
        tip = '<nobr>The energy at which the measurement starts.</nobr>'
        start = SettingsDialogOption('start_energy', widget, tooltip=tip,
                                     display_name=START_E_NAME,
                                     tags=SettingsTag.REGULAR)
        return (start,)

    @abstractmethod
    def ramp_finished(self):
        """Return whether the energy ramp has been finished."""

    @classmethod
    @abstractmethod
    def ramp_settings_ok(cls, energy_options):
        """Return whether the settings for the ramp are ok.

        Parameters
        ----------
        energy_options : tuple of SettingsDialogOptions
            The energy options to check. These are the widgets used to
            set the ramp settings, as returned by get_settings_widgets.

        Returns
        -------
        settings_ok : bool
            Whether the settings selected in the widgets are
            acceptable or not.
        reason : str
            A descriptive string elaborating why the settings
            are not acceptable.
        """
        return True, ''

    def set_ramp(self, settings):
        """Set the energy ramp from settings.

        This method can be extended in subclasses to collect
        the other remaining settings required to perform an
        energy ramp. super().set_ramp must be called at
        the start of reimplementations.

        Parameters
        ----------
        settings : ViPErLEEDSettings
            Settings containing the energy ramp.

        Returns
        -------
        None.
        """
        if not settings:
            self._min_energy = MINIMUM_ENERGY
            self._start_energy = max(self.min_energy, DEFAULT_START)            # TODO: without settings this will default to 0.0 right now. (Used to be 5.0 for energy calibration)
            self._step_profile = (ABRUPT,)
            return
        try:
            self._min_energy = settings.getfloat('energies', 'min_energy',
                                                 fallback=MINIMUM_ENERGY)
        except (TypeError, ValueError):
            # Not a float
            self._min_energy = MINIMUM_ENERGY
        try:
            start_e = settings.getfloat('energies', 'start_energy',
                                        fallback=DEFAULT_START)
        except (TypeError, ValueError):
            # Not a float
            self.emit_error(
                QObjectSettingsErrors.INVALID_SETTING_WITH_FALLBACK,
                '', 'energies/start_energy', DEFAULT_START
                )
            start_e = DEFAULT_START
        self._start_energy = max(self.min_energy, start_e)
        try:
            profile = settings.getsequence('energies', 'step_profile',
                                           fallback=(ABRUPT,))
        except NotASequenceError:
            profile = settings['energies']['step_profile']
        if isinstance(profile, str):
            profile = (profile,)
        self._set_step_profile(profile)

    def _get_linear_step(self, *params):
        """Return energies and times for a simple linear step."""
        n_steps, tot_time = (int(p) for p in params)
        delta_t = tot_time // n_steps
        delta_e = self.current_energy - self.previous_energy
        if not delta_t or abs(delta_e) < MINIMUM_DELTA:
            return tuple()

        # Make a line of the form f(t) = t/tot_time, with t == 0
        # the time at which self.current_energy is set, and choose
        # an (almost) equally-spaced time grid, with the first
        # interval compensating for non integer-divisibility
        times = [-tot_time, *(-(i-1)*delta_t for i in range(n_steps, 0, -1))]
        intervals = (tj - ti for ti, tj in zip(times, times[1:]))

        # The best way to approximate a function with a piecewise
        # constant signal is to have values fk equal to the average
        # of f over the k-th interval. For our line, this means
        # f[k] = (t[k] + t[k+1]) / (2*tot_time)
        slope = delta_e / (2 * tot_time)
        energies = (slope*(ti + tj) + self.current_energy
                    for ti, tj in zip(times, times[1:]))

        # Finally interleave energies and times
        return tuple(v for tup in zip(energies, intervals) for v in tup)

    def _set_custom_profile(self, profile):
        """Check and set custom step profile.

        Parameters
        ----------
        profile : sequence
            A sequence describing a step profile.

        Returns
        -------
        None.
        """
        try:
            n_items = len(profile)
        except TypeError:
            self.emit_error(QObjectSettingsErrors.INVALID_SETTINGS,
                            'energies/step_profile',
                            'Step profiles must be sequences.')
            self._step_profile = (ABRUPT,)
            return

        # Casting may break. This is caught in self._set_step_profile.
        for fraction in profile[::2]:
            float(fraction)
        for time_ in profile[1::2]:
            time_ = int(time_)
            if time_ < 0:
                self.emit_error(QObjectSettingsErrors.INVALID_SETTINGS,
                                'energies/step_profile',
                                'Time intervals must be non-negative.')
                self._step_profile = (ABRUPT,)
                return

        if n_items < 2 or n_items % 2 != 0:
            self.emit_error(
                QObjectSettingsErrors.INVALID_SETTINGS,
                'energies/step_profile',
                'Invalid custom step profile: expected an even number of '
                f'entries (fraction, time, ...), found {n_items}.'
            )
            self._step_profile = (ABRUPT,)
            return
        self._step_profile = profile

    def _set_linear_profile(self, profile):
        """Check and set linear step profile.

        Parameters
        ----------
        profile : sequence
            A sequence describing a step profile.

        Returns
        -------
        None.
        """
        _, *params = profile
        if len(params) != 2:
            self.emit_error(QObjectSettingsErrors.INVALID_SETTINGS,
                            'energies/step_profile',
                            'Too many/few parameters for linear profile. '
                            f'Expected 2, found {len(params)}')
            self._step_profile = (ABRUPT,)
            return
        try:
            n_steps, tot_time = (int(p) for p in params)
        except (TypeError, ValueError):
            self.emit_error(QObjectSettingsErrors.INVALID_SETTINGS,
                            'energies/step_profile',
                            'Could not convert to integer the '
                            'parameters for linear profile')
            self._step_profile = (ABRUPT,)
            return
        if n_steps <= 0 or tot_time <= 0:
            self.emit_error(QObjectSettingsErrors.INVALID_SETTINGS,
                            'energies/step_profile',
                            'Linear-step parameters should be '
                            'positive integers')
            self._step_profile = (ABRUPT,)
            return
        self._step_profile = profile

    def _set_step_profile(self, profile):
        """Set the step profile from settings.

        Parameters
        ----------
        profile : sequence
            A sequence describing a step profile.

        Returns
        -------
        None.
        """
        try:
            self._set_custom_profile(profile)
        except (ValueError, TypeError):
            # Not a custom profile.
            pass
        else:
            return

        shape, *_ = profile
        if shape.lower() == ABRUPT:
            pass
        elif shape.lower() == LINEAR:
            self._set_linear_profile(profile)
            return
        else:
            self.emit_error(QObjectSettingsErrors.INVALID_SETTINGS,
                            'energies/step_profile',
                            f'Unknown profile shape {shape}')
            profile = (ABRUPT,)
        self._step_profile = profile

    def _step_profile_from_strings(self, profile):
        """Return a tuple of energies and times from strings."""
        delta = self.current_energy - self.previous_energy
        if abs(delta) < MINIMUM_DELTA:
            return tuple()

        energies_times = [0]*len(profile)
        for i, fraction in enumerate(profile[::2]):
            # We shift by -1 here in order to display to the user that
            # the 'current_energy' (the energy before the energy step)
            # is equivalent to a step fraction of 0 and the next energy
            # is equivalent to 1. We have to do this as the
            # current_energy is already incremented to the next energy.
            this_delta = (float(fraction) - 1) * delta
            energies_times[2*i] = this_delta + self.current_energy
        for i, time_ in enumerate(profile[1::2]):
            time_ = int(time_)
            energies_times[2*i+1] = time_
        return tuple(energies_times)


class LinearEnergyRamp(EnergyRampABC):
    """Generic linear energy ramp."""

    display_name = 'Linear energy ramp'
    info_text = ('<nobr>Linearly increases or decreases the energy'
                 f'</nobr> until reaching the {END_E_NAME}.')

    def __init__(self, *args, **kwargs):
        """Initialize LinearEnergyRamp."""
        self._delta_energy = DEFAULT_DELTA
        self._end_energy = DEFAULT_END
        super().__init__(*args, **kwargs)

    @property
    def energy_range(self):
        """Return the energy range."""
        return abs(self._end_energy - self._start_energy)

    @property
    def n_steps(self):
        """Return the number of energy steps."""
        if abs(self._delta_energy) < MINIMUM_DELTA:
            return 0
        return 1 + math.floor(self.energy_range/abs(self._delta_energy))

    @classmethod
    def get_settings_widgets(cls):
        """Return the widgets necessary to set a linear ramp.

        Returns
        -------
        widgets : tuple
            A tuple of SettingsDialogOptions.
        """
        widgets = super().get_settings_widgets()
        delta_widget = CoercingDoubleSpinBox(decimals=1,
                                             soft_range=(-1000, 1000),
                                             suffix=' eV')
        delta_widget.setSingleStep(0.5)
        tip = ('<nobr>The energy difference between '
               'two measurement steps.</nobr>')
        delta = SettingsDialogOption('delta_energy', delta_widget, tooltip=tip,
                                     display_name=DELTA_E_NAME,
                                     tags=SettingsTag.REGULAR)
        end_widget = CoercingDoubleSpinBox(decimals=1, soft_range=(0, 1000),
                                           suffix=' eV')
        tip = ('<nobr>The energy value at which </nobr>'
               'the measurement will finish.')
        end = SettingsDialogOption('end_energy', end_widget, tooltip=tip,
                                   display_name=END_E_NAME,
                                   tags=SettingsTag.REGULAR)
        return (*widgets, delta, end)

    @classmethod
    def ramp_settings_ok(cls, energy_options):
        """Return whether the settings for the ramp are ok.

        Parameters
        ----------
        energy_options : tuple of SettingsDialogOptions
            The energy options to check. These are the widgets used to
            set the ramp settings, as returned by get_settings_widgets.

        Returns
        -------
        settings_ok : bool
            Whether the settings selected in the widgets are
            acceptable or not.
        reason : str
            A descriptive string elaborating why the settings
            are not acceptable.
        """
        energy_dict = {option.option_name: float(option.get_())
                       for option in energy_options}
        start_energy = energy_dict.get('start_energy', DEFAULT_START)
        delta_energy = energy_dict.get('delta_energy', DEFAULT_DELTA)
        end_energy = energy_dict.get('end_energy', DEFAULT_END)
        if abs(delta_energy) < MINIMUM_DELTA:
            return False, (f'{DELTA_E_NAME} cannot be zero.'
                           ' Use constant energy ramp.')
        if delta_energy > 0 and end_energy < start_energy:
            return False, (f'For positive {DELTA_E_NAME}, {END_E_NAME} '
                           f'should be greater than {START_E_NAME}.')
        if delta_energy < 0 and end_energy > start_energy:
            return False, (f'For negative {DELTA_E_NAME}, {END_E_NAME} '
                           f'should be less than {START_E_NAME}.')
        return True, ''

    def increment_energy(self):
        """Go to the next energy in the ramp."""
        self.current_energy += self._delta_energy

    def ramp_finished(self):
        """Return whether the energy ramp has been finished."""
        next_energy = self.current_energy + self._delta_energy
        if self._delta_energy > 0 and next_energy > self._end_energy:
            return True
        if self._delta_energy < 0 and next_energy < self._end_energy:
            return True
        return False

    def set_ramp(self, settings):
        """Set the energy ramp from settings.

        Parameters
        ----------
        settings : ViPErLEEDSettings
            Settings containing the energy ramp.

        Returns
        -------
        None.
        """
        super().set_ramp(settings)
        if not settings:
            self._delta_energy = DEFAULT_DELTA                                  # TODO: iv/time: fallback = 0.5, ecal: fallback = 5.0
            self._end_energy = DEFAULT_END                                      # TODO: iv/time: fallback = 0.0, ecal: fallback = 1000
            return
        try:
            self._delta_energy = settings.getfloat('energies', 'delta_energy')
        except (TypeError, ValueError):
            # Not a float
            self._delta_energy = DEFAULT_DELTA
            self.emit_error(QObjectSettingsErrors.INVALID_SETTINGS,
                            'energies/delta_energy', '')
        except (NoSectionError, NoOptionError):
            # Not present
            self._delta_energy = DEFAULT_DELTA
            self.emit_error(
                QObjectSettingsErrors.INVALID_SETTING_WITH_FALLBACK,
                '', 'energies/delta_energy', DEFAULT_DELTA
                )
        if abs(self._delta_energy) < MINIMUM_DELTA:
            self.emit_error(QObjectSettingsErrors.INVALID_SETTINGS,
                            'energies/delta_energy', f'{DELTA_E_NAME} was '
                            'set to 0. Use constant-energy mode instead.')
        try:
            self._end_energy = settings.getfloat('energies', 'end_energy')
        except (TypeError, ValueError):
            # Not a float
            self._end_energy = DEFAULT_END
            self.emit_error(QObjectSettingsErrors.INVALID_SETTINGS,
                            'energies/end_energy', '')
        except (NoSectionError, NoOptionError):
            # Not present
            self._end_energy = DEFAULT_END
            self.emit_error(
                QObjectSettingsErrors.INVALID_SETTING_WITH_FALLBACK,
                '', 'energies/end_energy', DEFAULT_END
                )


class ConstantEnergyRamp(EnergyRampABC):
    """Constant energy ramp."""

    display_name = 'Constant energy ramp'
    info_text = ('<nobr>Keeps the energy constant throughout the</nobr>'
                 ' measurement. Useful for time-resolved measurements '
                 'at fixed energy.')

    @property
    def n_steps(self):
        """Return the number of energy steps."""
        return 0

    @classmethod
    def get_settings_widgets(cls):
        """Return the settings widgets for the constant energy ramp."""
        return super().get_settings_widgets()

    @classmethod
    def ramp_settings_ok(cls, *_):
        """The settings for a constant energy ramp are always ok."""
        return True, ''

    def increment_energy(self):
        """Advance the ramp state without changing the energy."""
        # This is necessary to avoid recalling set_leed_energy.
        self._previous_energy = self.current_energy

    def ramp_finished(self):                                                    # TODO: we could add a counter here that allows us to repeat constant energy ramps only a limited amount of times
        """Return whether the energy ramp has been finished."""
        return False  # A constant energy ramp is never finished.


class SawtoothEnergyRamp(LinearEnergyRamp):
    """Sawtooth energy ramp."""

    display_name = 'Sawtooth energy ramp'
    info_text = ('<nobr>Linearly increases or decreases the energy'
                 f'</nobr> until reaching the {END_E_NAME}. Energy resets '
                 f'to {START_E_NAME} when reaching the {END_E_NAME}.')

    @property
    def n_steps(self):
        """Return the number of energy steps."""
        return 0

    def increment_energy(self):
        """Go to the next energy in the ramp."""
        if super().ramp_finished():
            self.current_energy = self.start_energy
            return
        super().increment_energy()

    def ramp_finished(self):
        """Return whether the energy ramp has been finished."""
        return False  # An endless energy ramp is never finished.
