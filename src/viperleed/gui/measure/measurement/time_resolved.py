"""Module time_resolved of viperleed.gui.measure.measurement.

This module contains the definition of the TimeResolved class
for measuring data in a time-based fashion.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    'Florian Dörr (@FlorianDoerr)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2021-07-19'
__license__ = 'GPLv3+'

from PyQt5 import QtCore as qtc

from viperleed.gui.measure import hardwarebase as base
from viperleed.gui.measure.controller.abc import ControllerABC
from viperleed.gui.measure.measurement.abc import MeasurementABC
from viperleed.gui.measure.measurement.abc import MeasurementErrors


_UNIQUE = qtc.Qt.UniqueConnection


class TimeResolved(MeasurementABC):
    """Time resolved measurement class."""

    display_name = 'Time resolved'

    _mandatory_settings = [*MeasurementABC._mandatory_settings,
                           ('measurement_settings', 'is_continuous'),
                           ('measurement_settings', 'energy_step_duration'),]

    __request_continuous_mode = qtc.pyqtSignal(bool)   # On/Off

    def __init__(self, measurement_settings):
        """Initialise measurement class."""
        super().__init__(measurement_settings)
        self.__end_energy = 0                                                   # --> generator; It can happen that there is no proper "end energy" (endless measurement, e.g., wiggle)
        self.__delta_energy = 1                                                 # --> generator; There may not be a delta_energy (e.g., I(t) at fixed energy)
        self.__endless = False                                                  # --> generator;
        self.__constant_energy = False                                          # --> generator;
        self.__n_digits = 0                                                     # Again, used for saving images.

        # One timer to trigger a change of electron energy: we will
        # sit at each energy for self.energy_step_duration millisecs
        self.__energy_step_timer = qtc.QTimer(parent=self)
        self.__energy_step_timer.setSingleShot(True)
        self.__energy_step_timer.timeout.connect(
            self._ready_for_next_measurement
            )

        # A second timer to trigger multiple measurements within
        # each energy step. This is used only if the measurement
        # is triggered (i.e., not self.is_continuous). It cannot
        # time out faster than the longest time it takes any
        # controller to return measurements. The time interval
        # is read from the settings.
        self.__trigger_one_measurement = qtc.QTimer(parent=self)
        for ctrl in self.controllers:
            self.__connect_trigger_timeout(ctrl)

        if self.settings:
            # pylint: disable=redefined-variable-type
            self.__delta_energy = self.settings.getfloat(
                'measurement_settings', 'delta_energy', fallback=10
                )
            self.__end_energy = self.settings.getfloat(
                'measurement_settings', 'end_energy', fallback=10
                )
            self.__endless = self.settings.getboolean(
                'measurement_settings', 'endless', fallback=False
                )
            self.__constant_energy = self.settings.getboolean(
                'measurement_settings', 'constant_energy', fallback=False
                )
        if not self.primary_controller:
            return

        num_meas = (1 + round((self.__end_energy - self.start_energy)
                              / self.__delta_energy))
        self.data_points.time_resolved = True
        if self.is_continuous:
            self.__prepare_continuous_mode()
            self.data_points.nr_steps_total = num_meas
        else:
            self.__n_digits = len(str(num_meas))

    @property
    def energy_step_duration(self):
        """Return the duration of one energy step (msec)."""
        if not self.controllers:
            return -1

        min_t = self.energy_step_duration_min
        try:
            interval = self.settings.getint('measurement_settings',
                                            'energy_step_duration')
        except (TypeError, ValueError):
            # Not an int
            base.emit_error(self, MeasurementErrors.INVALID_SETTINGS,
                            'measurement_settings/energy_step_duration', '')
            interval = min_t

        if interval < min_t:
            base.emit_error(
                self, MeasurementErrors.INVALID_SETTING_WITH_FALLBACK,
                f'{interval} (too short)',
                'measurement_settings/energy_step_duration', min_t
                )
            interval = min_t
            self.settings.set('measurement_settings', 'energy_step_duration',
                              str(interval))
            self.__energy_step_timer.setInterval(interval)
        return interval

    @energy_step_duration.setter
    def energy_step_duration(self, new_interval):
        """Set the duration of one energy step (msec)."""
        min_t = self.energy_step_duration_min
        if new_interval < min_t:
            raise ValueError(f"{new_interval} is too small for the "
                             f"controllers. Minimum is {min_t}")
        self.settings.set('measurement_settings', 'measurement_interval',
                          str(new_interval))
        self.__energy_step_timer.setInterval(new_interval)

    @property
    def energy_step_duration_min(self):
        """Return the smallest duration of an energy step."""
        min_interval = self.measurement_interval_min
        # The minimum duration of a step depends on whether we are
        # in continuous or in 'triggered' mode: in the latter case,
        # we are not measuring for the first __trigger_one_measurement
        # interval after setting the energy.
        if self.is_continuous:
            # Here the extra 5 msec ensure that energy_step_duration
            # is always a tiny bit longer than measurement_interval.
            return min_interval + 5
        return 2 * min_interval

    @property
    def is_continuous(self):
        """Return whether the measurement is continuous."""
        try:
            return self.settings.getboolean('measurement_settings',
                                            'is_continuous')
        except ValueError:
            # Not a valid boolean
            base.emit_error(self, MeasurementErrors.INVALID_SETTINGS,
                            'measurement_settings/is_continuous', '')
            return False

    @property
    def measurement_interval(self):
        """Return the time interval between 'triggered' measurements (msec)."""
        if self.is_continuous or not self.controllers:
            return -1
        # Does not make much sense to use intervals shorter than
        # 50 msec for a "triggered" mode. One should rather use
        # continuous mode.
        min_t = self.measurement_interval_min
        max_t = self.measurement_interval_max
        try:
            interval = self.settings.getint('measurement_settings',
                                            'measurement_interval',
                                            fallback=min_t)
        except (TypeError, ValueError):
            # Not an int
            base.emit_error(self, MeasurementErrors.INVALID_SETTINGS,
                            'measurement_settings/measurement_interval', '')
            interval = min_t

        if not min_t <= interval <= max_t:
            txt = f"{interval} (too "
            txt += "short)" if interval < min_t else "long)"
            _, interval, _ = sorted((min_t, interval, max_t))
            base.emit_error(
                self, MeasurementErrors.INVALID_SETTING_WITH_FALLBACK,
                txt, 'measurement_settings/measurement_interval', interval
                )
            self.settings.set('measurement_settings', 'measurement_interval',
                              str(interval))
            self.__trigger_one_measurement.setInterval(interval)
        return interval

    @measurement_interval.setter
    def measurement_interval(self, new_interval):
        """Set the time interval between 'triggered' measurements (msec)."""
        # Does not make much sense to use intervals shorter than
        # 50 msec for a "triggered" mode. One should rather use
        # continuous mode.
        min_t = self.measurement_interval_min
        max_t = self.measurement_interval_max
        if not min_t <= new_interval <= max_t:
            raise ValueError(f"{new_interval} out of bounds. "
                             f"Should be between {min_t} and {max_t}")
        self.settings.set('measurement_settings', 'measurement_interval',
                          str(new_interval))
        self.__trigger_one_measurement.setInterval(new_interval)

    @property
    def measurement_interval_max(self):
        """Return the largest measurement-time interval."""
        return (self.energy_step_duration
                # We have to wait for measurements to come back
                - self.measurement_interval_min - 5)

    @property
    def measurement_interval_min(self):                                         # TODO: depending on how we handle cameras, should also account for them
        """Return the smallest measurement-time interval."""
        _MIN = 50                                                               # TODO: should probably be OS-dependent? This estimate is actually pretty bad (too small) when there's a lot of data to be plotted.
        min_interval = _MIN
        if self.controllers:
            min_interval = 2 * max(c.initial_delay for c in self.controllers)
            # add a little, in case initial_delay == N*measurement_interval
            min_interval += 5
        return max(_MIN, round(min_interval))

    def abort(self):
        """Interrupt measurement, save partial data, set energy to zero."""
        try:
            timers = (self.__trigger_one_measurement, self.__energy_step_timer)
        except AttributeError:
            # .abort() happened during super().__init__
            timers = tuple()

        for timer in timers:
            timer.stop()
        super().abort()

    def energy_generator(self):     # TODO: move to parent; improve.
        """Determine next energy to set.

        Determine the next energy to set using parameters from
        the config and self.current_energy.

        Returns
        -------
        energy : float
            Next energy to set.
        """
        if self.__constant_energy:
            return self.start_energy
        energy = self.current_energy + self.__delta_energy
        if self.__endless:
            self.new_data_available.emit()
            if energy > self.__end_energy:
                energy = self.start_energy
        return energy

    def start_next_measurement(self):
        """Set energy and measure.

        Set energy via the primary controller. Once this is done
        the returned about_to_trigger signal will start the
        measurement on all secondary controllers.

        Returns
        -------
        None.
        """
        super().start_next_measurement()
        self.__energy_step_timer.setInterval(self.energy_step_duration)
        _continuous = self.is_continuous
        about_to_trigger = self.primary_controller.about_to_trigger

        # Each energy step will begin when the energy has settled,
        # i.e., when the primary_controller is done setting the
        # last energy (self.current_energy), i.e., when it emits
        # .about_to_trigger. This is necessary, since we may be
        # using a non-abrupt self.step_profile.
        base.safe_connect(about_to_trigger, self.__energy_step_timer.start,
                          type=_UNIQUE)

        if _continuous and self.current_step_nr == 2:
            # In continuous mode, the secondary controllers
            # return measurements undisturbed, and should not
            # be triggered, except for the first energy step.
            for ctrl in self.secondary_controllers:
                base.safe_disconnect(about_to_trigger, ctrl.measure_now)

        # Notice that it does not make sense to have an actual settle
        # time for measurements in a time-resolved: in continuous mode
        # we want to get measurements as quickly as possible (and look
        # at short-term time traces); "triggered" mode can be used for
        # either long-term stability measurements or I(t). Both would
        # normally be done at fixed energy.
        # For a "triggered" measurement type, let the triggering be             # TODO: couldn't we instead always measure, and rather connect the about_to_trigger to __trigger_one_measurement.start?
        # done by the __trigger_one_measurement timer rather than by
        # the energy setter. This means that we will start measuring
        # after .measurement_interval msec rather than right now.
        self.set_leed_energy(*self.step_profile, self.current_energy, 0,
                             trigger_meas=_continuous)

        if _continuous:
            # No camera images saved when we
            # spam measurements at max speed
            return

        self.__trigger_one_measurement.setInterval(self.measurement_interval)
        self.__trigger_one_measurement.start()

        # TODO: use a flag to decide if we want to save images
        # TODO: not really clear how to handle cameras in this case:
        # Right now we are acquiring a single camera image per each
        # "energy step" since the __trigger_one_measurement timer
        # is not connected to camera.trigger_now. Also, we are not
        # really waiting some time for starting measurements, while
        # we are waiting for the camera. So images and measurements
        # are not overlapping in time. This is less than ideal.
        image_name = (f"{self.current_step_nr:0>{self.__n_digits}}_"
                      f"{self.current_energy:.1f}eV.tiff")
        for camera in self.cameras:
            camera.process_info.filename = image_name
            self.data_points.add_image_names(image_name)
        self._camera_timer.start(self.hv_settle_time)

    def __check_is_finished(self, _):                                           # TODO: I'm not happy with the name. __continue_when_primary_stopped?
        """Check if the measurement is finished in continuous mode.

        This method is used only if self.is_continuous. It does the
        same checks as in super()._ready_for_next_measurement(), i.e.,
        decides whether we should go to the next step or if the whole
        loop is over. It is needed in continuous mode so we can
        wait the acknowledgment from the primary controller that it
        has been in fact stopped at the end of an energy step.

        Returns
        -------
        None.
        """
        base.safe_disconnect(self.primary_controller.controller_busy,
                             self.__check_is_finished)
        if self.aborted:
            # We entered this call after the measurement was aborted,
            # likely while processing an unprocessed timeout event
            return

        self.data_points.calculate_times(continuous=True)
        self.data_points.nr_steps_done += 1
        self.new_data_available.emit()

        if self._is_finished():
            self._prepare_finalization()
        else:
            self.start_next_measurement()

    def _connect_controller(self, ctrl):
        """Connect necessary controller signals."""
        super()._connect_controller(ctrl)
        base.safe_connect(self.__request_continuous_mode,
                          ctrl.set_continuous_mode, type=_UNIQUE)
        self.__connect_trigger_timeout(ctrl)

    def __connect_trigger_timeout(self, ctrl):
        """Connect self.__trigger_one_measurement to ctrl."""
        try:
            base.safe_connect(self.__trigger_one_measurement.timeout,
                              ctrl.measure_now, type=_UNIQUE)
        except AttributeError:
            # Not a MeasureController or called via
            # _connect_controller during super().__init__
            pass

    def _disconnect_controller(self, ctrl):
        """Disconnect necessary controller signals."""
        super()._disconnect_controller(ctrl)
        base.safe_disconnect(self.__request_continuous_mode,
                             ctrl.set_continuous_mode)
        try:
            base.safe_disconnect(self.__trigger_one_measurement.timeout,
                                 ctrl.measure_now)
        except AttributeError:
            # Not a MeasureController or called during super().__init__
            pass

    def _finalize(self, *_):
        """Finish the measurement cycle.

        Save data and set energy to zero.

        Parameters
        ----------
        busy : bool
            Needed if the _finalize is called by the controller
            busy state change signal. (continuous measurement mode)

        Returns
        -------
        None.
        """
        if self.is_continuous:
            self.data_points.recalculate_last_step_times()
        super()._finalize()

    def _is_finished(self):
        """Check if the full measurement cycle is done.

        If the desired number of steps/measurements has been reached
        evaluate data and return True. Otherwise generate next energy
        and return False.

        Returns
        -------
        bool
        """
        super()._is_finished()
        if self.current_energy > self.__end_energy:
            return True
        self.current_energy = self.energy_generator()
        return False

    def _make_cameras(self):
        """Make cameras from self.settings, none in continuous."""
        if self.is_continuous:
            # No cameras directly controlled in continuous mode
            self.settings.set('devices', 'cameras', '()')
        super()._make_cameras()

    # Not sure how to handle this appropriately: when continuous,               # TODO
    # we should never enter this at all, as the measurement never
    # has any camera to handle (see reimplemented _make_cameras).
    # In triggered: if we are in the middle of a step (timer active)
    # we should not do anything. If we're at the end of a step and
    # waiting for the last image we may end up with one more image
    # than the other data (if controllers return data after the image).
    # In practice this is not going to happen now, as we're acquiring
    # only a single image per energy step.
    def _on_camera_busy_changed(self, busy):
        """Do nothing when continuous."""
        if busy:
            return
        if self.is_continuous:
            # Should never happen as there is no camera!
            raise RuntimeError("SOMETHING WRONG: THIS SHOULD NEVER HAPPEN")
        if not self.__energy_step_timer.isActive():
            # We are at the end of an energy step, and
            # waiting for the last image to be acquired
            super()._on_camera_busy_changed(busy)

    def _on_controller_data_ready(self, data):
        """Receive measurement data from the controller.

        Append received data to the internal dictionary. Emit an
        error if received dictionary contains a section that does
        not exist in the internal dictionary.

        Parameters
        ----------
        data : dict
            A dictionary containing the measurements. Keys are
            QuantityInfo, values are Sequences (most commonly
            1-element long).

        Returns
        -------
        None.
        """
        controller = self.sender()
        if not isinstance(controller, ControllerABC):
            # This is a safguard, and should never happen,
            # although it did happen for me a couple of times
            # at random (i.e., not reproducibly.)
            base.emit_error(
                self, MeasurementErrors.RUNTIME_ERROR,
                "_on_controller_data_ready got an unexpected sender "
                f"{controller}. (?== self: {self == controller}). "
                "Was expecting a ControllerABC. Energy is "
                f"{self.current_energy}."
                )
            return
        self.data_points.add_data(data, controller)

    def __prepare_continuous_mode(self):
        """Adjust the preparations to fit continuous mode.

        The number of measurements to average over is always 1 in
        continuous measurements. At the end of their preparations
        all controllers have to turn on continuous mode.

        Returns
        -------
        None.
        """
        for ctrl in self.controllers:
            ctrl.settings.set(
                'measurement_settings', 'num_meas_to_average', '1'
                )
            ctrl.continue_prepare_todos['set_continuous_mode'] = True

    def _prepare_finalization(self):
        """Prepare for finalization.

        Connect controller busy to self._finalize, set controllers
        busy and tell them to switch continuous mode off.

        Returns
        -------
        None.

        Emits
        -----
        __request_continuous_mode(False)
            Tell the controller to turn continuous mode off.
        """
        # First disconnect controller_busy from slots that may be
        # inadvertently called in the super() call below
        for ctrl in self.controllers:
            base.safe_disconnect(ctrl.controller_busy,
                                 self.__check_is_finished)
        super()._prepare_finalization()

        # Note: here we are effectively asking to send two commands
        # one after the other. The "change continuous mode" will end
        # in the unsent messages from the controller, which cannot
        # turn "not busy" until all unsent messages are sent. In
        # practice this means that the controllers will turn "not busy"
        # only after the answer to the second command.
        self.__request_continuous_mode.emit(False)

    def _ready_for_next_measurement(self):
        """Start check if all measurements have been received.

        After all measurements have been received, a check if
        the loop is done will be called after the primary
        controller has been stopped.

        Returns
        -------
        None.
        """
        if not self.is_continuous:
            # Stop triggering for this step...
            self.__trigger_one_measurement.stop()
            # ...and check if we can intiate the next one
            super()._ready_for_next_measurement()
            return

        # In continuous mode, we have to explicitly stop the primary
        # controller from spamming us with measurements before we can
        # set a new energy (or conclude).  Force it to be busy (it
        # usually never is while spamming with measurements) to be
        # sure it propery turns "not busy" when it is done stopping.
        primary = self.primary_controller
        primary.busy = True
        primary.controller_busy.connect(self.__check_is_finished, type=_UNIQUE)
        self._request_stop_primary.emit()
