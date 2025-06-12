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

from math import ceil

from PyQt5 import QtCore as qtc

from viperleed.gui.measure import hardwarebase as base
from viperleed.gui.measure.classes.abc import QObjectSettingsErrors
from viperleed.gui.measure.classes.datapoints import QuantityInfo
from viperleed.gui.measure.measurement.abc import MeasurementABC


_INVOKE = qtc.QMetaObject.invokeMethod
_QUEUED = qtc.Qt.QueuedConnection
_UNIQUE = qtc.Qt.UniqueConnection


class TimeResolved(MeasurementABC):  # too-many-instance-attributes
    """Time resolved measurement class."""

    display_name = 'Time resolved'

    _mandatory_settings = (
        *MeasurementABC._mandatory_settings,
        ('measurement_settings', 'is_continuous'),
        ('measurement_settings', 'energy_step_duration'),
        )

    __request_continuous_mode = qtc.pyqtSignal(bool)   # On/Off                 # TODO: could be done with QMetaObject.invokeMethod

    def __init__(self, measurement_settings):
        """Initialise measurement class."""
        super().__init__(measurement_settings)
        self.__end_energy = 0                                                   # --> generator; It can happen that there is no proper "end energy" (endless measurement, e.g., wiggle)
        self.__delta_energy = 1                                                 # --> generator; There may not be a delta_energy (e.g., I(t) at fixed energy)
        self.__endless = False                                                  # --> generator;
        self.__constant_energy = False                                          # --> generator;

        # One timer to trigger a change of electron energy: we will
        # sit at each energy for self.energy_step_duration millisecs.
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
        trigger = self.__trigger_one_measurement = qtc.QTimer(parent=self)
        trigger.setSingleShot(True)
        trigger.timeout.connect(self.__on_one_measurement_triggered)

        # A third timer that attempts to trigger again if
        # __trigger_one_measurement failed to do so.
        self._retry_triggering = qtc.QTimer(parent=self)
        self._retry_triggering.setSingleShot(True)
        self._retry_triggering.setInterval(10)
        self._retry_triggering.timeout.connect(
            self.__on_one_measurement_triggered
            )

        # Finally, set the _camera_timer interval to zero, so
        # we can fire it at the same time as measurements are
        # acquired (in triggered mode).
        self._camera_timer.setInterval(0)

        # TODO: use a flag to decide if we want to save images.
        # If not, we can just disconnect _camera_timer, and we
        # should also not invoke in __on_one_measurement_triggered

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

        _continuous = self.is_continuous
        self.data_points.time_resolved = True
        self.data_points.continuous = _continuous

        num_meas = (1 + round((self.__end_energy - self.start_energy)
                              / self.__delta_energy))
        if _continuous:
            self.__prepare_continuous_mode()
            self.data_points.nr_steps_total = num_meas                          # TODO: incorrect for endless and constant energy. Also, unused now that we use different DataPoints for plotting and measurement.
        else:
            # With the next connection, triggering the first
            # measurement at the current energy also starts
            # __trigger_one_measurement: This triggers more
            # measurements at .measurement_interval intervals
            about_to_trigger = self.primary_controller.about_to_trigger
            about_to_trigger.connect(trigger.start)

            # The first measurement should also start the camera
            about_to_trigger.connect(self._camera_timer.start)

    @property
    def _n_digits(self):
        """Return the appropriate number of digits for padding image names."""
        # Used for zero-padding counter in image names.
        num_meas = (1 + round((self.__end_energy - self.start_energy)           # TODO: get it from generator when implemented correctly
                              / self.__delta_energy))
        num_meas *= ceil(self.energy_step_duration / self.measurement_interval)
        return len(str(num_meas))

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
            base.emit_error(self, QObjectSettingsErrors.INVALID_SETTINGS,
                            'measurement_settings/energy_step_duration', '')
            interval = min_t

        if interval < min_t:
            base.emit_error(
                self, QObjectSettingsErrors.INVALID_SETTING_WITH_FALLBACK,
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
        # The minimum time the slowest device needs to return data.
        return self.measurement_interval_min

    @property
    def is_continuous(self):
        """Return whether the measurement is continuous."""
        try:
            return self.settings.getboolean('measurement_settings',
                                            'is_continuous')
        except ValueError:
            # Not a valid boolean
            base.emit_error(self, QObjectSettingsErrors.INVALID_SETTINGS,
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
        try:
            interval = self.settings.getint('measurement_settings',
                                            'measurement_interval',
                                            fallback=min_t)
        except (TypeError, ValueError):
            # Not an int
            base.emit_error(self, QObjectSettingsErrors.INVALID_SETTINGS,
                            'measurement_settings/measurement_interval', '')
            interval = min_t

        if interval < min_t:
            txt = f"{interval} (too short)"
            interval = min_t
            base.emit_error(
                self, QObjectSettingsErrors.INVALID_SETTING_WITH_FALLBACK,
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
        if new_interval < min_t:
            raise ValueError(f"{new_interval} out of bounds. "
                             f"Should be at least {min_t}")
        self.settings.set('measurement_settings', 'measurement_interval',
                          str(new_interval))
        self.__trigger_one_measurement.setInterval(new_interval)

    @property
    def measurement_interval_min(self):
        """Return the smallest measurement-time interval."""
        _min = 50
        if not self.devices:
            return _min
        min_interval_cam = min_interval_ctrl = _min

        # Controllers:
        measuring_ctrls = [c for c in self.controllers if c.measures()]
        if any(measuring_ctrls):
            min_interval_ctrl = max(c.time_to_first_measurement
                                    for c in measuring_ctrls)
            # Add a little, in case delay % measurement_interval == 0
            min_interval_ctrl += 5

        # Cameras:
        if self.cameras:
            min_interval_cam = max(cam.time_to_image_ready
                                   for cam in self.cameras)
            min_interval_cam += 5  # Same as controller

        return max(_min, round(min_interval_ctrl), round(min_interval_cam))

    @qtc.pyqtSlot()
    def abort(self):
        """Interrupt measurement, save partial data, set energy to zero."""
        try:
            timers = (
                self.__trigger_one_measurement,
                self.__energy_step_timer,
                self._retry_triggering,
                )
        except AttributeError:
            # .abort() happened during super().__init__
            timers = tuple()

        for timer in timers:
            timer.stop()
        super().abort()

    def energy_generator(self):                                                 # TODO: move to parent; improve.
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
        ramp_finished = energy > self.__end_energy
        if self.__delta_energy < 0:
            ramp_finished = energy < self.__end_energy

        if self.__endless:
            self.new_data_available.emit(self.data_points[-1])
            if ramp_finished:
                energy = self.start_energy
        return energy

    @qtc.pyqtSlot(object)
    def set_settings(self, new_settings):
        """Change settings of the measurement."""
        settings_ok = super().set_settings(new_settings)

        # Since we should not wait for data to be stored to move
        # to the next step, invalidate the dictionary of super here
        self._data_stored = {}
        return settings_ok

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
        _continuous = self.is_continuous
        if _continuous:
            self._missing_data = dict.fromkeys(self._missing_data.keys(), 0)
        self.__energy_step_timer.setInterval(self.energy_step_duration)
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

        # Pick a different settle time for continuous and triggered:
        # in continuous mode we want to get measurements as quickly as
        # possible (and look at short-term time traces); in "triggered"
        # mode give the user freedom to choose by using hv_settle_time          # TODO: should we use something else, perhaps?
        settle_time = 0 if _continuous else self.hv_settle_time

        # Always trigger measurements. For a "triggered" measurement
        # type the about_to_trigger emitted when the first measurement
        # is requested also causes __trigger_one_measurement to start.
        # This will keep firing at .measurement_interval, calling
        # ctrl.measure_now() and camera.trigger_now() each time
        self.set_leed_energy(*self.step_profile,
                             self.current_energy, settle_time)

        if _continuous:
            return
        self.__trigger_one_measurement.setInterval(self.measurement_interval)

    @qtc.pyqtSlot(bool)
    def __check_is_finished(self, _):                                           # TODO: I'm not happy with the name. __continue_when_primary_stopped?
        """Check if the measurement is finished in continuous mode.

        This method is used only if self.is_continuous. It does the
        same checks as in super()._ready_for_next_measurement(), i.e.,
        decides whether we should go to the next step or if the whole
        loop is over. It is needed in continuous mode so we can
        wait for the acknowledgment from the primary controller that
        it has been in fact stopped at the end of an energy step.

        Returns
        -------
        None.
        """
        base.safe_disconnect(self.primary_controller.busy_changed,
                             self.__check_is_finished)
        if self.aborted:
            # We entered this call after the measurement was aborted,
            # likely while processing an unprocessed timeout event
            return
        super()._ready_for_next_measurement()

    def _connect_controller(self, ctrl):
        """Connect necessary controller signals."""
        super()._connect_controller(ctrl)
        try:
            base.safe_connect(self.__request_continuous_mode,
                              ctrl.set_continuous_mode, type=_UNIQUE)
        except AttributeError:
            # Not a MeasureControllerABC or
            # called during super().__init__
            pass

    def _disconnect_controller(self, ctrl):
        """Disconnect necessary controller signals."""
        super()._disconnect_controller(ctrl)
        # The next disconnect makes sense only for the primary
        # controller, but does not hurt to try it for all
        try:
            base.safe_disconnect(ctrl.about_to_trigger,
                                 self.__trigger_one_measurement.start)
        except AttributeError:
            # Called during super().__init__(). No need to try others
            return

        if not ctrl.measures():
            return

        base.safe_disconnect(self.__request_continuous_mode,
                             ctrl.set_continuous_mode)

    @qtc.pyqtSlot(bool)
    @qtc.pyqtSlot()
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
            # Recalculate the times in the last data point, in case
            # secondary controllers returned more data in the meantime
            # Notice that we do not allow any errors to be emitted at
            # this point, since we may end up in an infinite loop.
            self.data_points.calculate_times(complain=False)
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
        ramp_finished = self.energy_generator() > self.__end_energy
        if self.__delta_energy < 0:
            ramp_finished = self.current_energy < self.__end_energy
        if ramp_finished:
            return True
        self.current_energy = self.energy_generator()
        return False

    def _make_cameras(self):
        """Make cameras from self.settings, none in continuous."""
        if self.is_continuous:
            # No cameras directly controlled in continuous mode
            self.settings.set('devices', 'cameras', '()')
        super()._make_cameras()

    @qtc.pyqtSlot(bool)
    def _on_camera_busy_changed(self, busy):
        """Go to next energy if time is over (in triggered mode)."""
        if busy:
            # Just triggered
            return

        # Collected all frames, and will save a processed image later
        camera = self.sender()
        self.data_points.add_image(camera)
        camera.process_info.count += 1

        if self.is_continuous:
            # Does currently not happen as cameras are discarded
            # in continuous mode. If we ever support this we should
            # also add time stamps for camera images.
            return
        # We must make sure that we have already stored in the
        # datapoints all the images that we have ever triggered for.
        self._missing_data[camera] -= 1
        if not self.__energy_step_timer.isActive():
            # We are at the end of an energy step, and just finished
            # waiting for the last image. See if we can go on.
            self._ready_for_next_measurement()

    @qtc.pyqtSlot(dict)
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
        self.data_points.add_data(data, controller)
        if self.is_continuous:
            return
        self._missing_data[controller] -= 1
        if not self.__energy_step_timer.isActive():
            # We are at the end of an energy step, and just finished
            # waiting for the last data. See if we can go on.
            self._ready_for_next_measurement()

    @qtc.pyqtSlot()
    def __on_one_measurement_triggered(self):
        """Increment the number of missing data for all devices."""
        if self.is_continuous or self.aborted:
            # Also return if aborted because the _retry_triggering timer
            # might time out before the abort method stops it, which
            # would lead to this method being called after an abort.
            return
        if any(self._missing_data.values()):
            self._retry_triggering.start()
            return

        # Notice that we always invoke with a QueuedConnection (which is
        # the default for objects in another thread) even if the object
        # is in the same thread. This is to ensure that the next portion
        # of code runs as fast as possible, and does not stall in calling
        # methods in the same thread as the measurement
        remaining_time = self.__energy_step_timer.remainingTime()
        for ctrl in self.controllers:
            if not ctrl.measures():
                continue
            if remaining_time < ctrl.time_to_first_measurement:
                continue
            _INVOKE(ctrl, 'measure_now', _QUEUED)
            self._missing_data[ctrl] += 1
        for camera in self.cameras:
            if remaining_time < camera.time_to_image_ready:
                continue
            _INVOKE(camera, 'trigger_now', _QUEUED)
            self._missing_data[camera] += 1
        self.__trigger_one_measurement.start()

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
            if not ctrl.measures():
                continue
            ctrl.settings.set('measurement_settings', 'nr_samples', '1')
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
        # First disconnect controller busy_changed from slots that may
        # be inadvertently called in the super() call below.
        for ctrl in self.controllers:
            base.safe_disconnect(ctrl.busy_changed,
                                 self.__check_is_finished)

        # Disconnect other signals that we will not need any longer.
        if self.primary_controller:
            about_to_trigger = self.primary_controller.about_to_trigger
            base.safe_disconnect(about_to_trigger, self._camera_timer.start)
            try:
                trigger = self.__trigger_one_measurement
            except AttributeError:
                # Error during super().__init__
                pass
            else:
                base.safe_disconnect(about_to_trigger, trigger.start)
                base.safe_disconnect(trigger.timeout,
                                     self.__on_one_measurement_triggered)

        super()._prepare_finalization()

        # Note: here we are effectively asking to send two commands
        # one after the other: stop() in super() call, "change mode"
        # below. The second one will end in the unsent messages from
        # the controller, which cannot turn "not busy" till all unsent
        # messages are sent. This means that the controllers will turn
        # "not busy" only after the answer to the second command.
        self.__request_continuous_mode.emit(False)

    @qtc.pyqtSlot()
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
            # ...and check if we can initiate the next one
            super()._ready_for_next_measurement()
            return

        # In continuous mode, we have to explicitly stop the primary
        # controller from spamming us with measurements before we can
        # set a new energy (or conclude).  Force it to be busy (it
        # usually never is while spamming with measurements) to be
        # sure it properly turns "not busy" when it is done stopping.
        primary = self.primary_controller
        primary.busy = True
        primary.busy_changed.connect(self.__check_is_finished, type=_UNIQUE)
        primary.stop()
