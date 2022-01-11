"""Module time_resolved of viperleed.guilib.measure.measurement.

========================================
   ViPErLEED Graphical User Interface
========================================

Created: 2021-07-19
Author: Michele Riva
Author: Florian Doerr

This module contains the definition of the TimeResolved class
which gives commands to the controller classes.
"""
import ast

from PyQt5 import QtCore as qtc

# ViPErLEED modules
from viperleed.guilib.measure.measurement.abc import (MeasurementABC,
                                                      MeasurementErrors)
from viperleed.guilib.measure.datapoints import QuantityInfo


class TimeResolved(MeasurementABC):
    """Time resolved measurement class."""
    continuous_mode = qtc.pyqtSignal(bool)
    display_name = 'Time resolved'

    def __init__(self, measurement_settings):
        """Initialise measurement class."""
        super().__init__(measurement_settings)
        self.__settle_time = 0
        # Settle time has to be 0 for calibration
        self.__time_over = False
        self.__end_energy = 0
        self.__delta_energy = 1
        self.__endless = False
        self.__measurement_time = 0
        self.__constant_energy = False
        self.__limit_continuous = 0
        self.__cycle_time = 0
        if self.settings:
            self.__delta_energy = self.settings.getfloat(
                'measurement_settings', 'delta_energy', fallback=10
                )
            self.__end_energy = self.settings.getfloat(
                'measurement_settings', 'end_energy', fallback=10
                )
            self.__endless = self.settings.getboolean(
                'measurement_settings', 'endless', fallback=False
                )
            self.__measurement_time = self.settings.getint(
                'measurement_settings', 'measurement_time', fallback=100
                )
            self.__constant_energy = self.settings.getboolean(
                'measurement_settings', 'constant_energy', fallback=False
                )
            self.__limit_continuous = self.settings.getint(
                'measurement_settings', 'limit_continuous', fallback=10
                )
            self.__cycle_time = self.settings.getint(
                'measurement_settings', 'cycle_time', fallback=100
                )

        self.timer = qtc.QTimer(parent=self)
        self.timer.setSingleShot(True)
        self.timer.timeout.connect(self.ready_for_next_measurement)

        if self.__cycle_time > 0:
            self.cycle_timer = qtc.QTimer(parent=self)
            self.cycle_timer.setSingleShot(True)
            self.cycle_timer.timeout.connect(self.__set_time_over)

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
        if self.__measurement_time <= self.__limit_continuous:
            for controller in self.controllers:
                # Necessary to force secondaries into busy,
                # before the primary returns not busy anymore.
                controller.busy = True
            self.connect_continuous_mode_set()
            self.continuous_mode.emit(True)
        else:
            self.timer.start(self.__measurement_time)
            self.set_LEED_energy(self.current_energy, self.__settle_time)

    def continuous_mode_set(self, busy):
        """Check if continuous mode has been set on all controllers.

        Return if any controller has not been set to continuous
        mode yet. Otherwise set the LEED energy and start the
        measurement.

        Returns
        -------
        None.
        """
        if busy or any(controller.busy for controller in self.controllers):
            return
        self.disconnect_continuous_mode_set()
        self.timer.start(self.__measurement_time)
        self.set_LEED_energy(self.current_energy, self.__settle_time)

    def is_finished(self):
        """Check if the full measurement cycle is done.

        If the desired number of steps/measurements has been reached
        evaluate data and return True. Otherwise generate next energy
        and return False.

        Returns
        -------
        bool
        """
        self.new_data_available.emit()
        if self.__time_over or self.current_energy >= self.__end_energy:
            self.on_finished()
            return True
        self.current_energy = self.energy_generator()
        return False

    def on_finished(self):
        """Calculate settling time.

        Takes measured energies and calculates settling time
        from them. Writes settling time to primary controller
        configuration file afterwards.

        Returns
        -------
        None.
        """
        # TODO: This will either be moved to a subclass or a separate processor
        # TODO: currently using nominal energy on an uncalibrated energy measurement: offset might be larger than step height!!!
        quantity = self.settings.get(
            'measurement_settings', 'measure_this', fallback='None'
            )
        percentage = self.settings.getfloat('measurement_settings',
                                            'percentage', fallback=0.1)
        points = self.settings.getint(
            'measurement_settings', 'relevant_points', fallback=10
            )
        if quantity == QuantityInfo.HV.label:
            to_change = 'hv_settle_time'
        elif quantity == QuantityInfo.I0.label:
            to_change = 'i0_settle_time'
        else:
            print('not one of the expected quantities measured')
            # TODO: emit error
            return

        measured, _ = self.data_points.get_time_resolved_data(
            quantity, separate_steps=True, absolute_times=True
            )
        interval = self.primary_controller.measurement_interval

        if self.__measurement_time <= self.__limit_continuous:
            for j, step in enumerate(measured[0]):
                if j == 0:
                    previous_height = sum(step[-points:])/points
                    continue
                current_height = sum(step[-points:])/points
                max_deviation = (current_height - previous_height) * percentage
                points = 0
                for measurement in reversed(step):
                    if abs(measurement - current_height) < max_deviation:
                        points += 1
                    else:
                        break
                if points == 0:
                    points = 1
                previous_height = current_height
                settle_time = int(1000*points*interval)
                if settle_time > self.__settle_time:
                    self.__settle_time = settle_time

            self.primary_controller.settings.set(
                'measurement_settings', to_change, str(self.__settle_time))
            file_name = ast.literal_eval(
                            self.settings.get('devices', 'primary_controller')
                            )[0]
            with open(file_name, 'w') as configfile:
                self.primary_controller.settings.write(configfile)
        else:
            pass
            # TODO: save I0 (or I00) or only images for wiggle

    def is_preparation_finished(self):
        """Check if measurement preparation is done.

        Whenever a controller is done with its preparation
        this function will be called. The busy attribute is
        used to check if all controllers are done with their
        preparation. If this is true, the busy signal going to
        this function will be disconnected and the first
        measurement will immediately be started.

        Returns
        -------
        None.
        """
        if any(controller.busy for controller in self.controllers):
            return
        for controller in self.controllers:
            controller.controller_busy.disconnect()
        if self.__cycle_time > 0:
            self.cycle_timer.start(self.__cycle_time)
        self.prepared.emit()
        self.start_next_measurement()

    def connect_cameras(self):
        """Connect necessary camera signals."""
        for camera in self.cameras:
            # self.continuous_mode.connect(TODO: Function here.,
                                         # type=qtc.Qt.UniqueConnection)
            pass
        super().connect_cameras()

    def connect_secondary_controllers(self):
        """Connect necessary controller signals."""
        for controller in self.secondary_controllers:
            if not controller:
                continue
            self.continuous_mode.connect(controller.set_continuous_mode,
                                         type=qtc.Qt.UniqueConnection)
        super().connect_secondary_controllers()

    def connect_primary_controller(self):
        """Connect signals of the primary controller."""
        self.continuous_mode.connect(
            self.primary_controller.set_continuous_mode,
            type=qtc.Qt.UniqueConnection
            )
        super().connect_primary_controller()

    def disconnect_cameras(self):
        """Disconnect necessary camera signals."""
        super().disconnect_cameras()
        for camera in self.cameras:
            # self.continuous_mode.disconnect(TODO: Function here.)
            pass

    def disconnect_secondary_controllers(self):
        """Disconnect necessary controller signals."""
        super().disconnect_secondary_controllers()
        for controller in self.secondary_controllers:
            if not controller:
                continue
            try:
                self.continuous_mode.disconnect(controller.set_continuous_mode)
            except TypeError:
                pass

    def disconnect_primary_controller(self):
        """Disconnect signals of the primary controller."""
        super().disconnect_primary_controller()
        if self.primary_controller is not None:
            try:
                self.continuous_mode.disconnect(
                    self.primary_controller.set_continuous_mode)
            except TypeError:
                pass

    def abort(self):
        """Abort all current actions.

        Abort and reset all variables.

        Returns
        -------
        None.
        """
        self.timer.stop()
        super().abort()

    def receive_from_camera(self, busy):
        """Do nothing."""
        # TODO: We may want a live stream if we do the energy wiggle.
        pass

    def receive_from_controller(self, controller, receive):
        """Receive measurement data from the controller.

        Append received data to the internal dictionary. Emit an
        error if received dictionary contains a section that does
        not exist in the internal dictionary.

        Parameters
        ----------
        controller : ControllerABC
            The controller object that is sending data.
        receive : dictionary
            A dictionary containing the measurements.

        Returns
        -------
        None.

        Emits
        -----
        error_occurred
            If the received measurement contains a label that is
            not specified in the data_points[-1] dictionary.
        """
        if controller == self.primary_controller:
            self.data_points.add_data(receive, controller, self.primary_delay)
        else:
            self.data_points.add_data(receive, controller)

    def connect_continuous_mode_set(self):
        """Connect controller busy signal to continuous_mode_set.

        Use the controller_busy signal to see if all controllers
        have been set to continuous mode. Has to be disconnected
        once all controllers have been set to continuous mode,
        see disconnect_continuous_mode_set function.

        Returns
        -------
        None.
        """
        for controller in self.controllers:
            controller.controller_busy.connect(self.continuous_mode_set,
                                               type=qtc.Qt.UniqueConnection)

    def disconnect_continuous_mode_set(self):
        """Connect controller busy signal to continuous_mode_set.

        The controller_busy signal has to be disconnected once
        all controllers have been set to continuous mode,
        otherwise the signal from setting the voltage will start
        an infinite loop.

        Returns
        -------
        None.
        """
        for controller in self.controllers:
            controller.controller_busy.disconnect()

    def energy_generator(self):
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
            if energy >= self.__end_energy:
                energy = self.start_energy
        return energy

    def __set_time_over(self):
        """Set __time_over to True

        Serves to end a measurement after a specific time.

        Returns
        -------
        None.
        """
        self.__time_over = True

    def prepare_finalization(self):
        """Prepare for finalization.

        Connect controller busy to self.finalize, set controllers
        busy and tell them to switch continuous mode off.

        Returns
        -------
        None.

        Emits
        -----
        continuous_mode(False)
            Tell the controller to turn continuous mode off.
        """
        for controller in self.controllers:
            try:
                controller.controller_busy.disconnect()
            except TypeError:
                pass
            # Necessary to force secondaries into busy,
            # before the primary returns not busy anymore.
            controller.busy = True
            controller.controller_busy.connect(self.finalize,
                                               type=qtc.Qt.UniqueConnection)
        self.continuous_mode.emit(False)
        self.abort_action.emit()

    def ready_for_next_measurement(self):
        """Start check if all measurements have been received.

        After all measurements have been received, a check if
        the loop is done will be called after the continuous
        mode has been turned off on the primary controller.

        Returns
        -------
        None.
        """
        for controller in self.controllers:
            # Necessary to force secondaries into busy,
            # before the primary returns not busy anymore.
            controller.busy = True
            controller.controller_busy.connect(self.check_is_finished,
                                               type=qtc.Qt.UniqueConnection)
        if self.__measurement_time <= self.__limit_continuous:
            self.continuous_mode.emit(False)

    def check_is_finished(self):
        """Check if the full measurement is finished."""
        if any(device.busy for device in self.devices):
            return
        for controller in self.controllers:
            controller.controller_busy.disconnect()
        self.data_points.calculate_times()
        if self.is_finished():
            self.prepare_finalization()
        else:
            self.start_next_measurement()

