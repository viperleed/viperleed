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
from configparser import ConfigParser

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

        num_meas = (1 + round((self.__end_energy - self.start_energy)
                   / self.__delta_energy))
        if self.is_continuous_measurement:
            self.prepare_continuous_mode()
            self.data_points.num_measurements = num_meas
        else:
            self.__hv_settle_time = 0
            self.__hv_settle_time = self.primary_controller.settings.getint(
                'measurement_settings', 'hv_settle_time', fallback=0
                )
            self.__n_digits = len(str(num_meas))
        self.timer = qtc.QTimer(parent=self)
        self.timer.setSingleShot(True)
        self.timer.timeout.connect(self.ready_for_next_measurement)

        if self.__cycle_time > 0:
            self.cycle_timer = qtc.QTimer(parent=self)
            self.cycle_timer.setSingleShot(True)
            self.cycle_timer.timeout.connect(self.__set_time_over)

    @property
    def is_continuous_measurement(self):
        """Return whether the measurement is continuous."""
        return self.__measurement_time <= self.__limit_continuous

    def prepare_continuous_mode(self):
        """Adjust the preparations to fit continuous mode.

        The number of measurements to average over is always
        1 in continuous measurements. At the end of their
        preparations all controllers have to turn on
        continuous mode.
        """
        for ctrl in self.controllers:
            ctrl.settings.set(
                'measurement_settings', 'num_meas_to_average', '1'
                )
            ctrl.continue_prepare_todos['set_continuous_mode'] = True

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
        self.timer.start(self.__measurement_time)
        self.set_LEED_energy(self.current_energy, self.__settle_time)
        if not self.is_continuous_measurement:
            # TODO: decide if we really want to save images (flag)
            self.counter += 1
            image_name = (f"{self.counter:0>{self.__n_digits}}_"
                          f"{self.current_energy:.1f}eV_.tiff")
            for i, camera in enumerate(self.cameras):
                camera.process_info.filename = image_name
                self.data_points.add_image_names(image_name)
            self.camera_timer.start(self.__hv_settle_time)

    def is_finished(self):
        """Check if the full measurement cycle is done.

        If the desired number of steps/measurements has been reached
        evaluate data and return True. Otherwise generate next energy
        and return False.

        Returns
        -------
        bool
        """
        # self.new_data_available.emit()
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
        # TODO: calculation will be moved to datapoints class
        # TODO: currently using nominal energy on an uncalibrated energy measurement: offset might be larger than step height!!!
        # quantity = self.settings.get(
            # 'measurement_settings', 'measure_this', fallback='None'          measure_this not in config anymore
            # )
        return
        quantity_obj = QuantityInfo.from_label(quantity)
        percentage = self.settings.getfloat('measurement_settings',
                                            'percentage', fallback=0.1)
        points = self.settings.getint(
            'measurement_settings', 'relevant_points', fallback=10
            )
        if quantity_obj == QuantityInfo.HV:
            to_change = 'hv_settle_time'
        elif quantity_obj == QuantityInfo.I0:
            to_change = 'i0_settle_time'
        else:
            print('not one of the expected quantities measured')
            # TODO: emit error
            return
        # TODO:                             get_time_resolved_data currently throws an index error
        # measured, _ = self.data_points.get_time_resolved_data(
            # quantity_obj, separate_steps=True, absolute_times=True
            # )
        interval = self.primary_controller.measurement_interval

        if self.is_continuous_measurement:
            for ctrl in self.controllers:
                del ctrl.continue_prepare_todos['set_continuous_mode']

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
            # Create a dummy of the settings in order to keep the
            # original number of measurements to average over.
            orig_settings = ConfigParser(comment_prefixes='/',
                                         allow_no_value=True)
            orig_settings.read(file_name)
            orig_settings.set(
                'measurement_settings', to_change, str(self.__settle_time))
            with open(file_name, 'w') as configfile:
                orig_settings.write(configfile)
        else:
            pass
            # TODO: save I0 (or I00) or only images for wiggle

    def is_preparation_finished(self):
        """Check if measurement preparation is done.
        """
        # TODO: remove this once the energy generators have the ability to do the same energy multiple times.
        if any(device.busy for device in self.devices):
            return
        if self.__cycle_time > 0:
            self.cycle_timer.start(self.__cycle_time)
        super().is_preparation_finished()

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
        return

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
        self.request_stop_devices.emit()

    def ready_for_next_measurement(self):
        """Start check if all measurements have been received.

        After all measurements have been received, a check if
        the loop is done will be called after the primary
        controller has been stopped.

        Returns
        -------
        None.
        """
        if self.is_continuous_measurement:
            self.primary_controller.busy = True
            self.primary_controller.controller_busy.connect(
                self.check_is_finished, type=qtc.Qt.UniqueConnection
                )
            self.request_stop_primary.emit()
        else:
            self.check_is_finished()

    def check_is_finished(self):
        """Check if the full measurement is finished."""
        if self.is_continuous_measurement:
            if any(device.busy for device in self.devices):
                return
            self.primary_controller.controller_busy.disconnect()
            self.data_points.calculate_times(continuous=True)
        else:
            self.data_points.calculate_times()
        self.new_data_available.emit()
        if self.is_finished():
            self.prepare_finalization()
        else:
            self.start_next_measurement()

    def do_next_measurement(self):
        """Do the next measurement.

        Start measuring on all secondary controllers if the
        measurement is not continuous.
        """
        if self.is_continuous_measurement:
            self.primary_controller.about_to_trigger.disconnect(
                self.do_next_measurement
                )
        self.ready_for_measurement.emit()

    def finalize(self, busy=False):
        """Finish the measurement cycle.

        Save data and set energy to zero.

        Parameters
        ----------
        busy : bool
            Needed if the finalize is called by the controller
            busy state change signal. (continuous measurement mode)

        Returns
        -------
        None.
        """
        if self.is_continuous_measurement:
            self.data_points.recalculate_last_setp_times()
        super().finalize(busy=busy)

