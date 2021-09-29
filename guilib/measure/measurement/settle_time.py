"""Module settle_time of viperleed.
========================================
   ViPErLEED Graphical User Interface
========================================

Created: 2021-07-19
Author: Michele Riva
Author: Florian Doerr

This module contains the definition of the DetermineSettletime class
which gives commands to the controller classes.
"""
import ast

from PyQt5 import QtCore as qtc

# ViPErLEED modules
from measurementabc import MeasurementABC


class DetermineSettletime(MeasurementABC):
    """Settle time determination class."""
    continuous_mode = qtc.pyqtSignal(bool)

    def __init__(self, measurement_settings):
        """Initialise measurement class.

        This is an upgraded version of its parent class.
        __end_energy, __delta_energy and __points_to_take are
        read from the settings and made private properties.
        """
        super().__init__(measurement_settings)
        self.__end_energy = self.settings.getfloat('measurement_settings',
                                                   'end_energy')
        self.__delta_energy = self.settings.getfloat('measurement_settings',
                                                     'delta_energy')
        self.__settle_time = 0
        # Settle time has to be 0 for calibration
        self.__measurement_time = self.settings.getint('measurement_settings',
                                                       'measurement_time')
        self.timer = qtc.QTimer()
        self.timer.setSingleShot(True)
        self.timer.timeout.connect(self.ready_for_next_measurement)
        # TODO: think of something better than setting the most
        # important quantity as a string like below.
        self.thing_to_count = 'HV'

    def start_next_measurement(self):
        """Set energy and measure.

        Set energy via the primary controller. Once this is done
        the returned about_to_trigger signal will start the
        measurement on all secondary controllers.

        Returns
        -------
        None.
        """
        self.primary_controller.busy = True
        for controller in self.secondary_controllers:
            controller.busy = True
        self.connect_continuous_mode_set()
        self.continuous_mode.emit(True)

    def continuous_mode_set(self, busy):
        """Check if continuous mode has been set on all controllers.

        Return if any controller has not been set to continuous
        mode yet. Otherwise set the LEED energy and start the
        measurement.

        Returns
        -------
        None.
        """
        if busy:
            return
        if self.primary_controller.busy:
            return
        if any(controller.busy for controller in self.secondary_controllers):
            return
        self.disconnect_continuous_mode_set()
        self.timer.start(self.__measurement_time)
        self.set_LEED_energy(self.current_energy, self.__settle_time)

    def is_finished(self):
        """Check if the full measurement cycle is done.

        If the number of measurement points has been reached the
        cycle is completed. If not, one is added to the counter
        and the next measurement is started.

        Returns
        -------
        bool
        """
        self.continuous_mode.emit(False)
        if self.current_energy >= self.__end_energy:
            self.on_finished()
            return True
        self.current_energy += self.__delta_energy
        return False

    def on_finished(self):
        """Calculate settling time.

        Takes measured energies and calculates settling time
        from them. Writes settling time to primary controller
        configuration file afterwards."""
        if True:
            print(len(self.data_points['HV']), len(self.data_points['I0']))
            return
        # TODO: Need to do this for multiple steps and take the highest value
        data_points = 0
        step_height = 0.5
        # This step height will be the aimed for height used later in the LEED I(V) video
        int_update_rate = self.primary_controller.settings.getint(
                            'controller', 'update_rate')
        update_rate = self.primary_controller.settings.getint(
                            'adc_update_rate', int_update_rate)
        measured_energies = self.data_points['HV']

        length = len(measured_energies) - 1
        for i, energy in enumerate(measured_energies):
            if abs(measured_energies[length - i] - self.current_energy) < step_height:
                data_points += 1
            else:
                break
        self.__settle_time = 1000*data_points/update_rate
        self.primary_controller.settings.set(
            'measurement_settings', 'settle_time', self.__settle_time)

        file_name = ast.literal_eval(
                        self.settings.get('devices', 'primary_controller')
                        )[0]
        with open(file_name, 'w') as configfile:
            self.primary_controller.settings.write(configfile)
        # May want to run this class for each step and overwrite settle_time only if it is higher than the previous one
        # Need to use continuous mode --> one measurement ever 20 ms at update_rate 4

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
        if self.primary_controller.busy:
            return
        if any(controller.busy for controller in self.secondary_controllers):
            return
        for controller in self.secondary_controllers:
            if controller:
                controller.controller_busy.disconnect()
        self.primary_controller.controller_busy.disconnect()
        if self.__measurement_time <= 1000:
            # Above one second waiting time we will use a timer
            # self.continuous_mode.emit(True)
            pass
        else:
            # TODO: QTimer instantiated here which starts measurements.
            pass
        self.start_next_measurement()

    def connect_cameras(self, cameras=None):
        """Connect necessary camera signals."""
        super().connect_cameras(cameras)
        for camera in cameras:
            if camera:
                # self.continuous_mode.connect(TODO: Function here.,
                                             # type=qtc.Qt.UniqueConnection)
                pass

    def connect_controllers(self, controllers=None):
        """Connect necessary controller signals."""
        super().connect_controllers(controllers)
        for controller in controllers:
            if controller:
                self.continuous_mode.connect(controller.set_continuous_mode,
                                             type=qtc.Qt.UniqueConnection)

    def connect_primary_controller(self):
        """Connect signals of the primary controller."""
        super().connect_primary_controller()
        self.continuous_mode.connect(
            self.primary_controller.set_continuous_mode,
            type=qtc.Qt.UniqueConnection
            )

    def disconnect_cameras(self, cameras):
        """Disconnect necessary camera signals."""
        super().disconnect_cameras(cameras)
        for camera in cameras:
            if camera:
                # self.continuous_mode.disconnect(TODO: Function here.)
                pass

    def disconnect_controllers(self, controllers):
        """Disconnect necessary controller signals."""
        super().disconnect_controllers(controllers)
        for controller in controllers:
            if controller:
                self.continuous_mode.disconnect(controller.set_continuous_mode)

    def disconnect_primary_controller(self):
        """Disconnect signals of the primary controller."""
        super().disconnect_primary_controller()
        if self.primary_controller is not None:
            self.continuous_mode.disconnect(
                self.primary_controller.set_continuous_mode)

    def abort(self):
        """Abort all current actions.

        Abort and reset all variables.

        Returns
        -------
        None.
        """
        super().abort()
        # TODO: add stuff to reset

    def ready_for_next_measurement(self):
        """Check if continuous measurement is done.

        Returns
        -------
        None.
        """
        if self.is_finished():
            self.finalize()
        else:
            self.start_next_measurement()

    def receive_from_camera(self, busy):
        """Do nothing."""
        # TODO: We may want a live stream if we do the energy wiggle.
        pass

    def receive_from_controller(self, receive):
        """Receive measurement data from the controller.

        Append received data to the internal dictionary. Emit an
        error if received dictionary contains a section that does
        not exist in the internal dictionary. After appending all
        measurements check if all of the connected controllers are
        ready for the next measurement.

        Parameters
        ----------
        receive : dictionary
            A dictionary containing the measurements.

        Returns
        -------
        None.

        Emits
        -----
        error_occurred
            If the received measurement contains a label that is
            not specified in the data_points dictionary.
        """
        print(receive)
        for key in receive:
            if key not in self.data_points.keys():
                emit_error(self, MeasurementErrors.INVALID_MEASUREMENT)
            else:
                self.data_points[key].append(receive[key])
        # This has to be '==', do not change to 'is'
        # if any(key == self.thing_to_count for key in receive):
            # self.ready_for_next_measurement()

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
        for controller in self.secondary_controllers:
            if controller:
                controller.controller_busy.connect(
                    self.continuous_mode_set,
                    type=qtc.Qt.UniqueConnection
                    )
        self.primary_controller.controller_busy.connect(
            self.continuous_mode_set,
            type=qtc.Qt.UniqueConnection
            )

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
        for controller in self.secondary_controllers:
            if controller:
                controller.controller_busy.disconnect()
        self.primary_controller.controller_busy.disconnect()
