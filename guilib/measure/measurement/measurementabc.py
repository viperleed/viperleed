"""Module measurementabc of viperleed.
========================================
   ViPErLEED Graphical User Interface
========================================

Created: 2021-07-19
Author: Michele Riva
Author: Florian Doerr

This module contains the definition of the measurementABC class
which gives commands to the controller classes and its associated
ViPErLEEDErrorEnum class MeasurementErrors.
"""
from collections import defaultdict
from abc import abstractmethod
import csv
import ast
from configparser import ConfigParser
# from time import localtime, strftime

from PyQt5 import QtCore as qtc

# ViPErLEED modules
# from viperleed.guilib.measure.hardwarebase import (
    # emit_error, ViPErLEEDErrorEnum, QMetaABC,
    # config_has_sections_and_options, class_from_name
    # )

# test
import os, sys

p = os.path.abspath('C:/Users/Florian/Documents/Uni/Masterarbeit/ViperLEED/')
sys.path.insert(1, p)
from viperleed.guilib.measure.hardwarebase import (
    emit_error, ViPErLEEDErrorEnum, QMetaABC,
    config_has_sections_and_options, class_from_name
    )


class MeasurementErrors(ViPErLEEDErrorEnum):
    """Errors that might occur during a measurement cycle."""
    INVALID_MEASUREMENT = (300,
                           "The returned data dictionary contained a section "
                           "that was not specified in the measurement class.")
    MISSING_SETTINGS = (301,
                        "Measurements cannot be taken without settings. "
                        "Load an appropriate settings file before "
                        "proceeding.")
    INVALID_MEAS_SETTINGS = (302,
                             "Invalid measurement settings: Required "
                             "settings {!r} missing or wrong. Check "
                             "configuration file.")
    MISSING_CLASS_NAME = (303,
                          "{!r} is missing the name "
                          "of its related class object.")


class MeasurementABC(qtc.QObject, metaclass=QMetaABC):
    """Generic measurement class.

    The plot_info dictionary in this class may be reimplemented
    in subclasses. Each section (label) needs to contain
    the unit and the scaling ('lin' for linear and 'log' for
    logarithmic scaling.
    """
    # Is emitted if a measurement cycle has been completed.
    # Contains measurement data as a tuple of dictionaries.
    finished = qtc.pyqtSignal(tuple)
    # Is emitted when an error occurs
    error_occurred = qtc.pyqtSignal(tuple)
    # Start measurement on all controllers
    # (after receiving about_to_trigger from primary controller)
    ready_for_measurement = qtc.pyqtSignal()
    # Start first part of the preparation for measurements
    # Is emitted before starting energy is set.
    begin_preparation = qtc.pyqtSignal(tuple)
    # Start second part of the preparation for measurements
    # Is emitted after the starting energy is set.
    continue_preparation = qtc.pyqtSignal()
    # Abort current task on the controller side.
    # Is emitted if measurement is aborted.
    abort_action = qtc.pyqtSignal()

    # The reimplementation may introduce more/other keys.
    # See ViPErinoController for an example on how to do this.
    plot_info = defaultdict(list)
    plot_info['nominal_energy'] = ['eV', 'lin']
    plot_info['I0'] = ['uA', 'lin']
    plot_info['HV'] = ['eV', 'lin']
    plot_info['elapsed_time'] = ['ms', 'lin']
    plot_info['Isample'] = ['V', 'log']
    plot_info['temperature'] = ['°C', 'lin']
    plot_info['cold_junction'] = ['°C', 'lin']

    _mandatory_settings = [
        ('devices', 'primary_controller'),
        ('measurement_settings', 'start_energy')
    ]

    def __init__(self, measurement_settings):
        """Initialise measurement class"""

        super().__init__()

        self.current_energy = 0
        self.__settings = None
        self.__primary_controller = None
        self.__secondary_controllers = []
        self.__cameras = []
        self.__start_energy = 0
        self.thread = qtc.QThread()
        # The reimplementation may introduce more/other keys.
        self.data_points = defaultdict(list)
        for key in self.plot_info:
            self.data_points[key] = []

        self.set_settings(measurement_settings)

        self.__start_energy = self.settings.getfloat('measurement_settings',
                                                     'start_energy')
        self.__settle_time = self.primary_controller.settings.getint(
            'measurement_settings', 'settle_time')
        self.__long_settle_time = self.primary_controller.settings.getint(
            'measurement_settings', 'first_settle_time')

    @property
    def cameras(self):
        """Return the cameras used by this class."""
        return self.__cameras

    @cameras.setter
    def cameras(self, new_cameras):
        """Set the cameras which should be
        used and handle signals.

        Parameters
        ----------
        new_cameras : list
            List of camera class objects.
        """
        self.disconnect_cameras(self.__cameras)
        self.__cameras = new_cameras
        self.connect_cameras(self.__cameras)

    def __get_settings(self):
        """Return the current settings for the measurement.

        Returns
        -------
        settings : ConfigParser
            The ConfigParser containing the measurement settings.
        """
        return self.__settings

    def set_settings(self, new_settings):
        """Change settings of the measurement.

        Settings are loaded only if they are valid. Otherwise
        the previous settings stay in effect. If the settings
        have been accepted, controller and camera objects as
        specified in the settings will be instantiated, told
        what they will be measuring, moved to their respective
        properties and connected to all necessary signals.

        Parameters
        ----------
        new_settings : dict, ConfigParser, string or path
            Configuration of the measurement.

        Raises
        ------
        TypeError
            If new_settings is neither a dict, ConfigParser, string
            or path and if an element of the mandatory_settings is
            None or has a length greater than 3.

        Emits
        -----
        ExtraSerialErrors.MISSING_SETTINGS
            If new_settings is missing.
        ExtraSerialErrors.INVALID_MEAS_SETTINGS
            If any element of the new_settings does not fit the
            mandatory_settings.
        """
        if new_settings is None:
            emit_error(self, MeasurementErrors.MISSING_SETTINGS)
            return

        new_settings, invalid = config_has_sections_and_options(
            self,
            new_settings,
            self._mandatory_settings
            )
        for setting in invalid:
            emit_error(self, MeasurementErrors.INVALID_MEAS_SETTINGS, setting)

        if invalid:
            return

        self.__settings = new_settings

        # Instantiate primary controller class
        primary_config, primary_measures = ast.literal_eval(
            self.__settings.get('devices', 'primary_controller'))
        self.primary_controller = self.__make_controller(primary_config,
                                                         is_primary=True)
        self.primary_controller.what_to_measure(primary_measures)

        # Instantiate secondary controller classes
        secondary_controllers = []
        secondary_set = ast.literal_eval(
             self.__settings.get('devices', 'secondary_controllers'))
        for secondary_config, secondary_measures in secondary_set:
            ctrl = self.__make_controller(secondary_config, is_primary=False)
            ctrl.what_to_measure(secondary_measures)
            ctrl.moveToThread(self.thread)
            secondary_controllers.append(ctrl)
        self.secondary_controllers = secondary_controllers
        self.thread.start()

        # Instantiate camera classes
        cameras = []
        camera_set = ast.literal_eval(
             self.__settings.get('devices', 'cameras'))
        for camera_settings in camera_set:
            cameras.append(self.__make_camera(camera_settings))
        self.cameras = cameras

    settings = property(__get_settings, set_settings)

    @property
    def secondary_controllers(self):
        """Return the controllers used by this class."""
        return self.__secondary_controllers

    @secondary_controllers.setter
    def secondary_controllers(self, new_controllers):
        """Set the controllers which should be
        used and handle signals.

        Parameters
        ----------
        new_controllers : list
            List of controller class objects.
        """
        self.disconnect_controllers(self.__secondary_controllers)
        self.__secondary_controllers = new_controllers
        self.connect_controllers(self.__secondary_controllers)

    @property
    def primary_controller(self):
        """Return the primary controllers used by this class."""
        return self.__primary_controller

    @primary_controller.setter
    def primary_controller(self, new_controller):
        """Set the primary controller.

        Parameters
        ----------
        new_controller : controller class object
            Controller that sets the energy.
        """
        self.disconnect_primary_controller()
        self.__primary_controller = new_controller
        self.connect_primary_controller()

    @abstractmethod
    def is_finished(self):
        """Check if the full measurement cycle is done.

        This function must be reimplemented in subclasses. It
        should check if the measurement cycle is done via the
        settings.

        Returns
        -------
        bool
        """
        return True

    def finalize(self):
        """Finish the measurement cycle.

        Emit data and reset the class.

        Returns
        -------
        None.
        """
        try:
            self.save_data()
        finally:
            self.current_energy = 0
            # Set LEED energy to 0.
            self.__primary_controller.data_ready.disconnect()
            self.disconnect_controllers(self.secondary_controllers)
            self.disconnect_cameras(self.cameras)
            self.set_LEED_energy(self.current_energy, 1000)
            self.disconnect_primary_controller()
            self.finished.emit((self.plot_info, self.data_points))

    def save_data(self):
        """Save data.

        Returns
        -------
        None.
        """
        # class_name = self.__class__.__name__
        # clock = strftime("_%Y-%m-%d_%H-%M-%S", localtime())
        # csv_name = class_name + clock + ".csv"
        csv_name = "measurement_data.csv"
        length = len(self.data_points['nominal_energy'])
        # TODO: Where to save this? Need to add folder creation, this solution is only temporary. Images and measurement data from controllers will be saved with the configuration files in a zip folder.
        with open(csv_name, 'w', encoding='UTF8', newline='') as file_name:
            writer = csv.writer(file_name)
            writer.writerow(self.data_points.keys())
            for i in range(length-1):
                values = []
                for key in self.data_points.keys():
                    if self.data_points[key]:
                        values.append(self.data_points[key][i])
                writer.writerow(values)

    def set_LEED_energy(self, *message):
        """Set the electron energy used for LEED.

        In order to achieve quicker settling times for the
        LEED electronics one can do quick steps after each
        other whose effects cancel each other out.
        i.e.: A small energy overshoot with an immediate
        correction afterwards.

        Parameters
        ----------
        message : tuple
            Contains data necessary to set the energy.

        Returns
        -------
        None.
        """
        self.__primary_controller.set_energy(*message)

    @abstractmethod
    def abort(self):
        """Abort all current actions.

        This function needs to be reimplemented in subclasses.
        Implementation needs to reset all variables used in loop
        operation. Then call super.

        Returns
        -------
        None.
        """
        self.abort_action.emit()
        self.current_energy = 0
        # Set LEED energy to 0
        self.__primary_controller.data_ready.disconnect()
        self.disconnect_controllers(self.secondary_controllers)
        for controller in self.secondary_controllers:
            controller.busy = False
        self.disconnect_cameras(self.cameras)
        for camera in self.cameras:
            camera.busy = False
        self.set_LEED_energy(self.current_energy, 1000)
        self.disconnect_primary_controller()
        self.__primary_controller.busy = False

        for key in self.data_points:
            self.data_points[key] = []

    def connect_cameras(self, cameras=None):
        """Connect necessary camera signals."""
        if not cameras:
            cameras = self.__cameras
        for camera in cameras:
            if camera:
                camera.camera_busy.connect(
                    self.receive_from_camera, type=qtc.Qt.UniqueConnection
                    )
                self.ready_for_measurement.connect(
                    camera.trigger_now, type=qtc.Qt.UniqueConnection
                    )
        # camera.disconnect does not need to be hooked up to the
        # abort_action signal as it is called in the disconnecting
        # of the camera signals anyway.

    def connect_controllers(self, controllers=None):
        """Connect necessary controller signals."""

        if not controllers:
            controllers = self.__secondary_controllers
        for controller in controllers:
            if controller:
                controller.data_ready.connect(
                    self.receive_from_controller,
                    type=qtc.Qt.UniqueConnection
                    )
                self.ready_for_measurement.connect(
                    controller.measure_now,
                    type=qtc.Qt.UniqueConnection
                    )
                self.abort_action.connect(
                    controller.abort_and_reset,
                    type=qtc.Qt.UniqueConnection
                    )
                self.begin_preparation.connect(
                    controller.trigger_begin_preparation,
                    type=qtc.Qt.UniqueConnection
                    )
                self.continue_preparation.connect(
                    controller.trigger_continue_preparation,
                    type=qtc.Qt.UniqueConnection
                    )

    def connect_primary_controller(self):
        """Connect signals of the primary controller."""

        self.__primary_controller.data_ready.connect(
            self.receive_from_controller, type=qtc.Qt.UniqueConnection
            )
        self.abort_action.connect(
            self.__primary_controller.abort_and_reset,
            type=qtc.Qt.UniqueConnection
            )
        self.begin_preparation.connect(
            self.__primary_controller.trigger_begin_preparation,
            type=qtc.Qt.UniqueConnection
            )
        self.continue_preparation.connect(
            self.__primary_controller.trigger_continue_preparation,
            type=qtc.Qt.UniqueConnection
            )
        self.__primary_controller.about_to_trigger.connect(
            self.do_next_measurement, type=qtc.Qt.UniqueConnection
            )

    def disconnect_cameras(self, cameras):
        """Disconnect necessary camera signals."""
        for camera in cameras:
            if camera:
                camera.disconnect()
                camera.camera_busy.disconnect(self.receive_from_camera)
                self.ready_for_measurement.disconnect(camera.trigger_now)

    def disconnect_controllers(self, controllers):
        """Disconnect necessary controller signals."""
        for controller in controllers:
            if controller:
                controller.serial.serial_disconnect()
                controller.data_ready.disconnect(
                    self.receive_from_controller
                    )
                self.ready_for_measurement.disconnect(controller.measure_now)
                self.abort_action.disconnect(controller.abort_and_reset)
                self.begin_preparation.disconnect(
                    controller.trigger_begin_preparation
                    )
                self.continue_preparation.disconnect(
                    controller.trigger_continue_preparation
                    )
                try:
                    controller.controller_busy.disconnect()
                except TypeError:
                    # controller_busy was not connected at the time
                    pass

    def disconnect_primary_controller(self):
        """Disconnect signals of the primary controller."""
        if self.__primary_controller is not None:
            self.__primary_controller.serial.serial_disconnect()
            try:
                self.__primary_controller.data_ready.disconnect()
            except TypeError:
                # data_ready is already disconnected from all its slots
                pass
            self.abort_action.disconnect(self.__primary_controller.abort_and_reset)
            self.begin_preparation.disconnect(
                self.__primary_controller.trigger_begin_preparation
                )
            self.continue_preparation.disconnect(
                self.__primary_controller.trigger_continue_preparation
                )
            self.__primary_controller.about_to_trigger.disconnect()
            try:
                self.__primary_controller.controller_busy.disconnect()
            except TypeError:
                    # controller_busy was not connected at the time
                    pass

    def prepare_cameras(self):
        """Prepare cameras for a measurement.

        Can be a no-op if the connected cameras do not need
        to be prepared. Otherwise this function should do
        everything needed to prepare the cameras for a full
        measurement cycle.

        Returns
        -------
        None.
        """
        return
        # TODO: call it when controllers prepare, basic implementation?

    def switch_signals_for_preparation(self):
        """Switch signals for preparation.

        The about_to_trigger signal has to be disconnected
        during the first preparation step as the starting
        energy is set in it. Since the primary controller
        will return measurements after setting an energy,
        the data_ready signal has to be disconnected as well.

        In order to know if the first preparation step has
        been done, the controller_busy must be connected and
        returned as soon as all points in the begin_prepare_todos
        dictionary have been executed.

        Before the second step starts, this function is
        called again. It will reconnect all disconnected
        signals and disconnect controller_busy.

        If about_to_trigger is connected:
            Disconnect:
                about_to_trigger
                data_ready
            Connect:
                controller_busy

        If about_to_trigger is disconnected:
            Disconnect:
                controller_busy
            Connect:
                about_to_trigger
                data_ready

        Returns
        -------
        None.
        """
        try:
            self.__primary_controller.about_to_trigger.connect(
                self.do_next_measurement, type=qtc.Qt.UniqueConnection
                )
        except TypeError:
            self.__primary_controller.about_to_trigger.disconnect()
            self.__primary_controller.data_ready.disconnect()
            for controller in self.__secondary_controllers:
                if controller:
                    controller.data_ready.disconnect()
                    controller.controller_busy.connect(
                        self.continue_measurement_preparation
                        )
            self.__primary_controller.controller_busy.connect(
                self.continue_measurement_preparation
                )
        else:
            self.__primary_controller.controller_busy.disconnect()
            for controller in self.__secondary_controllers:
                if controller:
                    controller.controller_busy.disconnect()
                    controller.data_ready.connect(
                        self.receive_from_controller,
                        type=qtc.Qt.UniqueConnection
                        )
            self.__primary_controller.data_ready.connect(
                self.receive_from_controller, type=qtc.Qt.UniqueConnection
                )

    def begin_measurement_preparation(self):
        """Start preparation for measurements.

        Prepare the controllers for a measurement which starts
        the measurement cycle. This function only performs the
        first part. (Everything that is done before the starting
        energy is set.)

        Signals are switched to allow preparation to take place.

        Returns
        -------
        None.

        Emits
        -----
        begin_preparation
            Starts the measurement preparation and carries
            a tuple of energies and times with it.

        """
        self.switch_signals_for_preparation()
        self.current_energy = self.__start_energy
        self.begin_preparation.emit((self.__start_energy,
                                     self.__long_settle_time))

    def continue_measurement_preparation(self, busy):
        """Continue preparation for measurements.

        Prepare the controllers for a measurement which starts
        the measurement cycle. This function only performs the
        second part. (Everything that is done after the starting
        energy is set.)

        Signals are switched back and the busy signal is
        connected to a function which checks if all the
        controllers are done with the preparation.

        Parameters
        ----------
        busy : bool
            Busy state of the emitting controller.

        Returns
        -------
        None.

        Emits
        -----
        continue_preparation
            Starts the second part of the measurement
            preparation.
        """
        if busy:
            return
        if self.__primary_controller.busy:
            return
        if any(controller.busy for controller in self.__secondary_controllers):
            return
        self.switch_signals_for_preparation()
        for controller in self.secondary_controllers:
            if controller:
                controller.controller_busy.connect(
                    self.is_preparation_finished,
                    type=qtc.Qt.UniqueConnection
                    )
        self.__primary_controller.controller_busy.connect(
            self.is_preparation_finished, type=qtc.Qt.UniqueConnection
            )
        self.continue_preparation.emit()

    def ready_for_next_measurement(self):
        """Check if all measurements have been received.

        After all measurements have been received, a check if
        the loop is done will be called.

        Returns
        -------
        None.
        """
        if self.__primary_controller.busy:
            return
        if any(controller.busy for controller in self.__secondary_controllers):
            return
        if any(camera.busy for camera in self.__cameras):
            return
        if self.is_finished():
            self.finalize()
        else:
            self.start_next_measurement()

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
        if self.__primary_controller.busy:
            return
        if any(controller.busy for controller in self.__secondary_controllers):
            return
        for controller in self.secondary_controllers:
            if controller:
                controller.controller_busy.disconnect()
        self.__primary_controller.controller_busy.disconnect()
        self.start_next_measurement()

    def receive_from_camera(self, busy):
        """Receive not busy signal from camera.

        Receive the camera_busy signal and check if all cameras
        and controllers are ready if the camera that emitted the
        signal was not busy.

        Parameters
        ----------
        busy : bool
            busy state of the camera

        Returns
        -------
        None.
        """
        if not busy:
            self.ready_for_next_measurement()

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
        self.ready_for_next_measurement()

    def do_next_measurement(self):
        """Do the next measurement.

        Start measuring on all secondary controllers and cameras
        as soon as the about_to_trigger signal has been received
        from the primary controller.

        Returns
        -------
        None.
        """
        self.ready_for_measurement.emit()

    @abstractmethod
    def start_next_measurement(self):
        """Set energy (if required) and measure.

        If the primary controller sets an energy:
        Set energy via the primary controller. Once this is done
        the returned about_to_trigger signal should start the
        measurement on all secondary controllers.

        If the primary controller does not set an energy:
        Start the measurement on the primary and all secondary
        controllers immediately without waiting for an
        about_to_trigger signal.

        Returns
        -------
        None.
        """
        # Set every controller/camera busy at the beginning <-- DO THIS ALWAYS!!!!
        # And append self.current_energy to self.data_points['nominal_energy']
        # If the primary controller sets an energy:
        # set_LEED_energy(self.current_energy, self.__settle_time)
        # ^ set energy and wait, triggers measurement afterwards

        # If the primary controller does not set an energy:
        # self.__primary_controller.measure_now()
        # do_next_measurement()

    def __make_controller(self, controller_settings, is_primary=False):
        """Instantiate controller class object.

        Take controller settings and generate a controller object
        from it.

        Parameters
        ----------
        controller_settings : dict, ConfigParser, str, path
            Settings used for the instantiated controller.
            Has to contain the name of the controller class
            to be instantiated.
        is_primary : boolean
            True if the controller is the primary controller.

        Returns
        -------
        instance : controller class object
            The controller that is going to handle the hardware.

        Emits
        -----
        MeasurementErrors.MISSING_CLASS_NAME
            If the controller class name is missing.
        MeasurementErrors.INVALID_MEAS_SETTINGS
            If the controller could not be instantiated
            from the given name.
        """
        config = ConfigParser()
        config, invalid = config_has_sections_and_options(
            self,
            controller_settings,
            [('controller', 'controller_class'), ('controller', 'port_name')]
            )
        if invalid:
            emit_error(self, MeasurementErrors.MISSING_CLASS_NAME,
                       ('controller', 'controller_class'))
            return

        controller_cls_name = config.get('controller', 'controller_class')
        port_name = config.get('controller', 'port_name')
        try:
            controller_class = class_from_name('controller',
                                               controller_cls_name)
        except ValueError:
            emit_error(self, MeasurementErrors.INVALID_MEAS_SETTINGS,
                       'controller/controller_class')
            return
        instance = controller_class(config, port_name, sets_energy=is_primary)
        return instance

    def __make_camera(self, camera_settings):
        """Instantiate camera class object.

        Take camera settings and generate a camera object from it.

        Parameters
        ----------
        camera_settings : dict, ConfigParser, str, path
            Settings used for the instantiated camera.
            Has to contain the name of the camera class
            to be instantiated.

        Returns
        -------
        instance : camera class object
            The camara class used for the connected camera.

        Emits
        -----
        MeasurementErrors.MISSING_CLASS_NAME
            If the camera class name is missing.
        MeasurementErrors.INVALID_MEAS_SETTINGS
            If the camera could not be instantiated
            from the given name.
        """
        config = ConfigParser()

        config, invalid = config_has_sections_and_options(
            self,
            camera_settings,
            [('camera', 'camera_class'), ('camera', 'port_name')]
            )
        if invalid:
            emit_error(self, MeasurementErrors.MISSING_CLASS_NAME,
                       ('camera', 'camera_class'))
            return

        camera_cls_name = config.get('camera', 'camera_class')
        port_name = config.get('camera', 'port_name')
        try:
            camera_class = class_from_name('camera', camera_cls_name)
        except ValueError:
            emit_error(self, MeasurementErrors.INVALID_MEAS_SETTINGS,
                       'camera/camera_class')
            return
        instance = camera_class(config, port_name)
        return instance
