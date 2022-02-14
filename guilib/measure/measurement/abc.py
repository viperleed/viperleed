"""Module abc of viperleed.guilib.measure.measurement.

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
import os
import shutil
from abc import abstractmethod
import ast
from configparser import ConfigParser
from time import localtime, strftime
import inspect
from pathlib import Path

from PyQt5 import QtCore as qtc

# ViPErLEED modules
from viperleed.guilib.measure.hardwarebase import (
    emit_error, ViPErLEEDErrorEnum, QMetaABC,
    config_has_sections_and_options, class_from_name
    )
from viperleed.guilib.measure.datapoints import DataPoints
from viperleed.guilib import measure as vpr_measure


SYSTEM_CONFIG_PATH = (Path(inspect.getfile(vpr_measure)).parent
                       / '_defaults/_system_settings.ini')


class MeasurementErrors(ViPErLEEDErrorEnum):
    """Errors that might occur during a measurement cycle."""
    MISSING_SETTINGS = (300,
                        "Measurements cannot be taken without settings. "
                        "Load an appropriate settings file before "
                        "proceeding.")
    INVALID_MEAS_SETTINGS = (301,
                             "Invalid measurement settings: Required "
                             "settings {!r} missing or wrong. Check "
                             "configuration file.")
    MISSING_CLASS_NAME = (302,
                          "{!r} is missing the name "
                          "of its related class object.")
    RUNTIME_ERROR = (303, "Runtime error. Info: {}")


class MeasurementABC(qtc.QObject, metaclass=QMetaABC):
    """Generic measurement class."""
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
    # Abort current tasks on the connected devices.
    request_stop_devices = qtc.pyqtSignal()
    # Abort current tasks on the primary controller.
    request_stop_primary = qtc.pyqtSignal()
    # Make the GUI update the current display of collected data.
    new_data_available = qtc.pyqtSignal()
    # Let the GUI know that the preparation for the
    # measurement is finished.
    prepared = qtc.pyqtSignal()

    _mandatory_settings = [
        ('devices', 'primary_controller'),
        ('measurement_settings', 'start_energy'),
        ('measurement_settings', 'save_here'), # TODO: change name, move to global settings
        ]

    def __init__(self, measurement_settings):
        """Initialise measurement class"""
        super().__init__()
        self._other_mandatory_settings = [('measurement_settings',
                                           'measurement_class',
                                           (self.__class__.__name__,))]
        self.current_energy = 0
        self.__settings = None
        self.__primary_controller = None
        self.__secondary_controllers = []
        self.__cameras = []
        self.__aborted = False
        self.counter = 0
        self.start_energy = 0
        self.__long_settle_time = 0
        self.threads = []
        self.running = False
        self.primary_delay = 0
        self.data_points = DataPoints()
        self.data_points.error_occurred.connect(self.error_occurred)

        self.__init_errors = []  # Report these with a little delay
        self.__init_err_timer = qtc.QTimer(self)
        self.__init_err_timer.setSingleShot(True)

        self.error_occurred.connect(self.__on_init_errors)
        self.__init_err_timer.timeout.connect(self.__report_init_errors)

        self.camera_timer = qtc.QTimer()
        self.camera_timer.setSingleShot(True)
        self.camera_timer.setParent(self)

        self.set_settings(measurement_settings)

        if self.settings:
            self.start_energy = self.settings.getfloat(
                'measurement_settings', 'start_energy', fallback=0
                )
            self.__long_settle_time = self.primary_controller.long_settle_time

        self.force_return_timer = qtc.QTimer(parent=self)
        self.force_return_timer.setSingleShot(True)
        self.force_return_timer.timeout.connect(self.return_to_gui)


        if self.__init_errors:
            self.__init_err_timer.start(20)
        self.error_occurred.disconnect(self.__on_init_errors)

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
        self.disconnect_cameras()
        self.__cameras = new_cameras
        self.connect_cameras()

    @property
    def controllers(self):
        """Return a tuple of all controllers."""
        controllers = (self.primary_controller, *self.secondary_controllers)
        return tuple(c for c in controllers if c)

    @property
    def devices(self):
        """Return all controllers and cameras."""
        return *self.controllers, *self.cameras

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
            (*self._mandatory_settings, *self._other_mandatory_settings)
            )
        for setting in invalid:
            emit_error(self, MeasurementErrors.INVALID_MEAS_SETTINGS, setting)

        if invalid:
            return

        self.__settings = new_settings
        
        sys_config = ConfigParser()
        sys_config, invalid = config_has_sections_and_options(
            self,
            SYSTEM_CONFIG_PATH,
            [('settings', 'configuration_path'),]
            )
        device_config_path = sys_config.get('settings', 'configuration_path')

        # Instantiate primary controller class
        primary_config, primary_measures = ast.literal_eval(
            self.settings.get('devices', 'primary_controller')
            )
        primary_config = primary_config.replace('__CONFIG__',
                                                device_config_path)
        self.primary_controller = self.__make_controller(primary_config,
                                                         is_primary=True)
        self.primary_controller.what_to_measure(primary_measures)

        # Instantiate secondary controller classes
        secondary_controllers = []
        secondary_set = ast.literal_eval(
             self.settings.get(
                'devices', 'secondary_controllers', fallback='()'
                )
             )
        for secondary_config, secondary_measures in secondary_set:
            secondary_config = secondary_config.replace('__CONFIG__',
                                                        device_config_path)
            ctrl = self.__make_controller(secondary_config, is_primary=False)
            ctrl.what_to_measure(secondary_measures)
            self.threads.append(qtc.QThread())
            ctrl.moveToThread(self.threads[-1])
            secondary_controllers.append(ctrl)
        self.secondary_controllers = secondary_controllers
        for thread in self.threads:
            thread.start(priority=thread.TimeCriticalPriority)

        # Instantiate camera classes
        cameras = []
        camera_set = ast.literal_eval(
             self.settings.get(
                 'devices', 'cameras', fallback='()'
                 )
             )
        for camera_settings in camera_set:
            camera_settings = camera_settings.replace('__CONFIG__',
                                                      device_config_path)
            cameras.append(self.__make_camera(camera_settings))
        self.cameras = cameras
        path = self.settings.get('measurement_settings', 'save_here')
        try:
            os.mkdir(path + '__tmp__/')
        except FileExistsError:
            # Folder already exists.
            pass
        for camera in self.cameras:
            camera.process_info.base_path = path + '__tmp__/' + camera.name
            try:
                os.mkdir(camera.process_info.base_path)
            except FileExistsError:
                # Folder already exists.
                pass

        for device in self.devices:
            device.error_occurred.connect(self.__on_hardware_error)

        self.data_points.primary_controller = self.primary_controller

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
        self.disconnect_secondary_controllers()
        self.__secondary_controllers = new_controllers
        self.connect_secondary_controllers()

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
        primary = self.primary_controller
        if busy or any(controller.busy for controller in self.controllers):
            return
        try:
            self.save_data()
        finally:
            self.current_energy = 0
            # Set LEED energy to 0.
            try:
                primary.data_ready.disconnect()
            except TypeError:
                pass
            try:
                primary.controller_busy.disconnect()
            except TypeError:
                pass
            self.disconnect_secondary_controllers()
            self.disconnect_cameras()
            primary.controller_busy.connect(self.return_to_gui,
                                            type=qtc.Qt.UniqueConnection)
            self.set_LEED_energy(self.current_energy, 50, trigger_meas=False)

    def save_data(self):
        """Save data.

        Returns
        -------
        None.
        """
        path = self.settings.get('measurement_settings', 'save_here')
        clock = strftime("%Y-%m-%d_%H-%M-%S/", localtime())
        os.mkdir(path + clock)
        to_move_list = os.listdir(path + '__tmp__/')
        for to_move in to_move_list:
            shutil.move(path + '__tmp__/' + to_move, path + clock + to_move)
        os.rmdir(path + '__tmp__/')
        if not self.data_points:
            return
        csv_name = path + clock + 'measurement.csv'
        self.data_points.save_data(csv_name)

        for controller in self.controllers:
            file_name = (path + clock + 'controller_' +
                controller.name + '.ini')
            with open(file_name, 'w') as configfile:
                controller.settings.write(configfile)
        for camera in self.cameras:
            name = camera.name.replace(' ', '_')
            file_name = path + clock + 'camera_' + name + '.ini'
            with open(file_name, 'w') as configfile:
                camera.settings.write(configfile)

        file_name = path + clock + 'measurement.ini'
        with open(file_name, 'w') as configfile:
            self.settings.write(configfile)

    def set_LEED_energy(self, *message, trigger_meas=True, **kwargs):
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

        trigger_meas : bool, optional
            True if the controllers are supposed to take
            measurements after the energy has been set.
            Default is True.
        **kwargs : object
            Other unused keyword arguments.

        Returns
        -------
        None.
        """
        self.primary_controller.set_energy(*message, trigger_meas=trigger_meas)
        self.primary_delay = 0
        for i, value in enumerate(message):
            if i % 2 != 0:
                self.primary_delay += value

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
        if self.__aborted:
            return
        self.__aborted = True
        self.settings.set('measurement_settings', 'was_aborted', 'True')
        self.force_return_timer.start(4500)
        self.prepare_finalization()

    def connect_cameras(self):
        """Connect necessary camera signals."""
        for camera in self.cameras:
            camera.camera_busy.connect(self.receive_from_camera,
                                       type=qtc.Qt.UniqueConnection)
            self.camera_timer.timeout.connect(camera.trigger_now,
                                              type=qtc.Qt.UniqueConnection)
            self.begin_preparation.connect(camera.start,
                                           type=qtc.Qt.UniqueConnection)
        # camera.disconnect does not need to be hooked up to the
        # request_stop_devices signal as it is called in the disconnecting
        # of the camera signals anyway.

    def connect_secondary_controllers(self):
        """Connect necessary controller signals."""
        for controller in self.secondary_controllers:
            if not controller:
                continue
            controller.data_ready.connect(self.receive_from_controller,
                                          type=qtc.Qt.UniqueConnection)
            self.ready_for_measurement.connect(controller.measure_now,
                                               type=qtc.Qt.UniqueConnection)
            self.request_stop_devices.connect(controller.stop,
                                              type=qtc.Qt.UniqueConnection)
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
        primary = self.primary_controller
        primary.data_ready.connect(self.receive_from_controller,
                                   type=qtc.Qt.UniqueConnection)
        self.request_stop_devices.connect(primary.stop,
                                          type=qtc.Qt.UniqueConnection)
        self.request_stop_primary.connect(primary.stop,
                                          type=qtc.Qt.UniqueConnection)
        self.begin_preparation.connect(primary.trigger_begin_preparation,
                                       type=qtc.Qt.UniqueConnection)
        self.continue_preparation.connect(primary.trigger_continue_preparation,
                                          type=qtc.Qt.UniqueConnection)
        primary.about_to_trigger.connect(self.do_next_measurement,
                                         type=qtc.Qt.UniqueConnection)

    def disconnect_cameras(self):
        """Disconnect necessary camera signals."""
        for camera in self.cameras:
            camera.disconnect()
            try:
                camera.camera_busy.disconnect(self.receive_from_camera)
            except TypeError:
                pass
            try:
                self.camera_timer.timeout.disconnect(camera.trigger_now)
            except TypeError:
                pass
            try:
                self.begin_preparation.disconnect(camera.start)
            except TypeError:
                pass

    def disconnect_secondary_controllers(self):
        """Disconnect necessary controller signals."""
        for controller in self.secondary_controllers:
            if not controller:
                continue
            controller.disconnect_()
            try:
                controller.data_ready.disconnect(self.receive_from_controller)
            except TypeError:
                pass
            try:
                self.ready_for_measurement.disconnect(controller.measure_now)
            except TypeError:
                pass
            try:
                self.request_stop_devices.disconnect(controller.stop)
            except TypeError:
                pass
            try:
                self.begin_preparation.disconnect(
                    controller.trigger_begin_preparation
                    )
            except TypeError:
                pass
            try:
                self.continue_preparation.disconnect(
                    controller.trigger_continue_preparation
                    )
            except TypeError:
                pass
            try:
                controller.controller_busy.disconnect()
            except TypeError:
                pass

    def disconnect_primary_controller(self):
        """Disconnect signals of the primary controller."""
        primary = self.primary_controller
        if primary is None:
            return
        primary.disconnect_()
        try:
            primary.data_ready.disconnect()
        except TypeError:
            pass
        try:
            self.request_stop_devices.disconnect(primary.stop)
        except TypeError:
            pass
        try:
            self.request_stop_primary.disconnect(primary.stop)
        except TypeError:
            pass
        try:
            self.begin_preparation.disconnect(
                primary.trigger_begin_preparation
                )
        except TypeError:
            pass
        try:
            self.continue_preparation.disconnect(
                primary.trigger_continue_preparation
                )
        except TypeError:
            pass
        try:
            primary.about_to_trigger.disconnect()
        except TypeError:
            pass
        try:
            primary.controller_busy.disconnect()
        except TypeError:
            pass

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
            self.primary_controller.about_to_trigger.connect(
                self.do_next_measurement, type=qtc.Qt.UniqueConnection
                )
        except TypeError:
            self.primary_controller.about_to_trigger.disconnect()
            for controller in self.controllers:
                controller.data_ready.disconnect()
                controller.controller_busy.connect(
                    self.continue_measurement_preparation
                    )
        else:
            for controller in self.controllers:
                controller.controller_busy.disconnect()
                controller.data_ready.connect(self.receive_from_controller,
                                              type=qtc.Qt.UniqueConnection)

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
        self.__aborted = False
        self.settings.set('measurement_settings', 'was_aborted', 'False')
        self.running = True
        self.switch_signals_for_preparation()
        self.current_energy = self.start_energy
        self.begin_preparation.emit((self.start_energy,
                                     self.__long_settle_time))

    def continue_measurement_preparation(self, _):
        """Continue preparation for measurements.

        Prepare the controllers for a measurement which starts
        the measurement cycle. This function only performs the
        second part. (Everything that is done after the starting
        energy is set.)

        Signals are switched back and the busy signal is
        connected to a function which checks if all the
        controllers are done with the preparation.

        Returns
        -------
        None.

        Emits
        -----
        continue_preparation
            Starts the second part of the measurement
            preparation.
        """
        if any(controller.busy for controller in self.controllers):
            return
        self.switch_signals_for_preparation()
        for controller in self.controllers:
            controller.controller_busy.connect(self.is_preparation_finished,
                                               type=qtc.Qt.UniqueConnection)
        for camera in self.cameras:
            camera.camera_busy.connect(self.is_preparation_finished,
                                       type=qtc.Qt.UniqueConnection)
        self.continue_preparation.emit()

    def ready_for_next_measurement(self):
        """Check if all measurements have been received.

        After all measurements have been received, a check if
        the loop is done will be called.

        Returns
        -------
        None.
        """
        if any(device.busy for device in self.devices):
            return
        if self.is_finished():
            self.prepare_finalization()
            self.data_points.calculate_times()
            self.new_data_available.emit()
        else:
            self.data_points.calculate_times()
            self.new_data_available.emit()
            self.start_next_measurement()

    def is_preparation_finished(self):
        """Check if measurement preparation is done.

        Whenever a device is done with its preparation
        this function will be called. The busy attribute is
        used to check if all devices are done with their
        preparation. If this is true, the busy signal going to
        this function will be disconnected and the first
        measurement will immediately be started.

        Returns
        -------
        None.
        """
        if any(device.busy for device in self.devices):
            return
        for controller in self.controllers:
            controller.controller_busy.disconnect()
        for camera in self.cameras:
            camera.camera_busy.disconnect(self.is_preparation_finished)
        self.prepared.emit()
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
        self.ready_for_next_measurement()

    def receive_from_controller(self, controller, receive):
        """Receive measurement data from the controller.

        Append received data to the internal dictionary. Emit an
        error if received dictionary contains a section that does
        not exist in the internal dictionary. After appending all
        measurements check if all of the connected controllers are
        ready for the next measurement.

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
        # TODO: check if one controller can return data while
        # another controller has changed his busy state but
        # hasn't returned data yet. (race condition)
        if controller == self.primary_controller:
            self.data_points.add_data(receive, controller, self.primary_delay)
        else:
            self.data_points.add_data(receive, controller)
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

        super().start_next_measurement() must be called in
        subclasses in order to create a new data_point for
        the next measurement at the beginning of
        start_next_ measurement.

        Returns
        -------
        None.
        """
        self.data_points.new_data_point(self.current_energy, self.controllers,
                                        self.cameras)

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
            # TODO: we do this if the class name AND if the port name is missing, will need to be edited once we add the serial numbers
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
            [('camera_settings', 'class_name'), ]
            )
        if invalid:
            emit_error(self, MeasurementErrors.MISSING_CLASS_NAME,
                       ('camera_settings', 'class_name'))
            return

        camera_cls_name = config.get('camera_settings', 'class_name')
        try:
            camera_class = class_from_name('camera', camera_cls_name)
        except ValueError:
            emit_error(self, MeasurementErrors.INVALID_MEAS_SETTINGS,
                       'camera_settings/class_name')
            return
        instance = camera_class(settings=config)
        return instance

    def return_to_gui(self, *__args):
        """Return collected data to gui.

        Quit thread and emit data.

        Returns
        -------
        None.

        Emits
        -----
        finished
            A signal containing all collected data.
        """
        self.disconnect_primary_controller()
        self.force_return_timer.stop()
        for thread in self.threads:
            thread.quit()
        self.running = False
        # self.finished.emit(self.data_points)
        self.finished.emit(())

    def prepare_finalization(self):
        """Prepare for finalization.

        This function may need to be reimplemented in subclasses.
        Ensure that finalization is prepared properly and that
        self.finalize() gets called. Can be only a call to
        self.finalize() if no preparation is needed.

        Returns
        -------
        None.
        """
        for controller in self.secondary_controllers:
            if not controller:
                continue
            try:
                self.ready_for_measurement.disconnect(controller.measure_now)
            except TypeError:
                pass
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
        self.request_stop_devices.emit()

    def __on_hardware_error(self, *_):
        """Abort if a hardware error occurs."""
        self.abort()

    def __on_init_errors(self, err):
        """Collect initialization errors to report later."""
        self.__init_errors.append(err)

    def __report_init_errors(self):
        """Emit error_occurred for each initialization error."""
        for error in self.__init_errors:
            self.error_occurred.emit(error)
        self.__init_errors = []
