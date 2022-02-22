"""Module abc of viperleed.guilib.measure.measurement.

========================================
   ViPErLEED Graphical User Interface
========================================

Created: 2021-07-19
Author: Michele Riva
Author: Florian Doerr

This module contains the definition of the MeasurementABC class
which gives commands to the controller classes and its associated
ViPErLEEDErrorEnum class MeasurementErrors.
"""

import shutil
from abc import abstractmethod
from time import localtime, strftime
from pathlib import Path

from PyQt5 import QtCore as qtc

from viperleed.guilib.measure import hardwarebase as base
from viperleed.guilib.measure.datapoints import DataPoints
from viperleed.guilib.measure.classes.settings import (
    ViPErLEEDSettings, NoSettingsError, get_system_config, NotASequenceError
    )


class MeasurementErrors(base.ViPErLEEDErrorEnum):                                    # TODO: fix text
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
                          "{} is missing the name "
                          "of its related class object.")
    RUNTIME_ERROR = (303, "Runtime error. Info: {}")
    INVALID_SETTING_WITH_FALLBACK = (
        304,
        "Invalid/unreadable measurement settings value {} for setting {!r}. "
        "Using {} instead. Consider fixing your configuration file."
        )


class MeasurementABC(qtc.QObject, metaclass=base.QMetaABC):                          # TODO: doc about inner workings
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
        ]

    def __init__(self, measurement_settings):
        """Initialise measurement class"""
        super().__init__()
        self._other_mandatory_settings = [('measurement_settings',
                                           'measurement_class',
                                           (self.__class__.__name__,))]
        self.current_energy = 0
        self.__settings = ViPErLEEDSettings()
        self.__primary_controller = None
        self.__secondary_controllers = []
        self.__cameras = []
        self.__aborted = False

        self.threads = []
        self.running = False  # used for aborting from outside
        self.primary_delay = 0                                                  # TODO: bad name?
        self.data_points = DataPoints()
        self.data_points.error_occurred.connect(self.error_occurred)

        self.__init_errors = []  # Report these with a little delay
        self.__init_err_timer = qtc.QTimer(parent=self)
        self.__init_err_timer.setSingleShot(True)
        self.__init_err_timer.timeout.connect(self.__report_init_errors)

        self.error_occurred.connect(self.__on_init_errors)

        self.camera_timer = qtc.QTimer(parent=self)
        self.camera_timer.setSingleShot(True)

        self.force_return_timer = qtc.QTimer(parent=self)
        self.force_return_timer.setSingleShot(True)
        self.force_return_timer.timeout.connect(self.return_to_gui)

        self.set_settings(measurement_settings)

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
    def current_step_nr(self):
        """Return an incremental number for the current energy step."""
        return len(self.data_points)

    @property
    def devices(self):
        """Return all controllers and cameras."""
        return *self.controllers, *self.cameras

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
        try:
            new_settings = ViPErLEEDSettings.from_settings(new_settings)
        except (ValueError, NoSettingsError):
            base.emit_error(self, MeasurementErrors.MISSING_SETTINGS)
            return

        invalid = new_settings.has_settings(*self._mandatory_settings,
                                            *self._other_mandatory_settings)

        if invalid:
            base.emit_error(self, MeasurementErrors.INVALID_MEAS_SETTINGS,
                            ', '.join(invalid))
            return

        self.__settings = new_settings

        self.__make_primary_ctrl()
        if not self.primary_controller:
            # Something went wrong (already reported in __make_primary)
            return

        self.__make_secondary_ctrls()
        self.__make_cameras()
        self.__make_tmp_directory_tree()

        self.data_points.primary_controller = self.primary_controller

    settings = property(__get_settings, set_settings)

    @property
    def start_energy(self):
        """Return the first energy for the energy ramp."""
        if not self.settings:
            return 0
        try:
            return self.settings.getfloat('measurement_settings',
                                          'start_energy', fallback=0)
        except (TypeError, ValueError):
            # Not a float
            base.emit_error(self, MeasurementErrors.INVALID_MEAS_SETTINGS,
                            'measurement_settings/start_energy')
            return 0

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
                                     self.primary_controller.long_settle_time))

    def connect_cameras(self):
        """Connect necessary camera signals."""
        # camera.disconnect does not need to be hooked up to the
        # request_stop_devices signal as it is always called in
        # finalize() via .disconnect_cameras()
        for camera in self.cameras:
            camera.camera_busy.connect(self.receive_from_camera,
                                       type=qtc.Qt.UniqueConnection)
            self.camera_timer.timeout.connect(camera.trigger_now,
                                              type=qtc.Qt.UniqueConnection)
            self.begin_preparation.connect(camera.start,
                                           type=qtc.Qt.UniqueConnection)

    def connect_primary_controller(self):
        """Connect signals of the primary controller."""
        primary = self.primary_controller
        self.request_stop_primary.connect(primary.stop,
                                          type=qtc.Qt.UniqueConnection)
        primary.about_to_trigger.connect(self.do_next_measurement,
                                         type=qtc.Qt.UniqueConnection)
        primary.data_ready.connect(self.receive_from_controller,                # TODO: same for all
                                   type=qtc.Qt.UniqueConnection)
        self.request_stop_devices.connect(primary.stop,
                                          type=qtc.Qt.UniqueConnection)
        self.begin_preparation.connect(primary.trigger_begin_preparation,
                                       type=qtc.Qt.UniqueConnection)
        self.continue_preparation.connect(primary.trigger_continue_preparation,
                                          type=qtc.Qt.UniqueConnection)

    def connect_secondary_controllers(self):
        """Connect necessary controller signals."""
        for controller in self.secondary_controllers:
            self.ready_for_measurement.connect(controller.measure_now,
                                               type=qtc.Qt.UniqueConnection)
            controller.data_ready.connect(self.receive_from_controller,
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

    def disconnect_cameras(self):
        """Disconnect necessary camera signals."""
        for camera in self.cameras:
            camera.disconnect()                                                 # TODO: rename disconnect_!!
            base.safe_disconnect(camera.camera_busy, self.receive_from_camera)
            base.safe_disconnect(self.camera_timer.timeout, camera.trigger_now)
            base.safe_disconnect(self.begin_preparation, camera.start)

    def disconnect_primary_controller(self):
        """Disconnect signals of the primary controller."""
        primary = self.primary_controller
        if primary is None:
            return
        base.safe_disconnect(self.request_stop_primary, primary.stop)
        base.safe_disconnect(primary.about_to_trigger)
        primary.disconnect_()
        base.safe_disconnect(primary.data_ready)
        base.safe_disconnect(self.request_stop_devices, primary.stop)
        base.safe_disconnect(self.begin_preparation,
                             primary.trigger_begin_preparation)
        base.safe_disconnect(self.continue_preparation,
                             primary.trigger_continue_preparation)
        base.safe_disconnect(primary.controller_busy)

    def disconnect_secondary_controllers(self):
        """Disconnect necessary controller signals."""
        for ctrl in self.secondary_controllers:
            base.safe_disconnect(self.ready_for_measurement, ctrl.measure_now)
            ctrl.disconnect_()                                                  # TODO: same for all controllers
            base.safe_disconnect(ctrl.data_ready, self.receive_from_controller)
            base.safe_disconnect(self.request_stop_devices, ctrl.stop)
            base.safe_disconnect(self.begin_preparation,
                                 ctrl.trigger_begin_preparation)
            base.safe_disconnect(self.continue_preparation,
                                 ctrl.trigger_continue_preparation)
            base.safe_disconnect(ctrl.controller_busy)

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
            base.safe_disconnect(primary.data_ready)
            base.safe_disconnect(primary.controller_busy)
            self.disconnect_secondary_controllers()
            self.disconnect_cameras()
            primary.controller_busy.connect(self.return_to_gui,
                                            type=qtc.Qt.UniqueConnection)
            self.set_LEED_energy(self.current_energy, 50, trigger_meas=False)

    @abstractmethod
    def is_finished(self):
        """Check if the full measurement cycle is done.

        This function must be reimplemented in subclasses. It
        should check if the measurement cycle is done via the
        stings. super().is_finished() is supposed to be
        called in subclasses in order to increment the step
        counter properly.

        Returns
        -------
        bool
        """
        self.data_points.nr_steps_done += 1
        return True

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

    def save_data(self):
        """Save data."""
        sys_config = get_system_config()
        base_path = sys_config.get("PATHS", "measurements", fallback=None)
        if not base_path:
            # TODO: emit a SystemSettingsError of some kind
            raise RuntimeError("Invalid path in system config file")

        tmp_path = base_path + "__tmp__/"
        final_path = base_path + strftime("%Y-%m-%d_%H-%M-%S/", localtime())

        # Move all camera images
        shutil.move(tmp_path, final_path)

        if self.data_points:
            self.data_points.save_data(final_path + 'measurement.csv')

        ctrl_locations = []
        for ctrl in self.controllers:
            fname = "controller_" + ctrl.name.replace(' ', '_') + ".ini"
            with open(final_path + fname, 'w', encoding='utf-8') as fproxy:
                ctrl.settings.write(fproxy)
            ctrl_locations.append("./" + fname)

        cam_locations = []
        for camera in self.cameras:
            fname = "camera_" + camera.name.replace(' ', '_') + ".ini"
            with open(final_path + fname, 'w', encoding='utf-8') as fproxy:
                camera.settings.write(fproxy)
            cam_locations.append("./" + fname)

        # TODO: do the same for all controllers
        if cam_locations:
            self.settings.set("devices", "cameras", str(tuple(cam_locations)))

        file_name = final_path + "measurement.ini"
        with open(file_name, 'w', encoding='utf-8') as configfile:
            self.settings.write(configfile)

    def set_LEED_energy(self, energy, settle_time, *more_steps,
                        trigger_meas=True, **_):
        """Set the electron energy used for LEED.

        In order to achieve quicker settling times for the LEED
        electronics one can do quick steps after one another
        whose effects cancel each other out. E.g., a small energy
        overshoot quickly followed by a correction.

        Parameters
        ----------
        energy : float
            Nominal electron energy in electronvolts.
        settle_time : integer
            Interval in milliseconds that the controller will wait
            before deeming the LEED optics stable at set energy.
        *more_steps : Number
            If given, it should be an even number of elements.
            Odd elements are energies, even ones settle-time intervals.
            Multiple steps can be executed quickly after each other.
            The last step will be the final energy that is set and
            should ensure stabilization of the electronics.
        trigger_meas : bool, optional
            True if the controllers are supposed to take
            measurements after the energy has been set.
            Default is True.

        Returns
        -------
        None.
        """
        self.primary_controller.set_energy(energy, settle_time, *more_steps,
                                           trigger_meas=trigger_meas)
        self.primary_delay = settle_time
        for i, value in enumerate(more_steps):
            if i % 2 != 0:
                self.primary_delay += value

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

    def switch_signals_for_preparation(self):                                   # TODO: counterintitive. See if it can be made simpler to understand
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
            base.safe_disconnect(self.ready_for_measurement.disconnect,
                                 controller.measure_now)
        for controller in self.controllers:
            base.safe_disconnect(controller.controller_busy)
            # Necessary to force secondaries into busy,
            # before the primary returns not busy anymore.
            controller.busy = True
            controller.controller_busy.connect(self.finalize,
                                               type=qtc.Qt.UniqueConnection)
        self.request_stop_devices.emit()

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

        Raises
        ------
        RuntimeError
            If anything goes wrong with instantiation.

        Emits
        -----
        MeasurementErrors.MISSING_CLASS_NAME
            If the camera class name is missing or invalid.
        MeasurementErrors.INVALID_MEAS_SETTINGS
            If the camera could not be instantiated
            from the given name.
        """
        try:
            config = ViPErLEEDSettings.from_settings(camera_settings)
        except (ValueError, NoSettingsError):
            base.emit_error(self, MeasurementErrors.INVALID_MEAS_SETTINGS,
                            'devices/path to camera configuration')
            raise RuntimeError from None

        invalid = config.has_settings(('camera_settings', 'class_name'))
        if invalid:
            base.emit_error(self, MeasurementErrors.MISSING_CLASS_NAME,         # TODO: camera error?
                            'camera_settings/class_name')
            raise RuntimeError

        camera_cls_name = config.get('camera_settings', 'class_name')
        try:
            camera_class = base.class_from_name('camera', camera_cls_name)
        except ValueError:
            base.emit_error(self, MeasurementErrors.INVALID_MEAS_SETTINGS,      # TODO: camera error?
                            'camera_settings/class_name')
            raise RuntimeError

        instance = camera_class(settings=config)
        return instance

    def __make_cameras(self):
        """Make cameras from self.settings."""
        try:
            cam_settings = self.settings.getsequence('devices', 'cameras',
                                                     fallback=())
        except NotASequenceError:
            cam_settings = tuple()

        cameras = []
        for settings in cam_settings:
            try:
                cam = self.__make_camera(settings)
            except RuntimeError:
                continue
            cam.error_occurred.connect(self.__on_hardware_error)
            cameras.append(cam)
        self.cameras = cameras

    def __make_controller(self, configname, measurements, is_primary=False):
        """Instantiate controller class object.

        Take controller settings and generate a controller object
        from it. This method is used exclusively in .set_settings
        to generate controllers.

        Parameters
        ----------
        configname : str or path
            Path to the file of the settings to be given to the
            returned controller. In addition to the mandatory
            settings required by the controller, the file pointed
            by this path must contain the name of the controller
            class to be instantiated.
        measurements : Sequence of str
            The quantities measured by this controller. Expected to
            be compatible with the argument of .set_measurements in
            MeasureControllerABC.
        is_primary : boolean
            True if the controller is the primary controller.

        Returns
        -------
        controller : ControllerABC
            The controller that is going to handle the hardware.
            May be a subclass of ControllerABC.

        Raises
        ------
        RuntimeError
            If anything goes wrong with the creation of the controller

        Emits
        -----
        MeasurementErrors.MISSING_CLASS_NAME
            If the controller class name is missing.
        MeasurementErrors.INVALID_MEAS_SETTINGS
            If the controller could not be instantiated
            from the given name.
        """
        try:
            config = ViPErLEEDSettings.from_settings(controller_settings)
        except (ValueError, NoSettingsError):
            base.emit_error(self, MeasurementErrors.INVALID_MEAS_SETTINGS,
                            'devices/path to controller configuration')
            raise RuntimeError from None

        invalid = config.has_settings(('controller', 'controller_class'),
                                      ('controller', 'port_name'))
        if invalid:
            base.emit_error(self, MeasurementErrors.MISSING_CLASS_NAME,         # TODO: we do this if the class name AND if the port name is
                            config.last_file)                                   # missing, will need to be edited once we add the serial numbers.
            raise RuntimeError                                                  # Will need to get the port from list_devices.

        cls_name = config['controller']['controller_class']
        port_name = config.get('controller', 'port_name')                       # TODO: do we need this?
        try:
            cls = base.class_from_name('controller', cls_name)
        except ValueError:
            base.emit_error(self, MeasurementErrors.INVALID_MEAS_SETTINGS,      # TODO: isn't this a controller error?
                            'controller/controller_class')
            raise RuntimeError

        controller = cls(settings=config, port_name=port_name,
                         sets_energy=is_primary)
        controller.error_occurred.connect(self.__on_hardware_error)
        controller.set_measurements(measurements)
        return controller

    def __make_primary_ctrl(self):
        """Make primary controller from self.settings."""
        try:
            info = self.settings.getsequence('devices', 'primary_controller')
        except NotASequenceError:
            info = ()

        if len(info) != 2:
            base.emit_error(self, MeasurementErrors.INVALID_MEAS_SETTINGS,
                            'devices/primary_controller')
            return

        try:
            ctrl = self.__make_controller(*info, is_primary=True)
        except RuntimeError:
            # something went wrong with instantiating, and is
            # already reported by emitting in __make_controller
            return

        self.primary_controller = ctrl

    def __make_secondary_ctrls(self):
        """Make secondary controllers from self.settings."""
        # infos should be a sequence with elements
        # of the form (path, (stuff, to, measure))
        try:
            infos = self.settings.getsequence(
                'devices', 'secondary_controllers', fallback=()
                )
        except NotASequenceError:
            base.emit_error(self, MeasurementErrors.INVALID_MEAS_SETTINGS,
                            'devices/secondary_controllers')
            infos = tuple()

        secondary_controllers = []
        for info in infos:
            if not len(info) == 2:
                base.emit_error(self, MeasurementErrors.INVALID_MEAS_SETTINGS,
                                'devices/secondary_controllers')
                continue
            try:
                ctrl = self.__make_controller(*info, is_primary=False)
            except RuntimeError:
                continue

            self.threads.append(qtc.QThread())
            ctrl.moveToThread(self.threads[-1])
            secondary_controllers.append(ctrl)
        self.secondary_controllers = secondary_controllers
        for thread in self.threads:
            thread.start(priority=thread.TimeCriticalPriority)

    def __make_tmp_directory_tree(self):
        """Prepare temporary folder tree where data will be saved."""
        sys_config = get_system_config()
        base_path = sys_config.get("PATHS", "measurements", fallback=None)
        if not base_path:
            # TODO: emit a SystemSettingsError of some kind
            raise RuntimeError("Invalid path in system config file")

        base_path = Path(base_path) / "__tmp__"
        base_path.mkdir(parents=True, exist_ok=True)

        for camera in self.cameras:
            cam_dir = base_path / camera.name
            camera.process_info.base_path = str(cam_dir)
            cam_dir.mkdir(exist_ok=True)

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
