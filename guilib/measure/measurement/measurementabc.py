"""Module measurementabc of viperleed.?????.
========================================
   ViPErLEED Graphical User Interface
========================================

Created: 2021-07-19
Author: Michele Riva
Author: Florian Doerr

This module contains the definition of the measurementABC class
which gives commands to the controller classes and is inherited
by the other measurement classes.
"""
from collections import defaultdict
from abc import abstractmethod

from PyQt5 import QtCore as qtc

# ViPErLEED modules
from viperleed.guilib.measure.hardwarebase import (
    emit_error, ViPErLEEDErrorEnum, QMetaABC
    )


class MeasurementErrors(ViPErLEEDErrorEnum):
    """Errors that might occur during a measurement cycle."""
    INVALID_MEASUREMENT = (300,
                           "The returned data dictionary contained a section "
                           "that was not specified in the measurement class.")


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
    plot_info['measured_energy'] = ['eV', 'lin']
    plot_info['elapsed_time'] = ['ms', 'lin']
    plot_info['Isample'] = ['V', 'log']
    plot_info['temperature'] = ('°C', 'lin')
    plot_info['cold_junction'] = ('°C', 'lin')

    def __init__(self, measurement_settings=None, primary_controller=None,
                 controllers=None, cameras=None):
        """Initialise measurement class"""

        self.__settings = measurement_settings
        # Primary controller sets the energy.
        self.__primary_controller = primary_controller
        self.__primary_controller.sets_energy = True
        self.__secondary_controllers = controllers
        for controller in self.__secondary_controllers:
            controller.sets_energy = False
        self.__cameras = cameras

        # Connect signals to appropriate functions
        self.connect_primary_controller()
        self.connect_controllers(*self.__secondary_controllers)
        self.connect_cameras(self.__cameras)

        # The reimplementation may introduce more/other keys.
        self.data_points = defaultdict(list)
        for key in self.plot_info:
            self.data_points[key] = []

        # Attributes needed for loop operation
        self.current_energy = 0
        self.__start_energy = self.settings.getfloat('measurement_settings',
                                                     'start_energy')
        self.__end_energy = self.settings.getfloat('measurement_settings',
                                                   'end_energy')
        self.__delta_energy = self.settings.getfloat('measurement_settings',
                                                     'delta_energy')
        self.__settle_time = self.settings.getint('measurement_settings',
                                                  'dac_settle_time')
        self.__long_settle_time = self.settings.getint('measurement_settings',
                                                       'dac_first_settle_time')

    @property
    def cameras(self):
        """Return the cameras used by this class."""
        return self.__cameras

    @cameras.setter
    def cameras(self, *new_cameras):
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

    def add_cameras(self, *new_cameras):
        """Extend camera list by one or more
        cameras and connect signals.

        Parameters
        ----------
        new_cameras : list
            List of camera class objects.

        Returns
        -------
        None.
        """
        self.__cameras.extend(new_cameras)
        self.connect_cameras(new_cameras)

    def remove_cameras(self, *remove):
        """Remove unwanted cameras and disconnect signals.

        Parameters
        ----------
        remove : list
            List of camera class objects.

        Returns
        -------
        None.
        """
        for camera in remove:
            try:
                self.__cameras.remove(camera)
            except ValueError:
                pass
        self.disconnect_cameras(remove)

    @property
    def controllers(self):
        """Return the controllers used by this class."""
        return self.__secondary_controllers

    @controllers.setter
    def controllers(self, *new_controllers):
        """Set the controllers which should be
        used and handle signals.

        Parameters
        ----------
        new_controllers : list
            List of controller class objects.
        """
        self.disconnect_controllers(self.__secondary_controllers)
        self.__secondary_controllers = new_controllers
        for controller in self.__secondary_controllers:
            controller.sets_energy = False
        self.connect_controllers(self.__secondary_controllers)

    def add_controllers(self, *new_controllers):
        """Extend controller list by one or more
        controllers and connect signals.

        Parameters
        ----------
        new_controllers : list
            List of controller class objects.

        Returns
        -------
        None.
        """
        self.__secondary_controllers.extend(new_controllers)
        for controller in new_controllers:
            controller.sets_energy = False
        self.connect_controllers(new_controllers)

    def remove_controllers(self, *remove):
        """Remove unwanted controllers and disconnect signals.

        Parameters
        ----------
        remove : list
            List of controller class objects.

        Returns
        -------
        None.
        """
        for controller in remove:
            try:
                self.__secondary_controllers.remove(controller)
            except ValueError:
                pass
        self.disconnect_controllers(remove)

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
        self.__primary_controller.sets_energy = True
        self.connect_primary_controller()

    @abstractmethod
    def is_finished(self):
        """Check if the full measurement cycle is done.

        This function must be reimplemented in subclasses. It
        should check if the measurement cycle is done via the
        settings. If not it should call start_next_measurement.
        If the whole measurement is done it should emit a .finished()
        signal.

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
        self.finished.emit(self.plot_info, self.data_points)
        self.save_data()
        self.current_energy = 0
        # Set LEED energy to 0.
        self.set_LEED_energy(self.current_energy)
        for key in self.data_points:
            self.data_points[key] = []

    @abstractmethod
    def save_data(self):
        """Save data if required.

        Needs to reimplemented in subclasses. This function
        can be a no-op if saving data is not necessary.

        Returns
        -------
        None.
        """
        return

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
        self.__primary_controller.set_energy(message)

    def do_next_measurement(self):
        """Do the next measurement.

        Start measuring on all controllers and cameras as soon
        as the about_to_trigger signal has been received from the primary
        controller.

        Returns
        -------
        None.
        """
        self.ready_for_measurement.emit()

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
        # Set LEED energy to 0.
        self.set_LEED_energy(self.current_energy)
        for key in self.data_points:
            self.data_points[key] = []

    @abstractmethod
    def connect_cameras(self, cameras=None):
        """Connect necessary camera signals."""

        if not cameras:
            cameras = self.__cameras
        for camera in cameras:
            camera.camera_busy.connect(
                self.receive_from_camera, type=qtc.Qt.UniqueConnection
                )
            self.ready_for_measurement.connect(
                camera.TODO, type=qtc.Qt.UniqueConnection
                # TODO: ^ function?
                )
            self.abort_action.connect(
                camera.abort, type=qtc.Qt.UniqueConnection
                # TODO: ^ function?
                )
            # TODO: may need to connect preparation.

    @abstractmethod
    def connect_controllers(self, controllers=None):
        """Connect necessary controller signals."""

        if not controllers:
            controllers = self.__secondary_controllers
        for controller in controllers:
            controller.data_ready.connect(
                self.receive_from_controller, type=qtc.Qt.UniqueConnection
                )
            self.ready_for_measurement.connect(
                controller.measure_now, type=qtc.Qt.UniqueConnection
                )
            self.abort_action.connect(
                controller.abort, type=qtc.Qt.UniqueConnection
                )
            self.begin_preparation.connect(
                controller.trigger_begin_preparation, type=qtc.Qt.UniqueConnection
                )
            self.continue_preparation.connect(
                controller.trigger_continue_preparation, type=qtc.Qt.UniqueConnection
                )

    @abstractmethod
    def connect_primary_controller(self):
        """Connect signals of the primary controller."""

        self.__primary_controller.data_ready.connect(
            self.receive_from_controller, type=qtc.Qt.UniqueConnection
            )
        self.abort_action.connect(
            self.__primary_controller.abort, type=qtc.Qt.UniqueConnection
            )
        self.begin_preparation.connect(
            self.__primary_controller.trigger_begin_preparation,
            type=qtc.Qt.UniqueConnection
            )
        self.continue_preparation.connect(
            self.__primary_controller.trigger_continue_preparation,
            type=qtc.Qt.UniqueConnection
            )
        # Connections below depend on measurement type
        # and are exclusive to each other. Implemented in subclass.
        self.__primary_controller.about_to_trigger.connect(
            self.start_next_measurement, type=qtc.Qt.UniqueConnection
            )
        # TODO: ^ needed in preparation!!!
        # self.ready_for_measurement.connect(
            # self.__primary_controller.measure_now,
            # type=qtc.Qt.UniqueConnection
            # )

    @abstractmethod
    def disconnect_cameras(self, cameras):
        """Disconnect necessary camera signals."""

        for camera in cameras:
            camera.camera_busy.disconnect(
                self.receive_from_camera
                )
            self.ready_for_measurement.disconnect(camera.TODO)
                # TODO: ^ function?
            self.abort_action.connect(camera.abort)
            # TODO: ^ function?
            # TODO: may need to disconnect preparation.

    @abstractmethod
    def disconnect_controllers(self, controllers):
        """Disconnect necessary controller signals."""

        for controller in controllers:
            controller.data_ready.disconnect(
                self.receive_from_controller
                )
            self.ready_for_measurement.disconnect(controller.measure_now)
            self.abort_action.disconnect(controller.abort)
            self.begin_preparation.disconnect(controller.trigger_begin_preparation)
            self.continue_preparation.disconnect(
                controller.trigger_continue_preparation
                )

    @abstractmethod
    def disconnect_primary_controller(self):
        """Disconnect signals of the primary controller."""

        self.__primary_controller.data_ready.disconnect(
            self.receive_from_controller
            )
        self.abort_action.disconnect(self.__primary_controller.abort)
        self.begin_preparation.disconnect(
            self.__primary_controller.trigger_begin_preparation
            )
        self.continue_preparation.disconnect(
            self.__primary_controller.trigger_continue_preparation
            )
        # Connections below depend on measurement type
        # and are exclusive to each other. Implemented in subclass.
        self.__primary_controller.about_to_trigger.disconnect()
        # TODO: ^ needed in preparation!!!
        # self.ready_for_measurement.disconnect(
            # self.__primary_controller.measure_now
            # )

    @abstractmethod
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
        the data_ready signal has to be disconnected aswell.
        
        In order to know if the first preparation step has
        been done, the controller_busy must be connected and
        returned as soon as all points in the begin_prepare_todos
        dictionary have been executed.
        
        Before the second step starts, this functions is
        called again. It will reconnect all disconnected
        signals and disconnect controller_busy. The second
        preparation step has to emit a data_ready signal
        containing an empty dictionary at its end. This will 
        function like an ordinary measurement and as soon as 
        all controllers are done preparing this class will
        continue with the measurement cycle.
        
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
                self.start_next_measurement, type=qtc.Qt.UniqueConnection
                )
        except TypeError:
            self.__primary_controller.about_to_trigger.disconnect()
            self.__primary_controller.data_ready.disconnect()
            for controller in self.__secondary_controllers:
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
                controller.controller_busy.disconnect()
                controller.data_ready.connect(
                    self.receive_from_controller, type=qtc.Qt.UniqueConnection
                    )
            self.__primary_controller.data_ready.connect(
                self.receive_from_controller, type=qtc.Qt.UniqueConnection
                )

    def begin_measurement_preparation(self):
        """Start preparation for measurements.

        Prepare the controllers for a measurement which starts
        the measurement cycle. This function only performs the first part.
        (Everything that is done before the starting energy is set.)

        Returns
        -------
        None.
        """
        self.switch_signals_for_preparation()
        self.current_energy = self.start_energy
        self.begin_preparation.emit(self.start_energy,
                                    self.__long_settle_time)

    def continue_measurement_preparation(self, busy):
        """Continue preparation for measurements.

        Prepare the controllers for a measurement which starts
        the measurement cycle. This function only performs the second part.
        (Everything that is done after the starting energy is set.)
        
        Parameters
        ----------
        busy : bool
            Busy state of the emitting controller.

        Returns
        -------
        None.
        """
        if busy:
            return
        if self.__primary_controller.busy:
            return
        if any(controller.busy for controller in self.__secondary_controllers):
            return
        self.switch_signals_for_preparation()
        self.continue_preparation.emit()

    def ready_for_next_measurement(self):
        """Check if all measurements have been received.

        After all measurements have been received, a check if the
        loop is done will be called.

        Returns
        -------
        None.
        """

        if any(controller.busy for controller in self.__secondary_controllers):
            return
        if any(camera.busy for camera in self.__cameras):
            return
        if self.is_finished():
            self.finalize()
        else:
            self.start_next_measurement()
        # TODO: make this work

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
        """
        for key in receive:
            if key not in self.data_points.keys():
                emit_error(self, MeasurementErrors.INVALID_MEASUREMENT)
            else:
                self.data_points[key].append(receive[key])
        self.ready_for_next_measurement()
