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
    cycle_completed = qtc.pyqtSignal()
    # Is emitted when an error occurs
    error_occurred = qtc.pyqtSignal(tuple)

    # The reimplementation may introduce more/other keys.
    # See ViPErinoController for an example on how to do this.
    plot_info = defaultdict(list)
    plot_info['nominal_energies'] = ['eV', 'lin']
    plot_info['I0'] = ['uA', 'lin']
    plot_info['measured_energies'] = ['eV', 'lin']
    plot_info['elapsed_time'] = ['ms', 'lin']
    plot_info['aux0'] = ['V', 'log']

    def __init__(self, measurement_settings=None, controllers=None,
                 cameras=None):
        """Initialise measurement class"""
        
        self.__settings = measurement_settings
        self.__controllers = controllers
        self.__cameras = cameras
        
        # Connect signals to appropriate functions
        self.connect_controllers(self.__controllers)
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
        return self.__controllers

    @controllers.setter
    def controllers(self, *new_controllers):
        """Set the controllers which should be
        used and handle signals.

        Parameters
        ----------
        new_controllers : list
            List of controller class objects.
        """
        self.disconnect_controllers(self.__controllers)
        self.__controllers = new_controllers
        self.connect_controllers(self.__controllers)

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
        self.__controllers.extend(new_controllers)
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
                self.__controllers.remove(controller)
            except ValueError:
                pass
        self.disconnect_controllers(remove)

    @abstractmethod
    def do_next_measurement(self):
        """Do the next measurement.

        This method has to be reimplemented by subclasses.
        Depending on the measurement type conditions have
        to be introduced which decide if the measurement
        cycle is done or should be continued.

        Returns
        -------
        None.
        """
        return

    def abort_measurement(self):
        """Abort all current actions.

        This function needs to be reimplemented in subclasses.
        Implementation needs to reset all variables used in loop
        operation. Then call super.

        Returns
        -------
        None.
        """

        for controller in self.__controllers:
            controller.abort()
        for key in self.data_points:
            self.data_points[key] = []
            
    def connect_cameras(self, cameras):
        """Connect camera_busy signals to 
        receive_from_camera function.
        """
        for camera in cameras:
            camera.camera_busy.connect(
                self.receive_from_camera, type=qtc.Qt.UniqueConnection
                )    

    def connect_controllers(self, controllers):
        """Connect controller_ready signals to 
        receive_from_controller function.
        """
        for controller in controllers:
            controller.controller_ready.connect(
                self.receive_from_controller, type=qtc.Qt.UniqueConnection
                )

    def disconnect_cameras(self, cameras):
        """Disconnect camera_busy signals 
        from receive_from_camera function.
        """
        for camera in cameras:
            camera.camera_busy.disconnect(
                self.receive_from_camera
                )    

    def disconnect_controllers(self, controllers):
        """Disconnect controller_ready signals 
        from receive_from_controller function.
        """
        for controller in controllers:
            controller.controller_ready.disconnect(
                self.receive_from_controller
                )

    def start_cycle(self):
        """Start measurement cycle.

        Prepare the controllers for a measurement which starts
        the measurement cycle.

        Returns
        -------
        None.
        """
        self.current_energy = self.start_energy
        for controller in self.__controllers:
            controller.trigger_prepare(self.start_energy)

    def ready_for_next_measurement(self):
        """Check if all measurements have been received.

        The next measurement will be started as soon as both
        the image from the camera and the data from the controller
        have been received. If the controller class does not control
        a camera, it will not wait for an image to be returned. If
        the camera is in live mode, the controller class will not
        wait for measurement data from the controller.

        Returns
        -------
        None.
        """

        if any(controller.busy for controller in self.__controllers):
            return
        if any(camera.busy for camera in self.__cameras):
            return
        self.do_next_measurement()

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
