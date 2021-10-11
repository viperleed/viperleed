"""Module iv_video of viperleed.
========================================
   ViPErLEED Graphical User Interface
========================================

Created: 2021-07-19
Author: Michele Riva
Author: Florian Doerr

This module contains the definition of the IVVideo class
which gives commands to the controller classes.
"""

from PyQt5 import QtCore as qtc

# ViPErLEED modules
from viperleed.guilib.measure.measurement.measurementabc import MeasurementABC


class IVVideo(MeasurementABC):
    """Measurement class for LEED I(V) videos."""

    def __init__(self, measurement_settings):
        """Initialise measurement class."""
        super().__init__(measurement_settings)
        self.__end_energy = self.settings.getfloat('measurement_settings',
                                                   'end_energy')
        self.__delta_energy = self.settings.getfloat('measurement_settings',
                                                     'delta_energy')
        self.__hv_settle_time = self.primary_controller.settings.getint(
            'measurement_settings', 'hv_settle_time')
        self.__i0_settle_time = self.primary_controller.settings.getint(
            'measurement_settings', 'i0_settle_time')
        self.camera_timer = qtc.QTimer(parent=self)
        self.camera_timer.setSingleShot(True)

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
        self.data_points['nominal_energy'].append(self.current_energy)
        self.set_LEED_energy(self.current_energy, self.__i0_settle_time)
        self.camera_timer.start(self.__hv_settle_time)

    def is_finished(self):
        """Check if the full measurement cycle is done.

        If the energy is above the __end_energy the cycle is
        completed. If not, then the delta energy is added
        and the next measurement is started.

        Returns
        -------
        bool
        """
        if self.current_energy >= self.__end_energy:
            return True
        self.current_energy += self.__delta_energy
        return False

    def abort(self):
        """Abort all current actions.

        Abort and reset all variables.

        Returns
        -------
        None.
        """
        super().abort()

    def connect_cameras(self, cameras=None):
        """Connect necessary camera signals."""
        if not cameras:
            cameras = self.cameras
        for camera in cameras:
            if camera:
                camera.camera_busy.connect(
                    self.receive_from_camera, type=qtc.Qt.UniqueConnection
                    )
                self.camera_timer.timeout.connect(
                    camera.trigger_now,
                    type=qtc.Qt.UniqueConnection
                    )
        # camera.disconnect does not need to be hooked up to the
        # abort_action signal as it is called in the disconnecting
        # of the camera signals anyway.

