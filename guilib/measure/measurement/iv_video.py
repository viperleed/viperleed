"""Module iv_video of viperleed.guilib.measure.measurement.

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
from viperleed.guilib.measure.measurement.abc import MeasurementABC


class IVVideo(MeasurementABC):
    """Measurement class for LEED I(V) videos."""

    display_name = 'I(V) video'

    def __init__(self, measurement_settings):
        """Initialise measurement class."""
        self.camera_timer = qtc.QTimer()
        self.camera_timer.setSingleShot(True)

        super().__init__(measurement_settings)
        self.camera_timer.setParent(self)
        self.__end_energy = 0
        self.__delta_energy = 1
        self.__hv_settle_time = 0
        self.__i0_settle_time = 0
        if self.settings:
            self.__end_energy = self.settings.getfloat('measurement_settings',
                                                       'end_energy')
            self.__delta_energy = self.settings.getfloat(
                'measurement_settings', 'delta_energy'
                )
            self.__hv_settle_time = self.primary_controller.settings.getint(
                'measurement_settings', 'hv_settle_time'
                )
            self.__i0_settle_time = self.primary_controller.settings.getint(
                'measurement_settings', 'i0_settle_time'
                )
        num_meas = (1 + round((self.__end_energy - self.start_energy)
                              / self.__delta_energy))
        self.__n_digits = len(str(num_meas))

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
        for controller in self.controllers:
            controller.busy = True
        self.set_LEED_energy(self.current_energy, self.__i0_settle_time)
        self.counter += 1
        image_name = (f"{self.counter:0>{self.__n_digits}}_"
                      f"{self.current_energy:.1f}eV_.tiff")
        for i, camera in enumerate(self.cameras):
            camera.process_info.filename = image_name
            self.data_points.add_image_names(image_name)
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
        self.new_data_available.emit()
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

    def connect_cameras(self):
        """Connect necessary camera signals."""
        for camera in self.cameras:
            camera.camera_busy.connect(self.receive_from_camera,
                                       type=qtc.Qt.UniqueConnection)
            self.camera_timer.timeout.connect(camera.trigger_now,
                                              type=qtc.Qt.UniqueConnection)
            self.begin_preparation.connect(camera.start,
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
