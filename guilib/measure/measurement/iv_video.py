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
        super().__init__(measurement_settings)

        self.__end_energy = 0
        self.__delta_energy = 1
        self.__hv_settle_time = 0
        self.__i0_settle_time = 0
        if self.settings:
            self.__delta_energy = self.settings.getfloat(
                'measurement_settings', 'delta_energy', fallback=10
                )
            self.__end_energy = self.settings.getfloat(
                'measurement_settings', 'end_energy', fallback=10
                )
        self.__hv_settle_time = self.primary_controller.hv_settle_time
        self.__i0_settle_time = self.primary_controller.i0_settle_time
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
            # Necessary to force secondaries into busy,
            # before the primary returns not busy anymore.
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

