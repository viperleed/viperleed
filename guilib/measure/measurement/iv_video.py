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

from viperleed.guilib.measure import hardwarebase as base
from viperleed.guilib.measure.measurement.abc import (MeasurementABC,
                                                      MeasurementErrors)


# TODO: complain if started without a camera

class IVVideo(MeasurementABC):
    """Measurement class for LEED I(V) videos."""

    display_name = 'I(V) video'
    _mandatory_settings = [*MeasurementABC._mandatory_settings,
                           ('measurement_settings', 'end_energy'),
                           ('measurement_settings', 'delta_energy'),]

    def __init__(self, measurement_settings):
        """Initialise measurement instance."""
        super().__init__(measurement_settings)
        self.data_points.time_resolved = False

    @property
    def __delta_energy(self):
        """Return the amplitude of an energy step in eV."""
        # pylint: disable=redefined-variable-type
        # Seems a pylint bug

        # Eventually, this will be an attribute of an energy generator,
        # and it is unclear whether we will actualy need it.
        fallback = 0.5
        if not self.settings:
            return fallback
        try:
            delta = self.settings.getfloat('measurement_settings',
                                           'delta_energy')
        except (TypeError, ValueError):
            # Not a float
            delta = fallback
            base.emit_error(self, MeasurementErrors.INVALID_SETTINGS,
                            'measurement_settings/delta_energy')
        return delta

    @property
    def __end_energy(self):
        """Return the energy (in eV) at which the energy ramp ends."""
        # pylint: disable=redefined-variable-type
        # Seems a pylint bug

        # Eventually, this will be an attribute of an energy generator,
        # and it is unclear whether we will actualy need it.
        fallback = 0
        if not self.settings:
            return fallback
        try:
            egy = self.settings.getfloat('measurement_settings', 'end_energy')
        except (TypeError, ValueError):
            # Not a float
            egy = fallback
            base.emit_error(self, MeasurementErrors.INVALID_SETTINGS,
                            'measurement_settings/end_energy')
        return egy

    @property
    def __i0_settle_time(self):
        """Return the time interval for the settling of I0."""
        if not self.primary_controller:
            return 0
        return self.primary_controller.i0_settle_time

    @property
    def __n_digits(self):
        """Return the number of digts needed to represent each step."""
        # Used for zero-padding counter in image names.
        num_meas = (1 + round((self.__end_energy - self.start_energy)
                              / self.__delta_energy))
        return len(str(num_meas))

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
        self.set_leed_energy(self.current_energy, self.__i0_settle_time)
        image_name = (f"{self.current_step_nr:0>{self.__n_digits}}_"
                      f"{self.current_energy:.1f}eV_.tiff")
        for camera in self.cameras:
            camera.process_info.filename = image_name
            self.data_points.add_image_names(image_name)
        # TODO: here we should start the camera no earlier than
        # hv_settle_time, but such that image acquisition
        # overlaps as much as possible with the measurement time!
        self._camera_timer.start(self.hv_settle_time)

    def _is_finished(self):
        """Check if the full measurement cycle is done.

        If the energy is above the __end_energy the cycle is
        completed. If not, then the delta energy is added
        and the next measurement is started.

        Returns
        -------
        bool
        """
        super()._is_finished()
        if self.current_energy >= self.__end_energy:
            return True
        self.current_energy += self.__delta_energy
        return False

    # pylint: disable=useless-super-delegation
    # We don't have anything much to do in abort() that is not
    # already done in the ABC, but abort is abstract.
    def abort(self):
        """Abort all current actions."""
        super().abort()
    # pylint: enable=useless-super-delegation
