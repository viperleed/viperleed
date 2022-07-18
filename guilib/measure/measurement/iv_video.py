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

from viperleed.guilib.measure import hardwarebase as base
from viperleed.guilib.measure.measurement.abc import (MeasurementABC,
                                                      MeasurementErrors)


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
            # Force all controllers into busy, so we don't
            # risk the secondary controllers to be not yet
            # busy when the primary becomes not busy (which
            # would make us potentially move to the next step)
            controller.busy = True

        profile = self.step_profile
        self.set_leed_energy(*profile,
                             self.current_energy, self.__i0_settle_time)

        image_name = (f"{self.current_step_nr:0>{self.__n_digits}}_"
                      f"{self.current_energy:.1f}eV.tiff")
        for camera in self.cameras:
            camera.process_info.filename = image_name
            self.data_points.add_image_names(image_name)
            # Setting cameras busy prevents going to the next
            # energy before an image has been acquired: in fact,
            # triggering is delayed, using self._camera_timer
            camera.busy = True

        # TODO: here we should start the camera no earlier than
        # hv_settle_time, but such that image acquisition
        # overlaps as much as possible with the measurement time!
        # Frame delivery from the camera should take:
        # (exposure + 1000/fr_rate) + (n_frames - 1) * fr_interval
        # Also, we should probably have one timer per camera, as
        # cameras may potentially deliver frames at different rates!
        profile_duration = sum(profile[1::2])
        camera_delay = profile_duration + self.hv_settle_time
        self._camera_timer.start(camera_delay)

        # Let the user know how much time we are loosing, at
        # the second step, because the first one is a non-step
        if self.current_step_nr != 2:
            return

        for ctrl in self.controllers:
            if not ctrl.measured_quantities:
                continue
            ctrl_time = (ctrl.time_to_first_measurement + ctrl.time_to_trigger)
            txt = f"{ctrl.name} at {ctrl.port_name}:"
            print(txt, f"{ctrl_time:>{30-len(txt)}.2f} ms")
        for cam in self.cameras:
            cam_time = (cam.exposure + 1000/cam.get_frame_rate()
                        + (cam.n_frames - 1) * cam.frame_interval
                        + cam.extra_delay + camera_delay)
            txt = f"{cam.name}:"
            print(txt, f"{cam_time:>{30-len(txt)}.2f} ms")

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
    @qtc.pyqtSlot()
    def abort(self):
        """Abort all current actions."""
        super().abort()
    # pylint: enable=useless-super-delegation

    def set_settings(self, new_settings):
        """Change settings of the measurement.

        Settings are loaded only if they are valid. Otherwise
        the previous settings stay in effect. If the settings
        have been accepted, controller and camera objects as
        specified in the settings will be instantiated, told
        what they will be measuring, moved to their respective
        properties and connected to all necessary signals.

        This extension of the base-class method checks that
        at least one camera is available to the measurement.

        Parameters
        ----------
        new_settings : dict, ConfigParser, string or path
            Configuration of the measurement.

        Returns
        -------
        settings_valid : bool
            True if the new settings given were accepted.

        Raises
        ------
        TypeError
            If new_settings is neither a dict, ConfigParser, string
            or path and if an element of the mandatory_settings is
            None or has a length greater than 3.

        Emits
        -----
        MeasurementErrors.MISSING_SETTINGS
            If new_settings is missing.
        MeasurementErrors.INVALID_SETTINGS
            If any element of the new_settings does not fit the
            mandatory_settings.
        MeasurementErrors.MISSING_CAMERA
            If no camera is available to the IVVideo measurement
        """
        settings_valid = super().set_settings(new_settings)
        if settings_valid and not self.cameras:
            base.emit_error(self, MeasurementErrors.MISSING_CAMERA)
            settings_valid = False
        return settings_valid
