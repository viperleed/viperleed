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

from viperleed.guilib.measure.classes.abc import QObjectSettingsErrors
from viperleed.guilib.measure.measurement.abc import MeasurementABC
from viperleed.guilib.measure.measurement.abc import MeasurementErrors


class IVVideo(MeasurementABC):
    """Measurement class for LEED I(V) videos."""

    display_name = 'I(V) video'
    _mandatory_settings = (
        *MeasurementABC._mandatory_settings,
        ('measurement_settings', 'end_energy'),
        ('measurement_settings', 'delta_energy'),
        )

    @property
    def _delta_energy(self):
        """Return the amplitude of an energy step in eV."""
        # pylint: disable=redefined-variable-type
        # Seems a pylint bug

        # Eventually, this will be an attribute of an energy generator,
        # and it is unclear whether we will actually need it.
        fallback = 0.5
        if not self.settings:
            return fallback
        try:
            delta = self.settings.getfloat('measurement_settings',
                                           'delta_energy')
        except (TypeError, ValueError):
            # Not a float
            delta = fallback
            self.emit_error(QObjectSettingsErrors.INVALID_SETTINGS,
                            'measurement_settings/delta_energy', '')
        return delta

    @property
    def _end_energy(self):
        """Return the energy (in eV) at which the energy ramp ends."""
        # pylint: disable=redefined-variable-type
        # Seems a pylint bug

        # Eventually, this will be an attribute of an energy generator,
        # and it is unclear whether we will actually need it.
        fallback = 0
        if not self.settings:
            return fallback
        try:
            egy = self.settings.getfloat('measurement_settings', 'end_energy')
        except (TypeError, ValueError):
            # Not a float
            egy = fallback
            self.emit_error(QObjectSettingsErrors.INVALID_SETTINGS,
                            'measurement_settings/end_energy', '')
        return egy

    @property
    def _i0_settle_time(self):
        """Return the time interval for the settling of I0."""
        if not self.primary_controller:
            return 0
        return self.primary_controller.i0_settle_time

    @property
    def _n_digits(self):
        """Return the number of digits needed to represent each step."""
        # Used for zero-padding counter in image names.
        num_meas = (1 + round((self._end_energy - self.start_energy)
                              / self._delta_energy))
        return len(str(num_meas))

    # We don't have anything much to do in abort() that is not
    # already done in the ABC, but abort is abstract.
    # pylint: disable-next=useless-super-delegation
    @qtc.pyqtSlot()
    def abort(self):
        """Abort all current actions."""
        super().abort()

    def are_runtime_settings_ok(self):
        """Return whether runtime settings are ok.

        Check if at least one camera is available
        before starting a measurement.

        Returns
        -------
        settings_ok : bool
            True if the runtime settings are
            sufficient to start a measurement.

        Emits
        -----
        error_occurred
            If the primary controller or the camera is missing.
        """
        if not super().are_runtime_settings_ok():
            return False
        if not self.cameras:
            self.emit_error(MeasurementErrors.MISSING_CAMERA)
            return False
        return True

    def begin_next_energy_step(self):
        """Set energy and measure.

        Set energy via the primary controller. Once this is done
        the returned about_to_trigger signal will start the
        measurement on all secondary controllers.

        Returns
        -------
        None.
        """
        super().begin_next_energy_step()
        for device in self.devices:
            # Make all controllers and cameras busy, so we do not risk
            # going to the next energy step too early: the secondary
            # controllers may be not yet busy when the primary becomes
            # not busy; Same is true for cameras, as they may be
            # started later and we may go on without acquiring an image
            device.busy = True

        profile = self.step_profile
        self.set_leed_energy(*profile,
                             self.current_energy, self._i0_settle_time)

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

        if profile_duration:
            txt = "Setting each energy takes"
            print(txt, f"{profile_duration:>{30-len(txt)}.2f} ms")
        for ctrl in self.controllers:
            if not ctrl.measures():
                continue
            ctrl_time = ctrl.time_to_first_measurement + ctrl.time_to_trigger
            txt = f"{ctrl.name} at {ctrl.address}:"
            print(txt, f"{ctrl_time:>{30-len(txt)}.2f} ms")
        for cam in self.cameras:
            txt = f"{cam.name}:"
            cam_time = camera_delay + cam.time_to_image_ready
            print(txt, f"{cam_time:>{30-len(txt)}.2f} ms")

    def get_settings_handler(self):
        """Return a SettingsHandler object for displaying settings."""
        handler = super().get_settings_handler()
        return handler

    @qtc.pyqtSlot(object)
    def set_settings(self, new_settings):
        """Change settings of the measurement."""
        settings_ok = super().set_settings(new_settings)
        self.data_points.time_resolved = False
        return settings_ok

    def _is_finished(self):
        """Check if the full measurement cycle is done.

        If the energy is above the _end_energy the cycle is
        completed. If not, then the delta energy is added
        and the next measurement is started.

        Returns
        -------
        bool
        """
        super()._is_finished()
        if self.current_energy + self._delta_energy > self._end_energy:
            return True
        self.current_energy += self._delta_energy
        return False
