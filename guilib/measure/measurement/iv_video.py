"""Module iv_video of viperleed.?????.
========================================
   ViPErLEED Graphical User Interface
========================================

Created: 2021-07-19
Author: Michele Riva
Author: Florian Doerr

This module contains the definition of the measure_iv_video class
which gives commands to the controller classes.
"""

class measure_iv_video():
    """Measurement class for LEED I(V) videos."""





    def measure_iv_video(self):
        """Measure a ramp of energies.

        This method must be reimplemented in subclasses. This
        function should take measurements while doing an energy
        ramp. It should measure I0 (the current emitted by the
        filament used to create the electron beam) and take
        pictures of the LEED screen for every energy step.

        The energies at the beginning and at the end of the ramp
        have to be specified in the settings property derived
        from the configuration file. Furthermore, the delta
        energy between two energy steps (step height) has to be
        given as well.

        Parameters
        ----------
        current_energy

        Returns
        -------
        None.
        """
        # Edit docstring
        # Do stuff here
        self.set_energy(current_energy, settle_time)
        # # Need to wait for camera class
        # self.trigger_camera()

        self.current_energy += self.delta_energy
        if self.current_energy > self.end_energy:
            self.cycle_completed.emit()
        return