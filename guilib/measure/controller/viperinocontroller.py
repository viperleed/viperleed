"""ViPErino Controller

========================================
   ViPErLEED Graphical User Interface
========================================

Created: 2021-07-08
Author: Michele Riva
Author: Florian Doerr

This module contains the definition of the ViPErinoController class
which gives commands to the ViPErinoSerialWorker class.
"""
# ViPErLEED modules
from viperleed.guilib.measure.controller.controllerabc import ControllerABC


class ViPErinoController(ControllerABC):
    """Controller class for the ViPErLEED Arduino Micro."""

    def __init__(self, settings=None, port_name=''):
        """Initialise controller object."""

        super().init(settings=settings, port_name=port_name)

    def set_energy(self, energy, time, *more_steps):
        """Convert data from gui to usable values for the DAC.

        Take the energy (or energies), get setpoint energy (or
        energies) and send energy and time to

        Parameters
        ----------
        energy: float
        time: integer
        *more_steps: float and integer (alternating)
        """
