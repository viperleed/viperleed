"""Package measurement of viperleed.guilib.measure.

========================================
   ViPErLEED Graphical User Interface
========================================

Created: 2021-10-19
Author: Michele Riva
Author: Florian Doerr

This package contains classes to perform measurements.
"""
from viperleed.guilib.measure.measurement.energy_setpoint import (
    MeasureEnergySetpoint
    )
from viperleed.guilib.measure.measurement.time_resolved import TimeResolved
from viperleed.guilib.measure.measurement.iv_video import IVVideo

ALL_MEASUREMENTS = {cls.display_name: cls for cls in (MeasureEnergySetpoint,
                                                      TimeResolved,
                                                      IVVideo)}
