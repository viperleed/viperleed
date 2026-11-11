"""Package measurement of viperleed.gui.measure.

This package contains classes to perform measurements.

Modules
-------
abc
    Abstract base classes defining the interface for all measurement
    types.
energy_calibration
    A measurement for calibrating the LEED energy.
iv_video
    A measurement for acquiring LEED-I(V) data.
time_resolved
    Measurements for acquiring data in a not-(only-)energy-resolved
    manner.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    'Florian Dörr (@FlorianDoerr)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2021-10-19'
__license__ = 'GPLv3+'

from viperleed.gui.measure.measurement.energy_calibration import (
    MeasureEnergyCalibration,
    )
from viperleed.gui.measure.measurement.time_resolved import TimeResolved
from viperleed.gui.measure.measurement.iv_video import IVVideo

ALL_MEASUREMENTS = {cls.display_name: cls for cls in (MeasureEnergyCalibration,
                                                      TimeResolved,
                                                      IVVideo)}
