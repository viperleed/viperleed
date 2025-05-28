"""Package measurement of viperleed.gui.measure.

This package contains classes to perform measurements.
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
