"""Package classes of viperleed.gui.measure.

This package contains generic classes used throughout the measure
module, and that do not belong to other subpackages (i.e., are not
necessarily only camera-, controller-, measurement-, etc-related).

Modules
-------
abc
    Abstract base classes defining the interface of (most) other
    objects in the measure package.
calibrationtask
    Abstract base classes defining the interface for actions aimed
    at performing some calibration for a device.
datapoints
    Objects for storage of information about acquired data, as well
    as for reading/writing them from/to file.
ioarduinocli
    Functionality for installing, updating, and using the Arduino
    command-line interface for uploading firmware to Arduino boards.
settings
    Functionality for reading/writing configuration files.
thermocouple
    Functionality for reading temperatures from thermocouples.

Non-python data
---------------
thermocouple_coefficients.txt
    Information of calibration curves for thermocouples. Needed for
    the thermocouple module.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2022-02-16'
__license__ = 'GPLv3+'
