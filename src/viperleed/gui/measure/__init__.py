"""Package measure of viperleed.gui.

This package contains all the functionality necessary to perform
LEED-IV measurements with ViPErLEED.

Modules
-------
hardwarebase
    Functions and classes used in multiple spots within measure.
uimeasurement
    The main window of the measurement plug-in.

Packages
--------
camera
    Functionality for interacting with camera devices.
classes
    Python and Qt objects used in multiple spots with the measure
    package.
controller
    Functionality for handling the logical interface with non-camera
    devices.
dialogs
    Qt dialogs used within the measure package.
measurement
    Functionality for performing data acquisition.
serial
    Functionality for interacting with non-camera devices.
widgets
    Custom Qt widgets used throughout the measure package.

Non-python data
---------------
_defaults
    Collection of standard configuration files for known objects.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    'Florian Dörr (@FlorianDoerr)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2021-12-09'
__license__ = 'GPLv3+'
