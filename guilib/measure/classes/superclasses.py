"""Module superclasses of viperleed.guilib.measure.classes.

========================================
   ViPErLEED Graphical User Interface
========================================

Created: 2024-04-04
Author: Michele Riva
Author: Florian Doerr

This module contains QObjectWithErrorABC, QObjectWithSettingsABC, and
DeviceABC. QObjectWithErrorABC is the base class of
QObjectWithSettingsABC, Masure, CalibrationTask, and DataPoints.
QObjectWithSettingsABC is the base class of DeviceABC and
MeasurementABC. DeviceABC is the base class of ControllerABC, SerialABC,
and CameraABC.
"""

from PyQt5 import QtCore as qtc

from viperleed.guilib.measure import hardwarebase as base
from viperleed.guilib.measure.classes import settings


class QObjectWithErrorABC(qtc.QObject, metaclass=base.QMetaABC):
    """Metaclass of measurement objects with error detection."""

    # Emitted whenver an error has been detected. Contains
    # information about the error which occurred.
    error_occurred = qtc.pyqtSignal(tuple)


class QObjectWithSettingsABC(QObjectWithErrorABC):
    """Metaclass of measurement objects with settings."""

    def __init__(self, **kwargs):
        """Initialise instance."""
        self._settings = settings.ViPErLEEDSettings()
        super().__init__(**kwargs)


class DeviceABC(QObjectWithSettingsABC):
    """Metaclass of hardware device objects."""

    # Emitted whenever the busy state of the device changes.
    # Contains the busy state the device changed to.
    busy_changed = qtc.pyqtSignal(bool)

    def __init__(self, **kwargs):
        """Initialise instance."""
        # Busy state of the device.
        self._busy = False
        super().__init__(**kwargs)
