"""Module superclasses of viperleed.guilib.measure.classes.

========================================
   ViPErLEED Graphical User Interface
========================================

Created: 2024-04-04
Author: Michele Riva
Author: Florian Doerr

This module contains QObjectWithErrorABC, QObjectWithSettingsABC, and
HardwareABC. QObjectWithErrorABC is the base class of
QObjectWithSettingsABC, Masure, CalibrationTask, and DataPoints.
QObjectWithSettingsABC is the base class of HardwareABC and
MeasurementABC. HardwareABC is the base class of DeviceABC and
SerialABC. DeviceABC is the base class of ControllerABC and CameraABC.
"""

from PyQt5 import QtCore as qtc

from viperleed.guilib.measure import hardwarebase as base
from viperleed.guilib.measure.classes import settings


class QObjectWithErrorABC(qtc.QObject, metaclass=base.QMetaABC):
    """Abstract base class of measurement objects with error detection."""

    # Emitted whenver an error has been detected. Contains
    # information about the error which occurred.
    error_occurred = qtc.pyqtSignal(tuple)


class QObjectWithSettingsABC(QObjectWithErrorABC):
    """Abstract base class of measurement objects with settings."""

    def __init__(self, **kwargs):
        """Initialise instance."""
        self._settings = settings.ViPErLEEDSettings()
        super().__init__(**kwargs)
    # TODO: add set_settings and find_settings


class HardwareABC(QObjectWithSettingsABC):
    """Abstract base class of hardware related objects."""

    # Emitted whenever the busy state of the device changes.
    # Contains the busy state the device changed to.
    busy_changed = qtc.pyqtSignal(bool)

    def __init__(self, **kwargs):
        """Initialise instance."""
        # Busy state of the device.
        self._busy = False
        super().__init__(**kwargs)


class DeviceABC(HardwareABC):
    """Abstract base class of hardware device objects."""
    # TODO: add list_devices
