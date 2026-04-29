"""Script and functions for running device detection in a separate process.

Because ViPErLEED heavily interacts with OS driver calls (e.g. QSerialPort.open)
during detection, the Python GIL and C-level threads can stall the main GUI event loop
even when discovery is run locally within a QThread. Running it in a subprocess
solves all blocking issues seamlessly.
"""

import dataclasses
import json
import sys
from PyQt5.QtCore import QCoreApplication

from viperleed.gui.measure import hardwarebase as base
from viperleed.gui.measure.classes.settings import DefaultSettingsError


class JSONEncoderSafe(json.JSONEncoder):
    """Custom JSON encoder that gracefully handles datatypes like Version."""
    def default(self, o):
        try:
            return super().default(o)
        except TypeError:
            return str(o)

def run_device_detection():
    """Detect available, supported devices and serialize results as JSON data.

    Returns
    -------
    detected : dict
        A dictionary containing the detection results keyed by device type.
        Can contain either successfully serialized representations of SettingsInfo
        or structured error messages.
    """
    detected = {}

    # Qt networking and serial interfaces (e.g. QSerialPort, QThread, waitForReadyRead)
    # strictly require a QCoreApplication instance to process their underlying event loops.
    # Because we're executing headlessly in a new process, we must provide one.
    app = QCoreApplication.instance()
    if app is None:
        app = QCoreApplication(sys.argv)

    for device_type in ('camera', 'controller'):
        try:
            detected_devs = base.get_devices(device_type)
            serialized_devs = {}
            for name, (cls, device) in detected_devs.items():
                serialized_devs[name] = (
                    cls.__module__,
                    cls.__name__,
                    dataclasses.asdict(device),
                )
            detected[device_type] = {"success": True, "devices": serialized_devs}
        except DefaultSettingsError as exc:
            detected[device_type] = {
                "success": False,
                "error_type": "DEFAULT_SETTINGS_CORRUPTED",
                "error_msg": str(exc),
            }
        except Exception as exc:  # pylint: disable=broad-exception-caught
            detected[device_type] = {
                "success": False,
                "error_type": "RUNTIME_ERROR",
                "error_msg": str(exc),
            }
    return detected
