"""Script and functions for running device detection in a separate process.

Because ViPErLEED heavily interacts with OS driver calls (e.g.
QSerialPort.open) during detection, the Python GIL and C-level threads can
stall the main GUI event loop even when discovery is run locally within a
QThread. Running it in a subprocess solves all blocking issues seamlessly.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    'Florian Dörr (@FlorianDoerr)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-04-29'
__license__ = 'GPLv3+'

import dataclasses
from importlib import import_module
import json
import sys

import PyQt5.QtCore as qtc

from viperleed.gui.measure import hardwarebase as base
from viperleed.gui.measure.classes.abc import QObjectSettingsErrors
from viperleed.gui.measure.classes.abc import SettingsInfo
from viperleed.gui.measure.classes.settings import DefaultSettingsError


DEF_CORRUPTED = 'DEFAULT_SETTINGS_CORRUPTED'
ATTR_ERR = 'ATTRIBUTE_ERROR'
RUN_ERR = 'RUNTIME_ERROR'


class JSONEncoderSafe(json.JSONEncoder):
    """Custom JSON encoder that gracefully handles datatypes like Version."""
    def default(self, o):
        try:
            return super().default(o)
        except TypeError:
            return str(o)


class DeviceDetectionErrors(base.ViPErLEEDErrorEnum):
    """Class for errors occurring during device detection."""

    RUNTIME_ERROR = (1100, '{}')
    ATTRIBUTE_ERROR = (1101, '{}')


class DeviceDetectionWorker(qtc.QObject):
    """Worker object for detecting devices in a dedicated thread."""

    devices_detected = qtc.pyqtSignal(object)
    error_occurred = qtc.pyqtSignal(tuple)

    def __init__(self, parent=None):
        super().__init__(parent)
        self._process = None
        self._timeout = None

    @qtc.pyqtSlot()
    def detect_devices(self):
        """Detect all supported device types and emit the result via QProcess.

        Using a separate process guarantees the main event loop
        remains responsive during hardware connection checks.
        """
        # If a search is already running, avoid launching a new process.
        if self._process:
            return

        self._process = qtc.QProcess(self)
        self._process.finished.connect(self._on_process_finished)
        self._process.errorOccurred.connect(self._on_process_error)

        if getattr(sys, 'frozen', False):
            # Bundled executable
            exe = sys.executable
            args = ['--detect-devices']
        else:
            # Standard Python environment
            exe = sys.executable
            args = ['-m', 'viperleed.gui', '--detect-devices']

        self._timeout = qtc.QTimer(self)
        self._timeout.setSingleShot(True)
        self._timeout.timeout.connect(self._on_detection_timeout)
        self._timeout.start(20000)

        self._process.start(exe, args)

    @qtc.pyqtSlot()
    def stop(self):
        """Kill any in-flight detection subprocess and clean up."""
        proc = self._process
        if proc is None:
            return
        if self._timeout:
            self._timeout.stop()
            base.safe_disconnect(self._timeout.timeout,
                                 self._on_detection_timeout)
            self._timeout = None
        base.safe_disconnect(proc.finished, self._on_process_finished)
        base.safe_disconnect(proc.errorOccurred, self._on_process_error)
        if proc.state() != qtc.QProcess.NotRunning:
            proc.kill()
            proc.waitForFinished(1000)
        proc.deleteLater()
        self._process = None

    def _emit_detection_error(self, result):
        """Emit the appropriate detection error from detection results."""
        err_type = result.get('error_type')
        err_msg = result.get('error_msg', '')
        if err_type == DEF_CORRUPTED:
            err = QObjectSettingsErrors.DEFAULT_SETTINGS_CORRUPTED
        elif err_type == ATTR_ERR:
            err = DeviceDetectionErrors.ATTRIBUTE_ERROR
        else:
            err = DeviceDetectionErrors.RUNTIME_ERROR
        base.emit_error(self, err, f'Detection failed: {err_msg}')

    @qtc.pyqtSlot()
    def _on_detection_timeout(self):
        """Kill the subprocess if it takes too long."""
        if self._process is None:
            return
        base.emit_error(self, DeviceDetectionErrors.RUNTIME_ERROR,
                        'Device detection timed out.')
        self.stop()
        self.devices_detected.emit({'camera': {}, 'controller': {}})

    def _parse_detection_output(self):
        """Parse and return the output of the device detection."""
        out = self._process.readAllStandardOutput().data().decode().strip()
        json_str = out.splitlines()[-1] if out else '{}'
        return json.loads(json_str)

    @qtc.pyqtSlot(int, qtc.QProcess.ExitStatus)
    def _on_process_finished(self, exit_code, exit_status):
        detected_out = {'camera': {}, 'controller': {}}

        if exit_code or exit_status != qtc.QProcess.NormalExit:
            error_data = self._process.readAllStandardError().data().decode(
                            errors='replace'
                            )
            base.emit_error(self, DeviceDetectionErrors.RUNTIME_ERROR,
                            f'Detection failed: {error_data}')
            self.devices_detected.emit(detected_out)
            self.stop()
            return

        try:
            parsed = self._parse_detection_output()
        # pylint: disable-next=broad-exception-caught
        except Exception as exc:
            base.emit_error(self, DeviceDetectionErrors.RUNTIME_ERROR,
                            str(exc))
            self.devices_detected.emit(detected_out)
            self.stop()
            return
        # Device discovery prints connection warnings (e.g. Qt or
        # failed COMs) to stdout. We take only the last line, which
        # contains our dumped JSON.

        if not isinstance(parsed, dict):
            base.emit_error(
                self, DeviceDetectionErrors.RUNTIME_ERROR,
                f'Unexpected detection output: {parsed!r}'
                )
            self.devices_detected.emit(detected_out)
            self.stop()
            return

        # Restore objects
        for device_type, result in parsed.items():
            # pylint: disable-next=magic-value-comparison
            if not isinstance(result, dict) or 'success' not in result:
                base.emit_error(
                    self, DeviceDetectionErrors.RUNTIME_ERROR,
                    f'Invalid entry for {device_type!r}: {result!r}'
                    )
                continue
            if result.get('success'):
                try:
                    recreated = self._recreate_devices(result['devices'])
                # pylint: disable-next=broad-exception-caught
                except Exception as exc:
                    base.emit_error(self, DeviceDetectionErrors.RUNTIME_ERROR,
                                    str(exc))
                else:
                    detected_out[device_type] = recreated
            else:
                self._emit_detection_error(result)
        self.devices_detected.emit(detected_out)
        self.stop()

    @qtc.pyqtSlot(qtc.QProcess.ProcessError)
    def _on_process_error(self, err):
        if self._process is None:
            return
        err_str = self._process.errorString() if self._process else str(err)
        base.emit_error(self, DeviceDetectionErrors.RUNTIME_ERROR,
                        f'Subprocess launch error: {err_str}')
        self.devices_detected.emit({'camera': {}, 'controller': {}})
        self.stop()

    def _recreate_devices(self, devices):
        """Recreate and return devices from detection results."""
        recreated = {}
        for name, (mod_name, cls_name, kwargs) in devices.items():
            mod = import_module(mod_name)
            cls = getattr(mod, cls_name)
            recreated[name] = (cls, SettingsInfo(**kwargs))
        return recreated


def run_device_detection():
    """Detect available, supported devices and serialize results as JSON data.

    Returns
    -------
    detected : dict
        A dictionary containing the detection results keyed by device
        type. Can contain either successfully serialized representations
        of SettingsInfo or structured error messages.
    """
    detected = {}

    # Provide CoreApplication instance to process underlying event loops.
    app = qtc.QCoreApplication.instance()
    if app is None:
        app = qtc.QCoreApplication(sys.argv)

    for device_type in ('camera', 'controller'):
        try:
            detected_devs = base.get_devices(device_type)
        except DefaultSettingsError as exc:
            detected[device_type] = {
                'success': False,
                'error_type': DEF_CORRUPTED,
                'error_msg': str(exc),
            }
        except AttributeError as exc:
            detected[device_type] = {
                'success': False,
                'error_type': ATTR_ERR,
                'error_msg': str(exc),
            }
        # pylint: disable-next=broad-exception-caught
        except Exception as exc:
            detected[device_type] = {
                'success': False,
                'error_type': RUN_ERR,
                'error_msg': str(exc),
            }
        else:
            serialized_devs = {}
            for name, (cls, device) in detected_devs.items():
                serialized_devs[name] = (
                    cls.__module__,
                    cls.__name__,
                    dataclasses.asdict(device),
                )
            detected[device_type] = {'success': True,
                                     'devices': serialized_devs}
    return detected
