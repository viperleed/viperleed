"""Tests for module uimeasurement of viperleed.gui.measure."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    'Florian Dörr (@FlorianDoerr)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-03-31'
__license__ = 'GPLv3+'

import json

from types import SimpleNamespace
from pytest_cases import parametrize
from PyQt5 import QtCore as qtc

from viperleed.gui.measure.classes.abc import QObjectSettingsErrors
from viperleed.gui.measure.uimeasurement import _DeviceDetectionWorker
from viperleed.gui.measure.uimeasurement import Measure
from viperleed.gui.measure.uimeasurement import UIErrors


class _FakeSignal:  # pylint: disable=too-few-public-methods
    """A minimal signal-like object."""

    def __init__(self):
        self.connected = []
        self.emitted = 0
        self.disconnected = []

    def connect(self, slot):
        """Store connected slots."""
        self.connected.append(slot)

    def emit(self):
        """Count emissions."""
        self.emitted += 1

    def disconnect(self, slot):
        """Store disconnected slots."""
        self.disconnected.append(slot)


class _FakeAction:  # pylint: disable=too-few-public-methods
    """A minimal action-like object."""

    def __init__(self, text):
        self._text = text
        self._data = None
        self.triggered = _FakeSignal()

    def setData(self, data):
        """Store action data."""
        self._data = data

    @property
    def text(self):
        """Return action text."""
        return self._text

    @property
    def data(self):
        """Return action data."""
        return self._data


class _FakeSubMenu:
    """A minimal menu-like object."""

    def __init__(self):
        self._actions = []
        self.enabled = True

    def clear(self):
        """Clear existing actions."""
        self._actions = []

    def addAction(self, text):
        """Add and return an action."""
        action = _FakeAction(text)
        self._actions.append(action)
        return action

    def actions(self):
        """Return menu actions."""
        return self._actions

    def setEnabled(self, enabled):
        """Set enabled status."""
        self.enabled = enabled


class _FakeDevicesMenu:  # pylint: disable=too-few-public-methods
    """A devices menu exposing two sub-menus via actions()."""

    def __init__(self):
        self.cameras = _FakeSubMenu()
        self.controllers = _FakeSubMenu()
        self._actions = [
            SimpleNamespace(menu=lambda: self.cameras),
            SimpleNamespace(menu=lambda: self.controllers),
            ]

    def actions(self):
        """Return actions for submenus."""
        return self._actions


def test_device_search_allowed_states():
    """Check that active devices and running searches block new searches."""
    fake_measure = SimpleNamespace(running=False)
    camera_viewer = SimpleNamespace(isVisible=lambda: False)
    ctrl_dialog = SimpleNamespace(isVisible=lambda: False)
    fake = SimpleNamespace(
        _device_search_in_progress=False,
        measurement=fake_measure,
        _dialogs={
            'camera_viewers': [camera_viewer],
            'device_settings': {'ctrl': ctrl_dialog},
            },
        )

    assert Measure._device_search_allowed(fake)

    fake._device_search_in_progress = True
    assert not Measure._device_search_allowed(fake)
    fake._device_search_in_progress = False

    fake.measurement.running = True
    assert not Measure._device_search_allowed(fake)
    fake.measurement.running = False

    fake._dialogs['camera_viewers'] = [SimpleNamespace(isVisible=lambda: True)]
    assert not Measure._device_search_allowed(fake)
    fake._dialogs['camera_viewers'] = [camera_viewer]

    fake._dialogs['device_settings'] = {
        'ctrl': SimpleNamespace(isVisible=lambda: True)
        }
    assert not Measure._device_search_allowed(fake)


def test_update_device_lists_blocks_reentry():
    """Check that a second search is blocked while one is in progress."""
    signal = _FakeSignal()
    fake = SimpleNamespace(
        _device_search_in_progress=False,
        detect_devices_requested=signal,
        )
    fake._device_search_allowed = (
        lambda: not fake._device_search_in_progress
        )

    Measure.update_device_lists(fake)
    assert fake._device_search_in_progress
    assert signal.emitted == 1

    Measure.update_device_lists(fake)
    assert signal.emitted == 1


def test_on_devices_detected_updates_menu_and_unblocks_search():
    """Check menu updates with newly detected devices."""
    devices_menu = _FakeDevicesMenu()
    fake = SimpleNamespace(
        _ctrls={'menus': {'devices': devices_menu}},
        _on_camera_clicked=lambda: None,
        _on_controller_clicked=lambda: None,
        _device_search_in_progress=True,
        )
    detected_devices = {
        'camera': {'cam_a': ('camera_cls', 'camera_info')},
        'controller': {'ctrl_a': ('controller_cls', 'controller_info')},
        }

    Measure._on_devices_detected(fake, detected_devices)

    cameras = devices_menu.cameras.actions()
    controllers = devices_menu.controllers.actions()
    assert len(cameras) == 1
    assert len(controllers) == 1
    assert cameras[0].text == 'cam_a'
    assert controllers[0].text == 'ctrl_a'
    assert cameras[0].triggered.connected == [fake._on_camera_clicked]
    assert controllers[0].triggered.connected == [fake._on_controller_clicked]
    assert devices_menu.cameras.enabled
    assert devices_menu.controllers.enabled
    assert not fake._device_search_in_progress


def test_stop_device_search_triggers(mocker):
    """Check shutdown helper disables periodic and queued search triggers."""
    timer_signal = _FakeSignal()
    detect_signal = _FakeSignal()
    stop = mocker.Mock()
    refresh_timer = SimpleNamespace(stop=stop, timeout=timer_signal)
    fake_worker = SimpleNamespace(detect_devices=lambda: None)
    fake = SimpleNamespace(
        _timers={'refresh_devices': refresh_timer},
        detect_devices_requested=detect_signal,
        _device_detection_worker=fake_worker,
        update_device_lists=lambda: None,
        )

    # _stop_device_search_triggers also invokes QMetaObject.invokeMethod
    # to call "stop" on the worker thread. The fake worker is not a
    # QObject, so PyQt5 cannot resolve the overloaded call. Mock it out
    # since this is not what we are testing here.
    mocker.patch(
        'viperleed.gui.measure.uimeasurement.qtc.QMetaObject.invokeMethod'
        )

    Measure._stop_device_search_triggers(fake)

    stop.assert_called_once_with()
    assert timer_signal.disconnected == [fake.update_device_lists]
    assert detect_signal.disconnected == [fake_worker.detect_devices]

def test_detect_devices_does_not_restart_running_process(mocker):
    """A second detection request should be ignored while one is running."""
    proc = mocker.Mock()
    proc.state.return_value = qtc.QProcess.Running

    worker = _DeviceDetectionWorker()
    worker._process = proc

    qprocess = mocker.patch(
        "viperleed.gui.measure.uimeasurement.qtc.QProcess"
    )

    worker.detect_devices()

    qprocess.assert_not_called()


def test_process_finished_nonzero_exit(mocker):
    """A failing subprocess should emit an empty device list."""
    worker = _DeviceDetectionWorker()

    errors = []
    devices = []

    worker.error_occurred.connect(errors.append)
    worker.devices_detected.connect(devices.append)

    proc = mocker.Mock()
    proc.readAllStandardError.return_value.data.return_value = b"failure"
    worker._process = proc

    mocker.patch.object(worker, "stop")

    worker._on_process_finished(1, qtc.QProcess.NormalExit)

    assert devices == [{"camera": {}, "controller": {}}]
    assert errors
    worker.stop.assert_called_once()


def test_process_finished_invalid_json(mocker):
    """Invalid JSON should emit an error and empty device list."""
    worker = _DeviceDetectionWorker()

    errors = []
    devices = []

    worker.error_occurred.connect(errors.append)
    worker.devices_detected.connect(devices.append)

    proc = mocker.Mock()
    proc.readAllStandardOutput.return_value.data.return_value = b"not json"
    worker._process = proc

    mocker.patch.object(worker, "stop")

    worker._on_process_finished(0, qtc.QProcess.NormalExit)

    assert devices == [{"camera": {}, "controller": {}}]
    assert errors
    worker.stop.assert_called_once()


def test_process_finished_json_not_dict(mocker):
    """JSON output must be a dictionary."""
    worker = _DeviceDetectionWorker()

    errors = []
    devices = []

    worker.error_occurred.connect(errors.append)
    worker.devices_detected.connect(devices.append)

    proc = mocker.Mock()
    proc.readAllStandardOutput.return_value.data.return_value = (
        json.dumps([]).encode()
    )
    worker._process = proc

    mocker.patch.object(worker, "stop")

    worker._on_process_finished(0, qtc.QProcess.NormalExit)

    assert devices == [{"camera": {}, "controller": {}}]
    assert errors


def test_process_finished_invalid_device_entry(mocker):
    """Malformed device entries should be skipped."""
    worker = _DeviceDetectionWorker()

    errors = []
    devices = []

    worker.error_occurred.connect(errors.append)
    worker.devices_detected.connect(devices.append)

    proc = mocker.Mock()
    proc.readAllStandardOutput.return_value.data.return_value = json.dumps({
        "camera": [],
        "controller": {
            "success": True,
            "devices": {}
        }
    }).encode()

    worker._process = proc

    mocker.patch.object(worker, "stop")

    worker._on_process_finished(0, qtc.QProcess.NormalExit)

    assert devices == [{"camera": {}, "controller": {}}]
    assert errors


def test_process_finished_import_failure(mocker):
    """Import errors during object reconstruction should be reported."""
    worker = _DeviceDetectionWorker()

    errors = []
    devices = []

    worker.error_occurred.connect(errors.append)
    worker.devices_detected.connect(devices.append)

    proc = mocker.Mock()
    proc.readAllStandardOutput.return_value.data.return_value = json.dumps({
        "controller": {
            "success": True,
            "devices": {
                "ctrl": [
                    "does.not.exist",
                    "Missing",
                    {
                        "unique_name": "TEST",
                        "has_hardware_interface": True,
                        "more": {}
                    }
                ]
            }
        }
    }).encode()

    worker._process = proc

    mocker.patch.object(worker, "stop")

    worker._on_process_finished(0, qtc.QProcess.NormalExit)

    assert devices[0]["controller"] == {}
    assert errors


def test_detection_timeout(mocker):
    """Timeout should stop detection and emit an empty device list."""
    worker = _DeviceDetectionWorker()

    errors = []
    devices = []

    worker.error_occurred.connect(errors.append)
    worker.devices_detected.connect(devices.append)

    proc = mocker.Mock()
    proc.state.return_value = qtc.QProcess.Running

    worker._process = proc

    mocker.patch.object(worker, "stop")

    worker._on_detection_timeout()

    assert devices == [{"camera": {}, "controller": {}}]
    assert errors
    worker.stop.assert_called_once()


def test_process_error(mocker):
    """Process launch errors should emit an empty device list."""
    worker = _DeviceDetectionWorker()

    errors = []
    devices = []

    worker.error_occurred.connect(errors.append)
    worker.devices_detected.connect(devices.append)

    proc = mocker.Mock()
    proc.errorString.return_value = "boom"

    worker._process = proc

    mocker.patch.object(worker, "stop")

    worker._on_process_error(qtc.QProcess.FailedToStart)

    assert devices == [{"camera": {}, "controller": {}}]
    assert errors
    worker.stop.assert_called_once()


def test_stop_running_process(mocker):
    """Running processes should be terminated and cleaned up."""
    worker = _DeviceDetectionWorker()

    proc = mocker.Mock()
    proc.state.return_value = qtc.QProcess.Running

    worker._process = proc
    worker._timeout = None

    worker.stop()

    proc.kill.assert_called_once()
    proc.waitForFinished.assert_called_once_with(1000)
    proc.deleteLater.assert_called_once()
    assert worker._process is None


def test_stop_returns_early_when_no_process():
    """stop() should be a no-op when _process is None."""
    worker = _DeviceDetectionWorker()
    assert worker._process is None
    worker.stop()  # should not raise
    assert worker._process is None


def test_detection_timeout_no_process():
    """_on_detection_timeout should be a no-op when _process is None."""
    worker = _DeviceDetectionWorker()
    errors = []
    devices = []
    worker.error_occurred.connect(errors.append)
    worker.devices_detected.connect(devices.append)
    worker._on_detection_timeout()
    assert not errors
    assert not devices


def test_process_error_no_process():
    """_on_process_error should be a no-op when _process is None."""
    worker = _DeviceDetectionWorker()
    errors = []
    devices = []
    worker.error_occurred.connect(errors.append)
    worker.devices_detected.connect(devices.append)
    worker._on_process_error(qtc.QProcess.FailedToStart)
    assert not errors
    assert not devices


test_cases = (
    (
        {
            'camera': {
                'success': False,
                'error_type': 'DEFAULT_SETTINGS_CORRUPTED',
                'error_msg': 'bad defaults'
            },
            'controller': {
                'success': True,
                'devices': {
                    'test_ctrl': [
                         'types',
                         'SimpleNamespace',
                        {'unique_name': 'TEST',
                         'has_hardware_interface': True,
                         'more': {}}
                        ]
                }
            }
        },
        'test_ctrl',
        {},
        QObjectSettingsErrors.DEFAULT_SETTINGS_CORRUPTED.value[0]

    ),
    (
        {
            'camera': {'success': False, 'error_type': 'RUNTIME_ERROR',
                       'error_msg': 'something went wrong'},
            'controller': {'success': True, 'devices': {}}
        },
        {},
        {},
        UIErrors.RUNTIME_ERROR.value[0]
    )
)


@parametrize('devices, ctrl, camera, error', test_cases)
def test_device_detection_worker(mocker, devices, ctrl, camera, error):
    """Check device detection."""

    class FakeProcess(qtc.QObject):
        # Patch QProcess used by the worker to simulate subprocess.
        NormalExit = qtc.QProcess.NormalExit
        NotRunning = qtc.QProcess.NotRunning
        ExitStatus = qtc.QProcess.ExitStatus
        ProcessError = qtc.QProcess.ProcessError

        finished = qtc.pyqtSignal(int, qtc.QProcess.ExitStatus)
        errorOccurred = qtc.pyqtSignal(qtc.QProcess.ProcessError)

        def __init__(self, parent):
            super().__init__()
            self._finished_callbacks = []
            self._error_callbacks = []

        def start(self, exe, args):
            # Simulate immediate process completion with JSON on stdout.
            data = devices

            self._stdout = json.dumps(data).encode()
            self.finished.emit(0, self.NormalExit)

        def state(self):
            return self.NotRunning

        def readAllStandardOutput(self):
            class B:
                def __init__(self, b):
                    self._b = b
                def data(self):
                    return self._b
            return B(self._stdout)

        def readAllStandardError(self):
            class B:
                def data(self):
                    return b''
            return B()
    mocker.patch('viperleed.gui.measure.uimeasurement.qtc.QProcess',
                 FakeProcess)

    worker = _DeviceDetectionWorker()
    emitted_errors = []
    emitted_devices = []
    worker.error_occurred.connect(emitted_errors.append)
    worker.devices_detected.connect(emitted_devices.append)
    worker.detect_devices()

    if emitted_devices[0]['controller']:
        assert ctrl in emitted_devices[0]['controller']
    assert emitted_devices[0]['camera'] == camera
    assert emitted_errors[0][0] == error
