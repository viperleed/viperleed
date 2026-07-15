"""Tests for module uimeasurement of viperleed.gui.measure."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    'Florian Dörr (@FlorianDoerr)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-03-31'
__license__ = 'GPLv3+'

import json

from PyQt5 import QtCore as qtc
from pytest_cases import parametrize

from viperleed.gui.measure.classes.abc import QObjectSettingsErrors
from viperleed.gui.measure.devicedetection import DEF_CORRUPTED
from viperleed.gui.measure.devicedetection import RUN_ERR
from viperleed.gui.measure.devicedetection import DeviceDetectionErrors
from viperleed.gui.measure.devicedetection import DeviceDetectionWorker
from viperleed.gui.measure.uimeasurement import Measure

_QPROCESS = 'viperleed.gui.measure.uimeasurement.qtc.QProcess'

# pylint: disable=protected-access
class _FakeSignal:
    """A minimal signal-like object."""

    def __init__(self):
        """Initialize fake signal."""
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


class _FakeAction:
    """A minimal action-like object."""

    def __init__(self, text):
        """Initialize fake action."""
        self._text = text
        self._data = None
        self.triggered = _FakeSignal()

    def setData(self, data):    # pylint: disable=invalid-name
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
        """Initialize fake menu object."""
        self._actions = []
        self.enabled = True

    def clear(self):
        """Clear existing actions."""
        self._actions = []

    def addAction(self, text):      # pylint: disable=invalid-name
        """Add and return an action."""
        action = _FakeAction(text)
        self._actions.append(action)
        return action

    def actions(self):
        """Return menu actions."""
        return self._actions

    def setEnabled(self, enabled):  # pylint: disable=invalid-name
        """Set enabled status."""
        self.enabled = enabled


class _FakeDevicesMenu:  # pylint: disable=too-few-public-methods
    """A devices menu exposing two sub-menus via actions()."""

    def __init__(self, mocker):
        """Initialize fake device menu."""
        self.cameras = _FakeSubMenu()
        self.controllers = _FakeSubMenu()
        cam_action = mocker.Mock()
        cam_action.menu.return_value = self.cameras
        ctrl_action = mocker.Mock()
        ctrl_action.menu.return_value = self.controllers
        self._actions = [cam_action, ctrl_action]
    def actions(self):
        """Return actions for submenus."""
        return self._actions


def test_device_search_allowed_states(mocker):
    """Check that running searches and measurements block new searches."""
    fake_measure = mocker.MagicMock(running=False)
    camera_viewer = mocker.MagicMock()
    camera_viewer.isVisible.return_value = False
    ctrl_dialog = mocker.MagicMock()
    ctrl_dialog.isVisible.return_value = False
    fake = mocker.MagicMock(_device_search_in_progress=False,
                            measurement=fake_measure,
                            _dialogs={
                                'camera_viewers': [camera_viewer],
                                'device_settings': {'ctrl': ctrl_dialog}})

    assert Measure._device_search_allowed(fake)

    fake._device_search_in_progress = True
    assert not Measure._device_search_allowed(fake)
    fake._device_search_in_progress = False

    fake.measurement.running = True
    assert not Measure._device_search_allowed(fake)
    fake.measurement.running = False

    camera_viewer.isVisible.return_value = True
    fake._dialogs['camera_viewers'] = [camera_viewer]
    assert Measure._device_search_allowed(fake)
    fake._dialogs['camera_viewers'] = [camera_viewer]

    ctrl_dialog.isVisible.return_value = True
    fake._dialogs['device_settings'] = {'ctrl': ctrl_dialog}
    assert Measure._device_search_allowed(fake)


def test_update_device_lists_blocks_reentry(mocker):
    """Check that a second search is blocked while one is in progress."""
    signal = _FakeSignal()
    fake = mocker.Mock(_device_search_in_progress=False,
                       detect_devices_requested=signal)
    fake._device_search_allowed = mocker.Mock(side_effect=[True, False])

    Measure.update_device_lists(fake)
    assert fake._device_search_in_progress
    assert signal.emitted == 1

    Measure.update_device_lists(fake)
    assert signal.emitted == 1


def test_on_devices_detected_updates_menu_and_unblocks_search(mocker):
    """Check menu updates with newly detected devices."""
    devices_menu = _FakeDevicesMenu(mocker)
    fake = mocker.Mock(_ctrls={'menus': {'devices': devices_menu}},
                       _device_search_in_progress=True)
    fake._on_camera_clicked = mocker.Mock()
    fake._on_controller_clicked = mocker.Mock()
    detected_devices = {
        'camera': {'cam_a': ('camera_cls', 'camera_info')},
        'controller': {'ctrl_a': ('controller_cls', 'controller_info')},
        }

    Measure._on_devices_detected(fake, detected_devices)

    cameras = devices_menu.cameras.actions()
    controllers = devices_menu.controllers.actions()
    assert len(cameras) == 1
    assert len(controllers) == 1
    # pylint: disable-next=magic-value-comparison
    assert cameras[0].text == 'cam_a'
    # pylint: disable-next=magic-value-comparison
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
    timer_stop = mocker.Mock()
    worker_stop = mocker.Mock()
    delete_later = mocker.Mock()
    refresh_timer = mocker.Mock(stop=timer_stop, timeout=timer_signal)
    fake_worker = mocker.Mock(stop=worker_stop, deleteLater=delete_later)
    fake_worker.detect_devices = mocker.Mock()
    fake = mocker.Mock(_timers={'refresh_devices': refresh_timer},
                       detect_devices_requested=detect_signal,
                       _device_detection_worker=fake_worker)
    fake.update_device_lists = mocker.Mock()
    Measure._stop_device_search_triggers(fake)

    timer_stop.assert_called_once_with()
    worker_stop.assert_called_once_with()
    assert timer_signal.disconnected == [fake.update_device_lists]
    assert detect_signal.disconnected == [fake_worker.detect_devices]


def test_detect_devices_does_not_restart_running_process(mocker):
    """A second detection request should be ignored while one is running."""
    proc = mocker.Mock()
    proc.state.return_value = qtc.QProcess.Running

    worker = DeviceDetectionWorker()
    worker._process = proc
    qprocess = mocker.patch(_QPROCESS)

    worker.detect_devices()

    qprocess.assert_not_called()


@parametrize(
    'exit_code, stdout, stderr',
    [
        (1, b'', b'failure'),
        (0, b'not json', b''),
        (0, json.dumps([]).encode(), b''),
        (0, json.dumps({'camera': [],
                        'controller': {'success': True,'devices': {}}}
                       ).encode(), b''),
    ],
)
def test_process_finished_failure_cases(mocker, exit_code, stdout, stderr):
    """A failing subprocess should emit an empty device list."""
    worker = DeviceDetectionWorker()

    errors = []
    devices = []

    worker.error_occurred.connect(errors.append)
    worker.devices_detected.connect(devices.append)

    proc = mocker.Mock()
    proc.readAllStandardOutput.return_value.data.return_value = stdout
    proc.readAllStandardError.return_value.data.return_value = stderr
    worker._process = proc

    mocker.patch.object(worker, 'stop')

    worker._on_process_finished(exit_code, qtc.QProcess.NormalExit)

    assert devices == [{'camera': {}, 'controller': {}}]
    assert errors
    worker.stop.assert_called_once()


def test_process_finished_import_failure(mocker):
    """Import errors during object reconstruction should be reported."""
    worker = DeviceDetectionWorker()

    errors = []
    devices = []

    worker.error_occurred.connect(errors.append)
    worker.devices_detected.connect(devices.append)

    proc = mocker.Mock()
    proc.readAllStandardOutput.return_value.data.return_value = json.dumps({
        'controller': {
            'success': True,
            'devices': {
                'ctrl': [
                    'does.not.exist',
                    'Missing',
                    {
                        'unique_name': 'TEST',
                        'has_hardware_interface': True,
                        'more': {}
                    }
                ]
            }
        }
    }).encode()

    worker._process = proc

    mocker.patch.object(worker, 'stop')

    worker._on_process_finished(0, qtc.QProcess.NormalExit)

    assert devices[0]['controller'] == {}
    assert errors
    worker.stop.assert_called_once()


def test_detection_timeout(mocker):
    """Timeout should stop detection and emit an empty device list."""
    worker = DeviceDetectionWorker()

    errors = []
    devices = []

    worker.error_occurred.connect(errors.append)
    worker.devices_detected.connect(devices.append)

    proc = mocker.Mock()
    proc.state.return_value = qtc.QProcess.Running

    worker._process = proc

    mocker.patch.object(worker, 'stop')

    worker._on_detection_timeout()

    assert devices == [{'camera': {}, 'controller': {}}]
    assert errors
    worker.stop.assert_called_once()


def test_process_error(mocker):
    """Process launch errors should emit an empty device list."""
    worker = DeviceDetectionWorker()

    errors = []
    devices = []

    worker.error_occurred.connect(errors.append)
    worker.devices_detected.connect(devices.append)

    proc = mocker.Mock()
    proc.errorString.return_value = 'boom'

    worker._process = proc

    mocker.patch.object(worker, 'stop')

    worker._on_process_error(qtc.QProcess.FailedToStart)

    assert devices == [{'camera': {}, 'controller': {}}]
    assert errors
    worker.stop.assert_called_once()


def test_stop_running_process(mocker):
    """Running processes should be terminated and cleaned up."""
    worker = DeviceDetectionWorker()

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
    worker = DeviceDetectionWorker()
    assert worker._process is None
    worker.stop()  # should not raise
    assert worker._process is None


def test_detection_timeout_no_process():
    """_on_detection_timeout should be a no-op when _process is None."""
    worker = DeviceDetectionWorker()
    errors = []
    devices = []
    worker.error_occurred.connect(errors.append)
    worker.devices_detected.connect(devices.append)
    worker._on_detection_timeout()
    assert not errors
    assert not devices


def test_process_error_no_process():
    """_on_process_error should be a no-op when _process is None."""
    worker = DeviceDetectionWorker()
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
                'error_type': DEF_CORRUPTED,
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
            'camera': {
                'success': True,
                'devices': {
                    'test_cam': [
                         'types',
                         'SimpleNamespace',
                        {'unique_name': 'TEST1',
                         'has_hardware_interface': True,
                         'more': {}}
                        ]
                }
            },
            'controller': {
                'success': True,
                'devices': {
                    'test_ctrl': [
                         'types',
                         'SimpleNamespace',
                        {'unique_name': 'TEST2',
                         'has_hardware_interface': True,
                         'more': {}}
                        ]
                }
            }
        },
        'test_ctrl',
        'test_cam',
        []
    ),
    (
        {
            'camera': {'success': False, 'error_type': RUN_ERR,
                       'error_msg': 'something went wrong'},
            'controller': {'success': True, 'devices': {}}
        },
        {},
        {},
        DeviceDetectionErrors.RUNTIME_ERROR.value[0]
    )
)


@parametrize('devices, ctrl, camera, error', test_cases)
# pylint: disable-next=too-complex
def test_device_detection_worker(mocker, devices, ctrl, camera, error):
    """Check device detection."""

    class FakeProcess(qtc.QObject):
        """Fake process used by the worker to simulate detection subprocess."""
        # pylint: disable=invalid-name
        NormalExit = qtc.QProcess.NormalExit
        NotRunning = qtc.QProcess.NotRunning
        ExitStatus = qtc.QProcess.ExitStatus
        ProcessError = qtc.QProcess.ProcessError

        finished = qtc.pyqtSignal(int, qtc.QProcess.ExitStatus)
        errorOccurred = qtc.pyqtSignal(qtc.QProcess.ProcessError)
        # pylint: enable=invalid-name

        def __init__(self, _):
            """Initialize fake QProcess."""
            super().__init__()
            self._finished_callbacks = []
            self._error_callbacks = []
            self._stdout = None

        def start(self, *_):
            """Simulate immediate process completion with JSON on stdout."""
            data = devices
            self._stdout = json.dumps(data).encode()
            self.finished.emit(0, self.NormalExit)

        def state(self):
            """Return NotRunning as state."""
            return self.NotRunning

        # pylint: disable=invalid-name
        def readAllStandardOutput(self):
            """Return fake output."""
            # pylint: disable=missing-function-docstring
            class B:    # pylint: disable=too-few-public-methods
                """Dummy class to replicate call chain."""
                def __init__(self, b):
                    self._b = b
                def data(self):
                    return self._b
            # pylint: enable=missing-function-docstring
            return B(self._stdout)

        def readAllStandardError(self):
            """Return fake error."""
            class B:    # pylint: disable=too-few-public-methods
                """Dummy class to replicate call chain."""
                # pylint: disable-next=missing-function-docstring
                def data(self):
                    return b''
            return B()
        # pylint: enable=invalid-name
    mocker.patch(_QPROCESS, FakeProcess)

    worker = DeviceDetectionWorker()
    emitted_errors = []
    emitted_devices = []
    worker.error_occurred.connect(emitted_errors.append)
    worker.devices_detected.connect(emitted_devices.append)
    worker.detect_devices()

    if emitted_devices[0]['controller']:
        assert ctrl in emitted_devices[0]['controller']
    else:
        assert emitted_devices[0]['controller'] == ctrl
    if emitted_devices[0]['camera']:
        assert camera in emitted_devices[0]['camera']
    else:
        assert emitted_devices[0]['camera'] == camera
    if emitted_errors:
        assert emitted_errors[0][0] == error
    else:
        assert emitted_errors == error
