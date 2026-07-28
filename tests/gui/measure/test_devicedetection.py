"""Tests for module devicedetection of viperleed.gui.measure."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    'Florian Dörr (@FlorianDoerr)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-07-16'
__license__ = 'GPLv3+'

import json

from PyQt5 import QtCore as qtc
from pytest import fixture
from pytest_cases import parametrize

from viperleed.gui.measure.classes.abc import QObjectSettingsErrors
from viperleed.gui.measure.devicedetection import DEF_CORRUPTED
from viperleed.gui.measure.devicedetection import RUN_ERR
from viperleed.gui.measure.devicedetection import DeviceDetectionErrors
from viperleed.gui.measure.devicedetection import DeviceDetectionWorker


_QPROCESS = 'viperleed.gui.measure.devicedetection.qtc.QProcess'


@fixture
def _worker():
    """Create a DeviceDetectionWorker with error and device collectors."""
    worker = DeviceDetectionWorker()
    errors = []
    devices = []
    worker.error_occurred.connect(errors.append)
    worker.devices_detected.connect(devices.append)
    return worker, errors, devices


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

    def emit(self, *args, **kwargs):
        """Count emissions and call connected slots."""
        self.emitted += 1
        for slot in self.connected:
            slot(*args, **kwargs)

    def disconnect(self, slot):
        """Store disconnected slots."""
        self.disconnected.append(slot)


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
def test_process_finished_failure_cases(mocker, _worker, exit_code, stdout,
                                        stderr):
    """A failing subprocess should emit an empty device list."""
    worker, errors, devices = _worker

    proc = mocker.Mock()
    proc.readAllStandardOutput.return_value.data.return_value = stdout
    proc.readAllStandardError.return_value.data.return_value = stderr
    worker._process = proc

    mocker.patch.object(worker, 'stop')

    worker._on_process_finished(exit_code, qtc.QProcess.NormalExit)

    assert devices == [{'camera': {}, 'controller': {}}]
    assert errors
    worker.stop.assert_called_once()


def test_process_finished_import_failure(mocker, _worker):
    """Import errors during object reconstruction should be reported."""
    worker, errors, devices = _worker

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


def test_detection_timeout(mocker, _worker):
    """Timeout should stop detection and emit an empty device list."""
    worker, errors, devices = _worker

    proc = mocker.Mock()
    proc.state.return_value = qtc.QProcess.Running

    worker._process = proc

    mocker.patch.object(worker, 'stop')

    worker._on_detection_timeout()

    assert devices == [{'camera': {}, 'controller': {}}]
    assert errors
    worker.stop.assert_called_once()


def test_process_error(mocker, _worker):
    """Process launch errors should emit an empty device list."""
    worker, errors, devices = _worker

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


def test_detection_timeout_no_process(_worker):
    """_on_detection_timeout should be a no-op when _process is None."""
    worker, errors, devices = _worker
    worker._on_detection_timeout()
    assert not errors
    assert not devices


def test_process_error_no_process(_worker):
    """_on_process_error should be a no-op when _process is None."""
    worker, errors, devices = _worker
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


@parametrize('devices, ctrl, cam, error', test_cases)
def test_device_detection_worker(mocker, _worker, devices, ctrl, cam, error):
    """Check device detection."""

    finished_signal = _FakeSignal()
    error_signal = _FakeSignal()
    proc = mocker.Mock()
    proc.readAllStandardOutput.return_value.data.return_value = (
        json.dumps(devices).encode()
    )
    proc.readAllStandardError.return_value.data.return_value = b''
    proc.state.return_value = qtc.QProcess.NotRunning
    proc.finished = finished_signal
    proc.errorOccurred = error_signal

    def start_and_finish(*_):
        """Simulate process start and immediate completion."""
        finished_signal.emit(0, qtc.QProcess.NormalExit)

    proc.start.side_effect = start_and_finish

    mocker.patch(_QPROCESS, return_value=proc)
    mocker.patch('viperleed.gui.measure.devicedetection.qtc.QTimer')

    worker, emitted_errors, emitted_devices = _worker
    worker.detect_devices()

    # Check controllers: if a controller was found, verify the expected
    # controller key exists in the result. Otherwise verify the result
    # is empty.
    if ctrl:
        assert ctrl in emitted_devices[0]['controller']
    else:
        assert not emitted_devices[0]['controller']

    # Check cameras: same logic as controller.
    if cam:
        assert cam in emitted_devices[0]['camera']
    else:
        assert not emitted_devices[0]['camera']

    # Check errors: if errors were emitted, verify the first error matches
    # the expected error type. Else check that no errors were emitted.
    if error:
        assert emitted_errors[0][0] == error
    else:
        assert not emitted_errors
