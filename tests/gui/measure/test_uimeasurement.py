"""Tests for module uimeasurement of viperleed.gui.measure."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    'Florian Dörr (@FlorianDoerr)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-03-31'
__license__ = 'GPLv3+'

from viperleed.gui.measure.uimeasurement import Measure


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
