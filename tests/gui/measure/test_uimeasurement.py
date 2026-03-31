"""Tests for module uimeasurement of viperleed.gui.measure."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    'Florian Dörr (@FlorianDoerr)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-03-31'
__license__ = 'GPLv3+'

from types import SimpleNamespace

from viperleed.gui.measure.uimeasurement import Measure


class _FakeSignal:  # pylint: disable=too-few-public-methods
    """A minimal signal-like object."""

    def __init__(self):
        self.connected = []
        self.emitted = 0

    def connect(self, slot):
        """Store connected slots."""
        self.connected.append(slot)

    def emit(self):
        """Count emissions."""
        self.emitted += 1


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


def test_trigger_device_search_blocks_reentry():
    """Check that a second search is blocked while one is in progress."""
    signal = _FakeSignal()
    fake = SimpleNamespace(
        _device_search_in_progress=False,
        detect_devices_requested=signal,
        )
    fake._device_search_allowed = (
        lambda: not fake._device_search_in_progress
        )

    Measure._trigger_device_search(fake)
    assert fake._device_search_in_progress
    assert signal.emitted == 1

    Measure._trigger_device_search(fake)
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
