"""Tests for module energysetter of viperleed.gui.measure."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    'Florian Dörr (@FlorianDoerr)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-08-06'
__license__ = 'GPLv3+'

from pathlib import Path
import sys
import pytest

from PyQt5 import QtCore as qtc
from PyQt5 import QtWidgets as qtw

from viperleed.gui.measure.energysetter import EnergySetter

from .mock_qt import _FakeSignal


_ = qtw.QApplication(sys.argv)


# pylint: disable=protected-access
class _FakeController:
    """A minimal controller-like object for testing."""

    def __init__(self, connected=True):
        """Initialize fake controller."""
        self.connected = True
        self._connected = connected
        self.error_occurred = _FakeSignal()
        self._settings = _FakeControllerSettings()
        self._serial = _FakeSerial()

    @property
    def settings(self):
        """Return controller settings."""
        return self._settings

    @property
    def serial(self):
        """Return serial interface."""
        return self._serial

    def connect_(self):
        """Connect the controller."""
        self.connected = self._connected

    def disconnect_(self):
        """Disconnect the controller."""
        self.connected = False

    def set_energy(self, energy, ramp_time, trigger_meas):
        """Set energy on the controller."""

    def deleteLater(self):  # pylint: disable=invalid-name
        """Schedule for deletion."""


class _FakeControllerSettings:  # pylint: disable=too-few-public-methods
    """Fake controller settings object."""

    def __init__(self):
        """Initialize fake settings."""
        self.last_file = None


class _FakeSerial:
    """Fake serial interface."""

    def __init__(self):
        """Initialize fake serial."""
        self.busy_changed = _FakeSignal()
        self._busy = False

    @property
    def busy(self):
        """Return busy state."""
        return self._busy

    def set_busy(self, busy):
        """Set busy state and emit signal."""
        self._busy = busy
        self.busy_changed.emit(busy)


@pytest.fixture
def setter(tmp_path):
    """Create setter with primary path."""
    setter = EnergySetter()
    test_path = tmp_path / 'controller.ini'
    test_path.write_text('[controller]\n')
    setter.primary_path = test_path
    return setter


@pytest.fixture
def ctrl_setter(mocker, setter):
    """Create setter which create a fake primary.."""
    fake_ctrl = _FakeController()
    setter._get_primary_controller = mocker.Mock(return_value=fake_ctrl)
    return setter


def test_init_creates_widgets():
    """Check that initialization creates required widgets."""
    setter = EnergySetter()

    assert setter.energy_input is not None
    assert setter.set_energy is not None
    assert setter.primary_path is None
    assert setter._primary_controller is None
    assert setter._pending_energy is None
    assert not setter._operation_in_progress


def test_setting_energy_property(ctrl_setter):
    """Check setting_energy property reflects checkbox state."""
    ctrl_setter.set_energy.setCheckState(qtc.Qt.Unchecked)
    assert not ctrl_setter.setting_energy
    ctrl_setter.set_energy.setCheckState(qtc.Qt.Checked)
    assert ctrl_setter.setting_energy

def test_set_enabled_without_primary_path():
    """Check set_enabled disables when no primary path."""
    setter = EnergySetter()
    setter.primary_path = None

    setter.set_enabled(True)

    assert not setter.set_energy.isEnabled()
    assert not setter.energy_input.isEnabled()


def test_set_enabled_with_primary_path(setter):
    """Check set_enabled enables when primary path exists."""
    setter.set_enabled(True)

    assert setter.set_energy.isEnabled()
    assert not setter.energy_input.isEnabled()


def test_set_enabled_with_checked_checkbox(ctrl_setter):
    """Check energy input enabled when checkbox is checked."""
    ctrl_setter.set_energy.setChecked(True)

    ctrl_setter.set_enabled(True)

    assert ctrl_setter.set_energy.isEnabled()
    assert ctrl_setter.energy_input.isEnabled()


def test_on_set_energy_toggled_no_primary_path():
    """Check error emitted when toggling without primary path."""
    setter = EnergySetter()
    setter.primary_path = None
    setter.error_occurred = _FakeSignal()

    setter._on_set_energy_toggled(qtc.Qt.Checked)

    assert not setter.set_energy.isChecked()
    assert setter.error_occurred.emitted == 1


def test_on_set_energy_toggled_missing_file(tmp_path):
    """Check error emitted when primary path file is missing."""
    setter = EnergySetter()
    setter.primary_path = tmp_path / 'nonexistent.ini'
    setter.error_occurred = _FakeSignal()

    setter._on_set_energy_toggled(qtc.Qt.Checked)

    assert not setter.set_energy.isChecked()
    assert setter.error_occurred.emitted == 1


def test_on_set_energy_toggled_unchecked_sets_zero(mocker, ctrl_setter):
    """Check that unchecking sets energy to zero."""
    ctrl_setter._set_energy = mocker.Mock()
    ctrl_setter.set_energy.setChecked(True)

    ctrl_setter.energy_input.setValue(50.0)
    ctrl_setter._set_energy.assert_called_with(50.0)

    ctrl_setter.set_energy.setChecked(False)
    ctrl_setter._set_energy.assert_called_with(0.0)

    assert ctrl_setter.energy_input.value() == 0.0


def test_on_energy_changed_when_not_setting(mocker, ctrl_setter):
    """Check energy change ignored when not setting energy."""
    ctrl_setter.set_energy.setChecked(False)
    ctrl_setter._set_energy = mocker.Mock()

    ctrl_setter._on_energy_changed()

    ctrl_setter._set_energy.assert_not_called()


def test_on_energy_changed_operation_in_progress(ctrl_setter):
    """Check energy change queued when operation in progress."""
    ctrl_setter.set_energy.setChecked(True)
    ctrl_setter._operation_in_progress = True
    ctrl_setter.energy_input.setValue(75.0)

    ctrl_setter._on_energy_changed()

    assert ctrl_setter._pending_energy == 75.0


def test_on_energy_changed_sets_energy(mocker, ctrl_setter):
    """Check energy change triggers set_energy when idle."""
    ctrl_setter._set_energy = mocker.Mock()
    ctrl_setter.set_energy.setChecked(True)
    ctrl_setter._operation_in_progress = False
    ctrl_setter.energy_input.setValue(50.0)

    ctrl_setter._set_energy.assert_called_once_with(50.0)


def test_get_primary_controller_reuses_existing():
    """Check that existing controller is reused when path matches."""
    setter = EnergySetter()
    fake_ctrl = _FakeController()
    fake_ctrl.settings.last_file = Path('/test/path.ini')
    setter._primary_controller = fake_ctrl
    setter.primary_path = Path('/test/path.ini')

    result = setter._get_primary_controller()

    assert result is fake_ctrl
    assert fake_ctrl.connected


def test_get_primary_controller_creates_new(mocker):
    """Check that new controller is created when needed."""
    setter = EnergySetter()
    setter.primary_path = Path('/test/path.ini')
    setter._make_primary_controller = mocker.Mock(
        return_value=_FakeController()
        )

    result = setter._get_primary_controller()

    setter._make_primary_controller.assert_called_once()
    assert result is not None


def test_get_primary_controller_cleanup_on_different_path(mocker):
    """Check old controller is cleaned up when path changes."""
    setter = EnergySetter()
    old_ctrl = _FakeController()
    old_ctrl.settings.last_file = Path('/old/path.ini')
    setter._primary_controller = old_ctrl
    setter.primary_path = Path('/new/path.ini')
    setter.cleanup_primary_controller = mocker.Mock()
    setter._make_primary_controller = mocker.Mock(
        return_value=_FakeController()
        )

    setter._get_primary_controller()

    setter.cleanup_primary_controller.assert_called_once()


def test_get_primary_controller_load_failed(mocker):
    """Check error emitted when controller creation fails."""
    setter = EnergySetter()
    setter.primary_path = Path('/test/path.ini')
    setter.error_occurred = _FakeSignal()
    setter._make_primary_controller = mocker.Mock(
        side_effect=ValueError('test')
        )

    result = setter._get_primary_controller()

    assert result is None
    assert setter.error_occurred.emitted == 1


def test_get_primary_controller_connection_failed(mocker):
    """Check error emitted when controller connection fails."""
    setter = EnergySetter()
    setter.primary_path = Path('/test/path.ini')
    setter.error_occurred = _FakeSignal()
    fake_ctrl = _FakeController(connected=False)
    setter._make_primary_controller = mocker.Mock(return_value=fake_ctrl)

    result = setter._get_primary_controller()

    assert result is None
    assert setter.error_occurred.emitted == 1


def test_make_primary_controller(mocker, tmp_path):
    """Check controller creation from settings file."""
    setter = EnergySetter()
    test_path = tmp_path / 'controller.ini'
    test_path.write_text(
        '[controller]\ncontroller_class = TestCtrl\naddress = COM1\n'
        )
    setter.primary_path = test_path

    fake_settings = mocker.Mock()
    fake_cls = mocker.Mock(return_value=_FakeController())
    mocker.patch('viperleed.gui.measure.energysetter.ViPErLEEDSettings',
                 return_value=fake_settings)
    mocker.patch('viperleed.gui.measure.energysetter.base.class_from_name',
                 return_value=fake_cls)

    result = setter._make_primary_controller()

    fake_settings.read.assert_called_once_with(test_path)
    fake_cls.assert_called()
    assert result is not None


def test_cleanup_primary_controller(mocker):
    """Check controller cleanup disconnects signals."""
    setter = EnergySetter()
    fake_ctrl = _FakeController()
    setter._primary_controller = fake_ctrl
    mocker.patch('viperleed.gui.measure.energysetter.base.safe_disconnect')

    setter.cleanup_primary_controller()

    assert setter._primary_controller is None
    assert not fake_ctrl.connected


def test_cleanup_primary_controller_none():
    """Check cleanup handles None controller gracefully."""
    setter = EnergySetter()
    setter._primary_controller = None

    # Should not raise
    setter.cleanup_primary_controller()


def test_on_error(mocker):
    """Check error handler emits error and flushes."""
    setter = EnergySetter()
    setter.error_occurred = _FakeSignal()
    setter._flush = mocker.Mock()
    error_info = (2000, 'test error')

    setter._on_error(error_info)

    assert setter.error_occurred.emitted == 1
    setter._flush.assert_called_once()


def test_on_timeout(mocker):
    """Check timeout handler emits timeout error."""
    setter = EnergySetter()
    setter.error_occurred = _FakeSignal()
    setter._flush = mocker.Mock()

    setter._on_timeout()

    assert setter.error_occurred.emitted == 1
    setter._flush.assert_called_once()


def test_flush_resets_state(mocker):
    """Check flush resets all state."""
    setter = EnergySetter()
    setter.set_energy.setChecked(True)
    setter._operation_in_progress = True
    setter._pending_energy = 50.0
    setter._timeout_timer = mocker.Mock()
    setter.cleanup_primary_controller = mocker.Mock()

    setter._flush()

    assert not setter.set_energy.isChecked()
    assert not setter._operation_in_progress
    assert setter._pending_energy is None
    setter._timeout_timer.stop.assert_called_once()
    setter.cleanup_primary_controller.assert_called_once()


def test_set_energy_starts_timeout(mocker):
    """Check set_energy starts timeout timer."""
    setter = EnergySetter()
    setter._primary_controller = _FakeController()
    setter._timeout_timer = mocker.Mock()

    setter._set_energy(50.0)

    assert setter._operation_in_progress
    setter._timeout_timer.start.assert_called_once()


def test_set_energy_no_controller(mocker):
    """Check set_energy handles missing controller."""
    setter = EnergySetter()
    setter._primary_controller = None
    setter.set_energy = mocker.Mock()

    setter._set_energy(50.0)

    assert not setter._operation_in_progress
    setter.set_energy.setChecked.assert_called_once_with(False)


def test_on_ctrl_finished_busy(mocker):
    """Check ctrl finished ignores busy state."""
    setter = EnergySetter()
    setter._operation_in_progress = True
    setter._timeout_timer = mocker.Mock()

    setter._on_ctrl_finished(True)

    assert setter._operation_in_progress
    setter._timeout_timer.stop.assert_not_called()


def test_on_ctrl_finished_not_setting_cleanup(mocker):
    """Check ctrl finished cleans up when not setting."""
    setter = EnergySetter()
    setter._operation_in_progress = True
    setter.set_energy = mocker.Mock(isChecked=lambda: False)
    setter.cleanup_primary_controller = mocker.Mock()
    setter._timeout_timer = mocker.Mock()

    setter._on_ctrl_finished(False)

    assert not setter._operation_in_progress
    setter._timeout_timer.stop.assert_called_once()
    setter.cleanup_primary_controller.assert_called_once()


def test_on_ctrl_finished_pending_energy(mocker):
    """Check ctrl finished processes pending energy."""
    setter = EnergySetter()
    setter._operation_in_progress = True
    setter._pending_energy = 75.0
    setter.set_energy = mocker.Mock(isChecked=lambda: True)
    setter._set_energy = mocker.Mock()

    setter._on_ctrl_finished(False)

    assert setter._pending_energy is None
    setter._set_energy.assert_called_once_with(75.0)


def test_on_ctrl_finished_no_pending(mocker):
    """Check ctrl finished when no pending energy."""
    setter = EnergySetter()
    setter._operation_in_progress = True
    setter._pending_energy = None
    setter.set_energy = mocker.Mock(isChecked=lambda: True)
    setter._timeout_timer = mocker.Mock()
    setter._set_energy = mocker.Mock()

    setter._on_ctrl_finished(False)

    assert not setter._operation_in_progress
    setter._timeout_timer.stop.assert_called_once()
    setter._set_energy.assert_not_called()
