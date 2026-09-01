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

from PyQt5 import QtCore as qtc
from PyQt5 import QtWidgets as qtw
from pytest import fixture

from viperleed.gui.measure.energysetter import EnergySetter
from viperleed.gui.measure.classes.settings import SettingsError

from .mock_qt import _FakeSignal


_ = qtw.QApplication(sys.argv)


# pylint: disable=protected-access
# pylint: disable=redefined-outer-name
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

    def list_devices(self):
        """Return a list of available devices."""
        return []


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


@fixture
def setter(tmp_path):
    """Create setter with path."""
    setter = EnergySetter()
    test_path = tmp_path / 'controller.ini'
    test_path.write_text('[controller]\n')
    setter.path = test_path
    return setter


@fixture
def ctrl_setter(mocker, setter):
    """Create setter with mocked controller."""
    fake_ctrl = _FakeController()
    setter._get_controller = mocker.Mock(return_value=fake_ctrl)
    return setter


@fixture
def fake_controller(mocker, setter):
    """Return a controller whose settings are mocked."""
    fake_settings = mocker.Mock()
    fake_settings.last_file = setter.path
    fake_settings.read_again.return_value = True

    fake_ctrl = _FakeController(connected=True)
    fake_ctrl._settings = fake_settings
    # pylint: disable-next=attribute-defined-outside-init
    fake_ctrl.set_settings = mocker.Mock(return_value=True)
    return fake_ctrl


def test_cleanup_controller(mocker):
    """Check controller cleanup disconnects signals."""
    setter = EnergySetter()
    fake_ctrl = _FakeController()
    setter._controller = fake_ctrl
    mocker.patch('viperleed.gui.measure.energysetter.base.safe_disconnect')

    setter.cleanup_controller()

    assert setter._controller is None
    assert not fake_ctrl.connected


def test_cleanup_controller_none():
    """Check cleanup handles None controller gracefully."""
    setter = EnergySetter()
    setter._controller = None

    # Should not raise
    setter.cleanup_controller()


def test_flush_resets_state(mocker):
    """Check flush resets all state."""
    setter = EnergySetter()
    setter.set_energy.setChecked(True)
    setter._operation_in_progress = True
    setter._pending_energy = 50.0
    setter._timeout_timer = mocker.Mock()
    setter.cleanup_controller = mocker.Mock()

    setter._flush()

    assert not setter.set_energy.isChecked()
    assert not setter._operation_in_progress
    assert setter._pending_energy is None
    setter._timeout_timer.stop.assert_called_once()
    setter.cleanup_controller.assert_called_once()


def test_get_controller_cleanup_on_different_path(mocker, setter):
    """Check old controller is cleaned up when path changes."""
    old_ctrl = _FakeController()
    old_ctrl.settings.last_file = Path('/old/path.ini')
    setter._controller = old_ctrl
    setter.cleanup_controller = mocker.Mock()
    setter._make_controller = mocker.Mock(return_value=_FakeController())

    setter._get_controller()

    setter.cleanup_controller.assert_called_once()


def test_get_controller_connection_failed(mocker, setter):
    """Check error emitted when controller connection fails."""
    setter.error_occurred = _FakeSignal()
    fake_ctrl = _FakeController(connected=False)
    setter._make_controller = mocker.Mock(return_value=fake_ctrl)

    result = setter._get_controller()

    assert result is None
    assert setter.error_occurred.emitted == 1


def test_get_controller_creates_new(mocker, setter):
    """Check that new controller is created when needed."""
    setter._make_controller = mocker.Mock(
        return_value=_FakeController()
        )

    result = setter._get_controller()

    setter._make_controller.assert_called_once()
    assert result is not None


def test_get_controller_load_failed(mocker, setter):
    """Check error emitted when controller creation fails."""
    setter.error_occurred = _FakeSignal()
    setter._make_controller = mocker.Mock(side_effect=ValueError('test'))

    result = setter._get_controller()

    assert result is None
    assert setter.error_occurred.emitted == 1


def test_get_controller_read_again_failed(mocker, fake_controller, setter):
    """Check controller is cleaned up when re-reading settings fails."""
    fake_settings = fake_controller.settings
    fake_settings.read_again.return_value = False
    setter._controller = fake_controller
    setter.error_occurred = _FakeSignal()
    setter.cleanup_controller = mocker.Mock()

    result = setter._get_controller()

    assert result is None
    fake_settings.read_again.assert_called_once()
    assert setter.error_occurred.emitted == 1
    setter.cleanup_controller.assert_called_once()


def test_get_controller_read_again_raises(mocker, fake_controller, setter):
    """Check controller cleaned up when read_again raises an error."""
    fake_settings = fake_controller.settings
    fake_settings.read_again.side_effect = SettingsError('corrupted')
    setter._controller = fake_controller
    setter.error_occurred = _FakeSignal()
    setter.cleanup_controller = mocker.Mock()

    result = setter._get_controller()

    assert result is None
    fake_settings.read_again.assert_called_once()
    assert setter.error_occurred.emitted == 1
    setter.cleanup_controller.assert_called_once()


def test_get_controller_reuse_connection_failed(mocker, fake_controller,
                                                setter):
    """Check controller cleaned up when reconnection fails."""
    fake_settings = fake_controller.settings
    fake_controller._connected = False
    setter._controller = fake_controller
    setter.error_occurred = _FakeSignal()
    setter.cleanup_controller = mocker.Mock()

    result = setter._get_controller()

    assert result is None
    fake_settings.read_again.assert_called_once()
    fake_controller.set_settings.assert_called_once_with(fake_settings)
    assert setter.error_occurred.emitted == 1
    setter.cleanup_controller.assert_called_once()


def test_get_controller_reuses_existing(fake_controller, setter):
    """Check that existing controller is reused when path matches."""
    fake_settings = fake_controller.settings
    setter._controller = fake_controller

    result = setter._get_controller()

    assert result is fake_controller
    fake_settings.read_again.assert_called_once()
    fake_controller.set_settings.assert_called_once_with(fake_settings)
    assert fake_controller.connected


def test_get_controller_set_settings_failed(mocker, fake_controller, setter):
    """Check controller is cleaned up when set_settings fails."""
    fake_settings = fake_controller.settings
    fake_controller.set_settings.return_value = False
    setter._controller = fake_controller
    setter.error_occurred = _FakeSignal()
    setter.cleanup_controller = mocker.Mock()

    result = setter._get_controller()

    assert result is None
    fake_settings.read_again.assert_called_once()
    fake_controller.set_settings.assert_called_once_with(fake_settings)
    assert setter.error_occurred.emitted == 1
    setter.cleanup_controller.assert_called_once()


def test_init_creates_widgets():
    """Check that initialization creates required widgets."""
    setter = EnergySetter()

    assert setter.energy_input is not None
    assert setter.set_energy is not None
    assert setter.path is None
    assert setter._controller is None
    assert setter._pending_energy is None
    assert not setter._operation_in_progress


def test_is_busy_independent_of_checkbox(ctrl_setter):
    """Check is_busy is not tied to the checkbox state."""
    ctrl_setter._operation_in_progress = True
    ctrl_setter.set_energy.setChecked(False)
    assert not ctrl_setter.setting_energy
    assert ctrl_setter.is_busy


def test_is_busy_property(ctrl_setter):
    """Check is_busy reflects an in-flight operation."""
    assert not ctrl_setter.is_busy
    ctrl_setter._operation_in_progress = True
    assert ctrl_setter.is_busy
    ctrl_setter._operation_in_progress = False
    assert not ctrl_setter.is_busy


def test_make_controller(mocker, tmp_path):
    """Check controller creation from settings file."""
    setter = EnergySetter()
    test_path = tmp_path / 'controller.ini'
    test_path.write_text(
        '[controller]\ncontroller_class = TestCtrl\naddress = COM1\n'
        )
    setter.path = test_path

    fake_settings = mocker.Mock()
    fake_ctrl = _FakeController()
    fake_cls = mocker.Mock(return_value=fake_ctrl)
    mocker.patch('viperleed.gui.measure.energysetter.ViPErLEEDSettings',
                 return_value=fake_settings)
    mocker.patch('viperleed.gui.measure.energysetter.base.class_from_name',
                 return_value=fake_cls)

    result = setter._make_controller()

    fake_settings.read.assert_called_once_with(test_path)
    fake_cls.assert_called()
    assert result is not None


def test_on_ctrl_finished_busy(mocker):
    """Check ctrl finished ignores busy state."""
    setter = EnergySetter()
    setter._operation_in_progress = True
    setter._timeout_timer = mocker.Mock()

    setter._on_ctrl_finished(True)

    assert setter._operation_in_progress
    setter._timeout_timer.stop.assert_not_called()


def test_on_ctrl_finished_no_pending(mocker):
    """Check ctrl finished when no pending energy."""
    setter = EnergySetter()
    setter._operation_in_progress = True
    setter._pending_energy = None
    setter.set_energy = mocker.Mock()
    setter.set_energy.isChecked.return_value = True
    setter._timeout_timer = mocker.Mock()
    setter._set_energy = mocker.Mock()

    setter._on_ctrl_finished(False)

    assert not setter._operation_in_progress
    setter._timeout_timer.stop.assert_called_once()
    setter._set_energy.assert_not_called()


def test_on_ctrl_finished_not_setting_cleanup(mocker):
    """Check ctrl finished cleans up when not setting."""
    setter = EnergySetter()
    setter._operation_in_progress = True
    setter.set_energy = mocker.Mock()
    setter.set_energy.isChecked.return_value = False
    setter.cleanup_controller = mocker.Mock()
    setter._timeout_timer = mocker.Mock()
    setter._controller = mocker.Mock()

    setter._on_ctrl_finished(False)

    assert not setter._operation_in_progress
    setter._timeout_timer.stop.assert_called_once()
    setter._controller.disconnect_.assert_called_once()


def test_on_ctrl_finished_pending_energy(mocker):
    """Check ctrl finished processes pending energy."""
    setter = EnergySetter()
    setter._operation_in_progress = True
    setter._pending_energy = 75.0
    setter.set_energy = mocker.Mock()
    setter.set_energy.isChecked.return_value = True
    setter._set_energy = mocker.Mock()

    setter._on_ctrl_finished(False)

    assert setter._pending_energy is None
    setter._set_energy.assert_called_once_with(75.0)


def test_on_ctrl_finished_queued_zero_waits_for_completion(ctrl_setter):
    """Check that a deferred zeroing is sent and kept busy afterwards."""
    ctrl_setter._operation_in_progress = True
    ctrl_setter._pending_energy = 0.0
    ctrl_setter._controller = _FakeController()

    ctrl_setter._on_ctrl_finished(False)

    assert ctrl_setter._pending_energy is None
    assert ctrl_setter._operation_in_progress
    assert ctrl_setter._controller.connected


def test_on_ctrl_finished_zero_applied_then_disconnects(ctrl_setter):
    """Check the setter disconnects after the deferred zero was applied."""
    ctrl = _FakeController()
    ctrl_setter._controller = ctrl

    ctrl_setter._on_ctrl_finished(False)

    assert not ctrl_setter._operation_in_progress
    assert not ctrl.connected


def test_on_energy_changed_operation_in_progress(ctrl_setter):
    """Check energy change queued when operation in progress."""
    ctrl_setter.set_energy.setChecked(True)
    ctrl_setter._operation_in_progress = True
    ctrl_setter.energy_input.setValue(75.0)

    ctrl_setter._on_energy_changed()

    # pylint: disable-next=magic-value-comparison
    assert ctrl_setter._pending_energy == 75.0


def test_on_energy_changed_sets_energy(mocker, ctrl_setter):
    """Check energy change triggers set_energy when idle."""
    ctrl_setter._set_energy = mocker.Mock()
    ctrl_setter.set_energy.setChecked(True)
    ctrl_setter._operation_in_progress = False
    ctrl_setter.energy_input.setValue(50.0)
    ctrl_setter.energy_input.editingFinished.emit()

    ctrl_setter._set_energy.assert_called_once_with(50.0)


def test_on_energy_changed_when_not_setting(mocker, ctrl_setter):
    """Check energy change ignored when not setting energy."""
    ctrl_setter.set_energy.setChecked(False)
    ctrl_setter._set_energy = mocker.Mock()

    ctrl_setter._on_energy_changed()

    ctrl_setter._set_energy.assert_not_called()


def test_on_error(mocker):
    """Check error handler emits error and flushes."""
    setter = EnergySetter()
    setter.error_occurred = _FakeSignal()
    setter._flush = mocker.Mock()
    error_info = (2000, 'test error')

    setter._on_error(error_info)

    assert setter.error_occurred.emitted == 1
    setter._flush.assert_called_once()


def test_on_set_energy_toggled_missing_file(tmp_path):
    """Check error emitted when path file is missing."""
    setter = EnergySetter()
    setter.path = tmp_path / 'nonexistent.ini'
    setter.error_occurred = _FakeSignal()

    setter._on_set_energy_toggled(qtc.Qt.Checked)

    assert not setter.set_energy.isChecked()
    assert setter.error_occurred.emitted == 1


def test_on_set_energy_toggled_no_path():
    """Check error emitted when toggling without path."""
    setter = EnergySetter()
    setter.path = None
    setter.error_occurred = _FakeSignal()

    setter._on_set_energy_toggled(qtc.Qt.Checked)

    assert not setter.set_energy.isChecked()
    assert setter.error_occurred.emitted == 1


def test_on_set_energy_toggled_unchecked_in_flight_queues_zero(ctrl_setter):
    """Check that unchecking during an operation defers the zeroing."""
    ctrl_setter._operation_in_progress = True
    ctrl_setter.set_energy.setChecked(True)

    ctrl_setter.set_energy.setChecked(False)

    # pylint: disable-next=use-implicit-booleaness-not-comparison-to-zero
    assert ctrl_setter._pending_energy == 0.0
    # pylint: disable-next=use-implicit-booleaness-not-comparison-to-zero
    assert ctrl_setter.energy_input.value() == 0.0


def test_on_set_energy_toggled_unchecked_sets_zero(mocker, ctrl_setter):
    """Check that unchecking sets energy to zero."""
    ctrl_setter._set_energy = mocker.Mock()
    ctrl_setter.set_energy.setChecked(True)

    ctrl_setter.energy_input.setValue(50.0)
    ctrl_setter.energy_input.editingFinished.emit()
    ctrl_setter._set_energy.assert_called_with(50.0)

    ctrl_setter.set_energy.setChecked(False)
    ctrl_setter._set_energy.assert_called_with(0.0)

    # pylint: disable-next=use-implicit-booleaness-not-comparison-to-zero
    assert ctrl_setter.energy_input.value() == 0.0


def test_on_timeout(mocker):
    """Check timeout handler emits timeout error."""
    setter = EnergySetter()
    setter.error_occurred = _FakeSignal()
    setter._flush = mocker.Mock()

    setter._on_timeout()

    assert setter.error_occurred.emitted == 1
    setter._flush.assert_called_once()


def test_set_enabled_with_checked_checkbox(ctrl_setter):
    """Check energy input enabled when checkbox is checked."""
    ctrl_setter.set_energy.setChecked(True)

    ctrl_setter.set_enabled(True)

    assert ctrl_setter.set_energy.isEnabled()
    assert ctrl_setter.energy_input.isEnabled()


def test_set_enabled_with_path(setter):
    """Check set_enabled enables when path exists."""
    setter.set_enabled(True)

    assert setter.set_energy.isEnabled()
    assert not setter.energy_input.isEnabled()


def test_set_enabled_without_path():
    """Check set_enabled disables when no path."""
    setter = EnergySetter()
    setter.path = None

    setter.set_enabled(True)

    assert not setter.set_energy.isEnabled()
    assert not setter.energy_input.isEnabled()


def test_set_energy_no_controller(mocker):
    """Check set_energy handles missing controller."""
    setter = EnergySetter()
    setter._controller = None
    setter.set_energy = mocker.Mock()

    setter._set_energy(50.0)

    assert not setter._operation_in_progress
    setter.set_energy.setChecked.assert_called_once_with(False)


def test_set_energy_starts_timeout(mocker):
    """Check set_energy starts timeout timer."""
    setter = EnergySetter()
    setter._controller = _FakeController()
    setter._timeout_timer = mocker.Mock()

    setter._set_energy(50.0)

    assert setter._operation_in_progress
    setter._timeout_timer.start.assert_called_once()


def test_setting_energy_property(ctrl_setter):
    """Check setting_energy property reflects checkbox state."""
    ctrl_setter.set_energy.setCheckState(qtc.Qt.Unchecked)
    assert not ctrl_setter.setting_energy
    ctrl_setter.set_energy.setCheckState(qtc.Qt.Checked)
    assert ctrl_setter.setting_energy
