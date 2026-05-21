"""Tests for module _meassettings of viperleed.gui.measure.measurement."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    'Florian Dörr (@FlorianDoerr)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2026-03-31'
__license__ = 'GPLv3+'

import pytest

_meassettings = pytest.importorskip(
    'viperleed.gui.measure.measurement._meassettings',
    reason='PyQt5 GUI components are an optional dependency',
    )


def test_get_step_profile_limits_defaults_for_invalid_settings():
    """Return default controller limits for invalid primary settings."""
    limits = _meassettings.get_step_profile_limits('not a sequence')
    expected = (_meassettings.ControllerABC.MAX_NUM_STEPS,
                _meassettings.ControllerABC.MAX_DELAY)
    assert limits == expected


def test_get_step_profile_limits_uses_controller_class_limits(monkeypatch):
    """Return limits from the primary controller class."""
    class FakeConfig:  # pylint: disable=too-few-public-methods
        """Simple config object with get method."""
        @staticmethod
        def get(*_args, **_kwargs):
            return 'FakeController'

    class FakeController:  # pylint: disable=too-few-public-methods
        """Simple controller class with class limits."""
        MAX_NUM_STEPS = 3
        MAX_DELAY = 1200

    monkeypatch.setattr(_meassettings.ViPErLEEDSettings, 'from_settings',
                        lambda *_: FakeConfig())
    monkeypatch.setattr(_meassettings.base, 'class_from_name',
                        lambda *_: FakeController)

    limits = _meassettings.get_step_profile_limits(('controller.ini', ()))
    assert limits == (3, 1200)


def test_get_step_profile_limits_uses_defaults_when_class_missing_limits(monkeypatch):
    """Return defaults if controller class does not declare limits."""
    class FakeConfig:  # pylint: disable=too-few-public-methods
        """Simple config object with get method."""
        @staticmethod
        def get(*_args, **_kwargs):
            return 'FakeController'

    class FakeController:  # pylint: disable=too-few-public-methods
        """Simple controller class without class limits."""

    monkeypatch.setattr(_meassettings.ViPErLEEDSettings, 'from_settings',
                        lambda *_: FakeConfig())
    monkeypatch.setattr(_meassettings.base, 'class_from_name',
                        lambda *_: FakeController)

    limits = _meassettings.get_step_profile_limits(('controller.ini', ()))
    expected = (_meassettings.ControllerABC.MAX_NUM_STEPS,
                _meassettings.ControllerABC.MAX_DELAY)
    assert limits == expected


def test_get_step_profile_limits_falls_back_for_class_lookup_errors(monkeypatch):
    """Return defaults if class lookup raises non-ValueError."""
    class FakeConfig:  # pylint: disable=too-few-public-methods
        """Simple config object with get method."""
        @staticmethod
        def get(*_args, **_kwargs):
            return 'FakeController'

    monkeypatch.setattr(_meassettings.ViPErLEEDSettings, 'from_settings',
                        lambda *_: FakeConfig())
    def _raise_runtime_error(*_args, **_kwargs):
        raise RuntimeError

    monkeypatch.setattr(_meassettings.base, 'class_from_name',
                        _raise_runtime_error)

    limits = _meassettings.get_step_profile_limits(('controller.ini', ()))
    expected = (_meassettings.ControllerABC.MAX_NUM_STEPS,
                _meassettings.ControllerABC.MAX_DELAY)
    assert limits == expected


def test_get_step_profile_limits_falls_back_for_invalid_limit_values(monkeypatch):
    """Return defaults for non-numeric or negative class limits."""
    class FakeConfig:  # pylint: disable=too-few-public-methods
        """Simple config object with get method."""
        @staticmethod
        def get(*_args, **_kwargs):
            return 'FakeController'

    class FakeController:  # pylint: disable=too-few-public-methods
        """Controller class with invalid limits."""
        MAX_NUM_STEPS = 'abc'
        MAX_DELAY = -1

    monkeypatch.setattr(_meassettings.ViPErLEEDSettings, 'from_settings',
                        lambda *_: FakeConfig())
    monkeypatch.setattr(_meassettings.base, 'class_from_name',
                        lambda *_: FakeController)

    limits = _meassettings.get_step_profile_limits(('controller.ini', ()))
    expected = (_meassettings.ControllerABC.MAX_NUM_STEPS,
                _meassettings.ControllerABC.MAX_DELAY)
    assert limits == expected
