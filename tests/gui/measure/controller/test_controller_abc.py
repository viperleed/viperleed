"""Tests for controller data dispatch in controller abc module."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    'Florian Dörr (@FlorianDoerr)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2026-03-31'
__license__ = 'GPLv3+'

from viperleed.gui.measure.controller.abc import ControllerABC


def test_is_measurement_data_default_false():
    """Check that non-overridden data-type detection defaults to False."""
    controller = object.__new__(ControllerABC)
    assert not ControllerABC.is_measurement_data(controller, object())


def test_on_data_ready_dispatches_to_hardware_handler_by_default(mocker):
    """Check that data dispatch defaults to hardware-information handler."""
    controller = object.__new__(ControllerABC)
    process_hw = mocker.Mock()
    process_meas = mocker.Mock()
    mocker.patch.object(controller, 'process_hardware_information', process_hw)
    mocker.patch.object(controller, 'process_measurement_data', process_meas)

    payload = object()
    ControllerABC.on_data_ready(controller, payload)

    process_hw.assert_called_once_with(payload)
    process_meas.assert_not_called()


def test_on_data_ready_dispatches_to_measurement_handler(mocker):
    """Check that measurement data are dispatched to measurement handler."""
    controller = object.__new__(ControllerABC)
    process_hw = mocker.Mock()
    process_meas = mocker.Mock()
    mocker.patch.object(controller, 'process_hardware_information', process_hw)
    mocker.patch.object(controller, 'process_measurement_data', process_meas)
    mocker.patch.object(controller, 'is_measurement_data',
                        return_value=True)

    payload = object()
    ControllerABC.on_data_ready(controller, payload)

    process_meas.assert_called_once_with(payload)
    process_hw.assert_not_called()
