"""Tests for data-processing hooks in ViPErinoController."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    'Florian Dörr (@FlorianDoerr)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2026-03-31'
__license__ = 'GPLv3+'

from contextlib import nullcontext

from viperleed.gui.measure.classes.datapoints import QuantityInfo
from viperleed.gui.measure.controller.viperinocontroller import ViPErinoController


def test_is_measurement_data_distinguishes_list_from_dict():
    """Check data-type detection for measurement vs hardware data."""
    fake_ctrl = object()
    assert ViPErinoController.is_measurement_data(fake_ctrl, [1, 2, 3])
    assert not ViPErinoController.is_measurement_data(fake_ctrl, {'a': 1})


def test_process_hardware_information_updates_state_and_emits(mocker):
    """Check hardware-information processing side effects."""
    fake_ctrl = mocker.Mock()
    fake_ctrl.lock = nullcontext()
    fake_ctrl.hardware_info_arrived = mocker.Mock()

    payload = {'firmware': 'x.y.z'}
    ViPErinoController.process_hardware_information(fake_ctrl, payload)

    assert fake_ctrl.hardware == payload
    fake_ctrl._check_box_id.assert_called_once_with()
    fake_ctrl._ViPErinoController__check_measurements_possible.assert_called_once_with()
    fake_ctrl.hardware_info_arrived.emit.assert_called_once_with()


def test_process_measurement_data_stores_and_converts(mocker):
    """Check measurement-data processing and completion signaling."""
    fake_ctrl = mocker.Mock()
    fake_ctrl._ViPErinoController__adc_measurement_types = (
        QuantityInfo.I0,
        QuantityInfo.TEMPERATURE,
        None,
        )
    fake_ctrl.measurements = {}
    fake_ctrl._ViPErinoController__convert_i0_value.side_effect = (
        lambda value: value + 0.5
        )
    fake_ctrl.measures.side_effect = lambda quantity: (
        quantity is QuantityInfo.TEMPERATURE
        )

    ViPErinoController.process_measurement_data(fake_ctrl, [1.0, 2.0, 3.0])

    assert fake_ctrl.measurements[QuantityInfo.I0] == [1.5]
    assert fake_ctrl.measurements[QuantityInfo.TEMPERATURE] == [2.0]
    fake_ctrl._ViPErinoController__convert_thermocouple_voltages.assert_called_once_with()
    fake_ctrl.measurements_done.assert_called_once_with()
