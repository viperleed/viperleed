"""Module viperinocontroller of viperleed.guilib.measure.controller.

========================================
   ViPErLEED Graphical User Interface
========================================

Created: 2021-07-08
Author: Michele Riva
Author: Florian Doerr

This module contains the definition of the ViPErinoController
class and its associated ViPErLEEDErrorEnum class ViPErinoErrors
which gives commands to the ViPErinoSerialWorker class.
"""

import time
from collections import defaultdict

from PyQt5 import QtCore as qtc
from PyQt5 import QtSerialPort as qts

# ViPErLEED modules
from viperleed.guilib.measure.controller.abc import (MeasureControllerABC,
                                                     ControllerErrors)
from viperleed.guilib.measure.hardwarebase import (ViPErLEEDErrorEnum,
                                                   emit_error)
from viperleed.guilib.measure.datapoints import QuantityInfo
from viperleed.guilib.measure.classes.settings import NotASequenceError


_MANDATORY_CMD_NAMES = (
    "pc_configuration", "pc_set_up_adcs", "pc_ok", "pc_reset",
    "pc_set_voltage", "pc_autogain", "pc_error", "pc_calibration",
    "pc_measure_only", "pc_change_meas_mode", "pc_stop", "pc_set_voltage_only",
    "pc_set_serial_nr"
    )

class ViPErinoErrors(ViPErLEEDErrorEnum):
    """Errors specific to Arduino-based ViPErLEED controllers."""
    TOO_MANY_MEASUREMENT_TYPES = (150,
                                  "The ViPErinoController can only handle "
                                  "as many measurement requests as there are "
                                  "ADCs installed but more have been "
                                  "requested.")
    OVERLAPPING_MEASUREMENTS = (151,
                                "At least two requested measurements concern "
                                "the same ADC. Consider installing another "
                                "controller to do both measurements in "
                                "parallel.")
    INVALID_REQUEST = (152,                                                     # TODO: fix
                       "A measurement type has been requested that is not "
                       "implemented on the controller side. Maybe "
                       "implemented type has not been added in the "
                       "controller configuration file.")
    REQUESTED_ADC_OFFLINE = (153,
                             "Requested measurement type needs associated "
                             "ADC to be available for measurements but it is "
                             "not.")


class ViPErinoController(MeasureControllerABC):
    """Controller class for the ViPErLEED Arduino Micro."""

    __request_info = qtc.pyqtSignal()

    _mandatory_settings = [
        *MeasureControllerABC._mandatory_settings,
        ('available_commands',),
        ('controller', 'measurement_devices'),
        ('controller', 'FIRMWARE_VERSION'),  # also mandatory on serial
        ('measurement_settings', 'v_ref_dac'),
        ]

    def __init__(self, settings=None, port_name='', sets_energy=False):
        """Initialise ViPErino controller object.

        Initialise prepare_todos dictionaries. The key is a
        string that is used to call a function, the value is a
        boolean that is used to determine if the connected
        function has already been called.

        Parameters
        ----------
        settings : ConfigParser
            The controller settings
        port_name : str, optional
            Name of the serial port to be used to communicate with
            the controller. This parameter is optional only in case
            settings contains a 'controller'/'port_name' option. If
            this is given, it will also be stored in the settings
            file, overriding the value that may be there. Default is
            an empty string.
        sets_energy : bool, optional
            Used to determine whether this controller is responsible
            for setting the electron energy by communicating with the
            LEED optics. Only one controller may be setting the energy.
            Default is False.

        Raises
        ------
        TypeError
            If no port_name is given, and none was present in the
            settings file.
        """
        super().__init__(settings=settings, port_name=port_name,
                         sets_energy=sets_energy)
        # Initialise dictionaries for the measurement preparation.
        self.begin_prepare_todos['get_hardware'] = True
        self.begin_prepare_todos['calibrate_adcs'] = True
        self.begin_prepare_todos['set_up_adcs'] = True
        if sets_energy:
            self.begin_prepare_todos['set_energy'] = True
        self.continue_prepare_todos['start_autogain'] = True
        self.__adc_measurement_types = []
        self.__adc_channels = []
        self.hardware = {}

    @property
    def initial_delay(self):
        """Return the initial time delay of a measurement in seconds.

        Returns
        -------
        initial_delay : float
            The time interval between when a measurement was requested
            and when the measurement was actually acquired. If mutliple
            measurements are averaged over, the "time when measurements
            are actually acqiuired" is the middle time between the
            beginning and the end of the measurement.
        """
        try:
            nr_average = self.settings.getint(
                'measurement_settings', 'num_meas_to_average', fallback=1
                )
        except (TypeError, ValueError):
            nr_average = 1
            emit_error(self, ControllerErrors.INVALID_SETTING_WITH_FALLBACK,
                       '', 'measurement_settings/num_meas_to_average', 1)
        return (3 + (nr_average - 1) / 2)*self.measurement_interval

    @property
    def name(self):
        """Return a name for this controller.

        The name will be unique only if this controller has been
        requested (and responded) for a serial number.

        Returns
        -------
        name : str
            The name of this controller
        """
        serial_nr = self.hardware.get('serial_nr', 'UNKNOWN_SERIAL_NR')
        return f"ViPErLEED {serial_nr}"

    @property
    def firmware_version(self):
        """Return firware version of the hardware (or the settings).

        Returns
        -------
        firmware_version: str
            Firmware version of the form "<major>.<minor>".
        """
        version = self.hardware.get("firmware", None)
        if not version:
            # Get it from the settings. Notice that the are_settings_ok
            # reimplementation already checks that the firmware version
            # in the settings is present and valid.
            version = self.settings.get("controller", "firmware_version")
        return version

    def are_settings_ok(self, settings):
        """Return whether a ViPErLEEDSettings is compatible with self."""
        try:
            version = settings.getfloat("controller", "firmware_version",
                                        fallback=-1.0)
        except (TypeError, ValueError):
            version = -1.0
        if version < 0:
            emit_error(self, ControllerErrors.INVALID_CONTROLLER_SETTINGS,
                       "controller/firmware_version")
            return False

        mandatory_commands = [("available_commands", cmd)
                              for cmd in _MANDATORY_CMD_NAMES]
        # If new commands are added in newer versions, do:
        # if version >= something:
        #     mandatory_commands.extend([new_command_1, new_command_2, ...])

        self._mandatory_settings = [*self.__class__._mandatory_settings,
                                    *mandatory_commands]

        return super().are_settings_ok(settings)

    def set_energy(self, energy, settle_time, *more_steps,
                   trigger_meas=True, **_):
        """Set energy with associated settling time.

        Take the energy (or energies), get setpoint energy (or
        energies) and convert it to an integer value for the DAC.
        Afterwards send energy and time to the hardware via the
        serial. The controller will automatically trigger and
        start measuring after setting the voltage. Each energy
        will need a settling time associated with it which is the
        time the hardware will wait before setting the next energy.

        Parameters
        ----------
        energy : float
            Nominal electron energy in electronvolts.
        settle_time : integer
            Interval in milliseconds that the controller will wait
            before deeming the LEED optics stable at set energy.
        *more_steps : Number
            If given, it should be an even number of elements.
            Odd elements are energies, even ones settle-time intervals.
            Multiple steps can be executed quickly after each other.
            The last step will be the final energy that is set and
            should ensure stabilization of the electronics.
        trigger_meas : bool, optional
            True if the controllers are supposed to take
            measurements after the energy has been set.
            Default is True.

        Raises
        ------
        TypeError
            If an odd number of more_steps are given
        """
        _cmd_name = 'PC_SET_VOLTAGE' if trigger_meas else 'PC_SET_VOLTAGE_ONLY'
        cmd = self.settings.get('available_commands', _cmd_name)
        try:
            v_ref_dac = self.settings.getfloat('measurement_settings',
                                               'v_ref_dac')
        except (TypeError, ValueError):
            emit_error(self, ControllerErrors.INVALID_CONTROLLER_SETTINGS,
                       'measurement_settings/v_ref_dac')

        dac_out_vs_nominal_energy = 10/1000  # 10V / 1000 eV
        output_gain = 4  # Gain of the output stage on board
        conversion_factor = dac_out_vs_nominal_energy * 65536 / (v_ref_dac *
                                                                 output_gain)

        energies_and_times = [energy, settle_time, *more_steps]
        if len(more_steps) % 2:  # odd
            raise TypeError(f"{self.__class__.__name__}.set_energy: "
                            "Number of energy and time steps do not match. "
                            "Expected an even number of arguments, found "
                            f"{len(more_steps) + 2} arguments.")
        number_of_steps = int(len(energies_and_times)/2)

        for i in range(number_of_steps):
            tmp_energy = self.true_energy_to_setpoint(energies_and_times[2*i])
            tmp_energy = int(round(tmp_energy * conversion_factor))
            tmp_energy = max(min(tmp_energy, 65535), 0)
            energies_and_times[2*i] = tmp_energy
        self.send_message(cmd, energies_and_times)

    def start_autogain(self):
        """Determine starting gain.

        Determine gain after setting the starting energy.

        Returns
        -------
        None.
        """
        # TODO: Had issues in the beginning back then on omicron
        # (Gain too high)
        cmd = self.settings.get('available_commands', 'PC_AUTOGAIN')
        self.send_message(cmd)

    def set_up_adcs(self):
        """Set up ADCs.

        Set number of measurements to be taken and channels to
        be used. The used channels have to be calibrated beforehand,
        otherwise measurements will not be taken and an error will
        occur.

        Returns
        -------
        None.
        """
        cmd = self.settings.get('available_commands', 'PC_SET_UP_ADCS')
        try:
            num_meas_to_average = self.settings.getint(
                'measurement_settings', 'num_meas_to_average', fallback=1
                )
        except (TypeError, ValueError):
            num_meas_to_average = 1
            emit_error(self, ControllerErrors.INVALID_SETTING_WITH_FALLBACK,
                       '', 'measurement_settings/num_meas_to_average', 1)
        message = [num_meas_to_average, *self.__adc_channels[:2]]
        self.send_message(cmd, message)

    def get_hardware(self):
        """Get hardware connected to micro controller.

        This function determines the hardware controlled by
        the micro controller. It has to be run at least once
        after starting the controller up. It will tell both
        the micro controller and this controller class which
        hardware is present in the setup. Trying to do
        anything but getting the hardware or resetting the
        controller before executing this function will result
        in an error.

        Upon receiving this command the micro controller will
        first return the firmware version and then the hardware
        configuration. The conversion is already handled in the
        viperleed serial.

        Returns
        -------
        None.
        """
        cmd = self.settings.get('available_commands', 'PC_CONFIGURATION')
        self.send_message(cmd)

    def calibrate_adcs(self):
        """Calibrate ADCs.

        Set update rate, select channels for calibration and
        calibrate them. This function has to be run at least once
        after starting the micro controller. It will store all of
        the necessary calibration data on the micro controller
        itself. It is only done for the channels selected. If the
        controller is trying to measure channels that have not been
        calibrated yet an error will occur.

        Returns
        -------
        None.
        """
        cmd = self.settings.get('available_commands', 'PC_CALIBRATION')
        try:
            update_rate = self.settings.getint('controller', 'update_rate',
                                               fallback=4)
        except (TypeError, ValueError):
            # Cannot convert to int
            emit_error(self, ControllerErrors.INVALID_CONTROLLER_SETTINGS,
                       'controller/update_rate')
            return
        message = [update_rate, *self.__adc_channels[:-1]]
        self.send_message(cmd, message)

    def receive_measurements(self, receive):
        """Receive measurements from the serial.

        For measurements:
        Append measurements to the according section. Done via the
        settings. The chosen ADC channels will determine which
        value was measured.

        The receive parameter has to have the measurements listed
        in the same order as the measurement devices are listed in
        the controller configuration, otherwise the measurements
        will not be saved in the correct section.

        For hardware:
        Save received data into a dictionary to store the hardware
        configuration for future use. Hardware is only sent after
        the get_hardware function has been executed and before the
        calibration is done.

        Parameters
        ----------
        receive : list or dict
            Data received from the serial.
            list if a measurement has been received
            dict if the hardware configuration has been received

        Returns
        -------
        None.
        """
        if isinstance(receive, dict):
            # Got hardware info
            self.hardware = receive
            return

        # Otherwise it is data
        for i, measurement in enumerate(self.__adc_measurement_types):
            if measurement is not None:
                self.measurements[measurement] = [receive[i]]
        self.measurements_done()

    def measure_now(self):
        """Take a measurement.

        Measure without setting the energy. This function is
        supposed to be used in time resolved measurements and
        by secondary controllers which will not set the energy.

        Returns
        -------
        None.
        """
        cmd = self.settings.get('available_commands', 'PC_MEASURE_ONLY')
        self.send_message(cmd)

    def abort_and_reset(self):
        """Abort current task and reset the controller.

        Abort what the controller is doing right now, reset
        it and return to waiting for further instructions.

        Returns
        -------
        None.
        """
        pc_reset = self.settings.get('available_commands', 'PC_RESET')
        self.send_message(pc_reset)

        self.reset_preparation_todos()
        self.__adc_measurement_types = []
        self.__adc_channels = []
        self.hardware = {}
        self.measurements = {}
        self.__energies_and_times = []

    def set_measurements(self, quantities):
        """Decide what to measure.

        Receive requested measurement types from MeasurementABC
        class, check if request is valid and set channels
        accordingly.

        Parameters
        ----------
        requested : list of strings
            Contains all of the requested
            measurement types.

        Returns
        -------
        None.
        """
        if not self.settings:
            return
        try:
            measurement_devices = self.settings.getsequence(
                'controller', 'measurement_devices'
                )
        except NotASequenceError:
            emit_error(self, ControllerErrors.INVALID_CONTROLLER_SETTINGS,
                       'controller/measurement_devices')
            return

        n_devices = len(measurement_devices)
        self.__adc_measurement_types = [None]*n_devices
        self.__adc_channels = [0]*n_devices
        if len(quantities) > n_devices:
            emit_error(self, ViPErinoErrors.TOO_MANY_MEASUREMENT_TYPES)
            return
        for quantity in quantities:
            for i, measurement_device in enumerate(measurement_devices):
                if quantity not in measurement_device:
                    continue
                if self.__adc_measurement_types[i] is not None:
                    emit_error(self, ViPErinoErrors.OVERLAPPING_MEASUREMENTS)
                    return
                try:
                    channel = self.settings.getint('controller', quantity,
                                                   fallback=-1)
                except (TypeError, ValueError):
                    # Cannot convert to int
                    channel = -1
                if channel < 0:
                    emit_error(
                        self, ControllerErrors.INVALID_CONTROLLER_SETTINGS,
                        f'controller/{quantity}'
                        )
                    return
                self.__adc_channels[i] = channel
                self.__adc_measurement_types[i] = (
                    QuantityInfo.from_label(quantity)
                    )
                break
            else:
                emit_error(self, ViPErinoErrors.INVALID_REQUEST)
                return
        super().set_measurements(quantities)

    def set_continuous_mode(self, continuous):
        """Set continuous mode.

        If continuous is true the controller will continue
        returning measurements without further instructions.

        Parameters
        ----------
        continuous : bool
            Wether continuous mode should be on.

        Returns
        -------
        None.
        """
        super().set_continuous_mode(continuous)
        cmd = self.settings.get('available_commands', 'PC_CHANGE_MEAS_MODE')
        mode_on = int(bool(continuous))
        self.send_message(cmd, [mode_on, 0])

    def stop(self):
        """Stop.

        Stop whatever the controller is doing right now
        and return to idle state.

        Returns
        -------
        None.
        """
        if not self.settings:
            return
        try:
            super().stop()
        except AttributeError:
            # super().stop accesses serial_busy, which may
            # not exist if settings are somewhat funky.
            return
        stop = self.settings.get('available_commands', 'pc_stop')
        self.send_message(stop)

    @property
    def measurement_interval(self):
        """Return the time interval between measurements in seconds."""
        update_rate_raw = self.settings.get('controller', 'update_rate',
                                            fallback='4')
        try:
            meas_f = self.settings.getint('adc_update_rate', update_rate_raw,
                                          fallback=50)
        except (TypeError, ValueError):
            meas_f = 50
            emit_error(self, ControllerErrors.INVALID_CONTROLLER_SETTINGS,
                       f"adc_update_rate/{update_rate_raw}")
        return 1/meas_f

    def list_devices(self):
        """List Arduino Micro VipErLEED hardware.

        This function will take between 70 to slightly more
        than 100 ms to run.
        """
        ports = qts.QSerialPortInfo().availablePorts()
        device_list = []
        threads = []
        controllers = []
        for port in ports:
            ctrl = ViPErinoController(port_name=port.portName(),
                                      sets_energy=False)
            threads.append(qtc.QThread())
            ctrl.moveToThread(threads[-1])
            controllers.append(ctrl)
        for thread in threads:
            thread.start(priority=thread.TimeCriticalPriority)
        for ctrl in controllers:
            self.__request_info.connect(ctrl.get_hardware)
        self.__request_info.emit()
        if controllers:
            controllers[0].serial.port.waitForReadyRead(100)
        for ctrl in controllers:
            serial_nr = ctrl.hardware.get('serial_nr', None)
            ctrl.disconnect_()
            if serial_nr:
                device_list.append(f"{ctrl.name} ({ctrl.port_name})")
            else:
                print("Not a ViPErLEED controller at", ctrl.port_name,
                      flush=True)
        for thread in threads:
            thread.quit()
        time.sleep(0.001)
        return device_list
