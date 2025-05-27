"""Module viperinocontroller of viperleed.gui.measure.controller.

Defines the ViPErinoController class for interacting with a ViPErLEED
Data Acquisition box.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    'Florian Dörr (@FlorianDoerr)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2021-07-08'
__license__ = 'GPLv3+'

from PyQt5 import QtCore as qtc
from PyQt5 import QtSerialPort as qts

from viperleed.gui.measure import hardwarebase as base
from viperleed.gui.measure.controller.abc import ControllerErrors
from viperleed.gui.measure.controller.abc import MeasureControllerABC
from viperleed.gui.measure.datapoints import QuantityInfo
from viperleed.gui.measure.classes.settings import NotASequenceError


# TODO: it looks like the nominal ADC frequencies are not necessarily
#       quite correct. This is at the origin of a misalignment on the
#       time axis when using multiple controllers. It appears as if
#       the 50Hz is quite correct, but the 500Hz can deviate a few
#       percent from the nominal value. It may make sense to have a
#       tool for calibrating such frequencies.


_MANDATORY_CMD_NAMES = (
    "PC_AUTOGAIN", "PC_CONFIGURATION", "PC_SET_UP_ADCS", "PC_OK", "PC_RESET",
    "PC_SET_VOLTAGE", "PC_ERROR", "PC_CALIBRATION", "PC_MEASURE_ONLY",
    "PC_CHANGE_MEAS_MODE", "PC_STOP", "PC_SET_VOLTAGE_ONLY", "PC_SET_SERIAL_NR"
    )


class ViPErinoErrors(base.ViPErLEEDErrorEnum):
    """Errors specific to Arduino-based ViPErLEED controllers."""
    TOO_MANY_MEASUREMENT_TYPES = (
        150,
        "Can measure up to {} quantities (one per ADC), but {} were requested."
        )
    OVERLAPPING_MEASUREMENTS = (
        151,
        "Cannot measure {} and {} at the same time. "
        "Measurements are performed by the same ADC."
        )
    INVALID_REQUEST = (
        152,
        "Cannot measure {0}. If your hardware can measure {0}, "
        "update its corresponding configuration file (check also typos!)."
        )
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
        ('controller', 'firmware_version'),  # also mandatory on serial
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
        """Return the initial time delay of a measurement (msec).

        Returns
        -------
        initial_delay : float
            The time interval between when a measurement was triggered
            and when the measurement was actually acquired. If multiple
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
            base.emit_error(
                self, ControllerErrors.INVALID_SETTING_WITH_FALLBACK,
                '', 'measurement_settings/num_meas_to_average', nr_average
                )
        return (3 + (nr_average - 1) / 2) * self.measurement_interval

    @property
    def measurement_interval(self):
        """Return the time interval between measurements (msec)."""
        update_rate_raw = self.settings.get('controller', 'update_rate',
                                            fallback='4')
        try:
            meas_f = self.settings.getfloat('adc_update_rate', update_rate_raw,
                                            fallback=50)
        except (TypeError, ValueError):
            # pylint: disable=redefined-variable-type
            # Seems a pylint bug.
            meas_f = 50.0
            base.emit_error(self, ControllerErrors.INVALID_SETTINGS,
                            f"adc_update_rate/{update_rate_raw}", "")
        return 1000 / meas_f

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
        version = settings.get("controller", "firmware_version", fallback="-1")
        try:
            fversion = float(version)
        except (TypeError, ValueError):
            fversion = -1.0
        if fversion < 0:
            base.emit_error(self, ControllerErrors.INVALID_SETTINGS,
                            "controller/firmware_version", "")
            return False

        mandatory_commands = [("available_commands", cmd)
                              for cmd in _MANDATORY_CMD_NAMES]
        # If new commands are added in newer versions, do:
        # if version >= something:
        #     mandatory_commands.extend([new_command_1, new_command_2, ...])

        # pylint: disable=protected-access
        # Need to access the classe's _mandatory_settings rather
        # than the one of self to prevent mandatory settings from
        # newer firmware versions to stay if settings for earlier
        # versions are loaded.
        self._mandatory_settings = [*self.__class__._mandatory_settings,
                                    *mandatory_commands]

        return super().are_settings_ok(settings)

    def set_energy(self, energy, settle_time, *more_steps, trigger_meas=True):
        """Set energy with associated settling time.

        Take the energy (or energies), get setpoint energy (or
        energies) and convert it to an integer value for the DAC.
        Afterwards send energy and time to the hardware via the
        serial. The controller will automatically trigger and
        start measuring after setting the voltage if trigger_meas
        is True. Each energy will need a settling time associated
        with it which is the time the hardware will wait before
        setting the next energy.

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
            True if the controller is supposed to take measurements
            after the energy has been set. Default is True.

        Raises
        ------
        TypeError
            If an odd number of more_steps are given
        """
        # super() call needed to store time_to_trigger.
        super().set_energy(energy, settle_time, *more_steps,
                           trigger_meas=trigger_meas)
        # pylint: disable=too-many-locals
        # Count = 16. Quite a few, but makes it easier to understand
        # the nominal-energy-to-DAC-value conversion.
        # pylint: disable=redefined-variable-type
        # Seems a pylint bug.
        _cmd_name = 'PC_SET_VOLTAGE' if trigger_meas else 'PC_SET_VOLTAGE_ONLY'
        cmd = self.settings.get('available_commands', _cmd_name)
        try:
            v_ref_dac = self.settings.getfloat('measurement_settings',
                                               'v_ref_dac')
        except (TypeError, ValueError):
            base.emit_error(self, ControllerErrors.INVALID_SETTINGS,
                            'measurement_settings/v_ref_dac', "")
            v_ref_dac = 2.5

        dac_out_vs_nominal_energy = 10/1000  # 10 V / 1000 eV
        output_gain = 4  # Gain of the output stage on board
        conversion_factor = dac_out_vs_nominal_energy * 65536 / (v_ref_dac *
                                                                 output_gain)

        energies_and_times = [energy, settle_time, *more_steps]
        if len(more_steps) % 2:  # odd
            raise TypeError(f"{self.__class__.__name__}.set_energy: "
                            "Number of energy and time steps do not match. "
                            "Expected an even number of arguments, found "
                            f"{len(more_steps) + 2} arguments.")

        for i, energy in enumerate(energies_and_times[::2]):
            energy = self.true_energy_to_setpoint(energy)
            energy = int(round(energy * conversion_factor))
            energy = max(0, min(energy, 65535))
            energies_and_times[2*i] = energy

        # Since we may wait a potentially long time here, we
        # explicitly set the timeout keyword to be the timeout
        # from the settings plus the total waiting time.
        try:
            timeout = self.settings.getint("serial_port_settings", "timeout",
                                           fallback=0)
        except (TypeError, ValueError):
            base.emit_error(self, ControllerErrors.INVALID_SETTINGS,
                            'serial_port_settings/timeout', "")
            timeout = 0
        timeout = max(timeout, 0) + sum(energies_and_times[1::2])
        self.send_message(cmd, energies_and_times, timeout=timeout)

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
            # pylint: disable=redefined-variable-type
            # Seems a pylint bug
            num_meas_to_average = 1
            base.emit_error(
                self, ControllerErrors.INVALID_SETTING_WITH_FALLBACK,
                '', 'measurement_settings/num_meas_to_average', 1
                )
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
            base.emit_error(self, ControllerErrors.INVALID_SETTINGS,
                            'controller/update_rate', '')
            return
        message = [update_rate, *self.__adc_channels[:-1]]
        self.send_message(cmd, message)

    def on_data_ready(self, data):
        """Receive and store data from the serial."""
        # We may receive two types of data: hardware&firmware
        # information (a dictionary) and actual measurements
        if isinstance(data, dict):
            # Got hardware info
            self.hardware = data
            return

        # Otherwise it is a list of data: [ADC0, ADC1, LM35]
        # The order is the same as stored in the config file under
        # controller/measurement_devices. The ADC channels chosen
        # via set_measurements() determine which value was measured
        for value, quantity in zip(data, self.__adc_measurement_types):
            if quantity is None:
                # Was not requested (but measured nontheless)
                continue
            if quantity is QuantityInfo.I0:
                # The value returned by the Arduino is correctly
                # in microamps only under the assumption that
                #   (i) in 0--10V range 1uA produces 1V at the BNC
                #  (ii) in 0--2.5V range 1mA produces 1V at the BNC
                try:
                    conv_ = self.settings.getfloat('conversions', 'i0_gain',
                                                   fallback=1.0)
                except (TypeError, ValueError):
                    conv_ = 1.0
                value *= conv_
            self.measurements[quantity] = [value]

        # TODO: here we should convert the mV of the thermocouple
        # and the cold junction into a compensated temperature

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
        super().measure_now()
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
        super().abort_and_reset()
        pc_reset = self.settings.get('available_commands', 'PC_RESET')
        self.send_message(pc_reset)

        self.reset_preparation_todos()
        self.__adc_measurement_types = []
        self.__adc_channels = []
        self.hardware = {}
        self.measurements = {}
        self.first_energies_and_times = []

    def set_measurements(self, quantities):
        """Decide what to measure.

        Receive requested measurement types from MeasurementABC
        class, check if request is valid and set channels
        accordingly.

        Parameters
        ----------
        quantities : Sequence
            All the quantities to be measured. Each should be
            convertible to a QuantityInfo.

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
            base.emit_error(self, ControllerErrors.INVALID_SETTINGS,
                            'controller/measurement_devices', '')
            return

        n_devices = len(measurement_devices)
        self.__adc_measurement_types = [None]*n_devices
        self.__adc_channels = [0]*n_devices
        if len(quantities) > n_devices:
            base.emit_error(self, ViPErinoErrors.TOO_MANY_MEASUREMENT_TYPES,
                            n_devices, len(quantities))
            return
        for quantity in quantities:
            for i, measurement_device in enumerate(measurement_devices):
                if quantity not in measurement_device:
                    continue
                if self.__adc_measurement_types[i] is not None:
                    base.emit_error(
                        self, ViPErinoErrors.OVERLAPPING_MEASUREMENTS,
                        quantity, self.__adc_measurement_types[i].label
                        )
                    return

                try:
                    channel = self.settings.getint('controller', quantity,
                                                   fallback=-1)
                except (TypeError, ValueError):
                    # pylint: disable=redefined-variable-type
                    # Seems a pylint bug.
                    # Cannot convert to int
                    channel = -1
                if channel < 0:
                    base.emit_error(self, ControllerErrors.INVALID_SETTINGS,
                                    f'controller/{quantity}', '')
                    return

                self.__adc_channels[i] = channel
                self.__adc_measurement_types[i] = (
                    QuantityInfo.from_label(quantity)
                    )
                break
            else:
                base.emit_error(self, ViPErinoErrors.INVALID_REQUEST, quantity)
                return
        super().set_measurements(quantities)

    def set_continuous_mode(self, continuous=True):
        """Set continuous mode.

        If continuous is true the controller will continue
        returning measurements without further instructions.

        Parameters
        ----------
        continuous : bool, optional
            Wether continuous mode should be on.
            Default is True

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
        stop = self.settings.get('available_commands', 'PC_STOP')
        self.send_message(stop)

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
                      ctrl.hardware, flush=True)
        for thread in threads:
            thread.quit()
        for thread in threads:
            # wait max 100 ms for each thread to quit
            if not thread.wait(100):
                thread.terminate()
                thread.wait()
        return device_list
