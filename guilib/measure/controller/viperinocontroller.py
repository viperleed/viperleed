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
import re
import threading
import time

from PyQt5 import QtCore as qtc
from PyQt5 import QtSerialPort as qts

from viperleed.guilib.measure import hardwarebase as base
from viperleed.guilib.measure.classes.abc import QObjectSettingsErrors
from viperleed.guilib.measure.classes.abc import SettingsInfo
from viperleed.guilib.measure.classes.datapoints import QuantityInfo
from viperleed.guilib.measure.classes.settings import DefaultSettingsError
from viperleed.guilib.measure.classes.settings import NotASequenceError
from viperleed.guilib.measure.classes.thermocouple import Thermocouple
from viperleed.guilib.measure.controller import abc
from viperleed.guilib.measure.dialogs.settingsdialog import SettingsHandler
from viperleed.guilib.measure.dialogs.settingsdialog import SettingsTag

# For settings dialog:
from viperleed.guilib.measure.controller import _vprctrlsettings as _settings


_MANDATORY_CMD_NAMES = (
    'PC_AUTOGAIN',
    'PC_CALIBRATION',
    'PC_CHANGE_MEAS_MODE',
    'PC_CONFIGURATION',
    'PC_ERROR',
    'PC_MEASURE_ONLY',
    'PC_OK',
    'PC_RESET',
    'PC_SET_SERIAL_NR',
    'PC_SET_UP_ADCS',
    'PC_SET_VOLTAGE',
    'PC_SET_VOLTAGE_ONLY',
    'PC_STOP',
    )

_INVOKE = qtc.QMetaObject.invokeMethod


class ViPErinoErrors(base.ViPErLEEDErrorEnum):
    """Errors specific to Arduino-based ViPErLEED controllers."""

    TOO_MANY_MEASUREMENT_TYPES = (
        150,
        'Can measure up to {} quantities (one per ADC), but {} were requested.'
        )
    OVERLAPPING_MEASUREMENTS = (
        151,
        'Cannot measure {} and {} at the same time. '
        'Measurements are performed by the same ADC.'
        )
    UNSUPPORTED_QUANTITY = (
        152,
        'Cannot measure {0}. If your hardware can measure {0}, update '
        'its corresponding configuration file. Make sure to check for '
        'spelling mistakes.'
        )
    REQUESTED_ADC_OFFLINE = (
        153,
        'Cannot measure {} because its associated ADC was not detected.'
        )
    CANNOT_CONVERT_THERMOCOUPLE = (
        154,
        'Invalid/missing setting "conversions"/"thermocouple_type". '
        'Cannot convert temperature to °C. Temperatures will appear in '
        'millivolts (i.e., thermocouple voltage).{}'
        )
    HARDWARE_INFO_MISSING = (
        155,
        'No hardware information present. Cannot '
        '{} before get_hardware() is called.'
        )
    ERROR_WRONG_BOX_ID = (
        156,
        'The box ID {arduino_id} of the hardware does not match the ID '
        '{local_id} of the software. Please report this error to the '
        'ViPErLEED team. Make sure to include which hardware, which '
        'firmware version and which ViPErLEED version was used.'
        )


# too-many-instance-attributes, too-many-public-methods
class ViPErinoController(abc.MeasureControllerABC):
    """Controller class for the ViPErLEED Arduino Micro."""

    # The box ID is the identifier that differentiates viperleed
    # controller types from each other. The ID of the class must
    # match the box ID returned by the hardware controller. The
    # box ID must never be changed! The box ID is 1-based, 0
    # should never be used for any controller.
    box_id = 1

    cls_lock = threading.RLock()  # Thread safe access to class properties
    _mandatory_settings = (
        # pylint: disable=protected-access
        *abc.MeasureControllerABC._mandatory_settings,
        ('available_commands',),
        ('controller', 'measurement_devices'),
        ('controller', 'firmware_version'),  # also mandatory on serial
        ('energy_calibration', 'v_ref_dac'),
        )

    hardware_info_arrived = qtc.pyqtSignal()

    def __init__(self, parent=None, settings=None,
                 address='', sets_energy=False):
        """Initialise ViPErino controller object.

        Initialise prepare_todos dictionaries. The key is a
        string that is used to call a function, the value is a
        boolean that is used to determine if the connected
        function has already been called.

        Parameters
        ----------
        settings : ConfigParser
            The controller settings
        address : str, optional
            Address (the serial port) to be used to communicate with
            the controller. This parameter is optional only in case
            settings contains a 'controller'/'address' option. If
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
            If no address is given, and none was present in the
            settings file.
        """
        super().__init__(parent=parent, settings=settings,
                         address=address, sets_energy=sets_energy)
        # Initialise dictionaries for the measurement preparation.
        self.begin_prepare_todos['get_hardware'] = True
        self.begin_prepare_todos['calibrate_adcs'] = True
        self.begin_prepare_todos['set_up_adcs'] = True
        if sets_energy:
            self.begin_prepare_todos['set_energy'] = True
        self.continue_prepare_todos['start_autogain'] = True

        self.hardware = {}

        # One lock per instance to make the access to self.hardware
        # safe: when reading/writing instance.hardware, one can
        # .acquire() instance.lock (or use a context manager). The
        # class-level lock instance.cls_lock should be used if any
        # attribute common to all instances is to be modified instead.
        self.lock = threading.RLock()

        # __adc_measurement_types[i] contains the QuantityInfo
        # that should be measured by ADC[i] (ADC0, ADC1, ..., LM35)
        self.__adc_measurement_types = []

        # __adc_channels[i] contains which channel of ADC[i]
        # is needed to measure __adc_measurement_types[i]
        self.__adc_channels = []

        # __added_cold_junction keeps track of whether we would
        # like to also measure the cold-junction temperature,
        # even if the user did not ask for it, if possible.
        self.__added_cold_junction = False

        # The thermocouple.Thermocouple instance for this unit.
        # None if not present (or if TEMPERATURE was not measured)
        self._thermocouple = None

    @property
    def initial_delay(self):
        """Return the initial time delay of a measurement (msec).

        Returns
        -------
        initial_delay : float
            The time interval between when a measurement was triggered
            and when the measurement was actually acquired. If multiple
            measurements are averaged over, the "time when measurements
            are actually acquired" is the middle time between the
            beginning and the end of the measurement.
        """
        n_intervals = (3 + (self.nr_samples - 1) / 2)
        return n_intervals * self.measurement_interval

    @property
    def measurement_interval(self):
        """Return the time interval between measurements (msec)."""
        update_rate_raw = self.settings.get('measurement_settings',
                                            'adc_update_rate', fallback='4')
        try:
            meas_f = self.settings.getfloat('adc_update_rate', update_rate_raw,
                                            fallback=50)
        except (TypeError, ValueError):
            # pylint: disable=redefined-variable-type
            # Seems a pylint bug.
            meas_f = 50.0
            self.emit_error(QObjectSettingsErrors.INVALID_SETTINGS,
                            f'adc_update_rate/{update_rate_raw}', '')
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
        with self.lock:
            serial_nr = self.hardware.get('serial_nr', 'UNKNOWN_SERIAL_NR')
        return f"ViPErLEED {serial_nr}"

    @property
    def name_clean(self):
        """Return a version of .name suitable for file names."""
        # Fall back on the settings if we have no serial number info
        with self.lock:
            try:
                name = f"ViPErLEED {self.hardware['serial_nr']}"
            except KeyError:
                name = self.settings.get("controller", "device_name")
        return base.as_valid_filename(name)

    @property
    def firmware_version(self):
        """Return firmware version of the hardware (or the settings).

        Returns
        -------
        firmware_version : hardwarebase.Version
            Firmware version of the device, if the information is
            available, otherwise the one in self.settings.
        """
        with self.lock:
            version = self.hardware.get("firmware", None)
        if version is None:
            # Get it from the settings. Notice that the
            # are_settings_invalid extension already checks that the
            # firmware version in the settings is present and valid
            version = base.Version(
                self.settings.get("controller", "firmware_version")
                )
        return version

    @property
    def thermocouple(self):
        """Return the Thermocouple used for measuring temperatures.

        Returns
        -------
        thermocouple : Thermocouple or None
            Returns None in case it was not possible to derive
            the correct thermocouple from the settings.
        """
        if self._thermocouple is None:
            tc_type = self.settings.get('conversions', 'thermocouple_type',
                                        fallback=None)
            if tc_type is None:
                self.emit_error(ViPErinoErrors.CANNOT_CONVERT_THERMOCOUPLE,
                                '\nInfo: missing entry in settings.')
                return None
            try:
                self._thermocouple = Thermocouple(tc_type)
            except ValueError as err:
                # Unknown thermocouple type
                self.emit_error(ViPErinoErrors.CANNOT_CONVERT_THERMOCOUPLE,
                                f'\nInfo: {err}.')
                return None
        return self._thermocouple

    @thermocouple.setter
    def thermocouple(self, new_tc):
        """Set a new thermocouple. The type is stored in settings."""
        if new_tc is None or isinstance(new_tc, Thermocouple):
            self._thermocouple = new_tc
        else:
            try:
                self._thermocouple = Thermocouple(new_tc)
            except (ValueError, TypeError):
                return
        try:
            type_ = self._thermocouple.type_
        except AttributeError:  # None
            type_ = 'None'
        self.settings.set('conversions', 'thermocouple_type', type_)

    @property
    def time_to_first_measurement(self):
        """Return the interval between trigger and 1st measurement (msec)."""
        return (self.nr_samples + 2) * self.measurement_interval

    def are_settings_invalid(self, settings):
        """Check if there are any invalid settings.

        Parameters
        ----------
        settings : ViPErLEEDSettings
            The new settings.

        Returns
        -------
        invalid_settings : list of tuples
            Invalid _mandatory_settings of self as a list of tuples.
            The first entry in each tuple can be either '<section>',
            '<section>/<option>', or
            '<section>/<option> not one of <value1>, <value2>, ...'.
            Further optional entries may be added by subclasses. They
            specify additional information on what is wrong with each
            invalid setting.
        """
        invalid_settings = settings.has_settings(('controller',
                                                  'firmware_version'))
        if not invalid_settings:
            try:
                version = base.Version(
                    settings['controller']['firmware_version']
                    )
            except (TypeError, ValueError) as err:
                invalid_settings.append(('controller/firmware_version',
                                         f'Info: Value is invalid -- {err}'))
        else:
            invalid_settings = [(invalid,) for invalid in invalid_settings]
        if invalid_settings:
            return invalid_settings
        self._thermocouple = None  # In case it changed

        mandatory_cmd_names = list(_MANDATORY_CMD_NAMES)
        if version and version >= '0.7':
            mandatory_cmd_names.append('PC_DEBUG')
        mandatory_commands = (('available_commands', cmd)
                              for cmd in mandatory_cmd_names)

        # pylint: disable=protected-access
        # Need to access the class's _mandatory_settings rather
        # than the one of self to prevent mandatory settings from
        # newer firmware versions to stay if settings for earlier
        # versions are loaded.
        self._mandatory_settings = (*self.__class__._mandatory_settings,
                                    *mandatory_commands)

        invalid_settings.extend(super().are_settings_invalid(settings))
        return invalid_settings

    def available_adcs(self):
        """Return a list of available ADC measurements.

        Returns
        -------
        adcs : list of dict
            The list is taken from the settings, and is checked
            against which ADCs are actually available from the
            hardware configuration, if known. Each element is a
            {QuantityInfo: channel} dictionary, where channel is
            the channel to be used for measuring the quantity.

        Raises
        ------
        NotASequenceError
            If the "controller"/"measurement_devices" option in the
            settings cannot be converted to a sequence.
        KeyError
            In case the "controller" section is missing any of the
            quantities in "measurement_devices" (corresponding to
            which channel is used for measuring the quantity).
        TypeError, ValueError
            If any of the channels in the settings is not an integer
        ValueError
            If any of the quantities that are measurable from the
            ADCs is not a valid QuantityInfo
        RuntimeError
            If the settings is inconsistent: (1) ADC quantities
            are repeated (also accounting for aliases); (2) ADC
            quantities are measured on the same channel; (3) The
            number of ADCs in the settings is different from those
            available at the hardware level (may mean that power
            is disconnected)
        """
        # Start from the settings.
        # The next line raises NotASequenceError if something is wrong
        settings_adcs = self.settings.getsequence("controller",
                                                  "measurement_devices")

        qinfo_adcs = []
        for quantities in settings_adcs:
            # Next line raises KeyError if there is no channel for
            # a quantity, ValueError or TypeError if the channels
            # cannot be converted to integer, and ValueError if a
            # quantity is unknown
            channels_idx = {
                QuantityInfo.from_label(q): int(self.settings["controller"][q])
                for q in quantities
                }
            if len(quantities) != len(channels_idx):
                raise RuntimeError(f"Duplicate quantities in {quantities}")
            if len(channels_idx) != len(set(channels_idx.values())):
                raise RuntimeError(
                    f"Some of {quantities} is measured with the same channel"
                    )
            qinfo_adcs.append(channels_idx)

        # Now compare with the info in the hardware.
        with self.lock:
            hardware = self.hardware.copy()
        active_adcs = [k for k, present in hardware.items()
                       if self.__is_adc(k) and present]
        has_power = any(k for k in active_adcs if 'adc' in k)  # skip LM35
        if has_power and len(active_adcs) != len(qinfo_adcs):
            raise RuntimeError(
                f"Number of measurement_devices ({len(qinfo_adcs)}) "
                "inconsistent with number of devices detected on board"
                )
        return qinfo_adcs

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
            update_rate = self.settings.getint('measurement_settings',
                                               'adc_update_rate', fallback=4)
        except (TypeError, ValueError):
            # Cannot convert to int
            self.emit_error(QObjectSettingsErrors.INVALID_SETTINGS,
                            'measurement_settings/adc_update_rate', '')
            return
        with self.lock:
            hardware = self.hardware.copy()
        if not hardware:
            self.emit_error(ViPErinoErrors.HARDWARE_INFO_MISSING,
                            'calibrate the ADCs')
            return
        lm35_idx = list(hardware.keys()).index('lm35')
        message = [update_rate, *self.__adc_channels[:lm35_idx]]
        self.send_message(cmd, message)

    @qtc.pyqtSlot()
    @abc.ensure_connected
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

    def get_settings_handler(self):
        """Return a SettingsHandler object for displaying settings.

        The handler returned manages the following settings:
        'controller'/'firmware_version'
            Read-only text, with color indicating warnings
        'controller'/'device_name'
            Displays and allows changing the serial number.
        'controller' & 'conversions'
            'measurement_devices', their associated ADC channels,
            including input ranges, and all the conversion to be
            used to transform measurements into physical values.
            This whole complex section is "advanced", and will be
            shown only if there is something missing.
        'measurement_settings'/'adc_update_rate'
            Allows selecting the ADC update rate from given values.
        All the settings handled by MeasureControllerABC.

        Returns
        -------
        handler : SettingsHandler
            The handler used in a SettingsDialog to display the
            settings of this controller to users.
        """
        self.check_creating_settings_handler_is_possible()
        handler = SettingsHandler(self.settings, show_path_to_config=True)
        handler.add_option('controller', 'firmware_version',
                           handler_widget=_settings.FWVersionViewer(self),
                           display_name='Firmware version',
                           tags=SettingsTag.READ_ONLY | SettingsTag.REGULAR)
        handler.add_option('controller', 'device_name',
                           handler_widget=_settings.SerialNumberEditor(self),
                           display_name='Serial No.',
                           tags=SettingsTag.REGULAR)
        handler.add_from_handler(super().get_settings_handler())
        handler.add_complex_section(
            _settings.HardwareConfigurationEditor(controller=self)
            )
        widget = _settings.UpdateRateSelector(self)
        tooltip = ('<nobr>The frequency at which the </nobr>controller '
                  'acquires measurements. For regular measurements it is '
                  'preferable to set the measurement frequency to the '
                  'line frequency.')
        handler.add_option(
            'measurement_settings', 'adc_update_rate', handler_widget=widget,
            display_name='Measurement frequency', tooltip=tooltip
            )

        return handler

    @classmethod
    def is_matching_default_settings(cls, obj_info, config, match_exactly):
        """Determine if the default `config` file is for this controller.

        Parameters
        ----------
        obj_info : SettingsInfo or None
            The information that should be used to check `config`.
            If not None, then obj_info.more['firmware'], i.e., the
            firmware version stored in the box, is matched against
            the 'firmware_version' of `config`.
        config : ConfigParser
            The settings to check.
        match_exactly : bool
            Whether the firmware version should match exactly or not.

        Returns
        -------
        sorting_info : tuple
            A tuple that can be used to sort the detected settings.
            The firmware minor of matching settings is returned as an
            indicator of how well the firmware matches. If an exact
            match is required, the firmware version has to be exactly
            the same. If no exact match is required, the highest
            firmware minor will be preferred. An empty tuple signifies
            that `config` does not match the requirements.
        """
        ver = cls._get_version(config)
        firmware_matches = (obj_info.more.get('firmware') == ver if obj_info
                            else False)
        if not match_exactly or firmware_matches:                               # TODO: check major in the future (compare class major to settings file major)
            # If match_exactly is False we are looking for any default
            # that may satisfy the requirements of the controller while
            # likely not knowing the firmware version of the controller.
            # If match_exactly is True, then we are looking for a
            # default that does specifically match the firmware of
            # the device we are looking at.
            return (ver.minor,)
        return ()

    @classmethod
    def is_matching_user_settings(cls, obj_info, config, match_exactly):
        """Determine if a `config` file is for this controller.

        Parameters
        ----------
        obj_info : SettingsInfo
            The information that should be used to check `config`.
            .more must contain the 'name' of the controller and the
            installed 'firmware' version.
        config : ConfigParser
            The settings to check.
        match_exactly : bool
            Whether the firmware version should match exactly or not.

        Returns
        -------
        sorting_info : tuple
            A tuple that can be used to sort the detected settings.
            The firmware minor of matching settings is returned as an
            indicator of how well the firmware matches. If an exact
            match is required, the firmware version has to be exactly
            the same, so all matching files will return the same
            sorting value. If no exact match is required, the highest
            firmware minor with a matching major will be preferred.
            An empty tuple signifies that `config` does not match
            the requirements.
        """
        super().is_matching_user_settings(obj_info, config, match_exactly)
        controller_name = config.get('controller', 'device_name',
                                     fallback=None)
        ver = cls._get_version(config)
        if obj_info.more['name'] != controller_name:
            return ()
        if obj_info.more['firmware'] == ver:
            return (ver.minor,)
        if (not match_exactly
              and obj_info.more['firmware'].major == ver.major):
            # For a tolerant match checking the firmware major
            # is good enough, as settings with the same major
            # are bound to be at least partially compatible.
            return (ver.minor,)
        return ()

    @classmethod
    def is_settings_for_this_class(cls, config):
        """Determine if a `config` file is for this controller.

        Parameters
        ----------
        config : ConfigParser
            The settings to check.

        Returns
        -------
        is_suitable : bool
            True if the settings file is for this controller.
        """
        controller_class = config.get('controller', 'controller_class',
                                      fallback=None)
        return cls.__name__ == controller_class

    def list_devices(self):  # too-complex, too-many-branches
        """List Arduino Micro VipErLEED hardware -- can be slow.

        Returns
        -------
        device_list : list
            Each element is a SettingsInfo instance containing the unique
            name of a controller and .more information as a dict.
            The .more dict contains the following keys:
                'name' : str
                    The controller name. This name may not be unique!
                'address' : str
                    The COM port on which the controller is connected.
                'firmware' : Version
                    The firmware version that is installed
                    on the controller.

        Raises
        ------
        DefaultSettingsError
            If the default settings are corrupted.
        """
        # TODO: When detecting controllers we will have to run list_devices
        # for each firmware major once and after having detected controllers
        # we have to return the SettingsInfo for each controller from the
        # best matching firmware version (on serial side). This means
        # we will discard the SettingsInfo obtained with serials running
        # on a different major if .unique_name matches.
        ports = qts.QSerialPortInfo().availablePorts()
        port_names = [p.portName() for p in ports]

        device_list = []
        threads = []
        controllers = []
        for port in port_names:
            ctrl = ViPErinoController(address=port)
            if not ctrl.has_valid_settings:
                raise DefaultSettingsError(
                    'The default settings are insufficient '
                    'to make a controller.'
                    )
            if not ctrl.serial.is_open:
                # Port is already in use
                continue
            threads.append(qtc.QThread())
            ctrl.moveToThread(threads[-1])
            controllers.append(ctrl)
        for thread in threads:
            thread.start(priority=thread.TimeCriticalPriority)
        for ctrl in controllers:
            _INVOKE(ctrl, 'get_hardware', qtc.Qt.QueuedConnection)
        if controllers:
            controllers[0].serial.port.waitForReadyRead(100)
        for ctrl in controllers:
            with ctrl.lock:
                hardware = ctrl.hardware.copy()
            serial_nr = hardware.get('serial_nr', None)
            _INVOKE(ctrl, 'disconnect_', qtc.Qt.BlockingQueuedConnection)
            if serial_nr:
                txt = f"{ctrl.name} ({ctrl.address})"
                more_info = {k: hardware[k] for k in ('firmware', )}
                more_info['address'] = ctrl.address
                more_info['name'] = ctrl.name
                device_list.append(SettingsInfo(txt, more_info))
            else:
                print("Not a ViPErLEED controller at", ctrl.address,
                      hardware, flush=True)
        for thread in threads:
            thread.quit()
        for thread in threads:
            # wait max 100 ms for each thread to quit, then force
            if not thread.wait(100):
                thread.terminate()
                thread.wait()
        return device_list

    @qtc.pyqtSlot()
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
        self.time_stamp = time.perf_counter()

    @qtc.pyqtSlot(object)
    def on_data_ready(self, data):
        """Receive and store data from the serial."""
        # We may receive two types of data: hardware&firmware
        # information (a dictionary) and actual measurements
        if isinstance(data, dict):
            # Got hardware info
            with self.lock:
                self.hardware = data
            # Now that we have info, we can check the box ID and whether
            # the quantities that we should measure can be measured.
            self._check_box_id()
            self.__check_measurements_possible()
            self.hardware_info_arrived.emit()
            return

        # Otherwise it is a list of data: [ADC0, ADC1, LM35]
        # The order is the same as stored in the config file under
        # controller/measurement_devices. The ADC channels chosen
        # via set_measurements() determine which value was measured
        for value, quantity in zip(data, self.__adc_measurement_types):
            if quantity is None:
                # Was not requested (but measured nonetheless)
                continue
            if quantity is QuantityInfo.I0:
                value = self.__convert_i0_value(value)
            self.measurements[quantity] = [value]

        # See if we should convert the thermocouple voltages
        if self.measures(QuantityInfo.TEMPERATURE):
            self.__convert_thermocouple_voltages()
        self.measurements_done()

    @qtc.pyqtSlot()
    @abc.ensure_connected
    def prepare_to_show_settings(self):
        """Prepare the controller to present settings to the user.

        Attempt a connection. Complain if the connection cannot be
        established. Otherwise initiate a request for hardware
        information, and prepare to wait for the reply. Will emit
        .ready_to_show_settings once the reply comes.

        Returns
        -------
        None.
        """
        base.safe_connect(self.hardware_info_arrived,
                          self.__almost_ready_to_show_settings,
                          type=qtc.Qt.UniqueConnection)
        self.get_hardware()

    @qtc.pyqtSlot(bool)
    def set_continuous_mode(self, continuous=True):
        """Set continuous mode.

        If continuous is true the controller will continue
        returning measurements without further instructions.

        Parameters
        ----------
        continuous : bool, optional
            Whether continuous mode should be on.
            Default is True

        Returns
        -------
        None.
        """
        super().set_continuous_mode(continuous)
        cmd = self.settings.get('available_commands', 'PC_CHANGE_MEAS_MODE',
                                fallback=None)
        if cmd is None:
            # Probably entered this after an __init__ error
            return
        mode_on = int(bool(continuous))
        self.send_message(cmd, [mode_on, 0])

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
            v_ref_dac = self.settings.getfloat('energy_calibration',
                                               'v_ref_dac')
        except (TypeError, ValueError):
            self.emit_error(QObjectSettingsErrors.INVALID_SETTINGS,
                            'energy_calibration/v_ref_dac', '')
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

        for i, egy in enumerate(energies_and_times[::2]):
            egy = self.true_energy_to_setpoint(egy)
            egy = int(round(egy * conversion_factor))
            egy = max(0, min(egy, 65535))
            energies_and_times[2*i] = egy

        # Since we may wait a potentially long time here, we
        # explicitly set the timeout keyword to be the timeout
        # from the settings plus the total waiting time.
        try:
            timeout = self.settings.getint("serial_port_settings", "timeout",
                                           fallback=0)
        except (TypeError, ValueError):
            self.emit_error(QObjectSettingsErrors.INVALID_SETTINGS,
                            'serial_port_settings/timeout', '')
            timeout = 0
        timeout = max(timeout, 0) + sum(energies_and_times[1::2])
        self.send_message(cmd, energies_and_times, timeout=timeout)
        self.time_stamp = time.perf_counter()

    def set_measurements(self, quantities):  # too-complex
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
            available_adcs = self.available_adcs()
        except NotASequenceError:
            self.emit_error(QObjectSettingsErrors.INVALID_SETTINGS,
                            'controller/measurement_devices', '')
            return
        except (KeyError, ValueError, TypeError, RuntimeError) as err:
            self.emit_error(QObjectSettingsErrors.INVALID_SETTINGS,
                            'controller', err)
            return

        n_devices = len(available_adcs)
        self.__adc_measurement_types = [None]*n_devices
        self.__adc_channels = [0]*n_devices
        if len(quantities) > n_devices:
            self.emit_error(ViPErinoErrors.TOO_MANY_MEASUREMENT_TYPES,
                            n_devices, len(quantities))
            return
        for quantity in quantities:
            _quantity = QuantityInfo.from_label(quantity)
            for i, measurable_quantities in enumerate(available_adcs):
                if _quantity not in measurable_quantities:
                    continue
                if self.__adc_measurement_types[i] is not None:
                    self.emit_error(
                        ViPErinoErrors.OVERLAPPING_MEASUREMENTS,
                        quantity, self.__adc_measurement_types[i].label
                        )
                    return
                self.__adc_channels[i] = measurable_quantities[_quantity]
                self.__adc_measurement_types[i] = _quantity
                break
            else:
                self.emit_error(ViPErinoErrors.UNSUPPORTED_QUANTITY, quantity)
                return

        # Now see if we should also measure the cold-junction
        # temperature (if possible). Notice: this adds one extra
        # quantity (which may be saved, and plot-able) that the user
        # did not ask for. Make sure we don't raise errors if this
        # quantity cannot be measured (using __added_cold_junction)
        self.__added_cold_junction = False
        measurements = self.__adc_measurement_types
        if (QuantityInfo.TEMPERATURE in measurements
            and QuantityInfo.COLD_JUNCTION not in measurements
                and measurements[-1] is not None):
            quantities = (*quantities, QuantityInfo.COLD_JUNCTION)
            self.__added_cold_junction = True

        super().set_measurements(quantities)
        self.__check_measurements_possible()

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
        with self.lock:
            hardware = self.hardware.copy()
        if not hardware:
            self.emit_error(ViPErinoErrors.HARDWARE_INFO_MISSING,
                            'set up the ADCs')
            return
        lm35_idx = list(hardware.keys()).index('lm35')
        message = [self.nr_samples, *self.__adc_channels[:lm35_idx]]
        self.send_message(cmd, message)

    @qtc.pyqtSlot(str)
    @abc.ensure_connected
    def set_serial_number(self, new_serial_nr):
        """Set a 4-characters-long serial number in the electronics."""
        new_serial_nr = new_serial_nr.upper()
        with self.lock:
            old_serial_nr = self.hardware.get('serial_nr', '')
        if not old_serial_nr:
            (*_,
            old_serial_nr) = self.settings['controller']['device_name'].split()
        if new_serial_nr == old_serial_nr:
            return

        if not re.match("^[A-Z0-9]{4,4}$", new_serial_nr):
            raise ValueError(f"Invalid serial number {new_serial_nr}")

        # We will send two requests: (1) setting serial number, and
        # (2) retrieving hardware information. This way our internal
        # info is surely up to date. For this reason we force self to
        # be busy here: this way the two messages are queued.
        self.busy =  True
        cmd = self.settings.get('available_commands', 'PC_SET_SERIAL_NR')
        self.send_message(cmd, bytes(new_serial_nr, encoding='utf-8'))
        self.get_hardware()

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

    @qtc.pyqtSlot()
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
            # super().stop accesses serial.busy_changed, which may
            # not exist if settings are somewhat funky.
            return
        stop = self.settings.get('available_commands', 'PC_STOP')
        self.send_message(stop)

    @qtc.pyqtSlot()
    def __almost_ready_to_show_settings(self):
        """Check consistency of serial number, then ready_to_show_settings."""
        # The box replied with the necessary information
        base.safe_disconnect(self.hardware_info_arrived,
                             self.__almost_ready_to_show_settings)
        self.disconnect_()

        # Check that the serial no. in self.settings and the one in
        # self.hardware match. Update the serial port if they don't.
        settings_name = self.settings.get("controller", "device_name",
                                          fallback='')
        if settings_name != self.name:
            ports = {ctrl.more['name']: ctrl.more['address']
                     for ctrl in self.list_devices()}
            correct_port = ports.get(self.name, '')
            if correct_port:
                self.address = correct_port  # Also connects
                self.disconnect_()
        self.ready_to_show_settings.emit()

    def _check_box_id(self):
        """Check if the detected box ID matches the box ID of the class."""
        with self.lock:
            hardware = self.hardware.copy()
        if not hardware:
            return
        arduino_id = hardware.get('box_id')
        if arduino_id and arduino_id != self.box_id:
            self.emit_error(ViPErinoErrors.ERROR_WRONG_BOX_ID,
                            arduino_id=arduino_id, local_id=self.box_id)

    def __check_measurements_possible(self):
        """Check that it is possible to measure stuff, given the hardware."""
        with self.lock:
            hardware = self.hardware.copy()
        if not hardware:
            return

        for adc_present, qty in zip(hardware, self.__adc_measurement_types):
            if qty is None or adc_present:
                # Did not ask for anything, or ADC can measure it
                continue
            # qty was asked, but its ADC is not available
            if (qty == QuantityInfo.COLD_JUNCTION
                    and self.__added_cold_junction):
                # We would have liked to measure it even if user did
                # not ask, but it is not possible. Do not complain.
                continue
            self.emit_error(ViPErinoErrors.REQUESTED_ADC_OFFLINE,
                            qty.label)
        if (self.measures(QuantityInfo.TEMPERATURE)
                and not self.measures(QuantityInfo.COLD_JUNCTION)):
            print("No cold-junction temperature measured. "                     # TODO: non-critical warning
                  "Thermocouple temperature will be wrong.")

    def __convert_i0_value(self, value):
        """Convert i0 measurement to microamperes."""
        # The value returned by the Arduino is correctly
        # in microamps only under the assumption that
        #   (i) in 0--10V range 1uA produces 1V at the BNC
        #  (ii) in 0--2.5V range 1mA produces 1V at the BNC
        # If, e.g., one uses a larger gain resistor (say, 2kohm)
        # 1uA will produce 2 V, which the Arduino returns as
        # the value "2". Here we correct for this:
        try:
            gain = self.settings.getfloat('conversions', 'i0_gain',
                                           fallback=1.0)
        except (TypeError, ValueError):
            gain = 1.0
        return value / gain

    def __convert_thermocouple_voltages(self):
        """Convert TC voltages in measurements to degrees centigrade."""
        if self.thermocouple is None:
            return
        tc_voltages = self.measurements[QuantityInfo.TEMPERATURE]
        cjc_temperatures = self.measurements.get(QuantityInfo.COLD_JUNCTION,
                                                 [None]*len(tc_voltages))
        self.measurements[QuantityInfo.TEMPERATURE] = [
            self.thermocouple.temperature(v, t0)
            for v, t0 in zip(tc_voltages, cjc_temperatures)
            ]

    @classmethod
    def _get_version(cls, config):
        """Get the firmware version of the hardware from `config`."""
        ver = config.get('controller', 'firmware_version', fallback='0.0')
        return base.Version(ver)

    @staticmethod
    def __is_adc(name):
        return "adc" in name or 'lm35' in name
