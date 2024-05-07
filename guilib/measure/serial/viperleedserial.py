"""Module viperleedserial of viperleed.guilib.measure.serial.

========================================
   ViPErLEED Graphical User Interface
========================================

Created: 2021-07-01
Author: Michele Riva
Author: Florian Doerr

This module contains the definitions of the ViPErLEEDSerial class and
its associated ViPErLEEDErrorEnum class ViPErLEEDHardwareError
that are used for serial communication with the Arduino Micro controller
used by ViPErLEED. The serial communication happens in a separate thread
to prevent stalling the Graphical User Interface.
"""

import re
import struct

from PyQt5 import QtCore as qtc

from viperleed.guilib.measure import hardwarebase as base
from viperleed.guilib.measure.classes.abc import QObjectABCErrors
from viperleed.guilib.measure.serial.abc import ExtraSerialErrors
from viperleed.guilib.measure.serial.abc import SerialABC


class ViPErLEEDHardwareError(base.ViPErLEEDErrorEnum):
    """This class contains all errors related to the Arduino."""

    ERROR_NO_ERROR = (0, 'No error')
    ERROR_SERIAL_OVERFLOW = (1, 'Overflow of the hardware serial.')
    ERROR_MSG_TOO_LONG = (2, 'Sent message too long.')
    ERROR_MSG_SENT_INCONSISTENT = (3,
                                   'Sent message has not the length '
                                   'it is supposed to have.')
    ERROR_MSG_RCVD_INCONSISTENT = (13,
                                   'Received message has not the '
                                   'length it is supposed to have.')
    ERROR_MSG_UNKNOWN = (4, 'Unknown request from PC.')
    ERROR_MSG_DATA_INVALID = (5,
                              'Request from PC contains invalid information.')
    ERROR_NEVER_CALIBRATED = (6,
                              'Cannot perform operation before ADCs have '
                              'been calibrated. Send PC_CALIBRATION.')
    ERROR_TIMEOUT = (7, 'Controller timed out while waiting for something.')
    ERROR_ADC_SATURATED = (8,
                           'One of the ADC values reached saturation '
                           'and the gain cannot be decreased further.')
    ERROR_TOO_HOT = (9, 'The temperature read by the LM35 is too high.')
    ERROR_HARDWARE_UNKNOWN = (10,
                              'Cannot perform operation before the hardware '
                              'configuration is known. Send PC_CONFIGURATION.')
    ERROR_RUNTIME = (255,
                     'Some function has been called from an '
                     'inappropriate state. This is to flag '
                     'possible bugs for future development.')
    ERROR_UNKNOWN_ERROR = (256,
                           'The error code received is not a known '
                           'error type. Possibly a newly introduced '
                           'error code that has not been added to the '
                           'ViPErLEEDHardwareError enum.Enum yet.')
    ERROR_ERROR_SLIPPED_THROUGH = (11,
                                   'For some reason an error identifier '
                                   'ended up in the unprocessed messages.')
    ERROR_MSG_RCVD_INVALID = (12,
                              'Received message length does not fit '
                              'the length of any data expected.')
    ERROR_TOO_MANY_COMMANDS = (14,
                               'Only one request command '
                               'can be processed at a time')
    ERROR_NO_HARDWARE_DETECTED = (15,
                                  'No ADC detected. External power '
                                  'supply may be disconnected.')
    ERROR_VERSIONS_DO_NOT_MATCH = (
        16,
        'The version of the firmware installed on the ViPErLEED '
        'hardware (v{arduino_version}) does not match the one on '
        'the PC (v{local_version}). May be incompatible.'
        )
    ADC_POWER_FAULT = (
        17,
        'Inconsistent hardware configuration from the unit: ADC#1 '
        'has power, but ADC#0 was not detected. This likely indicates a '
        'hardware fault on the board. Check that ADC#0 has power.'
        )
    ERROR_MSG_SENT_TOO_LONG = (
        254,
        'The controller tried to send a message that was longer '
        'than MSG_SPECIAL_BYTE. This means the firmware has a bug.'
        )


class ViPErLEEDSerial(SerialABC):
    """Class for communication with Arduino Micro ViPErLEED controller."""

    debug_info_arrived = qtc.pyqtSignal(str)

    _mandatory_settings = [*SerialABC._mandatory_settings,
                           ('hardware_bits',), ('available_commands',),
                           ('arduino_states',), ('error_bytes',),
                           ('controller', 'FIRMWARE_VERSION'),
                           ('serial_port_settings', 'special_byte'),]

    def __init__(self, settings, port_name='', **kwargs):
        """Initialize serial worker object.

        Parameters
        ----------
        settings : ConfigParser
            The settings used to initialize the serial port. For the
            serial port itself one only needs a 'serial_port_settings'
            section. In practice, it is likely safer to pass the whole
            ConfigParser object used for the controller class. It can
            be changed using the self.settings property, or with the
            set_settings() setter method. See the documentation
            of set_settings() for more details on other mandatory
            content of the settings argument.
        port_name : str or QSerialPortInfo, optional
            The name (or info) of the serial port. If not given,
            it must be set via the .port_name attribute before any
            communication can be established. Default is an empty
            string.

        Raises
        ------
        TypeError
            If no settings are given.
        """
        self.__last_request_sent = ''
        self.__measurements = []
        self.__may_receive_stray_data = False
        self.__is_waiting_for_debug_msg = False
        self.__should_emit_debug_msg = False

        super().__init__(settings, port_name=port_name, **kwargs)

    @property
    def commands_requiring_data(self):
        """Return a list of the commands that require extra data.

        Commands that need data are those that require multiple
        messages to be sent: one message with the command itself,
        one message with the data.

        Returns
        -------
        commands : Sequence of strings
            Each element is the string representation of the byte
            that will be sent over the serial.
        """
        cmd_names = ['PC_CHANGE_MEAS_MODE', 'PC_SET_VOLTAGE',
                     'PC_CALIBRATION', 'PC_SET_UP_ADCS',
                     'PC_SET_VOLTAGE_ONLY', 'PC_SET_SERIAL_NR']
        if self.firmware_version >= "0.7":
            cmd_names.append("PC_DEBUG")

        _cmds = self.settings['available_commands']
        return [_cmds[name] for name in cmd_names]

    @property
    def firmware_version(self):
        """Return the Version of the firmware stored in settings."""
        return base.Version(
            self.settings.get("controller", "firmware_version")
            )

    def clear_errors(self):
        """Clear all errors and reset to the same state as a fresh instance."""
        super().clear_errors()
        self.__last_request_sent = ''
        self.__measurements = []
        self.__may_receive_stray_data = False
        self.__is_waiting_for_debug_msg = False
        self.__should_emit_debug_msg = False

    def decode(self, message):
        """Return a decoded version of message.

        Parameters
        ----------
        message : bytes or bytearray
            The message to be decoded. SPECIAL_BYTEs and those that
            follow them are decoded into single bytes.

        Returns
        -------
        decoded_message : bytearray
            The decoded message. decoded_message is an empty
            bytearray if not is_decoded_message_acceptable().
        """
        special_byte = self.settings.getint('serial_port_settings',
                                            'SPECIAL_BYTE')
        msg_length, *msg_data = message
        # Iterate through message and add special bytes up
        for index, value in enumerate(msg_data):
            if value == special_byte:
                msg_data[index] += msg_data.pop(index + 1)
        msg_to_check = bytearray((msg_length, *msg_data))
        if not self.is_decoded_message_acceptable(msg_to_check):
            return bytearray()
        return bytearray(msg_data)

    def encode(self, message):
        """Encode messages sent to the Arduino.

        If there are SPECIAL_BYTEs in the message encode them and make
        them two bytes. Afterwards, insert the length of the decoded
        message at the beginning of the bytearray.

        Parameters
        ----------
        message : bytearray
            The message to be encoded

        Returns
        -------
        message : bytearray
            The encoded message

        Raises
        ------
        ValueError
            If, after encoding, the message is longer than 62 bytes,
            as this would, together with start- and end-marker chars
            fill up the whole serial buffer of the Arduino board.
        """
        special_byte = self.settings.getint('serial_port_settings',
                                            'SPECIAL_BYTE')

        payload_length = len(message)
        # Iterate through message and split values that are larger
        # or equals special_byte into two different bytes. The second
        # byte gets inserted right after the first one.
        # Enumerate returns bytes as integers
        for i, value in enumerate(message):
            if value >= special_byte:
                message[i] = special_byte
                message.insert(i + 1, value - special_byte)
        # Insert the message length at the beginning of the message
        message.insert(0, payload_length)
        # Get the message length and check if the Arduino can store the
        # whole message with its 64 byte sized buffer. The calculation
        # factors in the START and END marker.
        if len(message) > 62:
            raise ValueError(f"Length of message {message} too long "
                             "to be received by the Arduino.")
        return message

    def identify_error(self, messages_since_error):
        """Identify the error sent back by the Arduino.

        Error data will be returned right after sending the PC_ERROR
        byte. Therefore messages_since_error[0] will contain the PC_ERROR
        byte and messages_since_error[1] will contain both the error
        state and the error type in this order.

        Parameters
        ----------
        messages_since_error : list
            All the full messages that came after the last error
            message was received, including the error message itself.
            messages_since_error[0] is the error message. Each message
            is a bytes object, containing the already decoded message.
            All messages_since_error will be discarded after the error
            has been identified.

        Emits
        -----
        error_occurred(ViPErLEEDErrorEnum)
            After the error has been identified
        """
        # Check if the error data has arrived yet.
        if len(messages_since_error) < 2:
            return
        # Reading error state and error code from the returned data
        error_state, error_code = messages_since_error[1]

        # Flipping the arduino_states section in order to access
        # the name of the state (key)
        arduino_states = {int(code): state.upper()
                          for state, code
                          in self.settings['arduino_states'].items()}
        # Preparing message which will be formatted and emitted.
        msg_to_format = ("ViPErLEED hardware {error_name} occurred while in"
                         " {state}. Reason: {err_details}")
        # Try to get the error enum.Enum via the error code. If this
        # fails the ERROR_UNKNOWN_ERROR will be taken.
        try:
            current_error = ViPErLEEDHardwareError.from_code(error_code)
        except AttributeError:
            current_error = ViPErLEEDHardwareError.ERROR_UNKNOWN_ERROR
        # Set values and format the string which will be emitted.
        try:
            fmt_data = {'error_name': current_error.name,
                        'state': arduino_states[error_state],
                        'err_details': current_error[1]}
        except KeyError:
            # error_state is not present in the config file.
            base.emit_error(self, QObjectABCErrors.INVALID_SETTINGS,
                            type(self).__name__,
                            f'arduino_states with code {error_state}')
            fmt_data = {'error_name': current_error.name,
                        'state': f"state with code {error_state}",
                        'err_details': current_error[1]}
        # Emit error and clear errors.
        base.emit_error(self, (error_code, msg_to_format), **fmt_data)
        self.clear_errors()

    def is_decoded_message_acceptable(self, message):
        """Check whether a decoded message is ok.

        Check if message length is consistent with the length given
        in the first byte and if message is one of the five possible
        lengths sent by the Arduino (1, 2, 4, 8, 9).

        Parameters
        ----------
        message : bytearray
            The message to be checked

        Returns
        -------
        acceptable : bool
            True if message is acceptable

        Emits
        -----
        error_occurred
            If the message length is not consistent with the length
            sent in the message or if the message length is not
            one of the acceptable lengths.
        """
        msg_length, *msg_data = message
        # Check if message length greater than zero
        if not msg_length:
            base.emit_error(self, ExtraSerialErrors.NO_MESSAGE_ERROR)
            return False
        # Check if message length consistent with sent length
        if msg_length != len(msg_data):
            base.emit_error(self,
                            ViPErLEEDHardwareError.ERROR_MSG_RCVD_INCONSISTENT)
            return False

        if self.firmware_version >= "0.7":
            _debug = self.settings.getint("available_commands", "pc_debug")
            if self.__is_waiting_for_debug_msg:
                self.__is_waiting_for_debug_msg = False
                return True
            if msg_length == 1 and msg_data[0] == _debug:
               self.__is_waiting_for_debug_msg = True
               return True

        # Check if message length is one of the expected lengths
        if msg_length not in (1, 2, 4, 8, 9):
            base.emit_error(self,
                            ViPErLEEDHardwareError.ERROR_MSG_RCVD_INVALID)
            return False
        return True

    def is_error_message(self, message):
        """Check if message is an error message.

        Check if the message length is equal 1 and has the value of
        PC_ERROR. If true return that the message is an error
        message.

        Parameters
        ----------
        message : bytes or bytearray
            A decoded message received through the serial line.

        Returns
        -------
        bool
            Whether message is an error message.
        """
        pc_error = self.settings.getint('available_commands', 'PC_ERROR')
        pc_error = pc_error.to_bytes(1, self.byte_order)
        return len(message) == 1 and message == pc_error

    def is_message_supported(self, message):
        """Check if message contains one allowed command.

        Check whether the message is a supported command or not and
        keep track of which request was sent last to the Arduino.
        If a command is unsupported, emit an error_occurred signal
        with ExtraSerialErrors.UNSUPPORTED_COMMAND_ERROR. If there
        are more than two commands in one set of messages emit an
        error as the current structure does not allow this. If
        a change measurement mode command has been sent, update the
        continuous measurement setting on the PC side. After clearing
        the message, remember the sent request and set the serial to
        busy.

        Parameters
        ----------
        message : Sequence
            The same message given to self.send_message(), right
            before encoding and sending to the device. message[0]
            is a string (the command to send). The other elements
            contain the data associated with the command.

        Returns
        -------
        message_acceptable : bool
            True if message is acceptable. Only acceptable
            messages will then be sent.

        Emits
        -----
        error_occurred
            If the first message is not a known command or if
            the message contains more than one command.

        Raises
        ------
        RuntimeError
            If the command to be sent requires data, but there is
            no data in message.
        RuntimeError
            If there is more than one block of data in message.
        """
        # Remember!: message[0] has to be of type Str
        # Check if message is valid data. Since the length of message
        # can vary, the only reasonable check is if single bytes are
        # indeed a command specified in the config.
        available_commands = self.settings['available_commands'].values()
        command, *data = message

        if command not in available_commands:
            # The first message must be a command
            base.emit_error(self, ExtraSerialErrors.UNSUPPORTED_COMMAND_ERROR,
                            command)
            return False
        if any(message in available_commands for message in data):
            # Only the fist message can be a command; the others (if
            # any) can only be data to accompany the command itself.
            base.emit_error(self,
                            ViPErLEEDHardwareError.ERROR_TOO_MANY_COMMANDS)
            return False

        _need_data = self.commands_requiring_data
        if command in _need_data and not data:
            # Too little data available
            # This is meant as a safeguard for future code changes
            raise RuntimeError(
                f"{self.__class__.__name__}.is_message_supported: "
                "No data message for one of the Arduino commands "
                "that requires data. Check implementation!"
                )
        if command not in _need_data and data:
            # Too much data available
            # This is meant as a safeguard for future code changes
            raise RuntimeError(
                f"{self.__class__.__name__}.is_message_supported: "
                "Got a data message for one of the Arduino commands "
                "that does not require data. Check implementation!"
                )
        if len(data) > 1:
            # Too much data available: we never send more than
            # two messages (one command, one data).
            # This should never happen, as long as we pass all the
            # data for a command as a single entity. It is meant as
            # a safeguard for future changes of the code.
            raise RuntimeError(
                f"{self.__class__.__name__}.is_message_supported: "
                "At most one data message should be given for each "
                "command to the Arduino. Wrap data for one command into "
                "a single Sequence."
                )
        change_mode = self.settings.get('available_commands',
                                        'PC_CHANGE_MEAS_MODE')

        if command != change_mode:
            self.__last_request_sent = command

        # We use this boolean while continuous mode is still active
        # and data may return while we are ordering the controller
        # to change its measurement mode.
        self.__may_receive_stray_data = command == change_mode

        return True

    def message_requires_response(self, *messages):
        """Return whether the messages to be sent require a response.

        Currently all commands return a PC_OK or other data.

        Parameters
        ----------
        *messages : tuple
            Same arguments passed to send_message

        Returns
        -------
        bool
        """
        return True

    def prepare_message_for_encoding(self, message, *other_messages):
        """Prepare a message to be encoded.

        This method is guaranteed to be called exactly once on each
        call to send_message(message, *other_messages), right before
        each of the messages are (separately) encoded. Since .encode()
        expects each message to be a bytes or bytearray, this method
        can be called to turn messages into bytearrays.

        Parameters
        ----------
        message : str
            The command to be sent to the ViPErLEED Arduino
        *other_messages : int, bytes, bytearray, Sequence
            At most one *other_messages should be given for
            the ViPErLEED Arduino. When a Sequence, each element
            must be an int. This is the data needed by the command
            in message.

        Returns
        -------
        messages_to_return : tuple of bytearrays
            Consisting of the command and the data necessary to
            execute it.
        """
        # Convert 'message' (i.e., the command) into a bytearray if
        # necessary.  'other_messages' is the data for the command.
        if not other_messages:
            return (bytearray(message.encode()),)
        messages_to_return = [bytearray(message.encode()),]
        # Now process the data depending on the command
        commands = self.settings['available_commands']
        for data in other_messages:
            # In the following, we treat the special cases of:
            # "some of the pieces of information passed to send_message
            # is an int that has to be turned into bytes with length >= 2"
            if message == commands['PC_SET_UP_ADCS']:
                # data is a 3-element Sequence, with
                # [n_averaging_points, adc0_channel, adc1_channel],
                # where n_averaging_points has to be 2 bytes
                data[:1] = data[0].to_bytes(2, self.byte_order)
            elif message in (commands['PC_SET_VOLTAGE'],
                             commands['PC_SET_VOLTAGE_ONLY']):
                # data is a 2N-long sequence of the format
                # voltage, waiting, voltage, waiting, ....
                # where each entry has to be turned into 2
                # bytes.
                # Run over a copy to avoid infinite loop.
                for i, elem in enumerate(data.copy()):
                    data[2*i:2*i+1] = elem.to_bytes(2, self.byte_order)
            elif message == commands['PC_CHANGE_MEAS_MODE']:
                # Here the data is [mode, time] where time                      TODO: time is not used at this point.
                # has to be turned into 2 bytes. We do not convert
                # mode as it is hardcoded to 0 (off) or 1 (on).
                data[1:] = data[1].to_bytes(2, self.byte_order)
            # If more commands are introduced in newer firmware
            # versions, one should use the firmware version of the
            # settings here to preserve backwards compatibility.
            messages_to_return.append(bytearray(data))
        # Note that PC_CALIBRATION also requires data, but
        # all the elements in data are 1-byte-long integers.
        return messages_to_return

    def process_received_messages(self):                                        # TODO: too-complex, too-many-branches, confusing-consecutive-elif: split off parts
        """Convert received data into human understandable information.

        Process data according to its type and emit data that will be
        processed further. Set serial busy to false if either a PC_OK
        arrives of if the hardware configuration gets returned by the
        Arduino.

        Emits
        -----
        data_received
            Any time the message received contains data that are
            worth processing.
        error_occurred
            If the message is an error message or if it could not
            be identified.
        """
        if not self.unprocessed_messages:
            # No messages came yet
            return

        pc_configuration = self.settings.get('available_commands',
                                             'PC_CONFIGURATION')

        # The following check catches data that is received
        # from setting up a connection to the microcontroller.
        if (self.__last_request_sent == pc_configuration                       # TODO: This may swallow a configuration response that arrives before we enter this check.
                and len(self.unprocessed_messages) != 1):
            self.unprocessed_messages = []
            return

        for message in self.unprocessed_messages:
            self._process_single_message(message)
        self.unprocessed_messages = []

    def _process_single_message(self, message):
        """Process a single one of the arrived messages."""
        pc_configuration = self.settings.get('available_commands',
                                             'PC_CONFIGURATION')
        pc_ok = self.settings.get('available_commands', 'PC_OK')
        pc_autogain = self.settings.get('available_commands', 'PC_AUTOGAIN')
        last_cmd = self.__last_request_sent

        pc_debug = None
        if self.firmware_version >= "0.7":
            pc_debug = self.settings.getint('available_commands', 'PC_DEBUG')
        if self.__should_emit_debug_msg:
            self._process_debug_message(message)
        elif len(message) == 1 and message == pc_ok.encode():
            # If the length of the message is 1, then it has to be a
            # PC_OK byte.
            self._process_ok_message()
        elif len(message) == 2 and last_cmd == pc_autogain:
            # If the message is 2 bytes long and the last command was
            # PC_AUTOGAIN, then the controller is done searching for
            # the optimal gain values and no longer busy.
            self.busy = False
        elif len(message) in (4, 8, 9) and last_cmd == pc_configuration:
            # Firmware versions <0.6 returned 4-byte-, <0.9 8-byte-,
            # more recent ones 9-byte-long hardware configurations.
            self._process_hardware_information(message)
        elif len(message) == 4:
            self._process_numerical_data(message)
        elif len(message) == 2:
            # If the message is 2 bytes long, then it is most likely
            # the identifier for an error and ended up here somehow
            base.emit_error(
                self, ViPErLEEDHardwareError.ERROR_ERROR_SLIPPED_THROUGH
                )
        elif len(message) == 1 and message[0] == pc_debug:
            # If the message we received is a PC_DEBUG, then
            # the next message will be a debug message.
            self.__should_emit_debug_msg = True
        else:
            # If the message does not fit one of the lengths above, it
            # is no known message type.
            base.emit_error(
                self, ViPErLEEDHardwareError.ERROR_MSG_RCVD_INVALID
                )

    def _process_debug_message(self, message):
        """Emit debug message."""
        self.__should_emit_debug_msg = False
        # The debug_info_arrived signal can be used to display
        # debug messages that have been sent by the hardware.
        self.debug_info_arrived.emit(message.decode('utf-8'))

    def _process_ok_message(self):
        """Process PC_OK message."""
        pc_set_voltage = self.settings.get('available_commands',
                                           'PC_SET_VOLTAGE')
        pc_autogain = self.settings.get('available_commands', 'PC_AUTOGAIN')
        # We are implicitly using __may_receive_stray_data in the elif
        # check if the last command was pc_set_voltage.
        if self.__may_receive_stray_data:
            self.__may_receive_stray_data = False
        elif self.__last_request_sent == pc_set_voltage:
            self.about_to_trigger.emit()
        elif (self.firmware_version >= '0.9'
              and self.__last_request_sent == pc_autogain):
            # Autogain just started. Will return gains when done.
            return
        self.busy = False

    def _process_numerical_data(self, message):
        """Process one ADC measurement. Emit data when all arrived."""
        pc_set_voltage = self.settings.get('available_commands',
                                           'PC_SET_VOLTAGE')
        pc_measure_only = self.settings.get('available_commands',
                                            'PC_MEASURE_ONLY')
        # If continuous mode is on and data got returned after a
        # new request has been sent, nothing is done with it.
        if self.__last_request_sent in (pc_set_voltage, pc_measure_only):
            self.__measurements.append(self.__bytes_to_float(message))
            if len(self.__measurements) < 3:
                # Not enough data yet.
                return
            self.data_received.emit(self.__measurements.copy())

        # Always clear the data. This takes care of deleting partial
        # (<3) readings from ADCs and other stray 4-long messages.
        self.__measurements = []

    def _process_hardware_information(self, message):
        """Process hardware configuration."""
        info = self.__firmware_and_hardware(message)
        self.data_received.emit(info)
        self.busy = False

    def __bytes_to_float(self, bytes_in):
        """Convert four bytes to float.

        Parameters
        ----------
        bytes_in : bytes or bytearray
            Should have length 4.

        Returns
        -------
        float

        Raises
        ------
        TypeError
            If bytes_in is not bytes ot bytearray
        ValueError
            If bytes_in does not contain exactly 4 bytes.
        """
        if (not hasattr(bytes_in, '__len__')
                or not isinstance(bytes_in, (bytes, bytearray))):
            raise TypeError("__bytes_to_float: invalid type "
                            f"{type(bytes_in).__name__}. Expected 'bytes',"
                            "or 'bytearray'")
        if len(bytes_in) != 4:
            raise ValueError("__bytes_to_float: Invalid number of bytes "
                             f"({len(bytes_in)}). Expected 4.")
        float_fmt = '>f' if self.byte_order == 'big' else '<f'
        return struct.unpack(float_fmt, bytes_in)[0]

    def __firmware_and_hardware(self, message):
        """Return Arduino firmware version and hardware configuration.

        Decode the message containing the configuration and return the
        hardware. If the firmware_version does not fit the version
        specified in the configuration files emit an error. If no ADC
        has been detected emit an error. The hardware configuration
        contains data about the internal setting of the controller
        and the serial number stored in its EEPROM.

        Parameters
        ----------
        message : bytearray
            Should have length 8 or 9
            Contains box ID as byte 0, firmware as bytes 1 and 2,
            hardware configuration as bytes 3 and 4, and serial
            number as bytes 5 to 8. The box ID is not present for
            versions earlier than 0.9. Then the message is 8-bytes
            long. All other bytes are shifted back by one.

        Returns
        -------
        hardware_config : dict
            keys : {'adc_0', 'adc_1', 'lm35', 'relay', 'i0_range',
                    'aux_range', 'serial_nr', 'firmware'}
            values : bool or str or int or None
                Values are True/False for 'adc_0', 'adc_1', 'lm35',
                and 'relay', corresponding to the hardware having
                access to the devices; Values for 'i0_range' and
                'aux_range' are the strings '0 -- 2.5 V' or
                '0 -- 10 V'; a human-readable (numbers and letters)
                serial number as str for 'serial_nr', and a Version
                of the form '<major>.<minor>' for 'firmware'.
                'box_id' is an int. box_id may be None if the
                installed firmware version does not have the
                box ID yet.

        Emits
        -----
        error_occurred
            If the firmware version of the controller class and
            the hardware do not match up or if the hardware did
            not detect any ADCs to take measurements with.
        """
        info = {'adc_0': False,
                'adc_1': False,
                'lm35': False,
                'relay': False,
                'i0_range': '0 \u2013 10 V',
                'aux_range': '0 \u2013 10 V',
                'serial_nr': 'NO_SERIAL_NR',
                'firmware': None,
                'box_id': None}

        # TODO: from version 1.0 onwards we do not want to pop the
        # first byte of the message. We only do this for now to
        # ensure backwards compatibility.
        if len(message) == 9:
            info['box_id'] = message.pop(0)

        local_version = self.firmware_version
        major, minor, *hardware = message[:4]
        firmware_version = base.Version(major, minor)
        if firmware_version < local_version:
            base.emit_error(self,
                            ViPErLEEDHardwareError.ERROR_VERSIONS_DO_NOT_MATCH,
                            arduino_version=firmware_version,
                            local_version=local_version)
        # TODO: here we may want to report a (non fatal) warning in
        # case the firmware version in the hardware is newer than the
        # one of the settings.

        hardware_bits = self.settings['hardware_bits']
        hardware = int.from_bytes(hardware, self.byte_order)
        for key, value in hardware_bits.items():
            present_or_closed = bool(int(value) & hardware)
            key = key.lower().replace('_present', '')
            if 'closed' in key:
                if present_or_closed:
                    present_or_closed = '0 \u2013 2.5 V'
                else:
                    present_or_closed = '0 \u2013 10 V'
                key = 'i0_range' if 'i0' in key else 'aux_range'
            info[key] = present_or_closed
        if not (info['adc_0'] or info['adc_1']):
            base.emit_error(self,
                            ViPErLEEDHardwareError.ERROR_NO_HARDWARE_DETECTED)
        if info['adc_1'] and not info['adc_0']:
            # It cannot be that ADC1 is active while zero isn't.
            # Probably a hardware fault.
            base.emit_error(self, ViPErLEEDHardwareError.ADC_POWER_FAULT)

        try:
            serial_nr = message[4:].decode('utf-8')
        except UnicodeDecodeError:
            serial_nr = ''
        if re.match("^[A-Z0-9]{4,4}$", serial_nr):
            info['serial_nr'] = serial_nr
        info['firmware'] = firmware_version
        return info

