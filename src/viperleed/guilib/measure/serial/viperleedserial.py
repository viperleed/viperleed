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

import struct

from viperleed.guilib.measure.serial.abc import ExtraSerialErrors, SerialABC
from viperleed.guilib.measure.hardwarebase import (ViPErLEEDErrorEnum,
                                                   emit_error)


class ViPErLEEDHardwareError(ViPErLEEDErrorEnum):
    """This class contains all errors related to the Arduino."""
    ERROR_NO_ERROR = (0, "No error")
    ERROR_SERIAL_OVERFLOW = (1, "Overflow of the hardware serial.")
    ERROR_MSG_TOO_LONG = (2, "Sent message too long.")
    ERROR_MSG_SENT_INCONSISTENT = (3,
                                   "Sent message has not the length "
                                   "it is supposed to have.")
    ERROR_MSG_RCVD_INCONSISTENT = (13,
                                   "Received message has not the "
                                   "length it is supposed to have.")
    ERROR_MSG_UNKNOWN = (4, "Unknown request from PC.")
    ERROR_MSG_DATA_INVALID = (5,
                              "Request from PC contains invalid information.")
    ERROR_NEVER_CALIBRATED = (6,
                              "Cannot perform operation before ADCs have "
                              "been calibrated. Send PC_CALIBRATION.")
    ERROR_TIMEOUT = (7, "Controller timed out while waiting for something.")
    ERROR_ADC_SATURATED = (8,
                           "One of the ADC values reached saturation "
                           "and the gain cannot be decreased further.")
    ERROR_TOO_HOT = (9, "The temperature read by the LM35 is too high.")
    ERROR_HARDWARE_UNKNOWN = (10,
                              "Cannot perform operation before the hardware "
                              "configuration is known. Send PC_CONFIGURATION.")
    ERROR_RUNTIME = (255,
                     "Some function has been called from an "
                     "inappropriate state. This is to flag "
                     "possible bugs for future development.")
    ERROR_UNKNOWN_ERROR = (256,
                           "The error code received is not a known "
                           "error type. Possibly a newly introduced "
                           "error code that has not been added to the "
                           "ViPErLEEDHardwareError enum.Enum yet.")
    ERROR_ERROR_SLIPPED_THROUGH = (11,
                                   "For some reason and error identifier "
                                   "ended up in the unprocessed messages.")
    ERROR_MSG_RCVD_INVALID = (12,
                              "Received message length does not fit "
                              "the length of any data expected.")
    ERROR_TOO_MANY_COMMANDS = (14,
                               "Only one request command "
                               "can be processed at a time")
    ERROR_NO_HARDWARE_DETECTED = (15,
                                  "No ADC detected. External power "
                                  "supply may be disconnected.")
    ERROR_VERSIONS_DO_NOT_MATCH = (
        16,
        "The version of the firmware installed on the ViPErLEED "
        "hardware (v{arduino_version}) does not match the one on "
        "the PC (v{local_version}). May be incompatible."
        )


class ViPErLEEDSerial(SerialABC):
    """Class for communication with Arduino Micro ViPErLEED controller."""

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
            be changed using the .port_settings property, or with the
            set_port_settings() setter method. See the documentation
            of set_port_settings() for more details on other mandatory
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
        self.__changed_mode = False

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
        cmd_names = ('PC_CHANGE_MEAS_MODE', 'PC_SET_VOLTAGE',
                     'PC_CALIBRATION', 'PC_SET_UP_ADCS',
                     'PC_SET_VOLTAGE_ONLY', 'PC_SET_SERIAL_NR')
        # Here, for newer firmware versions, one would:
        # version = self.port_settings.getfloat('controller', 'firmware_version',
        #                                       fallback=1.0)
        # if version >= SOMETHING:
        #     cmd_names += (NEW_COMMAND_1, NEW_COMMAND_2, ...)

        _cmds = self.port_settings['available_commands']
        return [_cmds[name] for name in cmd_names]

    def decode(self, message):
        """Decodes messages received from the Arduino.

        If there are SPECIAL_BYTEs in the message, decode them and turn
        them into single bytes. Afterwards, check if message fits the
        length given in msg_length and check if the function has the
        length of any data expected. If true return bytearray of data,
        if false return empty bytearray.

        Parameters
        ----------
        message : bytes or bytearray
            The message to be decoded

        Returns
        -------
        message : bytearray
            The decoded message
        """
        special_byte = self.port_settings.getint('serial_port_settings',
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
        special_byte = self.port_settings.getint('serial_port_settings',
                                                 'SPECIAL_BYTE')

        # print(f"before special byte: {message=}")
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
                          in self.port_settings['arduino_states'].items()}
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
            emit_error(self, ExtraSerialErrors.INVALID_PORT_SETTINGS,
                       f'arduino_states with code {error_state}')
            fmt_data = {'error_name': current_error.name,
                        'state': f"state with code {error_state}",
                        'err_details': current_error[1]}
        # Emit error and clear errors.
        emit_error(self, (error_code, msg_to_format), **fmt_data)
        self.clear_errors()

    def is_decoded_message_acceptable(self, message):
        """Check whether a decoded message is ok.

        Check if message length is consistent with the length given
        in the first byte and if message is one of the three possible
        lengths sent by the Arduino (1, 2, 4).

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
            emit_error(self, ExtraSerialErrors.NO_MESSAGE_ERROR)
            return False
        # Check if message length consistent with sent length
        if msg_length != len(msg_data):
            emit_error(self,
                       ViPErLEEDHardwareError.ERROR_MSG_RCVD_INCONSISTENT)
            return False
        # Check if message length is one of the expected lengths
        if msg_length not in (1, 2, 4, 8):
            emit_error(self, ViPErLEEDHardwareError.ERROR_MSG_RCVD_INVALID)
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
        pc_error = self.port_settings.getint('available_commands', 'PC_ERROR')
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
        available_commands = self.port_settings['available_commands'].values()
        command, *data = message

        if command not in available_commands:
            # The first message must be a command
            emit_error(self, ExtraSerialErrors.UNSUPPORTED_COMMAND_ERROR,
                       command)
            return False
        if any(message in available_commands for message in data):
            # Only the fist message can be a command; the others (if
            # any) can only be data to accompany the command itself.
            emit_error(self, ViPErLEEDHardwareError.ERROR_TOO_MANY_COMMANDS)
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
        change_mode = self.port_settings.get('available_commands',
                                             'PC_CHANGE_MEAS_MODE')

        if command != change_mode:
            self.__last_request_sent = command
            self.__changed_mode = False
        else:
            self.__changed_mode = True

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
        commands = self.port_settings['available_commands']
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
                # Here the data is [mode, time] where time
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

    def process_received_messages(self):                                        # TODO: too-complex: split off parts
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

        pc_configuration = self.port_settings.get('available_commands',
                                                  'PC_CONFIGURATION')
        pc_set_voltage = self.port_settings.get('available_commands',
                                                'PC_SET_VOLTAGE')
        pc_measure_only = self.port_settings.get('available_commands',
                                                 'PC_MEASURE_ONLY')
        pc_ok = self.port_settings.get('available_commands', 'PC_OK')

        last_cmd = self.__last_request_sent

        # The following check catches data that is received
        # from setting up a connection to the micro controller.
        if (last_cmd == pc_configuration
                and len(self.unprocessed_messages) != 1):
            self.unprocessed_messages = []
            return

        for message in self.unprocessed_messages:
            # If the length of the message is 1, then it has to be a
            # PC_OK byte.
            if len(message) == 1 and message == pc_ok.encode():
                if self.__changed_mode:
                    self.__changed_mode = False
                elif last_cmd == pc_set_voltage:
                    self.about_to_trigger.emit()
                self.busy = False
            # Hardware config and measurement values are both 4 bytes long.
            # Both commands are differentiated by the __last_request_sent
            # attribute and get processed accordingly. If continuous mode
            # is on and data got returned after a new request has been
            # sent, nothing is done with it.
            elif len(message) == 4:
                if last_cmd == pc_configuration:
                    # Firmware already checked on serial side.
                    info = self.__firmware_and_hardware(message)
                    self.data_received.emit(info)
                    self.busy = False
                elif last_cmd in (pc_set_voltage, pc_measure_only):
                    self.__measurements.append(self.__bytes_to_float(message))
                    if len(self.__measurements) < 3:
                        # Not enough data yet.
                        continue
                    self.data_received.emit(self.__measurements.copy())
                    self.__measurements = []
                else:
                    # We may have read only a part of the three
                    # ADC measurements. Throw them all away. This
                    # may swallow stray 4-long messages.
                    self.__measurements = []
            elif len(message) == 8:
                if last_cmd == pc_configuration:
                    info = self.__firmware_and_hardware(message)
                    self.data_received.emit(info)
                    self.busy = False
            # If the message is 2 bytes long, then it is most likely the
            # identifier for an error and ended up here somehow
            elif len(message) == 2:
                emit_error(self,
                           ViPErLEEDHardwareError.ERROR_ERROR_SLIPPED_THROUGH)
            # If the message does not fit one of the lengths above, it
            # is no known message type.
            else:
                emit_error(self, ViPErLEEDHardwareError.ERROR_MSG_RCVD_INVALID)
        self.unprocessed_messages = []

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
            raise TypeError("bytes_to_float: invalid type "
                            f"{type(bytes_in).__name__}. Expected 'bytes',"
                            "or 'bytearray'")
        if len(bytes_in) != 4:
            raise ValueError("bytes_to_float: Invalid number of bytes "
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
            Should have length 8
            Contains firmware as bytes 0 and 1 and hardware
            configuration as bytes 2 and 3
            serial number as bytes 4 to 7

        Returns
        -------
        hardware_config : dict
            keys : {'adc_0', 'adc_1', 'lm35', 'relay', 'i0_range',
                    'aux_range', 'serial_nr', 'firmware'}
            values : bool or str
                Values are True/False for 'adc_0', 'adc_1', 'lm35',
                and 'relay', corresponding to the hardware having
                access to the devices; Values for 'i0_range' and
                'aux_range' are the strings '0 -- 2.5 V' or
                '0 -- 10 V'; a human readable (numbers and letters)
                serial number as str for 'serial_nr', and a string
                of the form '<major>.<minor>' for 'firmware'.

        Emits
        -----
        error_occurred
            If the firmware version of the controller class and
            the hardware do not match up or if the hardware did
            not detect any ADCs to take measurements with.
        """
        local_version = self.port_settings['controller']['FIRMWARE_VERSION']
        local_major, local_minor = (int(m) for m in local_version.split("."))
        major, minor, *hardware = message[:4]
        firmware_version = f"{major}.{minor}"
        if (major < local_major
                or (major == local_major and minor < local_minor)):
            emit_error(self,
                       ViPErLEEDHardwareError.ERROR_VERSIONS_DO_NOT_MATCH,
                       arduino_version=firmware_version,
                       local_version=local_version)
        # TODO: here we may want to report a (non fatal) warning in
        # case the firmware version in the hardware is newer than the
        # one of the settings.
        hardware_bits = self.port_settings['hardware_bits']
        info = {'adc_0': False,
                'adc_1': False,
                'lm35': False,
                'relay': False,
                'i0_range': '0 -- 10 V',
                'aux_range': '0 -- 10 V',
                'serial_nr': 'NO_SERIAL_NR',
                'firmware': None}

        hardware = int.from_bytes(hardware, self.byte_order)
        for key, value in hardware_bits.items():
            present_or_closed = bool(int(value) & hardware)
            key = key.lower().replace('_present', '')
            if 'closed' in key:
                if present_or_closed:
                    present_or_closed = '0 -- 2.5 V'
                else:
                    present_or_closed = '0 -- 10 V'
                key = 'i0_range' if 'i0' in key else 'aux_range'
            info[key] = present_or_closed
        if not (info['adc_0'] or info['adc_1']):
            emit_error(self, ViPErLEEDHardwareError.ERROR_NO_HARDWARE_DETECTED)
        info['serial_nr'] = ''.join([chr(v) for v in message[4:]])
        info['firmware'] = firmware_version
        return info

    def is_measure_command(self, command):
        """Returns true if the command sent is a
        command that will return a measurement."""
        pc_set_voltage = self.port_settings.get('available_commands',
                                                'PC_SET_VOLTAGE')
        pc_measure_only = self.port_settings.get('available_commands',
                                                 'PC_MEASURE_ONLY')
        if command in (pc_set_voltage, pc_measure_only):
            return True
        return False
