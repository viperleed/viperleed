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
# Python standard modules
import struct
from collections.abc import Sequence

# ViPErLEED modules
from viperleed.guilib.measure.serial.abc import ExtraSerialErrors, SerialABC
from viperleed.guilib.measure.hardwarebase import (
    ViPErLEEDErrorEnum,
    config_has_sections_and_options,
    emit_error
    )


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
    ERROR_VERSIONS_DO_NOT_MATCH = (16,
                                   "The version of the firmware installed on "
                                   "the ViPErLEED hardware (v{arduino_version})"
                                   " does not match the one on the PC "
                                   "(v{local_version}). May be incompatible.")


class ViPErLEEDSerial(SerialABC):
    """
    Class for serial communication with
    the Arduino Micro ViPErLEED controller.
    """

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
            The name (or info) of the serial port. If not given, it
            must be set via the .port_name attribute, or the setter
            set_port_name(), before any communication can be
            established. Default is an empty string.

        Raises
        ------
        TypeError
            If no settings are given.
        """

        self.__last_request_sent = ''
        self.__measurements = []
        self.__is_continuous_mode = False

        self._mandatory_settings.extend((
            ('hardware_bits',),
            ('available_commands',),
            ('arduino_states',),
            ('error_bytes',),
            ('controller', 'FIRMWARE_VERSION')
            ))

        super().__init__(settings, port_name=port_name, **kwargs)

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
                msg_data[index] += msg_data[index+1]
                msg_data.pop(index + 1)
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

        Returns
        -------
        None.

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
        arduino_states = {int(code): state
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
        fmt_data = {'error_name': current_error.name,
                    'state': arduino_states[error_state],
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
        if msg_length == 0:
            emit_error(self, ExtraSerialErrors.NO_MESSAGE_ERROR)
            return False
        # Check if message length consistent with sent length
        if msg_length != len(msg_data):
            emit_error(self,
                       ViPErLEEDHardwareError.ERROR_MSG_RCVD_INCONSISTENT)
            return False
        # Check if message length is one of the expected lengths
        if msg_length not in (1, 2, 4):
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

        return len(message) == 1 and message == pc_error

    def is_message_supported(self, messages):
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
        messages : str, object
            The same message given to self.send_message(), right
            before encoding and sending to the device

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
        """
        # Remember!: messages[0] has to be of type Str
        # Check if message is valid data. Since the length of messages
        # can vary, the only reasonable check is if single bytes are
        # indeed a command specified in the config.
        available_commands = self.port_settings['available_commands'].values()
        command, *data = messages

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

        need_data = [self.port_settings['available_commands'][value]
                        for value in ('PC_CHANGE_MEAS_MODE', 'PC_SET_VOLTAGE',
                                      'PC_CALIBRATION', 'PC_SET_UP_ADCS')]
        if command in need_data:
            if not data:
                # Too little data available
                # This is meant as a safeguard for future code changes
                raise RuntimeError(
                    "{self.__class.__name__}.is_message_supported: "
                    "No data message for one of the Arduino commands "
                    "that requires data. Check implementation!"
                    )
        if len(data) > 1:
            # Too much data available: we never send more than
            # two messages (one command, one data).
            # This should never happen, as long as we pass all the
            # data for a command as a single entity. It is meant as
            # a safeguard for future changes of the code.
            raise RuntimeError(
                "{self.__class.__name__}.is_message_supported: "
                "At most one data message should be given for each "
                "command to the Arduino. Wrap data for one command into "
                "a single Sequence."
                )
        change_mode = self.port_settings.get('available_commands',
                                             'PC_CHANGE_MEAS_MODE')
        if command == change_mode:
            # Continuous measurement is the only thing we have to check
            # as it will potentially spam us with a lot of measurements
            # so we should react appropriately
            # print(f"{data[0][0]=} and {data[0]=}")
            self.__is_continuous_mode = bool(data[0][0])

        self.__last_request_sent = command

        return True

    def message_requires_response(self, command, *__args):
        """Return whether the messages to be sent require a response.

        Currently all commands return a PC_OK or other data.

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
        message = message.encode()
        if not other_messages:
            return (bytearray(message),)
        messages_to_return = []
        messages_to_return.append(bytearray(message))
        # Now process the data depending on the command
        commands = self.port_settings['available_commands']
        for data in other_messages:
            # In the following, we treat the special cases of:
            # "some of the pieces of information passed to send_message
            # is an int that has to be turned into bytes with length >= 2"
            if message == commands['PC_SET_UP_ADCS'].encode():
                # data is a 3-element Sequence, with
                # [n_averaging_points, adc0_channel, adc1_channel],
                # where n_averaging_points has to be 2 bytes
                data[:1] = data[0].to_bytes(2, self.byte_order)
            elif message == commands['PC_SET_VOLTAGE'].encode():
                # data is a 2N-long sequence of the format
                # voltage, waiting, voltage, waiting, ....
                # where each entry has to be turned into 2
                # bytes.
                # Run over a copy to avoid infinite loop.
                # print(f"Before byte splitting: {data=}")
                for i, elem in enumerate(data.copy()):
                    data[2*i:2*i+1] = elem.to_bytes(2, self.byte_order)
                # print(f"After byte splitting: {data=}")
            elif message == commands['PC_CHANGE_MEAS_MODE'].encode():
                # Here the data is [mode, time] where time
                # has to be turned into 2 bytes.
                data[1:] = data[-1].to_bytes(2, self.byte_order)
            messages_to_return.append(bytearray(data))
        # Note that PC_CALIBRATION also requires data, but
        # all the elements in data are 1-byte-long integers.
        return messages_to_return

    def process_received_messages(self):
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
        # The following check catches data that is received
        # from setting up a connection to the micro controller.
        if self.__last_request_sent == pc_configuration:
            if len(self.unprocessed_messages) != 1:
                self.unprocessed_messages = []
                return
        # print(self.port_name, self.unprocessed_messages)
        for message in self.unprocessed_messages:
            # If the length of the message is 1, then it has to be a
            # PC_OK byte.
            if len(message) == 1:
                if message == pc_ok.encode():
                    if self.__last_request_sent == pc_set_voltage:
                        self.about_to_trigger.emit()
                    self.busy = False
            # Hardware config and measurement values are both 4 bytes long.
            # Both commands are differentiated by the __last_request_sent
            # attribute and get processed accordingly. If continuous mode
            # is on and data got returned after a new request has been
            # sent, nothing is done with it.
            elif len(message) == 4:
                if self.__last_request_sent == pc_configuration:
                    firmware, hardware = self.__firmware_and_hardware(message)
                    # self.data_received.emit((firmware, hardware))
                    # Firmware already checked on serial side.
                    self.data_received.emit(hardware)
                    self.busy = False
                elif self.__last_request_sent in (pc_set_voltage,
                                                  pc_measure_only):
                    self.__measurements.append(self.__bytes_to_float(message))
                    if len(self.__measurements) < 3:
                        # Not enough data yet
                        continue
                    self.data_received.emit(self.__measurements)
                    self.__measurements = []
                elif self.__is_continuous_mode:
                    pass
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
        return

    def __bytes_to_float(self, bytes_in):
        """Convert four bytes to float.

        Parameters
        ----------
        bytes_in : bytes or bytearray
            Should have length 4.

        Returns
        -------
        float
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
        has been detected emit an error.

        Parameters
        ----------
        message : bytearray
            Should have length 4
            Contains firmware as bytes 0 and 1 and hardware
            configuration as bytes 2 and 3

        Returns
        -------
        firmware_version : str
            Firmware version as "<major>.<minor>"
        hardware_config : dict
            keys : {'adc_0', 'adc_1', 'lm35', 'relay',
                    'i0_range', 'aux_range'}
            values : bool or str
                Values are True/False for 'adc_0', 'adc_1',
                'lm35', and 'relay', corresponding to the
                hardware having access to the devices; Values
                for 'i0_range' and 'aux_range' are the strings
                '0 -- 2.5 V' or '0 -- 10 V'.

        Emits
        -----
        error_occurred
            If the firmware version of the controller class and
            the hardware do not match up or if the hardware did
            not detect any ADCs to take measurements with.
        """
        local_version = self.port_settings['controller']['FIRMWARE_VERSION']
        major, minor, *hardware = message
        firmware_version = f"{major}.{minor}"
        if firmware_version != local_version:
            emit_error(self,
                       ViPErLEEDHardwareError.ERROR_VERSIONS_DO_NOT_MATCH,
                       arduino_version=firmware_version,
                       local_version=local_version)
        hardware_bits = self.port_settings['hardware_bits']
        hardware_config = {'adc_0': False,
                           'adc_1': False,
                           'lm35': False,
                           'relay': False,
                           'i0_range': '0 -- 10 V',
                           'aux_range': '0 -- 10 V'}

        hardware = int.from_bytes(hardware, self.byte_order)
        for key, value in hardware_bits.items():
            if key.startswith('#'):  # line is a comment
                continue
            present_or_closed = bool(int(value) & hardware)
            key = key.lower().replace('_present', '')
            if 'closed' in key:
                if present_or_closed:
                    present_or_closed = '0 -- 2.5 V'
                else:
                    present_or_closed = '0 -- 10 V'
                key = 'i0_range' if 'i0' in key else 'aux_range'
            hardware_config[key] = present_or_closed
        if not (hardware_config['adc_0'] or hardware_config['adc_1']):
            emit_error(self, ViPErLEEDHardwareError.ERROR_NO_HARDWARE_DETECTED)
        return firmware_version, hardware_config
