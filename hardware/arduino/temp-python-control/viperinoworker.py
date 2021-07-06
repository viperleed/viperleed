"""ViPErino Serialwoker

========================================
   ViPErLEED Graphical User Interface
========================================

Created: 2021-07-01
Author: Michele Riva
Author: Florian Doerr

This module contains the definition of the ViPErinoSerialWorker class
that is used for serial communication with the Arduino Micro controller
used by ViPErLEED. The serial communication happens in a separate thread
to prevent stalling the Graphical User Interface.
"""
# Python standard modules
import struct
from configparser import ConfigParser

# ViPErLEED modules
from serialworker import ExtraSerialErrors, SerialWorkerABC, ViPErLEEDErrorEnum


# TODO: update send_message doc string
# .super() take full doc string only change acceptable data types
class ViPErLEEDHardwareError(ViPErLEEDErrorEnum):
    """This class contains all errors related to the Arduino."""
    ERROR_NO_ERROR = (0, "No error")
    ERROR_SERIAL_OVERFLOW = (1, "Hardware overflow of Arduino serial.")
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
                              "The ADCs have never been "
                              "calibrated since bootup")
    ERROR_TIMEOUT = (7, "Arduino timed out while waiting for something")
    ERROR_ADC_SATURATED = (8,
                           "One of the ADC values reached saturation "
                           "and the gain cannot be decreased further.")
    ERROR_TOO_HOT = (9, "The temperature read by the LM35 is too high.")
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
                              "Received message length does not fit the "
                              "length of any data expected.")
    ERROR_TOO_MANY_COMMANDS = (14,
                               "Only one request command can be processed "
                               "at a time")


class ViPErinoSerialWorker(SerialWorkerABC):
    """
    Class for serial communication with
    the Arduino Micro ViPErLEED controller.
    """
    def __init__(self, settings=None, port_name=''):
        """Initialize serial worker object."""
        self.__last_request_sent = ''
        self.__measurements = []
        self.__is_continuous_mode = False
        super().__init__(port_name=port_name, settings=settings)

    def decode(self, message):
        """Decodes message received from the Arduino.

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
        message : object
            The message to be encoded

        Returns
        -------
        message : bytearray
            The encoded message
        """
        special_byte_int = self.port_settings.getint('serial_port_settings',
                                                     'SPECIAL_BYTE')
        special_byte = special_byte_int.to_bytes(1, self.byte_order)
        # Check for data type and convert message into a bytearray if
        # necessary
        if isinstance(message, int):
            n_bytes = (message.bit_length()+7)//8
            message = message.to_bytes(n_bytes, self.byte_order)
        elif isinstance(message, str):
            message = message.encode()
        if isinstance(message, bytes):
            message = bytearray(message)
        else:
            raise TypeError("Encoding of message to Arduino failed. "
                            f"Cannot encode {type(message)}.")
        # TODO: Arduino can not receive a message this long (possibly 513 bytes)
        # Get the message length and check if it can be sent with 1 byte
        payload_length = len(message)
        if payload_length > 255:
            raise ValueError(f"Length of message {message} too long "
                             "to be sent with a 1 byte long length.")
        # Iterate through message and split values that are larger
        # or equals special_byte into two different bytes. The second
        # byte gets inserted right after the first one.
        # Enumerate returns bytes as integers
        for i, value in enumerate(message):
            if value >= special_byte_int:
                message[i] = special_byte
                message.insert(
                    i + 1,
                    (value - special_byte_int).to_bytes(1, self.byte_order)
                    )
        # Insert the message length at the beginning of the message
        print(f"encode, before length: {message=}")
        message.insert(0, payload_length)
        print(f"encode, after length: {message=}")
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
        error_msg = msg_to_format.format(**fmt_data)
        # Emit error and clear errors.
        self.error_occurred.emit((error_code, error_msg))
        self.clear_errors()

    def is_decoded_message_acceptable(self, message):
        """Check whether a decoded message is ok.

        Check if meesage length is consistent with the length given
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
        """
        msg_length, *msg_data = message
        # Check if message length greater than zero
        if msg_length == 0:
            self.error_occurred.emit(ExtraSerialErrors.NO_MESSAGE_ERROR)
            return False
        if msg_length != len(msg_data):
            self.error_occurred.emit(
                ViPErLEEDHardwareError.ERROR_MSG_RCVD_INCONSISTENT
                )
            return False
        if msg_length not in (1, 2, 4):
            self.error_occurred.emit(
                ViPErLEEDHardwareError.ERROR_MSG_RCVD_INVALID
                )
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

        Returns
        -------
        bool
            Whether message is an error message
        """
        pc_error = self.port_settings.getint('available_commands', 'PC_ERROR')

        return len(message) == 1 and message[0] == pc_error

    def is_message_supported(self, messages):
        """Check if message contains one allowed command.

        Check whether the message is a supported command or not and
        keep track of which request was sent last to the Arduino.
        If a command is unsupported, emit an error_occurred signal
        with ExtraSerialErrors.UNSUPPORTED_COMMAND_ERROR. If there
        are more than two commands in one set of messages emit an
        error as the current structure does not allow this. If
        a change measurement mode command has been sent, update the
        countinuous measurement setting on the PC side. After clearing
        the message, remember the sent request and set the serial to
        busy.

        Parameters
        ----------
        message : str, object
            The same message given to self.send_message(), right
            before encoding and sending to the device

        Returns
        -------
        message_acceptable : bool
            True if message is acceptable. Only acceptable
            messages will then be sent.
        """
        # Remember!: messages[0] has to be of type Str
        # Check if message is valid data. Since the length of messages
        # can vary, the only reasonable check is if single bytes are
        # indeed an order specified in the config.
        available_commands = self.port_settings['available_commands'].values()
        command, *data = messages
        if command not in available_commands:
            # The first message must be a command
            self.error_occurred.emit(
                ExtraSerialErrors.UNSUPPORTED_COMMAND_ERROR
                )
            return False
        if any(message in available_commands for message in data):
            # Only the fist message can be a command; the others (if
            # any) can only be data to accompany the command itself.
            self.error_occurred.emit(
                ViPErLEEDHardwareError.ERROR_TOO_MANY_COMMANDS
                )
            return False

        pc_change_meas_mode = self.port_settings.get('available_commands',
                                                     'PC_CHANGE_MEAS_MODE')
        if command == pc_change_meas_mode:
            # Continuous measurement is the only thing we have to check
            # as it will potentially spam us with a lot of measurements
            # so we should react appropriately
            if not data:
                raise  # TODO: actually emit an error message
                return False
            self.__is_continuous_mode = bool(int.from_bytes(data[0],
                                             self.byte_order))

        self.__last_request_sent = command
        self.busy = True
        return True

    def process_messages(self):
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

        for message in self.unprocessed_messages:
            # If the length of the message is 1, then it has to be a
            # PC_OK byte
            if message.decode() == pc_ok:
                self.busy = False
            # Hardware config and measure values are both 4 bytes long.
            # Both commands are differentiated by the __last_request_sent
            # attribute and get processed accordingly. If continuous mode
            # is on and data got returned after a new request has been
            # sent, nothing is done with it.
            elif len(message) == 4:
                if self.__last_request_sent == pc_configuration:
                    self.busy = False
                    firmware, hardware = self.__firmware_and_hardware(message)
                    self.data_received.emit((firmware, hardware))
                elif (self.__last_request_sent == pc_set_voltage or
                        self.__last_request_sent == pc_measure_only):
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
                self.error_occurred.emit(
                    ViPErLEEDHardwareError.ERROR_ERROR_SLIPPED_THROUGH
                    )
            # If the message does not fit one of the lengths above, it
            # is no known message type.
            else:
                self.error_occurred.emit(
                    ViPErLEEDHardwareError.ERROR_MSG_RCVD_INVALID
                    )
        # TODO: I think we can empty unprocessed_messages every time after
        # checking the messages.
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
        """
        major, minor, *hardware = message
        firmware_version = f"{major}.{minor}"

        hardware_bits = self.port_settings['hardware_bits']                     # TODO: extend set_port_settings to check for all the extra sections that we are using in this implementation
        hardware_config = {'adc_0': False,
                           'adc_1': False,
                           'lm35': False,
                           'relay': False,
                           'i0_range': '0 -- 10 V',
                           'aux_range': '0 -- 10 V'}

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
            hardware_config[key] = present_or_closed

        return firmware_version, hardware_config


