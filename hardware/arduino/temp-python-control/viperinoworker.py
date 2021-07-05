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
    ERROR_MSG_SENT_INCONSISTENT = (3, "Sent message has not the length "
                                      "it is supposed to have.")
    ERROR_MSG_RCVD_INCONSISTENT = (13, "Received message has not the "
                                       "length it is supposed to have.")
    ERROR_MSG_UNKNOWN = (4, "Unknown request from PC.")
    ERROR_MSG_DATA_INVALID = (5, "Request from PC contains invalid "
                                 "information.")
    ERROR_NEVER_CALIBRATED = (6, "The ADCs have never been "
                                 "calibrated since bootup")
    ERROR_TIMEOUT = (7, "Arduino timed out while waiting for something")
    ERROR_ADC_SATURATED = (8, "One of the ADC values reached saturation "
                              "and the gain cannot be decreased further.")
    ERROR_TOO_HOT = (9, "The temperature read by the LM35 is too high.")
    ERROR_RUNTIME = (255, "Some function has been called from an "
                          "inappropriate state. This is to flag "
                          "possible bugs for future development.")
    ERROR_UNKNOWN_ERROR = (256, "The error code received is not a known "
                                "error type. Possibly a newly introduced "
                                "error code that has not been added to the "
                                "ViPErLEEDHardwareError enum.Enum yet.")
    ERROR_ERROR_SLIPPED_THROUGH = (11, "For some reason and error identifier "
                                       "ended up in the unprocessed messages.")
    ERROR_MSG_RCVD_INVALID = (12, "Received message length does not fit the "
                                  "length of any data expected.")


class ViPErinoSerialWorker(SerialWorkerABC):
    """
    Class for serial communication with
    the Arduino Micro ViPErLEED controller.
    """
    def __init__(self, settings=None, port_name=''):
        """Initialize serial worker object."""
        self.__last_request_sent = ''
        super().__init__(port_name=port_name, settings=settings)

    def decode(self, message):
        """
        Function decodes messages received from the Arduino. If there
        are SPECIAL_BYTEs in the message, it decodes them and turns them
        into single bytes. Afterwards it checks if the message length is
        consistent with the supposed length sent by the Arduino. If true
        it will return the message without the message length.

        Parameters
        ----------
        message : bytes or bytearray
            The message to be decoded

        Returns
        -------
        message : bytes or bytearray
            The decoded message
        """
        special_byte = self.port_settings.getint('serial_port_settings',
                                                 'SPECIAL_BYTE')
        # Check if message length greater than zero
        if int.from_bytes(message[0], self.byte_order) == 0:
            self.error_occurred.emit(ExtraSerialErrors.NO_MESSAGE_ERROR)
        # Iterate through message and add special bytes up
        for index, value in enumerate(message[1:]):
            if value == special_byte:
                message[index] += message[index+1]
                message.pop(index+1)
        # Check if the message length is consistent with the supposed
        # length sent by the Arduino and remove message length from
        # message.
        if message[0] != len(message[1:]):
            self.error_occurred.emit(
                ViPErLEEDHardwareError.ERROR_MSG_RCVD_INCONSISTENT
                )
        message.pop(0)

        return message

    def encode(self, message):
        """
        Function encodes messages sent to the Arduino. If there are
        SPECIAL_BYTEs in the message, it decodes them and turns them
        into single bytes. Afterwards the original message length gets
        added at the beginning of the message.

        This method is called on every message right before
        writing it to the serial line.

        Parameters
        ----------
        message : object
            The message to be encoded

        Returns
        -------
        message : bytes or bytearray
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
        elif isinstance(message, bytes):
            pass
        else:
            raise TypeError("Encoding of message to Arduino failed. "
                            f"Cannot encode {type(message)}.")
        # Get the message length and check if it can be sent with 1 byte
        payload_length = len(message)
        if payload_length > 255:
            raise ValueError(f"Length of message {message} too long "
                             "to be sent with a 1 byte long length.")
        # Iterate through message and split values that are larger
        # or equals special_byte into two different bytes. The second
        # byte gets inserted right after the first one.
        for i, value in enumerate(message):
            value_int = int.from_bytes(value, self.byte_order)
            if value_int >= special_byte_int:
                message[i] = special_byte
                message.insert(
                    i + 1,
                    (value_int - special_byte_int).to_bytes(1, self.byte_order)
                    )
        # Insert the message length at the beginning of the message
        message.insert(0, payload_length)
        return message

    def identify_error(self, messages_since_error):
        """
        Identifies the error sent back by the Arduino. Error data will
        be returned right after sending the PC_ERROR byte. Therefore
        messages_since_error[0] will contain the PC_ERROR byte and
        messages_since_error[1] will contain both the error state and
        the error type in this order.


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
        error_state = int.from_bytes(error_state, self.byte_order)
        error_code = int.from_bytes(error_code, self.byte_order)
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
        state = arduino_states[error_state]
        error_name = current_error.name
        err_details = current_error[1]
        msg = msg_to_format.format(error_name, state, err_details)
        # Emit error and clear errors.
        ViPErLEEDHardwareError.temp_error = (error_code, msg)
        self.error_occurred.emit(ViPErLEEDHardwareError.temp_error)
        del ViPErLEEDHardwareError.temp_error
        self.clear_errors()

        # Old solution before the .from_code(error_code) function
        # for code, error_name in ViPErLEEDHardwareError.as_dict().items():
            # if code == error_code:
                # _, err_details = getattr(ViPErLEEDHardwareError, error_name)
                # msg = msg_to_format.format(error_name,
                                           # arduino_states[error_state],
                                           # err_details)
                # self.error_occurred.emit((code, msg))
                # self.clear_errors()
                # return

    def is_error_message(self, message):
        """
        Checks if the message length is equal 1 and has the value of
        PC_ERROR. If true it will return that the message is an error
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
        value = int.from_bytes(message[0], self.byte_order)

        return len(message) == 1 and value == pc_error

    def is_message_supported(self, message):
        """
        Checks whether message is a supported command or not and
        keeps track of which request was sent last to the Arduino.
        If a command is unsupported, it emits an error_occurred
        signal with ExtraSerialErrors.UNSUPPORTED_COMMAND_ERROR.

        Parameters
        ----------
        message : object
            The same message given to self.send_message(), right
            before encoding and sending to the device

        Returns
        -------
        message_acceptable : bool
            True if message is acceptable. Only acceptable
            messages will then be sent.
        """
        # Check if message is valid data. Since the length of messages
        # can vary, the only reasonable check is if single bytes are
        # indeed an order specified in the config.
        if not isinstance(message, (str, bytes, bytearray)):
            return True
        if len(message) > 1:
            return True
        available_commands = self.port_settings['available_commands'].keys()
        if message[0].decode() in available_commands:
            self.__last_request_sent = message[0].decode()
            return True
        # Emit error message
        self.error_occurred.emit(ExtraSerialErrors.UNSUPPORTED_COMMAND_ERROR)
        return False
        
    def process_messages(self):
        """Process data received into human-understandable information.

        This function is called every time one (or more) messages
        arrive on the serial line, unless there currently is an
        error message that has not yet been handled. Reimplement
        identify_error() to handle errors.

        This method should be reimplemented by any concrete subclass.
        The reimplementation should look into the list
        self.unprocessed_messages and .pop() from there the messages
        it wants to process. Then it should emit the data_received
        signal for data that are worth processing (e.g., measurements),
        as this is caught by the controller class.

        When reimplementing this method, it is a good idea to either
        not process any message (if not enough messages arrived yet)
        or to process all .unprocessed_messages. This prevents data
        loss, should an error message be received while there are
        still unprocessed messages.

        Emits
        -----
        data_received
            Any time the message received contains data that are
            worth processing.
        """
        for message in self.unprocessed_messages:
            # Only bytes and bytearrays will be sent by the Arduino
            if isinstance(message, (bytes, bytearray)):
                # If the length of the message is 1, then it has to be a
                # PC_OK byte
                if len(message) == 1:
                    self.data_received.emit(message.decode())
                # Hardware config and measure values are both 4 bytes long
                # How are we going to check if we convert to float or not?
                elif len(message) == 4:
                    self.data_received.emit(message)
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
        
def bytes_to_float(bytes_in):
    """Convert four bytes to float.

    Parameters
    ----------
    bytes_in : bytes or bytearray
        Should have length 4. It is assumed that the bytes are in
        big-endian format.

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
    return struct.unpack(">f", bytes_in)[0]
