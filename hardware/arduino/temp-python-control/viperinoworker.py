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
from configparser import ConfigParser

# ViPErLEED modules
from serialworker import ExtraSerialErrors, SerialWorkerABC


class ViPErLEEDHardwareError(ViPErLEEDErrorEnum):
    ERROR_NO_ERROR = (0, "No error")
    ERROR_SERIAL_OVERFLOW = (1, "Hardware overflow of Arduino serial.")
    ERROR_MSG_TOO_LONG = (2, "Sent message too long.")
    ERROR_MSG_SENT_INCONSITENT = (3, "Sent message has not the length "
                                     "it is supposed to have.")
    ERROR_MSG_RCVD_INCONSITENT = (13,"Received message has not the "
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


class ViPErinoSerialWorker(SerialWorkerABC):
    """Class for serial communication with
    the Arduino Micro ViPErLEED controller.
    """
    def __init__(self, settings=None, port_name=''):
        """Initialize serial worker object."""
        super().__init__(port_name=port_name, settings=settings)

    def decode(self, message):
        """
        Function decodes messages received from the Arduino. If there
        are SPECIAL_BYTEs in the message, it decodes them and turns them
        into single bytes.

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

        if int.from_bytes(message[0], self.byte_order) == 0:
            self.error_occurred.emit(ExtraSerialErrors.NO_MESSAGE_ERROR)
        for index, value in enumerate(message[1:]):
            if value == special_byte:
                message[index] += message[index+1]
                message.pop(index+1)
        if message[0] != len(message[1:]):
            self.error_occurred.emit(
                ViPErLEEDHardwareError.ERROR_MSG_RCVD_INCONSITENT
                )
        message.pop(0)

        return message

    def encode(self, message):
        """
        Function encodes messages sent to the Arduino. If there are
        SPECIAL_BYTEs in the message, it decodes them and turns them
        into single bytes.

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
        payload_length = len(message)
        if payload_length > 255:
            raise ValueError(f"Length of message {message} too long "
                             "to be sent with a 1 byte long length.")
        for i, value in enumerate(message):
            value_int = int.from_bytes(value, self.byte_order)
            if value_int >= special_byte_int:
                message[i] = special_byte
                message.insert(
                    i + 1,
                    (value_int - special_byte_int).to_bytes(1, self.byte_order)
                    )
        message.insert(0, payload_length)
        return message

        def identify_error(self, messages_since_error):
        """
        Indentifies the error sent back by the Arduino. Error data will 
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
        if len(self.messages_since_error) < 2:
            return
        error_state, error_code = messages_since_error[1]
        error_state = int.from_bytes(error_state, self.byte_order)
        #
        arduino_states = {int(code): state
                          for state, code
                          in self.port_settings['arduino_states'].items()}
        msg_to_format = ("ViPErLEED hardware {error_name} occurred while in"
                         " {state}. Reason: {err_details}")

        for code, error_name in ViPErLEEDHardwareError.as_dict().items():
            if code == int.from_bytes(error_code, self.byte_order):
                _, err_details = getattr(ViPErLEEDHardwareError, error_name)
                msg = msg_to_format.format(error_name,
                                           arduino_states[error_state],
                                           err_details)
                self.error_occurred.emit((code, msg))
                self.clear_errors()
                return




