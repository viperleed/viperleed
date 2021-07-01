"""Module serialwoker of viperleed.?????.

========================================
   ViPErLEED Graphical User Interface
========================================
*** module guilib.leedsim.leedparameters ***

Created: 2021-06-230
Author: Michele Riva
Author: Florian Doerr

This module contains the definition of the BaseSerialWorker class
and ExtraSerialErrors enum that are used as base classes for serial
communcation with controllers supported by ViPErLEED. The serial
communication happens in a separate thread to prevent stalling the
Graphical User Interface.
"""

import enum
from configparser import ConfigParser

from PyQt5 import (QtCore as qtc,
                   QtSerialPort as qts)


SERIAL_ERROR_MESSAGES = {
    qts.QSerialPort.DeviceNotFoundError:
        ("No device on port {}. Please check the communication "
         "cables and/or update the list of ports."),
    qts.QSerialPort.PermissionError:
        ("Permission error while opening port {}. The port may be "
         "alredy in use, or you may not have sufficient privileges.")
    qts.QSerialPort.OpenError:
        ("Cannot open again a port on the same object. "
         "Close port {} before by calling disconnect().")
    qts.QSerialPort.NotOpenError:
        ("Cannot perform the requested operation on a "
         "closed port. Open {} by calling .connect().")
    qts.QSerialPort.WriteError:
        ("Writing data to port {} failed.")
    qts.QSerialPort.ReadError:
        ("Reading data from port {} failed.")
    qts.QSerialPort.ResourceError:
        ("Port {} became unavailable to system. Device may be "
         "disconnected from the system. Check connection cables.")
    qts.QSerialPort.UnsupportedOperationError:
        ("Cannot perform the requested operation on port {}: operation is "
         "either not supported or prohibited by the operating system.")
    qts.QSerialPort.TimeoutError:
        ("Serial timeout error on port {}. This should not normally occur. "
         "It means someone incorrectly implemented a subclass of "
         "BaseSerialWorker, using waitForBytesWritten or waitForReadyRead "
         "instead of asynchronous behavior.")
    qts.QSerialPort.UnknownError:
        ("An unknown error occurred while acessing port {}.")
}


# TODO: perhaps nicer to have this just an Enum subclass and have tuples
#       where the first element is the code and the second one the message.
#       This way we can have consistent descriptions all in one place.
#       Perhaps we should also explicitly have (tuple, Enum) inheritance.
#       Actually even better, have a base ViPErLEEDErrorEnum class with
#       no members that inherits tuple, and that can be used as a parent
#       class for all ViPErLEED errors!
class ExtraSerialErrors(enum.Enum):
    """Data class for basic serial errors not available in QSerialPort."""
    NO_MESSAGE_ERROR = 50
    NO_START_MARKER_ERROR = 51
    TIMEOUT_ERROR = 52
    UNSUPPORTED_COMMAND_ERROR = 53


class BaseSerialWorker:
    """Base class for serial communication for a ViPErLEED controller."""

    error_occurred = qtc.pyqtSignal(int, str)
    data_received = qtc.pyqtSignal(object)

    def __init__(self, port_name='', settings=None):
        """Initialize window."""
        self.__port = qts.QSerialPort(port_name)

        if settings is None:
            raise ValueError(f"{self.__class__.__name__}: missing "
                             "mandatory serial port settings")

        # The next three are set via the call to the .port_settings
        # property for extra checks and preprocessing
        self.__serial_settings = None
        self.__byte_order = ''
        self.__msg_markers = {'START': None,
                              'END': b''}
        self.port_settings = settings
        self.unprocessed_messages = []

        self.__last_partial_message = b''
        self.__messages_since_error = []

        self.__timeout = qtc.QTimer()
        self.__timeout.setSingleShot(True)
        self.__timeout.timeout.connect(self.__on_serial_timeout)

    @property
    def port_name(self):
        """Return the name of the current port."""
        return self.__port.portName()

    @port_name.setter
    def port_name(self, new_port_name):
        """Set a new port name.

        This will disconnect an already-open port, but will not
        connect to the port with the newly set name. Explicitly
        use .connect() after changing .port_settings, if needed.

        Parameters
        ----------
        new_port_name : str or QtSerialPortInfo
            Name of port to be opened. Use the availablePorts()
            method of QSerialPortInfo to determine appropriate
            port names.

        Returns
        -------
        None.

        Raises
        ------
        TypeError
            If new_port_name is neither 'str' nor 'QSerialPortInfo'
        """
        self.__port.close()
        if isinstance(new_port_name, str):
            self.__port.setPortName(new_port_name)
        elif isinstance(new_port_name, qts.QSerialPortInfo):
            self.__port.setPort(new_port_name)
        else:
            raise TypeError(
                f"{self.__class__.__name__}: Invalid port "
                f"name type {type(new_port_name).__name__!r}. "
                "Should be 'str' or 'QSerialPortInfo'"
                )

    @property
    def port_settings(self):
        """Return the current settings for the port."""
        return self.__serial_settings

    @port_settings.setter
    def port_settings(self, new_settings):
        """Change settings of the port.

        Parameters
        ----------
        new_settings : dict or ConfigParser
            Configuration of port. It must have a 'serial_port_settings'
            section.
            The following fields in new_settings['serial_port_settings']
            will be used (no quotes if using a ConfigParser!):
            BYTE_ORDER : str, {'big', 'little'}
                Order in which bytes are transmitted. 'big' means
                most-significant-byte first.
            MSG_END : bytes
                Byte to be used as a marker for the end of a message
            MSG_START : bytes, optional
                Byte to be used as a marker for the beginning of a
                message. Default is no start marker.
            'baud_rate': number
                Baud rate for the transmission. Will be interpreted
                as an integer.
            'break_enabled' : bool, optional
                Whether the serial device makes use of the break
                condition for communication. Default is False.
            'data_bits' : int, optional
                Number of data bits per data packet. Default is 8.
            'data_terminal_ready' : bool, optional
                Whether the Data Terminal Ready (DTR) line should
                be held high or not for the device to communicate.
                Default is False.
            'flow_control' : {'NoFlowControl', 'HardwareControl',
                              'SoftwareControl'}, optional
                The type of flow control used by the device. Default
                is NoFlowControl.
            'parity' : {'NoParity', 'EvenParity', 'OddParity',
                        'SpaceParity', 'MarkParity'}, optional
                Which parity bits are used for each data packet.
                Default is NoParity.
            'request_to_send' : bool, optional
                Whether the Request To Send (RTS) line should be held
                high or not for communication. This option has no effect
                if 'flow_control' is set to 'HardwareControl'. Default
                is False.
            'stop_bits' : int, optional
                Number of bits used as stop bits for each data packet.
                Default is 1.

        Raises
        ------
        TypeError
            If new_settings is neither a dict nor a ConfigParser
        ValueError
            If new_settings does not contain a 'serial_port_settings'
            section, of if the section does not contain the mandatory
            options MSG_END and BYTE_ORDER.
        """
        if not isinstance(new_settings, (dict, ConfigParser)):
            raise TypeError(
                f"{self.__class__.__name__}: invalid type "
                f"{type(new_settings).__name__} for settings. "
                "Should be 'dict' or 'ConfigParser'."
                )
        if 'serial_port_settings' not in new_settings:
            raise ValueError(
                f"{self.__class__.__name__}: settings must "
                "contain a 'serial_port_settings' section"
                )
        if any(key not in new_settings['serial_port_settings']
               for key in ('MSG_END', 'BYTE_ORDER')):
            raise ValueError(
                f"{self.__class__.__name__}: section 'serial_port_settings' "
                "of the settings should contain both 'MSG_END' and "
                "'BYTE_ORDER' fields."
                )

        self.__serial_settings = new_settings
        self.__byte_order = new_settings.get('serial_port_settings',
                                             'BYTE_ORDER')
        start_marker = new_settings.getint('serial_port_settings',
                                           'MSG_START', fallback=None)
        if start_marker is not None:
            start_marker.to_bytes(1, self.__byte_order)

        self.__msg_markers = {
            'START':  start_marker,
            'END': new_settings.getint(
                'serial_port_settings', 'MSG_END'
                ).to_bytes(1, self.__byte_order)
            }

    def clear_errors(self):
        """Clear all errors.

        This method must be called in the reimplementation of
        identify_error after the error has been correctly identified.
        """
        self.__port.clearErrors()
        self._timeout.stop()
        self.__messages_since_error = []

    def decode(self, message):
        """Decode a serial message.

        The base implementation of this method is a no-op, i.e.,
        it returns the same message it got. Subclasses of the
        BaseSerialWorker class should reimplement this function
        to actually decode possibly encoded messages.

        Parameters
        ----------
        message : bytes or bytearray

        Returns
        -------
        decoded_message : bytes or bytearray
        """
        return message

    def identify_error(self, messages_since_error):
        """Identify which error occurred.

        This function is called whenever an error occurred. Concrete
        subclasses should reimplement this function and emit a signal
        self.error_occurred(code, text) where code is the error code
        and text is an informative message that will be displayed as
        a message box in the GUI. It is advisable to define new error
        codes in a specific enum.Enum subclass.

        Error identification means picking the correct message(s)
        that contain information about the error. This information
        may be contained in the error message itself, or may come
        as data messages following the error message. All messages
        that came after an error message are discarded after
        identify_error() successfully identifies error information.
        New messages can be processed only if identify_error() is
        successfully completed, i.e., is calls clear_errors().

        The function should be a no-op if identification of the error
        should fail because the data did not arrive yet.

        The reimplementation must call the clear_errors() method to
        be considered successfully completed. It makes sense to call
        clear_errors() only after an error has been successfully
        identified (or if no error information is expected).

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

        Raises
        ------
        NotImplementedError
            If this method is not reimplemented by concrete subclasses.

        Emits
        -----
        error_occurred(error_code, descriptive_text)
            After the error has been identified
        """
        raise NotImplementedError(
            f"{self.__class__.__name__} must reimplement "
            "identify_error(self, messages_since_error)!"
            )

    def is_error_message(self, message):
        """Check if a message corresponds to an error.

        This method must be reimplemented by concrete classes.
        It should return True if message is an error message,
        False otherwise.

        Parameters
        ----------
        message : bytes or bytearray

        Returns
        -------
        bool
            Whether message is an error message
        """
        raise NotImplementedError(
            f"{self.__class__.__name__} must reimplement "
            "is_error_message(self, message)!"
            )

    def is_message_supported(self, message):
        """Check whether message is a supported command.

        This method is guaranteed to be called exactly once
        during send_message() right before actually sending
        the message, and before encoding(). The message will
        be encoded and sent only if this method returns True.

        This method can be reimplemented to (i) keep track of
        the last 'request' command sent, right before sending
        it, and (ii) check whether the command requested is one
        of those supported by the controller, and/or (iii) the
        syntax of the command is acceptable.

        The base implementation always returns True.

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
        return True

    def process_messages(self):
        """Process data received into human-understandable information.

        This function is called every time one (or more) messages
        arrive on the serial line, provided that there is no unhandled
        error that came before message. Reimplement identify_error()
        to handle errors.

        This method should be reimplemented by any concrete subclass.
        The reimplementation should look into the list
        self.unprocessed_messages and .pop() from there the messages
        it wants to process. Then it should emit the data_received
        signal for data that are worth processing (e.g., measurements),
        as this is caught by the controller class.

        Emits
        -----
        data_received
            Any time the message received contains data that are
            worth processing.
        """
        raise NotImplementedError(
            f"{self.__class__.__name__} must reimplement "
            "process_message(self, message)!"
            )

    def send_message(self, message, timeout=-1):
        """Send message to controller via serial port.

        Parameters
        ----------
        message : object
            The message to send. This message will be encoded as
            per the self.encode() method. Typically will be either
            a bytes/bytearray, a number, or a string.
        timeout : int, optional
            Maximum amount of time in milliseconds that the worker
            will wait before emitting a timeout error. If negative,
            it will never time out. A timeout of zero will cause
            an immediate time out. Default is -1.

        Returns
        -------
        None.

        Emits       # TODO: probably nothing actually!
        -----
        error_occurred
            If message could not be sent
        """
        timeout = int(timeout)
        if timeout >= 0:
            self.__timoout.start(timeout)

        if not self.is_message_supported(message):
            return
        self.__port.write(self.encode(message))
        # TODO: the error that we're handling below should already
        # be taken care of by __on_serial_error, but needs checking!
        # if self.__port.write(self.encode(message)) < 0:
            # self.error_occurred.emit(
                # "Error sending message. Check communication and port settings?"
                # )

    def serial_connect(self, *__args):
        """Connect to currently selected port."""
        if not self.__port.open(self.__port.ReadWrite):
            print('Not open', flush=True)
            self.__port.clearError()
            return

        self.print_port_config()  # TEMP
        self.__set_up_serial_port()
        self.print_port_config()  # TEMP

        self.__port.readyRead.connect(self.__on_bytes_ready_to_read)
        self.__port.errorOccurred.connect(self.__on_serial_error)

    def serial_disconnect(self, *__args):
        """Disconnect from connected serial port."""
        self.__port.close()
        self.__port.readyRead.disconnect(self.on_bytes_ready_to_read)
        self.__port.errorOccurred.disconnect(self.on_serial_error)

    def __check_and_preprocess_message(self, message):
        """Check integrity of message.

        Parameter
        ---------
        message : bytes or bytearray
            The message to be checked. Integrity checks only involve
            making sure that the message is not empty before and after
            removal of a start marker if used. If a start marker is
            used (i.e., a MSG_START is present in self.port_settings),
            the method considers the message valid only if the first
            character is a MSG_START.

        Returns
        -------
        processed_message : bytes or bytearray or None
            Message stripped of the MSG_START marker, if used, otherwise
            the message is returned unchanged. Returns None if message
            is invalid.

        Emits
        -----
        error_occurred(ExtraSerialErrors.NO_MESSAGE_ERROR, msg)
            If message is empty (before or after removal of MSG_START
        error_occurred(ExtraSerialErrors.NO_START_MARKER_ERROR)
            If the first byte of message is not MSG_START, if there
            is a MSG_START in self.port_settings.
        """
        if not message:
            self.error_occurred.emit(
                ExtraSerialErrors.NO_MESSAGE_ERROR,
                "Empty message received from controller"
                )
            return

        start_marker = self.__msg_markers['START']
        if start_marker is not None:
            if message[0] != start_marker:
                self.error_occurred.emit(
                    ExtraSerialErrors.NO_START_MARKER_ERROR,
                    "Inconsistent message received from "
                    "controller (missing start marker). "
                    "Probably a communication error."
                    )
                return
            message = message[1:]

        if not message:
            self.error_occurred.emit(
                ExtraSerialErrors.NO_MESSAGE_ERROR,
                "Empty message received from controller"
                )
            return
        return message

    def __on_bytes_ready_to_read(self):
        """Read the message received."""
        self.__timeout.stop()

        msg = bytes(self.__port.readAll())
        head, *messages, tail = msg.split(self.__msg_markers['END'])
        messages.insert(0, self.__last_partial_message + head)
        self.__last_partial_message = tail

        for i, message in enumerate(messages):
            message = self.__check_and_preprocess_message(message)
            if message is None:
                return
            messages[i] = self.decode(message)

        if self.__messages_since_error:
            # We got an error sometime before, but we did not yet
            # process information about the error because it did
            # not come yet.
            self.__messages_since_error.extend(messages)
            self.identify_error(self.__messages_since_error)
            return

        for i, message in enumerate(messages):
            if self.is_error_message(message):
                self.__messages_since_error = messages[i:]
                self.identify_error(self.__messages_since_error)
                return
        self.unprocessed_messages.extend(messages)
        self.process_messages()

    def __on_serial_error(self, error_code):
        """React to a serial-port error.

        Parameters
        ----------
        error_code : one of QSerialPort.SerialPortError
            The error code of the serial port

        Emits
        -----
        error_occurred(error_code, msg)
            Always
        """
        error_msg = SERIAL_ERROR_MESSAGES[error_code].format(self.port_name)
        self.error_occurred.emit(error_code, error_msg)

        print(error_msg)

    def __on_serial_timeout(self):
        """React on a serial timeout.

        Returns
        -------
        None.

        Emits
        -----
        error_occurred(ExtraSerialErrors.TIMEOUT_ERROR, msg)
            Always
        """
        self.error_occurred.emit(
            ExtraSerialErrors.TIMEOUT_ERROR,
            "Serial Timeout: No message received in the last"
            f"{self.__timeout.interval()/1000:.2f} s. Check the "
            "communication cable."
            )

    def __set_up_serial_port(self):
        """Load settings to serial port."""
        settings = self.__serial_settings
        self.__port.setBaudRate(
            int(settings.getfloat('serial_port_settings',
                                  'baud_rate',
                                  self.__port.baudRate()))
            )
        self.__port.setBreakEnabled(
            settings.getbool('serial_port_settings',
                             'break_enabled',
                             self.__port.isBreakEnabled())
            )
        self.__port.setDataBits(
            settings.getint('serial_port_settings',
                            'data_bits',
                            self.__port.dataBits())
            )
        self.__port.setDataTerminalReady(
            settings.getbool('serial_port_settings',
                             'data_terminal_ready',
                             self.__port.isDataTerminalReady())
            )

        flow_control = getattr(
            self.__port, settings.get('serial_port_settings',
                                      'flow_control',
                                      self.__port.flowControl())
            )
        self.__port.setFlowControl(flow_control)

        parity = getattr(
            self.__port, settings.get('serial_port_settings',
                                      'parity',
                                      self.__port.parity())
            )
        self.__port.setParity(parity)

        if flow_control is not self.__port.HardwareControl:
            self.__port.setRequestToSend(
                settings.getbool('serial_port_settings',
                                 'request_to_send',
                                 self.__port.isRequestToSend())
                )
        self.__port.setStopBits(
            settings.getint('serial_port_settings',
                            'stop_bits',
                            self.__port.stopBits())
            )

    def print_port_config(self):
        """Print the configuration of an open port."""
        print(
            f"### Serial port settings of {self.__port.portName()} ###",
            f"baudRate: {self.__port.baudRate()}",
            f"breakEnabled: {self.__port.isBreakEnabled()}",
            f"dataBits: {self.__port.dataBits()}",
            f"dataTerminalReady: {self.__port.isDataTerminalReady()}",
            f"flowControl: {self.__port.flowControl()}",
            f"parity: {self.__port.parity()}",
            f"requestToSend: {self.__port.isRequestToSend()}",
            f"stopBits: {self.__port.stopBits()}",
            sep='\n', end='\n\n'
            )
