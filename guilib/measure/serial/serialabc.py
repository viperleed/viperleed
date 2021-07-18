"""Module serialworkerabc of viperleed.?????.

========================================
   ViPErLEED Graphical User Interface
========================================

Created: 2021-06-23
Author: Michele Riva
Author: Florian Doerr

This module contains the definition of the SerialWorkerABC abstract
base class, ViPErLEEDErrorEnum base Enum class and ExtraSerialErrors
derived Enum that are used as base classes for serial communication
with controllers supported by ViPErLEED. The serial communication
can also happens in a separate QThread to prevent stalling the
Graphical User Interface.
"""

from abc import ABCMeta, abstractmethod
from configparser import ConfigParser

from PyQt5 import (QtCore as qtc,
                   QtSerialPort as qts)

from viperleed.guilib.measure.hardwarebase import (
    ViPErLEEDErrorEnum,
    config_has_sections_and_options,
    emit_error
    )

SERIAL_ERROR_MESSAGES = {
    qts.QSerialPort.NoError: "",
    qts.QSerialPort.DeviceNotFoundError:
        ("No device on port {}. Please check the communication "
         "cables and/or update the list of ports."),
    qts.QSerialPort.PermissionError:
        ("Permission error while opening port {}. The port may be "
         "already in use, or you may not have sufficient privileges."),
    qts.QSerialPort.OpenError:
        ("Cannot open again a port on the same object. "
         "Close port {} by calling disconnect(), then try again."),
    qts.QSerialPort.NotOpenError:
        ("Cannot perform the requested operation on a "
         "closed port. Open {} by calling .connect()."),
    qts.QSerialPort.WriteError:
        "Writing data to port {} failed.",
    qts.QSerialPort.ReadError:
        "Reading data from port {} failed.",
    qts.QSerialPort.ResourceError:
        ("Port {} became unavailable to the system. Device may be "
         "disconnected from the system. Check connection cables."),
    qts.QSerialPort.UnsupportedOperationError:
        ("Cannot perform the requested operation on port {}: operation is "
         "either not supported or prohibited by the operating system."),
    qts.QSerialPort.TimeoutError:
        ("Serial timeout error on port {}. This should not normally occur. "
         "It means someone incorrectly implemented a subclass of "
         "SerialWorkerABC, using waitForBytesWritten or waitForReadyRead "
         "instead of asynchronous behavior."),
    qts.QSerialPort.UnknownError:
        "An unknown error occurred while accessing port {}."
}


class ExtraSerialErrors(ViPErLEEDErrorEnum):
    """Data class for basic serial errors not available in QSerialPort."""
    NO_MESSAGE_ERROR = (50,
                        "Empty message received from controller.")
    NO_START_MARKER_ERROR = (51,
                             "Inconsistent message received from "
                             "controller (missing start marker). "
                             "Probably a communication error.")
    TIMEOUT_ERROR = (52,
                     "Serial Timeout: No message received in the "
                     "last {} sec. Check the communication cable.")
    UNSUPPORTED_COMMAND_ERROR = (53,
                                 "Command {} is not supported "
                                 "by the controller. Check implementation "
                                 "and/or your configuration file.")
    # The following two are fatal errors, and should make the GUI
    # essentially unusable, apart from loading appropriate settings
    INVALID_PORT_SETTINGS = (54,
                             "Invalid serial port settings: Required "
                             "settings {!r} missing or values "
                             "inappropriate. Check configuration file.")
    MISSING_SETTINGS = (55,
                        "Serial port cannot operate without settings. "
                        "Load an appropriate settings file before "
                        "proceeding.")


# TODO: Do we really need it to be a QObject? the only reason would be
# running in a separate thread, which seems we will not do
class SerialABC(qtc.QObject):
    """Base class for serial communication for a ViPErLEED controller."""

    error_occurred = qtc.pyqtSignal(tuple)
    data_received = qtc.pyqtSignal(object)
    serial_busy = qtc.pyqtSignal(bool)
    about_to_trigger = qtc.pyqtSignal()
    
    _mandatory_settings = [
            ('serial_port_settings', 'MSG_END'),
            ('serial_port_settings', 'BYTE_ORDER', ('big', 'little'))
            ]

    def __init__(self, settings=None, port_name=''):
        """Initialize serial worker object.

        Parameters
        ----------
        settings : dict or ConfigParser
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
        super().__init__()
        self.__port = qts.QSerialPort(port_name)

        # .__serial_settings is set via the following call to
        # set_port_settings() for extra checks and preprocessing
        self.__serial_settings = None
        self.set_port_settings(settings)

        # .unprocessed_messages is a list of all the messages
        # that came on the serial line and that have not been
        # processed yet. They should be processed in the call
        # to process_received_messages() by .pop()ping one at
        # a time. Each element is a bytearray.
        self.unprocessed_messages = []

        # __last_partial_message contains the 'head' of the last
        # message that was not yet fully read from the serial.
        # Used to prevent loosing bytes in asynchronous read.
        self.__last_partial_message = bytearray()

        # __messages_since_error is a list of all the serial
        # messages that came since the last unidentified error.
        # All these messages will be discarded automatically
        # after the error has been identified with identify_error()
        self.__messages_since_error = []

        # __busy keeps track of whether the serial line is
        # currently busy, e.g., a message was sent and we
        # should wait for the reply to come before we send
        # another message. Accessed via the .busy property.
        # Default behavior is not to consider the serial busy.
        # A reimplementation of is_message_supported can set
        # the .busy right before returning True.
        self.__busy = False

        # Keep track of whether we got an unacceptable message
        # after a .send_message()
        self.__got_unacceptable_response = False

        self.__timeout = qtc.QTimer()
        self.__timeout.setSingleShot(True)
        self.__timeout.timeout.connect(self.__on_serial_timeout)

    @property
    def byte_order(self):
        """Return the byte order used for serial messages.

        Returns
        -------
        byte_order : str
            Returns either 'big' or 'little'. The return value can
            be used in functions converting bytes to human-readable
            data types [e.g., int.from_bytes(bytes_obj, byte_order)].
        """
        return self.port_settings.get('serial_port_settings',
                                      'BYTE_ORDER')

    @property
    def busy(self):
        """Return whether the serial port is busy."""
        return self.__busy

    @busy.setter
    def busy(self, is_busy):
        """Set the serial to busy True/False.

        Parameters
        ----------
        is_busy : bool
            True if serial is busy
        """
        was_busy = self.busy
        is_busy = bool(is_busy)
        if was_busy is not is_busy:
            self.__busy = is_busy
            self.serial_busy.emit(self.busy)

    def get_port_name(self):
        """Return the name of the current port as a string."""
        return self.__port.portName()

    def set_port_name(self, new_port_name):
        """Set a new port name.

        This will disconnect an already-open port, but will not
        connect to the port with the newly set name. Explicitly
        use .connect() (after changing .port_settings, if needed).

        This method can be used as a slot for a signal carrying a
        string or a QSerialPortInfo.

        Parameters
        ----------
        new_port_name : str or QtSerialPortInfo
            Name of port to be opened. Use the availablePorts()
            method of QSerialPortInfo to determine appropriate
            port names.

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

    port_name = property(get_port_name, set_port_name)

    @property
    def msg_markers(self):
        """Return the markers signaling beginning/end of messages.

        Returns
        -------
        msg_markers : dict
            keys : {'START', 'END'}
            values : bytes or None
                Only the 'START' marker may be None, in case no
                start marker is used.
        """
        start_marker = self.port_settings.getint('serial_port_settings',
                                                 'MSG_START', fallback=None)
        if start_marker is not None:
            start_marker = start_marker.to_bytes(1, self.byte_order)

        return {'START':  start_marker,
                'END': self.port_settings.getint(
                    'serial_port_settings', 'MSG_END'
                    ).to_bytes(1, self.byte_order)}

    def get_port_settings(self):
        """Return the current settings for the port.

        Returns
        -------
        port_settings : ConfigParser
            The ConfigParser containing all the settings, among
            which the port settings that are available in section
            'serial_port_settings'.
        """
        return self.__serial_settings

    def set_port_settings(self, new_settings):
        """Change settings of the port.

        This will disconnect an already-open port, but will not
        connect to the port with the new_settings. Explicitly
        use .connect() after changing port settings.

        This method can be used as a slot for a signal carrying
        a dict or a ConfigParser.

        Parameters
        ----------
        new_settings : dict or ConfigParser
            Configuration of port. It must have a 'serial_port_settings'
            section.
            The following fields in new_settings['serial_port_settings']
            will be used:
            'BYTE_ORDER' : str, {'big', 'little'}
                Order in which bytes are transmitted. 'big' means
                most-significant-byte is transmitted earlier in time.
            'MSG_END' : bytes
                Byte to be used as a marker for the end of a message
            'MSG_START' : bytes, optional
                Byte to be used as a marker for the beginning of a
                message. Default is no start marker.
            'baud_rate': number
                Baud rate for the transmission. Will be interpreted
                as an integer.
            'break_enabled' : bool, optional
                Whether the transmission line should be placed in
                the break state. Default is False.
            'data_bits' : int, optional
                Number of data bits per data packet. Default is 8.
            'data_terminal_ready' : bool, optional
                Whether the Data Terminal Ready (DTR) line should
                be held high or not for the device to communicate.
                Default is False.
            'flow_control' : {'NoFlowControl', 'HardwareControl',
                              'SoftwareControl'}, optional
                The type of flow control used by the device.
                'HardwareControl' means the device uses RTS/CTS,
                'SoftwareControl' means the device uses XON/XOFF.
                Default is NoFlowControl.
            'parity' : {'NoParity', 'EvenParity', 'OddParity',
                        'SpaceParity', 'MarkParity'}, optional
                Which parity bits are used for each data packet.
                Default is NoParity.
            'request_to_send' : bool, optional
                Whether the Request To Send (RTS) line should be held
                high or not for communication. This option has no effect
                if 'flow_control' is set to 'HardwareControl', as the
                serial driver controls the RTS line. Default is False.
            'stop_bits' : int, optional
                Number of bits used as stop bits for each data packet.
                Default is 1.

        Raises
        ------
        TypeError
            If new_settings is neither a dict nor a ConfigParser
        KeyError
            If new_settings does not contain a 'serial_port_settings'
            section, of if the section does not contain the mandatory
            options MSG_END and BYTE_ORDER
        ValueError
            If 'BYTE_ORDER' is neither 'big' nor 'little'.
        """
        if new_settings is None:
            emit_error(self, ExtraSerialErrors.MISSING_SETTINGS)
            
        self.__serial_settings, invalid = config_has_sections_and_options(
            self,
            new_settings,
            self._mandatory_settings
            )
        for setting in invalid:
            emit_error(self, ExtraSerialErrors.INVALID_PORT_SETTINGS, setting)

        self.__serial_settings = new_settings
        self.__port.close()

    port_settings = property(get_port_settings, set_port_settings)

    def clear_errors(self):
        """Clear all errors.

        This method must be called in the reimplementation of
        identify_error after the error has been correctly identified.
        """
        self.__port.clearError()
        self.__timeout.stop()
        self.__messages_since_error = []
        self.__got_unacceptable_response = False

    # Disable pylint check as this is supposed
    # to be the signature for subclasses
    # pylint: disable=no-self-use
    def decode(self, message):
        """Decode a message received from the serial.

        The base implementation of this method is a no-op, i.e.,
        it returns the same message it got. Subclasses of the
        SerialWorkerABC class should reimplement this function
        to actually decode possibly encoded messages.

        Parameters
        ----------
        message : bytes or bytearray
            The message to be decoded

        Returns
        -------
        decoded_message : bytearray
            The decoded message. Should return an empty message
            if not is_decoded_message_acceptable().
        """
        if not self.is_decoded_message_acceptable(message):
            return bytearray()
        return bytearray(message)
    # pylint: enable=no-self-use

    # Disable pylint check as this is supposed
    # to be the signature for subclasses
    # pylint: disable=no-self-use
    def encode(self, message):
        """Encode a message to be sent via the serial port.

        The base implementation of this method is a no-op, i.e.,
        it returns the same message it got. Subclasses of the
        SerialWorkerABC class should reimplement this function
        to actually encode messages before sending them to the
        device via the serial line. The reimplemented method
        should care about appropriately converting a supported
        object to an actual message in bytes.

        This method is called on every message right before
        writing it to the serial line.

        The reimplemented encode() should NOT add start and end
        markers, as these are already handled by send_message().

        Parameters
        ----------
        message : object
            The message to be encoded

        Returns
        -------
        encoded_message : bytearray
            The encoded message
        """
        return message
    # pylint: enable=no-self-use

    @abstractmethod
    def identify_error(self, messages_since_error):
        """Identify which error occurred.

        This function is called whenever an error occurred. Concrete
        subclasses must reimplement this function and emit a signal
        self.error_occurred(ViPErLEEDErrorEnum). New error codes
        should be defined in a specific ViPErLEEDErrorEnum subclass.

        Error identification means picking the correct message(s)
        that contain information about the error. This information
        may be contained in the error message itself, or may come
        as data messages following the error message. All messages
        that are received after an error message and before the error
        is identified are discarded. New messages can be processed
        only if identify_error() is successfully completed, i.e.,
        it calls clear_errors().

        The function should be a no-op if identification of the error
        should fail because the data did not arrive yet.

        The reimplementation must call the clear_errors() method to
        be considered successfully completed. It only makes sense to
        call clear_errors() only after an error has been successfully
        identified (or if no error information is expected).

        Parameters
        ----------
        messages_since_error : list
            All the full messages that were received after the last
            error message, including the error message itself.
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
        return

    # Disable pylint check as this is supposed
    # to be the signature for subclasses
    # pylint: disable=no-self-use,unused-argument
    def is_decoded_message_acceptable(self, message):
        """Check whether a decoded message is ok.

        This method can be called inside .decode(). If a message
        is unacceptable, decode() should return an empty bytearray.
        Only non-empty decoded messages will be stored. Moreover,
        all messages that come after an unacceptable message are
        discarded, until the next time the method .send_message()
        runs, or after calling .clear_errors().

        Reimplement this method in subclasses. The base
        implementation always returns True.

        Parameters
        ----------
        message : bytearray
            The message to be checked

        Returns
        -------
        acceptable : bool
            True if message is acceptable
        """
        return True
    # pylint: enable=no-self-use,unused-argument

    @abstractmethod
    def is_error_message(self, message):
        """Check if a message corresponds to an error.

        This method must be reimplemented by concrete classes.
        It should return True if message is an error message,
        False otherwise. If an error message is received at any
        point, any message that is still unprocessed is discarded
        immediately.

        Parameters
        ----------
        message : bytes or bytearray
            A decoded message received through the serial line

        Returns
        -------
        bool
            Whether message is an error message
        """
        return

    # Disable pylint check as this is supposed
    # to be the signature for subclasses
    # pylint: disable=no-self-use,unused-argument
    def is_message_supported(self, message):
        """Check whether message is a supported command.

        This method is guaranteed to be called exactly once
        during send_message() right before encoding and sending
        the message. The message will be encoded and sent only
        if this method returns True.

        This method can be reimplemented to (i) keep track of
        the last 'request' command sent, right before sending
        it, and (ii) check whether the command requested is one
        of those supported by the controller, and/or (iii) the
        syntax of the command is acceptable.

        The base implementation always returns True.

        If a command is unsupported, emit an error_occurred
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
        return True
    # pylint: enable=no-self-use,unused-argument

    def prepare_message_for_encoding(self, message, *other_messages):
        """Prepare a message to be encoded.

        This method is guaranteed to be called exactly once on each
        call to send_message(message, *other_messages), right before
        each of the messages are (separately) encoded. Since .encode()
        expects each message to be a bytes or bytearray, this method
        can be called to turn messages into bytearrays.

        For example, assume that a command expects data as an integer,
        expressed as two bytes; one can pass send_message() the command
        (first argument) and the integer (second argument), and use this
        method to turn the integer to two bytes before encoding.

        The basic implementation returns messages as they were passed.

        Parameters
        ----------
        message : object
            The message to be sent. This is the same as the first
            argument given to send_message().
        *other_messages : object, optional
            Other messages to be sent. These are the same as the
            optional *other_messages given to send_message().

        Returns
        -------
        message_out : object
            Should have one of the types acceptable for encode()
        *other_messages_out : object
            should have one of the types acceptable for encode()
        """
        return message, *other_messages

    @abstractmethod
    def process_received_messages(self):
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

        Each element in .unprocessed_messages is a bytearray.

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
        return

    def send_message(self, message, *other_messages, timeout=None):
        """Send message to controller via serial port.

        Parameters
        ----------
        message : object
            The message to send. This message will be encoded as
            per the self.encode() method. Typically will be either
            a bytes/bytearray, a number, or a string.
        timeout : int or None, optional
            Maximum amount of time (in milliseconds) that the worker
            will spend while waiting to receive any messages following
            this message. If negative, it will never time out.
            A timeout of zero will cause an immediate time out. When
            the worker times out an ExtraSerialErrors.TIMEOUT_ERROR is
            emitted. If not given or None, the timeout will be taken
            from the option 'serial_port_settings'/'timeout' of
            self.port_settings. Should this option not exist, timeout
            will fall back to -1 (i.e., no timeout). Default is None

        Returns
        -------
        None.
        """
        all_messages = (message, *other_messages)
        if not self.is_message_supported(all_messages):
            return

        self.__got_unacceptable_response = False

        if timeout is None:
            timeout = self.port_settings.getint('serial_port_settings',
                                                'timeout',
                                                fallback=-1)

        timeout = int(timeout)
        if timeout >= 0:
            self.__timeout.start(timeout)
        all_messages = self.prepare_message_for_encoding(*all_messages)
        for message in all_messages:
            encoded = self.encode(message)
            if self.msg_markers['START'] is not None:
                encoded[:0] = self.msg_markers['START']
            encoded.extend(self.msg_markers['END'])
            print(f"{encoded=}")
            self.__port.write(encoded)

    def serial_connect(self, *__args):
        """Connect to currently selected port."""
        if not self.__port.open(self.__port.ReadWrite):
            # Could not open port. self.__on_serial_error()
            # will emit an error_occurred(ExtraSerialErrors)
            print('Not open', flush=True)  # TEMP
            print(self.port_name)
            self.print_port_config()
            return

        self.print_port_config()  # TEMP
        self.__set_up_serial_port()
        self.print_port_config()  # TEMP

        self.__port.readyRead.connect(self.__on_bytes_ready_to_read)
        self.__port.errorOccurred.connect(self.__on_serial_error)

    def serial_disconnect(self, *__args):
        """Disconnect from connected serial port."""
        self.clear_errors()
        self.__port.close()
        try:
            self.__port.readyRead.disconnect(self.__on_bytes_ready_to_read)
            self.__port.errorOccurred.disconnect(self.__on_serial_error)
        except TypeError:
            # port is already disconnected
            pass

    def __check_and_preprocess_message(self, message):
        """Check integrity of message.

        This function is called exactly once on each complete serial
        message received. Only messages passing this check will be
        decoded (and, later, processed).

        Parameters
        ----------
        message : bytes
            The message to be checked. Integrity checks only involve
            making sure that the message is not empty before and after
            removal of a start marker, if used. If a start marker is
            used (i.e., a MSG_START is present in self.port_settings),
            the method considers the message valid only if the first
            character is a MSG_START.

        Returns
        -------
        processed_message : bytearray
            Message stripped of the MSG_START marker, if used,
            otherwise the message is returned unchanged, unless
            it is invalid. An empty message is returned if the
            original message was invalid.

        Emits
        -----
        error_occurred(ExtraSerialErrors.NO_MESSAGE_ERROR)
            If message is empty (before or after removal of MSG_START)
        error_occurred(ExtraSerialErrors.NO_START_MARKER_ERROR)
            If the first byte of message is not MSG_START, if there
            is a MSG_START in self.port_settings.
        """
        if not message:
            emit_error(self, ExtraSerialErrors.NO_MESSAGE_ERROR)
            return bytearray()

        start_marker = self.msg_markers['START']
        if start_marker is not None:
            # Protocol uses a start marker
            if message[:1] != start_marker:
                emit_error(self, ExtraSerialErrors.NO_START_MARKER_ERROR)
                return bytearray()
            message = message[1:]

        if not message:
            emit_error(self, ExtraSerialErrors.NO_MESSAGE_ERROR)
            return bytearray()
        return bytearray(message)

    def __on_bytes_ready_to_read(self):
        """Read the message(s) received."""
        self.__timeout.stop()

        if self.__got_unacceptable_response:
            return

        msg = bytes(self.__port.readAll())

        if self.msg_markers['END'] not in msg:
            # Only a fraction of a message has been read.
            self.__last_partial_message.extend(msg)
            return

        head, *messages, tail = msg.split(self.msg_markers['END'])
        messages.insert(0, self.__last_partial_message + head)
        self.__last_partial_message = bytearray(tail)

        for i, message in enumerate(messages):
            # Make sure the messages are not empty, and get rid of
            # MSG_START if used. error_occurred is emitted in case
            # the message is empty (before of after stripping start)
            message = self.__check_and_preprocess_message(message)
            if not message:
                return
            decoded = self.decode(message)
            if not decoded:
                self.__got_unacceptable_response = True
                return
            messages[i] = decoded

        if self.__messages_since_error:
            # We got an error sometime before, but we did not yet
            # process information about the error because it did
            # not come yet.
            self.__messages_since_error.extend(messages)
            self.identify_error(self.__messages_since_error)
            return

        # Check if any of the messages we received is an error.
        # If yes, all unprocessed messages are discarded, and
        # all messages that arrived after the error are kept
        # until the error information has been identified
        for i, message in enumerate(messages):
            if self.is_error_message(message):
                self.unprocessed_messages = []
                self.__messages_since_error = messages[i:]
                self.identify_error(self.__messages_since_error)
                return
        self.unprocessed_messages.extend(messages)
        self.process_received_messages()

    def __on_serial_error(self, error_code):
        """React to a serial-port error.

        Parameters
        ----------
        error_code : one of QSerialPort.SerialPortError
            The error code of the serial port

        Emits
        -----
        error_occurred((error_code, msg))
            Always
        """
        if error_code == qts.QSerialPort.NoError:
            return
        emit_error(self, (error_code, SERIAL_ERROR_MESSAGES[error_code]),
                   self.port_name)
        self.clear_errors()

        print(error_msg)  # TEMP

    def __on_serial_timeout(self):
        """React to a serial timeout.

        Returns
        -------
        None.

        Emits
        -----
        error_occurred(ExtraSerialErrors.TIMEOUT_ERROR)
            Always
        """
        emit_error(self, ExtraSerialErrors.TIMEOUT_ERROR)

    def __set_up_serial_port(self):
        """Load settings to serial port."""
        settings = self.port_settings
        self.__port.setBaudRate(
            int(settings.getfloat('serial_port_settings',
                                  'baud_rate',
                                  fallback=self.__port.baudRate()))
            )
        self.__port.setBreakEnabled(
            settings.getboolean('serial_port_settings',
                                'break_enabled',
                                fallback=self.__port.isBreakEnabled())
            )
        self.__port.setDataBits(
            settings.getint('serial_port_settings',
                            'data_bits',
                            fallback=self.__port.dataBits())
            )
        self.__port.setDataTerminalReady(
            settings.getboolean('serial_port_settings',
                                'data_terminal_ready',
                                fallback=self.__port.isDataTerminalReady())
            )

        flow_control = getattr(
            self.__port, settings.get('serial_port_settings',
                                      'flow_control',
                                      fallback=self.__port.flowControl())
            )
        self.__port.setFlowControl(flow_control)

        parity = getattr(
            self.__port, settings.get('serial_port_settings',
                                      'parity',
                                      fallback=self.__port.parity())
            )
        self.__port.setParity(parity)

        if flow_control is not self.__port.HardwareControl:
            # Do the following only if not under HardwareControl,
            # otherwise the following causes UnsupportedOperationError
            self.__port.setRequestToSend(
                settings.getboolean('serial_port_settings',
                                    'request_to_send',
                                    fallback=self.__port.isRequestToSend())
                )
        self.__port.setStopBits(
            settings.getint('serial_port_settings',
                            'stop_bits',
                            fallback=self.__port.stopBits())
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
