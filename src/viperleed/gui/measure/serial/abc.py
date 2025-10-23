"""Module abc of viperleed.gui.measure.serial.

This module contains the definition of the SerialABC abstract base
class and the ExtraSerialErrors from ViPErLEEDErrorEnum derived
Enum that are used as base classes for serial communication with
hardware controllers supported by ViPErLEED. The serial communication
can also happen in a separate QThread to prevent stalling the
Graphical User Interface.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    'Florian Dörr (@FlorianDoerr)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2021-06-23'
__license__ = 'GPLv3+'

from abc import abstractmethod
from collections import deque
from configparser import NoOptionError
from configparser import NoSectionError

from PyQt5 import QtCore as qtc
from PyQt5 import QtSerialPort as qts

from viperleed.gui.measure.classes.abc import HardwareABC
from viperleed.gui.measure.dialogs.settingsdialog import SettingsHandler
from viperleed.gui.measure.hardwarebase import ViPErLEEDErrorEnum


_QU = qtc.Qt.QueuedConnection | qtc.Qt.UniqueConnection


SERIAL_ERROR_MESSAGES = {
    qts.QSerialPort.NoError: '',
    qts.QSerialPort.DeviceNotFoundError:
        ('No device on port {}. Please check the communication '
         'cables and/or update the list of ports.'),
    qts.QSerialPort.PermissionError:
        ('Permission error while opening port {}. The port may be '
         'already in use, or you may not have sufficient privileges.'),
    qts.QSerialPort.OpenError:
        ('Cannot open again a port on the same object. '
         'Close port {} by calling disconnect(), then try again.'),
    qts.QSerialPort.NotOpenError:
        ('Cannot perform the requested operation on a '
         'closed port. Open {} by calling .connect().'),
    qts.QSerialPort.WriteError:
        'Writing data to port {} failed.',
    qts.QSerialPort.ReadError:
        'Reading data from port {} failed.',
    qts.QSerialPort.ResourceError:
        ('Port {} became unavailable to the system. Device may be '
         'disconnected from the system. Check connection cables.'),
    qts.QSerialPort.UnsupportedOperationError:
        ('Cannot perform the requested operation on port {}: operation is '
         'either not supported or prohibited by the operating system.'),
    qts.QSerialPort.TimeoutError:
        ('Serial timeout error on port {}. This should not normally occur. '
         'It means someone incorrectly implemented a subclass of '
         'SerialABC, using waitForBytesWritten or waitForReadyRead '
         'instead of asynchronous behavior.'),
    qts.QSerialPort.UnknownError:
        'An unknown error occurred while accessing port {}.'
}


class ExtraSerialErrors(ViPErLEEDErrorEnum):
    """Data class for basic serial errors not available in QSerialPort."""

    NO_MESSAGE_ERROR = (50,
                        'Empty message received from controller.')
    NO_START_MARKER_ERROR = (51,
                             'Inconsistent message received from '
                             'controller (missing start marker). '
                             'Probably a communication error.')
    TIMEOUT_ERROR = (52,
                     'Serial Timeout: No message received in the '
                     'last {} sec. Check the communication cable '
                     'and/or serial-port settings in .ini file.')
    UNSUPPORTED_COMMAND_ERROR = (53,
                                 'Command {} is not supported '
                                 'by the controller. Check implementation '
                                 'and/or your configuration file.')
    PORT_NOT_OPEN = (54,
                     'Serial port could not be opened.')


# too-many-public-methods, too-many-instance-attributes, too-many-lines
class SerialABC(HardwareABC):
    """Base class for serial communication for a ViPErLEED controller."""

    data_received = qtc.pyqtSignal(object)
    about_to_trigger = qtc.pyqtSignal()

    __move_to_thread_requested = qtc.pyqtSignal(bool)  # True==connect          # TODO: Can be done with QMetaObject.invokeMethod

    _mandatory_settings = (
            ('serial_port_settings', 'MSG_END'),
            ('serial_port_settings', 'BYTE_ORDER', ('big', 'little'))
            )

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
        # Note that looking for default settings for the serial
        # will always fail, but the serial is guaranteed to be given
        # settings from the controller that instantiates it.
        super().__init__(settings=settings, **kwargs)

        self.__port = qts.QSerialPort(port_name, parent=self)

        self.__timeout = qtc.QTimer(parent=self)
        self.__timeout.setSingleShot(True)
        self.__timeout.timeout.connect(self.__on_serial_timeout)

        # Keeps track of whether the serial port is open or not.
        self._open = False

        with self.errors_delayed():
            self.set_settings(self._settings_to_load)

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

        # Keep track of whether we got an unacceptable message
        # after a .send_message()
        self.__got_unacceptable_response = False

        # unsent_messages is a deque of messages that have been
        # stored by the serial because it was not yet possible
        # to send them to the hardware controller. New messages
        # will automatically be appended to unsent_messages if
        # the serial is still waiting for a response from the
        # hardware. All messages in unsent_messages are processed
        # in the order they arrive at. When the serial receives
        # a valid answer from the hardware it will automatically
        # trigger an attempt to send the next message in
        # unsent_messages. Each element of unsent_messages is
        # a tuple whose first element is the command and its associated
        # data, and the second element is the timeout parameter.
        self.unsent_messages = deque()

        self.__move_to_thread_requested.connect(self.__on_moved_to_thread)
        self.busy_changed.connect(self.send_unsent_messages, type=_QU)

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
        return self.settings.get('serial_port_settings', 'BYTE_ORDER')

    @property
    def is_open(self):
        """Return whether this port is currently open."""
        return self._open

    def _set_is_open(self, open_status):
        """Set whether this port is currently open.

        Emits
        -----
        connection_changed
            If the connection status of the port has changed.
        """
        was_open = self._open
        self._open = open_status
        if open_status != was_open:
            self.connection_changed.emit(open_status)

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
        start_marker = self.settings.getint('serial_port_settings',
                                            'MSG_START', fallback=None)
        if start_marker is not None:
            start_marker = start_marker.to_bytes(1, self.byte_order)

        return {'START':  start_marker,
                'END': self.settings.getint(
                    'serial_port_settings', 'MSG_END'
                    ).to_bytes(1, self.byte_order)}

    @property
    def port_name(self):                                                        # TODO: change port_name to address
        """Return the name of the current port as a string."""
        return self.port.portName()

    @port_name.setter
    def port_name(self, new_port_name):
        """Set a new port name.

        This will disconnect an already-open port, but will not
        connect to the port with the newly set name. Explicitly
        use .connect() (after changing self.settings, if needed).

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
        if self.port:
            self.disconnect_()
        if isinstance(new_port_name, str):
            self.port.setPortName(new_port_name)
        elif isinstance(new_port_name, qts.QSerialPortInfo):
            self.port.setPort(new_port_name)
        else:
            raise TypeError(
                f"{self.__class__.__name__}: Invalid port "
                f"name type {type(new_port_name).__name__!r}. "
                "Should be 'str' or 'QSerialPortInfo'"
                )

    @property
    def port(self):
        """Return the underlying QSerialPort."""
        return self.__port

    @port.setter
    def port(self, port_name):
        """Create a new QSerialPort with name port_name.

        Port settings should already exist.

        Parameters
        ----------
        port_name : str or QSerialPortInfo
            Information needed to create and open a new port
        """
        if self.port:
            self.disconnect_()
        self.__port = qts.QSerialPort(port_name, parent=self)
        self.connect_()

    @classmethod
    def _cannot_check_settings(cls, *_):
        """Raise TypeError as SerialABC settings come from controllers."""
        # Generally speaking controllers have to find the approriate
        # settings for their serial. Therefore this method should never
        # be called.
        raise TypeError(f'{cls.__name__} was asked to perform a '
                        'settings check. Serials should never determine '
                        'appropriate settings on their own.')

    is_matching_default_settings = _cannot_check_settings
    is_matching_user_settings = _cannot_check_settings
    is_settings_for_this_class = _cannot_check_settings

    @qtc.pyqtSlot(object)
    def set_settings(self, new_settings):
        """Change settings of the port.

        This will disconnect an already-open port, but will not
        connect to the port with the new_settings. Explicitly
        use .connect() after changing port settings.

        Settings are loaded only if they are valid. Otherwise
        the previous settings stay in effect.

        This method can be used as a slot for a signal carrying
        a dict or a ConfigParser.

        Parameters
        ----------
        new_settings : dict or ConfigParser or str or Path or ViPErLEEDSettings
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
            'baud_rate': number, optional
                Baud rate for the transmission. Will be interpreted
                as an integer. Default is 9600.
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
            If new_settings is neither a dict, ConfigParser, string
            or path and if an element of the mandatory_settings is
            None or has a length greater than 3.

        Returns
        -------
        settings_valid : bool
            True if the new settings given were accepted.

        Emits
        -----
        QObjectSettingsErrors.MISSING_SETTINGS
            If new_settings is missing.
        QObjectSettingsErrors.INVALID_SETTINGS
            If any element of the new_settings does not fit the
            mandatory_settings.
        """
        if not super().set_settings(new_settings):
            return False
        self.disconnect_()
        return True

    def clear_errors(self):
        """Clear all errors.

        This method must be called in the overridden identify_error
        after the error has been correctly identified.

        Returns
        -------
        None.
        """
        if self.port.error() != qts.QSerialPort.NoError:
            self.port.clearError()
        self.__timeout.stop()
        self.__messages_since_error = []
        self.__got_unacceptable_response = False

    def decode(self, message):
        """Decode a message received from the serial.

        The base implementation of this method is a no-op, i.e., it
        returns the same message it got, unless the message is not
        acceptable (via is_decoded_message_acceptable). Subclasses
        can extend this method to decode possibly encoded messages.
        They should then "return super().decode(message)".

        Parameters
        ----------
        message : bytes or bytearray
            The message to be decoded.

        Returns
        -------
        decoded_message : bytearray
            The decoded message. Should return an empty message
            if not is_decoded_message_acceptable().
        """
        if not self.is_decoded_message_acceptable(message):
            return bytearray()
        return bytearray(message)

    def encode(self, message):
        """Encode a message to be sent via the serial port.

        The base implementation of this method is a no-op, i.e.,
        it returns the same message it got. Subclasses of the
        SerialABC class should override this method to actually
        encode messages before sending them to the device via the
        serial line. The overridden method should care about
        appropriately converting a supported object to an actual
        message in bytes.

        This method is called on every message right before
        writing it to the serial line.

        The overridden encode() should NOT add start and end
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

    def get_settings_handler(self):
        """Return a SettingsHandler object for displaying settings."""
        # super().get_settings_handler() is only peformed to check
        # whether a settings handler can be created from the current
        # settings.
        super().get_settings_handler()
        handler = SettingsHandler(self.settings, show_path_to_config=True)
        return handler

    @abstractmethod
    def identify_error(self, messages_since_error):
        """Identify which error occurred.

        This method is called whenever an error occurred. Concrete
        subclasses must override this method and emit a signal
        self.error_occurred(ViPErLEEDErrorEnum). New error codes
        should be defined in a specific ViPErLEEDErrorEnum subclass.

        Error identification means picking the correct message(s)
        that contain information about the error. This information
        may be contained in the error message itself, or may come
        as data messages following the error message. All messages
        that are received after an error message and before the error
        is identified can be used for identification, and will then
        be discarded. New messages will be processed only after
        identify_error() is successfully completed, i.e., it calls
        clear_errors().

        The method should be a no-op if identification of the error
        should fail because the data did not arrive yet.

        The reimplementation must call the clear_errors() method to
        be considered successfully completed. It makes sense to call
        clear_errors() only after an error has been successfully
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
    # pylint: disable=unused-argument
    def is_decoded_message_acceptable(self, message):
        """Check whether a decoded message is ok.

        This method can be called inside .decode(). If a message
        is unacceptable, decode() should return an empty bytearray.
        Only non-empty decoded messages will be stored. Moreover,
        all messages that come after an unacceptable message are
        discarded, until the next time the method .send_message()
        runs, or after calling .clear_errors().

        Override this method in subclasses. The base implementation
        always returns True.

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
    # pylint: enable=unused-argument

    @abstractmethod
    def is_error_message(self, message):
        """Check if a message corresponds to an error.

        This method must be overridden by concrete subclasses.
        It should return True if the message is an error message,
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
    # pylint: disable=unused-argument
    def is_message_supported(self, message):
        """Check whether message is a supported command.

        This method is guaranteed to be called exactly once
        during send_message() right before encoding and sending
        the message. The message will be encoded and sent only
        if this method returns True.

        This method can be overridden to (i) keep track of the
        last 'request' command sent, right before sending it,
        and (ii) check whether the command requested is one of
        those supported by the controller, and/or (iii) the
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
    # pylint: enable=unused-argument

    @abstractmethod
    def message_requires_response(self, *messages):
        """Return whether the messages to be sent require a response.

        Needs to be overridden in subclasses. This method should
        return True if the sent message requires a response from
        the connected hardware.

        Parameters
        ----------
        *messages : tuple
            Same arguments passed to send_message

        Returns
        -------
        bool
        """
        return True

    def moveToThread(self, thread):  # pylint: disable=invalid-name
        """Move self to a different thread by recreating self.port."""
        was_open = self.is_open
        if was_open:
            self.disconnect_()
        super().moveToThread(thread)
        self.__move_to_thread_requested.emit(was_open)

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
        return (message, *other_messages)

    @abstractmethod
    def process_received_messages(self):
        """Process data received into human-understandable information.

        This method is called every time one (or more) messages
        arrive at the serial line, unless there currently is an
        error message that has not yet been handled. Override
        identify_error() to handle errors.

        This method should be overridden by concrete subclass. The
        reimplementation should look into .unprocessed_messages and
        .pop() from there the messages it wants to process. Then it
        should emit the data_received signal for data that are worth
        processing (e.g., measurements).

        Each element in .unprocessed_messages is a bytearray.

        When overriding this method, it is a good idea to either not
        process any message (if not enough messages arrived yet) or
        to process all .unprocessed_messages. This prevents data loss,
        should an error message be received while there are still
        unprocessed messages.

        Emits
        -----
        data_received
            Any time the message received contains data that is
            worth processing.
        """
        return

    def send_message(self, message, *other_messages, timeout=None):
        """Send message to hardware via serial port.

        If the serial is busy when a message is supposed to be sent,
        the message is stored in .unsent_messages and sent at a later
        time when the serial is no longer busy. If a message is not
        supported send_message will return and not send the message.
        Afterwards a timeout timer is started and the message is
        encoded. START and END markers are added to the message. In the
        end the message is sent to the port specified in the settings
        via self.port_name.

        Parameters
        ----------
        message : object
            The message to send. This message will be encoded as
            per the self.encode() method. Typically will be either
            a bytes/bytearray, a number, or a string.
        *other_messages : objects
            Additional messages if the hardware needs to receive
            multiple consecutive messages.
        timeout : int or None, optional
            Maximum amount of time (in milliseconds) that the worker
            will spend while waiting to receive any messages following
            this message. If negative, it will never time out.
            A timeout of zero will cause an immediate time out. When
            the worker times out an ExtraSerialErrors.TIMEOUT_ERROR is
            emitted. If not given or None, the timeout will be taken
            from the option 'serial_port_settings'/'timeout' of
            self.settings. Should this option not exist, timeout
            will fall back to -1 (i.e., no timeout). Default is None

        Returns
        -------
        None.
        """
        if not self.is_open:
            self.emit_error(ExtraSerialErrors.PORT_NOT_OPEN)
            return

        all_messages = (message, *other_messages)

        if self.busy:
            self.unsent_messages.append((all_messages, timeout))
            return

        if not self.is_message_supported(all_messages):
            return

        self.__got_unacceptable_response = False

        if timeout is None:
            timeout = self.settings.getint('serial_port_settings', 'timeout',
                                           fallback=-1)

        timeout = int(timeout)
        if timeout >= 0:
            self.__timeout.start(timeout)
        all_messages = self.prepare_message_for_encoding(*all_messages)

        _requires_response = self.message_requires_response(*all_messages)
        self.busy = True
        encoded_messages = bytearray()
        for msg in all_messages:
            encoded = self.encode(msg)
            if self.msg_markers['START'] is not None:
                encoded[:0] = self.msg_markers['START']
            encoded.extend(self.msg_markers['END'])
            encoded_messages += encoded
        self.port.write(encoded_messages)

        if not _requires_response:
            self.busy = False

    @qtc.pyqtSlot()
    @qtc.pyqtSlot(bool)
    def send_unsent_messages(self, *_):
        """Send messages that have been stored."""
        if self.busy:
            return
        if self.unsent_messages:
            data, timeout = self.unsent_messages.popleft()
            self.send_message(*data, timeout=timeout)

    def set_private_busy(self, is_busy):
        """Set the _busy attribute without emitting busy_changed."""
        self._busy = is_busy

    def connect_(self, *__args):
        """Connect to currently selected port."""
        if self.is_open:
            return
        if not self.port.open(self.port.ReadWrite):
            self.emit_error(ExtraSerialErrors.PORT_NOT_OPEN)
            self.__print_port_config()
            self._set_is_open(False)
            return

        self._set_is_open(True)
        self.__set_up_serial_port()

        self.port.readyRead.connect(self.__on_bytes_ready_to_read)
        self.port.errorOccurred.connect(self.__on_serial_error)
        _ = self.port.readAll()

    def disconnect_(self, *__args):
        """Disconnect from connected serial port."""
        self.clear_errors()
        self.port.close()
        try:
            self.port.readyRead.disconnect(self.__on_bytes_ready_to_read)
        except TypeError:
            # port is already disconnected
            pass
        else:
            self.port.errorOccurred.disconnect(self.__on_serial_error)
        self._set_is_open(False)

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
            used (i.e., a MSG_START is present in self.settings),
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
            is a MSG_START in self.settings.
        """
        if not message:
            self.emit_error(ExtraSerialErrors.NO_MESSAGE_ERROR)
            return bytearray()

        start_marker = self.msg_markers['START']
        if start_marker is not None:
            # Protocol uses a start marker
            if message[:1] != start_marker:
                self.emit_error(ExtraSerialErrors.NO_START_MARKER_ERROR)
                return bytearray()
            message = message[1:]

        if not message:
            self.emit_error(ExtraSerialErrors.NO_MESSAGE_ERROR)
            return bytearray()
        return bytearray(message)

    @qtc.pyqtSlot()
    def __on_bytes_ready_to_read(self):
        """Read the message(s) received."""
        if self.__got_unacceptable_response:
            self.__timeout.stop()
            return

        msg = bytes(self.port.readAll())
        if self.msg_markers['END'] not in msg:
            # Only a fraction of a message has been read.
            self.__last_partial_message.extend(msg)
            return

        self.__timeout.stop()

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

    @qtc.pyqtSlot(bool)
    def __on_moved_to_thread(self, was_connected):
        """Remake QSerialPort in the new thread; connect if needed."""
        self.port = self.port_name
        if was_connected:
            self.connect_()

    @qtc.pyqtSlot(qts.QSerialPort.SerialPortError)
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
        self.emit_error((error_code, SERIAL_ERROR_MESSAGES[error_code]),
                        self.port_name)
        self.clear_errors()

    @qtc.pyqtSlot()
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
        timeout = self.settings.getint('serial_port_settings', 'timeout')
        self.emit_error(ExtraSerialErrors.TIMEOUT_ERROR,
                   round(timeout/1000, 1))

    def __set_up_serial_port(self):
        """Load settings to serial port."""
        settings = self.settings
        self.port.setBaudRate(
            int(settings.getfloat('serial_port_settings',
                                  'baud_rate',
                                  fallback=self.port.baudRate()))
            )
        self.port.setBreakEnabled(
            settings.getboolean('serial_port_settings',
                                'break_enabled',
                                fallback=self.port.isBreakEnabled())
            )
        self.port.setDataBits(
            settings.getint('serial_port_settings',
                            'data_bits',
                            fallback=self.port.dataBits())
            )
        self.port.setDataTerminalReady(
            settings.getboolean('serial_port_settings',
                                'data_terminal_ready',
                                fallback=self.port.isDataTerminalReady())
            )

        try:
            setting = settings.get('serial_port_settings', 'flow_control')
        except (NoSectionError, NoOptionError):
            flow_control = self.port.flowControl()
        else:
            flow_control = getattr(self.port, setting)
        self.port.setFlowControl(flow_control)

        try:
            setting = settings.get('serial_port_settings', 'parity')
        except (NoSectionError, NoOptionError):
            parity = self.port.parity()
        else:
            parity = getattr(self.port, setting)
        self.port.setParity(parity)

        if flow_control is not self.port.HardwareControl:
            # Do the following only if not under HardwareControl,
            # otherwise the following causes UnsupportedOperationError
            self.port.setRequestToSend(
                settings.getboolean('serial_port_settings',
                                    'request_to_send',
                                    fallback=self.port.isRequestToSend())
                )
        self.port.setStopBits(
            settings.getint('serial_port_settings',
                            'stop_bits',
                            fallback=self.port.stopBits())
            )

    def __print_port_config(self):  # For debug
        """Print the configuration of an open port."""
        print(
            f"### Serial port settings of {self.port.portName()} ###",
            f"baudRate: {self.port.baudRate()}",
            f"breakEnabled: {self.port.isBreakEnabled()}",
            f"dataBits: {self.port.dataBits()}",
            f"dataTerminalReady: {self.port.isDataTerminalReady()}",
            f"flowControl: {self.port.flowControl()}",
            f"parity: {self.port.parity()}",
            f"requestToSend: {self.port.isRequestToSend()}",
            f"stopBits: {self.port.stopBits()}",
            sep='\n', end='\n\n'
            )
