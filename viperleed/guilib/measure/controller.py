"""
Definition of controller classes for viperleed.measure

Created: 2020-01-10
Author: Michele Riva

Code is inspired by the project thesis work of Bernhard Mayr
"""

import PyQt5.QtCore as qtc
import PyQt5.QtSerialPort as qserial

class Controller:
    read_timeout = qtc.QTimer()
    read_timeout.setSingleShot(True)
    read_timeout.timout.connect(__on_read_timeout__)
    
    msg_read = qtc.pyqtSignal(str)

    def __init__(self):
        # self.name will be looked for among the available ports
        self.name = None

        # Communication port. The port object should have methods
        # port.write(bytes)
        # port.read()
        # port.inWaiting()  # returns the no. of bytes in the reading buffer that have not been processed yet  -> QSerialPort.bytesAvailable()
        self.port = qserial.QSerialPort()  # perhaps some other QIODevice?
        
        # self.__read is set every time one wants to read from the serial line.
        # If data is ready to be read and self.__read is True, then the data
        # will be read (if it conforms with the communication protocol)
        self.__read = False
        self.port.readyRead.connect(self.on_data_ready)
        
        # connect to the method handling a new message received via self.port
        self.msg_read.connect(self.on_new_message)
    
    @classmethod
    def __on_read_timeout__():
        raise IOError("Failed reading from controller. Connection timed out. "
                      "Disconnecting.")
        self.port.close()
    
    def on_new_message(self, msg):
        """
        This method defines the state of the controller after a new message is 
        received via self.port
        """
        pass ## TODO
    
    def send(self, bytes_to_send):
        """
        Send bytes bytes_to_send via self.port after encoding into the format
        required by the communication protocol
        """
        self.port.write(self.encode(bytes_to_send))
    
    def receive(self):
        """
        Handles the reading and decoding of messages
        """
        # The processing is actually done asynchronously, and only if the timer
        # does not time out, and it is done in the event handler conected to the
        # data ready signal
        self.read_timeout.start()
        self.__read = True
        self.__bytes_read = qtc.QByteArray()
    
    def on_data_ready(self):
        # stop the timeout timer: If there is new data, the port is active.
        self.read_timeout.stop()

        if self.__read:
            n = self.port.bytesAvailable()
            new_bytes = self.port.readAll()
            if new_bytes.length() < n:
                raise IOError("Not enough bytes read")
            self.__bytes_read.append(new_bytes)
            success, msg, rest = self.decode(self.__bytes_read)
            if success >= 0:
                # message decoded successfully. The remainder of the message
                # can be kept. WILL SEE IF NEED TO BE MORE PRECISE.
                self.__bytes_read = rest
                self.__read = False
                self.msg_read.emit(msg)
            elif success == -1:
                # message was partial. Still need to wait
                pass
            else:
                # some nonsense message. Empty the buffer and start over.
                self.__bytes_read = qtc.QByteArray()
    
    def connect(self, port=None):
        """
        Try opening port "port". If no port is given, try opening self.port. If unsuccessful find a port that has self.name in it, and try opening it.
        """
        if port is None:
            # here have to try opening self.port first, perhaps? or just
            # checking self.port.name (or whatever it's called)
            
            if self.name is None:
                raise RuntimeError("Controller has no name, and no port id "
                                   "given. Cannot connect.")
            # KEEP GOING FROM HERE
    
    def decode(self, bytes_in):
        """
        This function should be reimplemented to reflect the correct decoding
        of the bytes. Should return a tuple (exit_code, msg, rest). The exit
        codes are
        - code >= 0: message successfully decoded and it is meaningful
        - code = -1: message was partial (e.g., missing end character)
        - other values: message is rubbish and should be discarded.
        
        In all cases except code >= 0 "msg" should be an empty string.
        If code >=0, "msg" is the meaningful part of the bytes; all the bytes 
        that came after the meaningful message and have been read should be 
        returned unchanged in "rest".
        """
        raise NotImplementedError()
    
    def encode(self, bytes_in):
        """
        This method should be reimplemented to reflect the correct encoding 
        needed by the communication protocol, including whether the hardware
        expects MSB-first or LSB-first messages
        """
        return bytes_in