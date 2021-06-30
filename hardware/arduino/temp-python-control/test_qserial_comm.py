"""Module test_qserial_comm.py

Used for simple testing of the communication with Arduino.

Created: 2021-06-29
Author: Michele Riva
"""

# For a simple example of a serial on a different thread see
# https://www.programmersought.com/article/57762641596/

# Some good info: https://stackoverflow.com/questions/42576537/

import time
import sys

from PyQt5 import (QtWidgets as qtw,
                   QtCore as qtc,
                   QtSerialPort as qts)


TIMEOUT = 5000  # milliseconds
MSG_START = b'\x00'
MSG_END = b'\xff'


class MainWindow(qtw.QWidget):
    """Simple window for basic tests."""

    def __init__(self, parent=None):
        """Initialize window."""
        super().__init__(parent)

        self._ctrls = {'select_port': qtw.QComboBox(),
                       'update_ports': qtw.QPushButton('Update port list'),
                       'connect': qtw.QPushButton('Connect'),
                       'disconnect': qtw.QPushButton('Disconnect'),
                       'send': qtw.QPushButton('Send'),
                       'msg_to_send': qtw.QLineEdit(),
                       'msg_received': qtw.QLineEdit()}

        self.__port = qts.QSerialPort()

        self.__compose()
        self.__connect()

        self.update_port_list()

    def __compose(self):
        """Place children widgets."""
        layout = qtw.QGridLayout()

        layout.addWidget(self._ctrls['select_port'], 0, 0)
        layout.addWidget(self._ctrls['update_ports'], 0, 1)
        layout.addWidget(self._ctrls['connect'], 1, 0)
        layout.addWidget(self._ctrls['disconnect'], 1, 1)
        layout.addWidget(self._ctrls['msg_to_send'], 2, 0)
        layout.addWidget(self._ctrls['send'], 2, 1)
        layout.addWidget(self._ctrls['msg_received'], 3, 0)

        self._ctrls['msg_received'].setReadOnly(True)
        self._ctrls['msg_to_send'].setEnabled(False)
        self._ctrls['disconnect'].setEnabled(False)
        self._ctrls['send'].setEnabled(False)

        self.setLayout(layout)

    def __connect(self):
        """Connect relevant signals."""
        self._ctrls['connect'].clicked.connect(self.serial_connect)
        self._ctrls['disconnect'].clicked.connect(self.serial_disconnect)
        self._ctrls['send'].clicked.connect(self.send_message)
        self._ctrls['update_ports'].clicked.connect(self.update_port_list)

    def update_port_list(self, *__args):
        """Update list of available ports."""
        ports = qts.QSerialPortInfo().availablePorts()

        port_list = [f"{port.portName()}:{port.description()}"
                     for port in ports]

        self._ctrls['select_port'].clear()
        self._ctrls['select_port'].addItems(port_list)

    def serial_connect(self, *__args):
        """Connect to currently selected port."""
        name = self._ctrls['select_port'].currentText().split(':')[0]
        self.__port = qts.QSerialPort(name)
        if not self.__port.open(self.__port.ReadWrite):
            print('Not open', flush=True)
            ret = qtw.QMessageBox(qtw.QMessageBox.Critical,
                                 'Error connecting',
                                 'An error occurred while trying to '
                                 'open the serial port. Error code: '
                                 f'{self.__port.error()}.',
                            qtw.QMessageBox.Ok, self)
            self.__port.clearError()
            return

        self.print_port_config()

        # TODO: here we need to set up the port. At least:
        self.__port.setBaudRate(self.__port.Baud115200)
        # self.__port.setDataBits(self.__port.Data8)  # Default is 8 --> OK
        self.__port.setDataTerminalReady(True)
        # setParity()     # Default is QSerialPort::NoParity --> OK
        # setStopBits()   # Default is 1 --> OK

        self.print_port_config()

        self.__port.readyRead.connect(self.on_bytes_ready_to_read)
        self.__port.errorOccurred.connect(self.on_serial_error)

        for key in ('select_port', 'update_ports', 'connect'):
            self._ctrls[key].setEnabled(False)
        for key in ('disconnect', 'msg_to_send', 'send'):
            self._ctrls[key].setEnabled(True)
        self._ctrls['msg_received'].setText('')

    def serial_disconnect(self, *__args):
        """Disconnect from connected serial port."""
        self.__port.close()

        for key in ('select_port', 'update_ports', 'connect'):
            self._ctrls[key].setEnabled(True)
        for key in ('disconnect', 'msg_to_send', 'send'):
            self._ctrls[key].setEnabled(False)
        
        self.__port.readyRead.disconnect(self.on_bytes_ready_to_read)
        self.__port.errorOccurred.disconnect(self.on_serial_error)

    def send_message(self, *__args):
        """Send message to port, and wait for reply."""
        msg_out = self._ctrls['msg_to_send'].text()
        if not msg_out:
            print('No message', flush=True)
            ret = qtw.QMessageBox(qtw.QMessageBox.Information,
                                 'No message to send',
                                 'No message to send. Please provide one.',
                                 qtw.QMessageBox.Ok, self)
            return

        msg_out = MSG_START + bytes(msg_out, 'utf-8') + MSG_END
        print("Sending", msg_out)
        if self.__port.write(msg_out) < 0:
            print('Could not send', flush=True)
            ret = qtw.QMessageBox(qtw.QMessageBox.Critical,
                                  'Error sending message',
                                  'Could not send message. Try again.',
                                  qtw.QMessageBox.Ok, self)

    def on_bytes_ready_to_read(self):
        """Read the message received."""
        msg = bytes(self.__port.readAll())
        msg = msg.replace(MSG_START, b'').replace(MSG_END, b'')
        
        self._ctrls['msg_received'].setText(str(msg))
    
    def on_serial_error(self, error_code):
        """React to a serial-port error."""
        print("Serial error occurred. Code:", error_code)

    def closeEvent(self, event):
        """Reimplement closeEvent to also close open ports."""
        self.__port.close()
        super().closeEvent(event)
    
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


if __name__ == '__main__':
    app = qtw.QApplication(sys.argv)
    window = MainWindow()
    window.show()

    sys.exit(app.exec_())
