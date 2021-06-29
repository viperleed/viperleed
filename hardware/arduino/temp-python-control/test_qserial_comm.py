"""Module test_qserial_comm.py

Used for simple testing of the communication with Arduino.

Created: 2021-06-29
Author: Michele Riva
"""

# For a simple example of a serial on a different thread see
# https://www.programmersought.com/article/57762641596/

import sys

from PyQt5 import (QtWidgets as qtw,
                   QtCore as qtc,
                   QtSerialPort as qts)


TIMEOUT = 5000  # milliseconds


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

        self.__port.readyRead.connect(self.on_bytes_ready_to_read)
        self.__port.errorOccurred.connect(self.on_serial_error)

        for key in ('select_port', 'update_ports', 'connect'):
            self._ctrls[key].setEnabled(False)
        for key in ('disconnect', 'msg_to_send', 'send'):
            self._ctrls[key].setEnabled(True)

    def serial_disconnect(self, *__args):
        """Disconnect from connected serial port."""
        self.__port.close()

        for key in ('select_port', 'update_ports', 'connect'):
            self._ctrls[key].setEnabled(True)
        for key in ('disconnect', 'msg_to_send', 'send'):
            self._ctrls[key].setEnabled(False)

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

        if self.__port.write(bytes(msg_out, 'utf-8')) < 0:
            print('Could not send', flush=True)
            ret = qtw.QMessageBox(qtw.QMessageBox.Critical,
                                  'Error sending message',
                                  'Could not send message. Try again.',
                                  qtw.QMessageBox.Ok, self)
            return

        # if not self.__port.waitForBytesWritten(TIMEOUT):
            # print('Timeout while sending', flush=True)
            # ret = qtw.QMessageBox(qtw.QMessageBox.Critical,
                                 # 'Error sending message',
                                 # f"Timed out ({TIMEOUT/1000} s) while "
                                 # "sending message.",
                                 # qtw.QMessageBox.Ok, self)
            # return

        # self._ctrls['msg_received'].setText('')
        # if not self.__port.waitForReadyRead(TIMEOUT):
            # print('Timeout while reading', flush=True)
            # ret = qtw.QMessageBox(qtw.QMessageBox.Critical,
                                 # 'No message received',
                                 # f"Timed out ({TIMEOUT/1000} s) while "
                                 # "waiting for response.",
                                 # qtw.QMessageBox.Ok, self)

    def on_bytes_ready_to_read(self):
        """Read the message received."""
        msg = self.__port.readAll()
        self._ctrls['msg_received'].setText(bytes(msg))
    
    def on_serial_error(self, error_code):
        """React to a serial-port error."""
        print("Serial error occurred. Code:", error_code)

    def closeEvent(self, event):
        """Reimplement closeEvent to also close open ports."""
        self.__port.close()
        super().closeEvent(event)


if __name__ == '__main__':
    app = qtw.QApplication(sys.argv)
    window = MainWindow()
    window.show()

    sys.exit(app.exec_())
