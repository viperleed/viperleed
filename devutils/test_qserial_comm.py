"""Module test_qserial_comm.py

Used for simple testing of the communication with Arduino.

Created: 2021-06-29
Author: Michele Riva
"""

# For a simple example of a serial on a different thread see
# https://www.programmersought.com/article/57762641596/

# Some good info: https://stackoverflow.com/questions/42576537/

from configparser import ConfigParser
import os
from pathlib import Path
import random
import sys
import time

from PyQt5 import (QtWidgets as qtw,
                   QtCore as qtc,
                   QtSerialPort as qts)

from viperleed.gui.measure.serial.viperleedserial import ViPErLEEDSerial

TIMEOUT = 30000  # milliseconds


CONFIG_FILE = Path(__file__).resolve().parent / 'viperleed_hardware.ini'
assert CONFIG_FILE.exists(), f'No {CONFIG_FILE.name} at {CONFIG_FILE.parent}'
CONFIG = ConfigParser()
CONFIG.read(CONFIG_FILE)

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
                       'msg_received': qtw.QLineEdit(),
                       'serial_nr': qtw.QPushButton('Set serial nr')}

        self.__port = ViPErLEEDSerial(settings=CONFIG)

        self.__compose()
        self.__connect()

        self.update_port_list()

    def __compose(self):
        """Place children widgets."""
        layout = qtw.QGridLayout()

        self._ctrls['msg_received'].setSizePolicy(qtw.QSizePolicy.Preferred,
                                                  qtw.QSizePolicy.Expanding)

        layout.addWidget(self._ctrls['select_port'], 0, 0)
        layout.addWidget(self._ctrls['update_ports'], 0, 1)
        layout.addWidget(self._ctrls['connect'], 1, 0)
        layout.addWidget(self._ctrls['disconnect'], 1, 1)
        layout.addWidget(self._ctrls['msg_to_send'], 2, 0)
        layout.addWidget(self._ctrls['send'], 2, 1)
        layout.addWidget(self._ctrls['msg_received'], 3, 0)
        layout.addWidget(self._ctrls['serial_nr'], 3, 1)

        self._ctrls['msg_received'].setReadOnly(True)
        self._ctrls['msg_to_send'].setEnabled(False)
        self._ctrls['disconnect'].setEnabled(False)
        self._ctrls['send'].setEnabled(False)

        self.setLayout(layout)
        self.adjustSize()
        self.resize(2*self.size())

    def __connect(self):
        """Connect relevant signals."""
        self._ctrls['connect'].clicked.connect(self.serial_connect)
        self._ctrls['disconnect'].clicked.connect(self.serial_disconnect)
        self._ctrls['send'].clicked.connect(self.send_message)
        self._ctrls['update_ports'].clicked.connect(self.update_port_list)
        self._ctrls['msg_to_send'].editingFinished.connect(self.send_message)
        self._ctrls['serial_nr'].clicked.connect(self.set_serial_number)

    def update_port_list(self, *__args):
        """Update list of available ports."""
        ports = qts.QSerialPortInfo().availablePorts()

        port_list = [f"{port.portName()}:{port.description()} VID={hex(port.vendorIdentifier())} PID={hex(port.productIdentifier())}"
                     for port in ports]

        self._ctrls['select_port'].clear()
        self._ctrls['select_port'].addItems(port_list)

    def serial_connect(self, *__args):
        """Connect to currently selected port."""
        name = self._ctrls['select_port'].currentText().split(':')[0]
        self.__port.port_name = name
        self.__port.connect_()

        self.__port.data_received.connect(self.on_data_received)
        self.__port.error_occurred.connect(self.on_error_occurred)

        for key in ('select_port', 'update_ports', 'connect'):
            self._ctrls[key].setEnabled(False)
        for key in ('disconnect', 'msg_to_send', 'send'):
            self._ctrls[key].setEnabled(True)
        self._ctrls['msg_received'].setText('')

    def serial_disconnect(self, *__args):
        """Disconnect from connected serial port."""
        self.__port.disconnect_()

        for key in ('select_port', 'update_ports', 'connect'):
            self._ctrls[key].setEnabled(True)
        for key in ('disconnect', 'msg_to_send', 'send'):
            self._ctrls[key].setEnabled(False)

        self.__port.error_occurred.disconnect(self.on_error_occurred)

    def send_message(self, *__args, msg_to_send=None):
        """Send message to port, and wait for reply."""
        if not msg_to_send:
            msg_to_send = self._ctrls['msg_to_send'].text()
        if not msg_to_send:
            print('No message', flush=True)
            return

        msg_command, *data = msg_to_send.split(';')
        if data:
            data = data[0]
            try:
                msg_data = [int(d) for d in data.split(',')]
            except ValueError:
                print("Invalid message")
                return
            msg = (msg_command, msg_data)
        else:
            msg = (msg_command,)

        print(f"msg={msg}")
        self.__port.send_message(*msg, timeout=TIMEOUT)

    def on_data_received(self, data):
        """Read the message received."""
        # msg = bytes(self.__port.readAll())
        # msg = msg.replace(MSG_START, b'').replace(MSG_END, b'')

        self._ctrls['msg_received'].setText(str(data))

    def on_error_occurred(self, error_data):
        """React to a serial-port error."""
        error_code, error_msg = error_data
        print(error_msg, "\n", f"    (Error code: {error_code})")

    def closeEvent(self, event):
        """Reimplement closeEvent to also close open ports."""
        try:
            self.__port.disconnect_()
        except TypeError:
            # port is already disconnected
            pass
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

    def set_serial_number(self):
        """Generate and write serial number to device."""
        sr_nr = self.generate_serial_number()
        self.send_message(msg_to_send=f"s;{sr_nr}")
        print(sr_nr)
        return

    def generate_serial_number(self):
        """Generate a valid serial number."""
        sr_nr = 4*[None]
        for i, _ in enumerate(sr_nr):
            sr_nr[i] = self.get_random_sr_byte()
        sr_nr_string = ",".join(str(v) for v in sr_nr)
        return sr_nr_string

    def get_random_sr_byte(self):
        """Get random value."""
        unacceptable_values = [58 + i for i in range(7)]
        random_value = random.randint(48, 90)
        if random_value in unacceptable_values:
            random_value = self.get_random_sr_byte()
        return random_value


if __name__ == '__main__':
    app = qtw.QApplication(sys.argv)
    window = MainWindow()
    window.show()

    sys.exit(app.exec_())
