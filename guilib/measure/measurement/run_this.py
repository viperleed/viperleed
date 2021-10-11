import configparser
import sys
from pathlib import Path

from PyQt5 import (QtWidgets as qtw, QtCore as qtc)

from energy_setpoint import MeasureEnergySetpoint
from time_resolved import TimeResolved
from iv_video import IVVideo


class WindowDoesSomething(qtw.QWidget):

    error_occurred = qtc.pyqtSignal(tuple)

    def __init__(self):
        super().__init__()

        self.btn_do = qtw.QPushButton("Do stuff")
        self.btn_do.clicked.connect(self.do_stuff)

        self.btn_end = qtw.QPushButton("Close")
        self.btn_end.clicked.connect(self.__on_close_pressed)

        layout = qtw.QVBoxLayout()
        layout.addWidget(self.btn_do)
        layout.addWidget(self.btn_end)

        self.setLayout(layout)

        self.error_occurred.connect(self.__on_error_occurred)

    def do_stuff(self, checked):
        config = configparser.ConfigParser(comment_prefixes='/', allow_no_value=True)
        file_name = Path('C:/Users/Florian/Documents/Uni/Masterarbeit/ViperLEED/viperleed/guilib/measure/configuration/viperleed_config.ini')
        # TODO: path here not nice
        try:
            f = open(file_name, 'r')
            f.close()
            config.read(file_name)
        except OSError:
            raise FileNotFoundError(f"Could not open/read file: {file_name}. "
                                    "Check if this file exists and is in the "
                                    "correct folder.")
            return
        self.do_this = TimeResolved(config)

        self.do_this.error_occurred.connect(self.error_occurred)
        self.do_this.finished.connect(self.__on_close_pressed)
        self.do_this.finished.connect(self.__on_stuff_done)
        self.do_this.primary_controller.error_occurred.connect(self.error_occurred)
        for controller in self.do_this.secondary_controllers:
            if controller:
                controller.error_occurred.connect(self.error_occurred)

        self.do_this.begin_measurement_preparation()
        self.btn_do.setEnabled(False)

    def __on_error_occurred(self, error_code, error_message):
        print("Error occurred in", self.sender())
        print(f"code={error_code}")
        print(f"message={error_message}")

    def __on_close_pressed(self, *_):
        # After the measurement is done, close the serial ports.
        for controller in (self.do_this.primary_controller,
                           *self.do_this.secondary_controllers):
            if controller:
                controller.serial.serial_disconnect()
        self.btn_do.setEnabled(True)

    def __on_stuff_done(self, *_):
        print("\n#### DONE! ####")

    def closeEvent(self, event):
        self.__on_close_pressed()
        super().closeEvent(event)


if __name__ == '__main__':
    app = qtw.QApplication(sys.argv)
    window = WindowDoesSomething()
    window.show()

    sys.exit(app.exec_())
