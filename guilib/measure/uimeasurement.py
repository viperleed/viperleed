"""Module uimeasurement of viperleed.guilib.measure

===============================================
      ViPErLEED Graphical User Interface
===============================================

Created: 2021-10-12
Author: Michele Riva
Author: Florian Doerr

Defines the Measure class.
"""

import configparser
from pathlib import Path

import PyQt5.QtCore as qtc
import PyQt5.QtWidgets as qtw
import PyQt5.QtGui as qtg

import time

# ViPErLEED modules
from viperleed import guilib as gl
from viperleed.guilib.measure.measurement import ALL_MEASUREMENTS
from viperleed.guilib.measure.hardwarebase import class_from_name
from viperleed.guilib.basewidgets import MeasurementFigureCanvas as Figure
from viperleed.guilib.basewidgets import QDoubleValidatorNoDot

from viperleed.guilib.measure.camera.abc import CameraABC
from viperleed.guilib.measure.controller.abc import ControllerABC
from viperleed.guilib.measure.measurement.abc import MeasurementABC

TITLE = 'Measurement UI'


# @gl.broadcast_mouse
class Measure(gl.ViPErLEEDPluginBase):
    """A class that allows simulating LEED Patterns."""

    error_occurred = qtc.pyqtSignal(tuple)

    def __init__(self, parent=None):
        """Initialize window."""
        super().__init__(parent, name=TITLE)
        # Keep references to controls, dialogs, and some globals
        self._ctrls = {
            'measure': qtw.QPushButton("Start Measurement"),
            'abort': qtw.QPushButton("Abort"),
            'select': qtw.QComboBox(),
            'Start energy': qtw.QLineEdit(''),
            'End energy': qtw.QLineEdit(''),
            'Delta energy': qtw.QLineEdit(''),
            'save': qtw.QPushButton("Store settings"),
            'plots': [],
            'energy_input': qtw.QLineEdit(''),
            'set_energy': qtw.QPushButton("Set energy"),
            # TODO: add times/ make separate class/ read stuff from config and put it in the box upon opening extra window/add save button
            }

        self._dialogs = {}
        self._glob = {}

        # Set window properties
        self.setWindowTitle(TITLE)
        self.setAcceptDrops(True)

        self.__compose()
        self.do_this = None
        self.error_occurred.connect(self.__on_error_occurred)
        self.__errors = []
        self.__error_report_timer = qtc.QTimer(parent=self)
        self.__error_report_timer.setSingleShot(True)
        self.__error_report_timer.timeout.connect(self.__report_errors)

    def __compose(self):
        """Prepare menu."""
        self.setStatusBar(qtw.QStatusBar())
        self.setCentralWidget(qtw.QWidget())
        self.centralWidget().setLayout(qtw.QGridLayout())

        self._ctrls['measure'].setFont(gl.AllGUIFonts().buttonFont)
        self._ctrls['measure'].ensurePolished()
        self._ctrls['measure'].clicked.connect(self.__on_start_pressed)
        self._ctrls['measure'].setEnabled(True)

        self._ctrls['abort'].setFont(gl.AllGUIFonts().buttonFont)
        self._ctrls['abort'].ensurePolished()
        self._ctrls['abort'].setEnabled(False)

        self._ctrls['save'].setFont(gl.AllGUIFonts().buttonFont)
        self._ctrls['save'].ensurePolished()
        # self._ctrls['save'].clicked.connect(TODO: this)
        self._ctrls['save'].setEnabled(False)

        self._ctrls['select'].addItems(ALL_MEASUREMENTS.keys())
        self._ctrls['select'].setFont(gl.AllGUIFonts().buttonFont)
        self._ctrls['select'].ensurePolished()

        self._ctrls['set_energy'].setFont(gl.AllGUIFonts().buttonFont)
        self._ctrls['set_energy'].ensurePolished()
        self._ctrls['set_energy'].setEnabled(False)

        layout = self.centralWidget().layout()

        for key in ('Start energy', 'End energy', 'Delta energy', 'energy_input'):
            self._ctrls[key].setFont(gl.AllGUIFonts().labelFont)
            self._ctrls[key].ensurePolished()
            self._ctrls[key].setValidator(QDoubleValidatorNoDot())
            self._ctrls[key].validator().setLocale(qtc.QLocale.c())

        layout.addWidget(self._ctrls['measure'], 1, 1, 1, 1)
        layout.addWidget(self._ctrls['abort'], 1, 2, 1, 1)
        layout.addWidget(self._ctrls['select'], 2, 1, 1, 2)
        layout.addWidget(qtw.QLabel('Start energy ='), 3, 1, 1, 1)
        layout.addWidget(self._ctrls['Start energy'], 3, 2, 1, 1)
        layout.addWidget(qtw.QLabel('End energy ='), 4, 1, 1, 1)
        layout.addWidget(self._ctrls['End energy'], 4, 2, 1, 1)
        layout.addWidget(qtw.QLabel('Delta energy ='), 5, 1, 1, 1)
        layout.addWidget(self._ctrls['Delta energy'], 5, 2, 1, 1)
        layout.addWidget(self._ctrls['save'], 6, 1, 1, 2)
        layout.addWidget(self._ctrls['set_energy'], 7, 1, 1, 1)
        layout.addWidget(self._ctrls['energy_input'], 7, 2, 1, 1)
        self.statusBar().showMessage('Ready')

        fig = Figure()
        self._ctrls['plots'].append(fig)
        layout.addWidget(fig, 8, 1, 1, 2)

    def do_stuff(self):
        self.do_this.begin_measurement_preparation()
        self._ctrls['measure'].setEnabled(False)
        self._ctrls['select'].setEnabled(False)
        self.statusBar().showMessage('Busy')

    def __on_finished(self, *_):
        # After the measurement is done, close the serial ports.
        for controller in self.do_this.controllers:
            controller.serial.serial_disconnect()
        self._ctrls['measure'].setEnabled(True)
        self._ctrls['select'].setEnabled(True)
        self._ctrls['abort'].setEnabled(False)
        self.statusBar().showMessage('Ready')

    def __on_stuff_done(self, *_):
        print("\n#### DONE! ####")

    def __on_start_pressed(self):
        text = self._ctrls['select'].currentText()
        if self.do_this:
            self.do_this.thread.quit()
            self.__on_finished()
        config = configparser.ConfigParser(comment_prefixes='/',
                                           allow_no_value=True)
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
        measurement_cls = ALL_MEASUREMENTS[text]
        self.do_this = measurement_cls(config)

        if not isinstance(self.do_this, ALL_MEASUREMENTS['Time resolved']):  # TEMP. TODO: Handle data structure of time resolved
            self.do_this.new_data_available.connect(self.__on_new_data)

        self._ctrls['abort'].clicked.connect(self.do_this.abort)
        self._ctrls['set_energy'].clicked.connect(self.__on_set_energy)
        self.do_this.error_occurred.connect(self.error_occurred)
        for controller in self.do_this.controllers:
            controller.error_occurred.connect(self.error_occurred)
        for camera in self.do_this.cameras:
            camera.error_occurred.connect(self.error_occurred)
        self.do_this.finished.connect(self.__on_finished)
        self.do_this.finished.connect(self.__on_stuff_done)
        self.do_this.prepared.connect(self.__on_controllers_prepared)
        self._ctrls['abort'].setEnabled(True)                                   #tmp

        self.do_stuff()

    def __on_new_data(self):
        """Replot measured data."""
        fig = self._ctrls['plots'][0]
        fig.ax.cla()  # Clear old stuff
        meas = self.sender()

        fig.ax.plot(meas.data_points['nominal_energy'],
                    meas.data_points['I0'], '.')
        fig.ax.figure.canvas.draw_idle()

    def __on_controllers_prepared(self):
        """Enable abort button."""
        self._ctrls['abort'].setEnabled(True)

    def __on_set_energy(self):
        """Set energy on primary controller."""
        energy = float(self._ctrls['select'].currentText())
        self.do_this.set_LEED_energy((energy,0))

    def __on_error_occurred(self, error_info):
        """React to an error."""
        sender = self.sender()
        self.__errors.append((sender, *error_info))
        self.__error_report_timer.start(50)

    def __report_errors(self):
        if not self.__errors:
            return

        err_text = []
        for sender, error_code, error_message in self.__errors:
            if isinstance(sender, CameraABC):
                source = f"camera {sender.name}"
            elif isinstance(sender, ControllerABC):
                source = f"controller at {sender.serial.port_name}"
            elif isinstance(sender, MeasurementABC):
                source = f"measurement {sender.__class__.__name__}"
            else:
                source = "system or unknown"
            
            err_text.append(f"ERROR from {source}\n  "
                            f"error code: {error_code}"
                            f"\n{error_message}")
        
        _ = qtw.QMessageBox.critical(self, "Error", 
                                     "\n\n".join(err_text),
                                     qtw.QMessageBox.Ok)
        self.__errors = []