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
import inspect

import PyQt5.QtCore as qtc
import PyQt5.QtWidgets as qtw

# ViPErLEED modules
from viperleed import guilib as gl
from viperleed.guilib.measure.measurement import ALL_MEASUREMENTS
from viperleed.guilib.measure.hardwarebase import get_devices
from viperleed.guilib.basewidgets import MeasurementFigureCanvas as Figure
from viperleed.guilib.basewidgets import QDoubleValidatorNoDot

from viperleed.guilib.measure.camera.abc import CameraABC
from viperleed.guilib.measure.controller.abc import ControllerABC
from viperleed.guilib.measure.measurement.abc import MeasurementABC
from viperleed.guilib.measure.widgets.camerawidgets import CameraViewer
from viperleed.guilib.measure.uimeasurementsettings import SettingsEditor
from viperleed.guilib.measure.datapoints import DataPoints, QuantityInfo
from viperleed.guilib.measure.widgets.measurement_plot import MeasurementPlot
from viperleed.guilib.measure import dialogs

# temporary solution till we have a system config file
from viperleed.guilib import measure as vpr_measure


DEFAULT_CONFIG_PATH = (Path(inspect.getfile(vpr_measure)).parent
                       / 'configuration')

TITLE = 'Measurement UI'


class Measure(gl.ViPErLEEDPluginBase):
    """A class that allows to take measurements."""

    error_occurred = qtc.pyqtSignal(tuple)

    def __init__(self, parent=None):
        """Initialize window."""
        super().__init__(parent, name=TITLE)
        # Keep references to controls, dialogs, and some globals
        self._ctrls = {
            'measure': qtw.QPushButton("Start Measurement"),
            'abort': qtw.QPushButton("Abort"),
            'select': qtw.QComboBox(),
            'energy_input': qtw.QLineEdit(''),
            'set_energy': qtw.QPushButton("Set energy"),
            'settings_editor': qtw.QPushButton("Open settings editor"),
            'menus': {
                'file': qtw.QMenu("&File"),
                'devices': qtw.QMenu("&Devices"),
                'tools': qtw.QMenu("&Tools"),
                },
            }

        self._dialogs = {
            'change_settings': SettingsEditor(),
            'bad_px_finder': dialogs.BadPixelsFinderDialog(parent=self),
            }
        self._glob = {
            'plot': MeasurementPlot()
            }

        # Set window properties
        self.setWindowTitle(TITLE)
        self.setAcceptDrops(True)

        self.__compose()
        self.measurement = None
        self.error_occurred.connect(self.__on_error_occurred)
        self.__errors = []
        self.__error_report_timer = qtc.QTimer(parent=self)
        self.__error_report_timer.setSingleShot(True)
        self.__error_report_timer.timeout.connect(self.__report_errors)

        self.__retry_close = qtc.QTimer(parent=self)
        self.__retry_close.timeout.connect(self.close)

        self.__delayed_start = qtc.QTimer(parent=self)
        self.__delayed_start.setSingleShot(True)
        self.__delayed_start.timeout.connect(self.__run_measurement)
        self.__camera_viewers = []

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

        self._ctrls['select'].addItems(ALL_MEASUREMENTS.keys())
        self._ctrls['select'].setFont(gl.AllGUIFonts().buttonFont)
        self._ctrls['select'].ensurePolished()

        self._ctrls['set_energy'].setFont(gl.AllGUIFonts().buttonFont)
        self._ctrls['set_energy'].ensurePolished()
        self._ctrls['set_energy'].setEnabled(False)

        self._ctrls['settings_editor'].setFont(gl.AllGUIFonts().buttonFont)
        self._ctrls['settings_editor'].ensurePolished()
        self._ctrls['settings_editor'].clicked.connect(
            self.__on_change_settings_pressed
            )
        self._ctrls['settings_editor'].setEnabled(True)

        self._ctrls['energy_input'].setFont(gl.AllGUIFonts().labelFont)
        self._ctrls['energy_input'].ensurePolished()
        self._ctrls['energy_input'].setValidator(QDoubleValidatorNoDot())
        self._ctrls['energy_input'].validator().setLocale(qtc.QLocale.c())

        layout = self.centralWidget().layout()

        layout.addWidget(self._ctrls['measure'], 1, 1, 1, 1)
        layout.addWidget(self._ctrls['abort'], 1, 2, 1, 1)
        layout.addWidget(self._ctrls['select'], 2, 1, 1, 2)
        layout.addWidget(self._ctrls['set_energy'], 3, 1, 1, 1)
        layout.addWidget(self._ctrls['energy_input'], 3, 2, 1, 1)
        layout.addWidget(self._ctrls['settings_editor'], 4, 1, 1, 2)

        self._glob['plot'].show()
        self.statusBar().showMessage('Ready')

        # Menus
        menu = self.menuBar()
        file_menu = self._ctrls['menus']['file']
        menu.insertMenu(self.about_action, file_menu)
        act = file_menu.addAction("&Open...")
        act.setShortcut("Ctrl+O")
        act.triggered.connect(self.__on_read_pressed)

        devices_menu = self._ctrls['menus']['devices']
        menu.insertMenu(self.about_action, devices_menu)
        cam_devices = devices_menu.addMenu("Cameras")
        controller_devices = devices_menu.addMenu("Controllers")
        for camera in get_devices("camera"):
            cam_devices.addAction(camera)

        tools_menu = self._ctrls['menus']['tools']
        menu.insertMenu(self.about_action, tools_menu)
        act = tools_menu.addAction("Find bad pixels...")
        act.triggered.connect(self._dialogs['bad_px_finder'].show)
        self._dialogs['bad_px_finder'].setModal(True)

    def __run_measurement(self):
        self.measurement.begin_measurement_preparation()
        self.__switch_enabled(False)
        self.statusBar().showMessage('Busy')

    def __on_finished(self, *_):
        # After the measurement is done, close the serial ports.
        for controller in self.measurement.controllers:
            controller.disconnect_()
        self.__switch_enabled(True)
        self.statusBar().showMessage('Ready')

    def __switch_enabled(self, idle):
        """Switch enabled status of buttons.

        Parameters
        ----------
        idle : boolean
            False when a measurement is running.
        """
        self._ctrls['measure'].setEnabled(idle)
        self._ctrls['select'].setEnabled(idle)
        self._ctrls['abort'].setEnabled(not idle)
        self._ctrls['settings_editor'].setEnabled(idle)
        self.menuBar().setEnabled(idle)

    def __on_measurement_finished(self, *_):
        print("\n#### DONE! ####")

    def __on_start_pressed(self):
        text = self._ctrls['select'].currentText()
        if self.measurement:
            for thread in self.measurement.threads:
                thread.quit()
            self.__on_finished()
        for viewer in self.__camera_viewers:
            viewer.close()
        self.__camera_viewers = []
        config = configparser.ConfigParser(comment_prefixes='/',
                                           allow_no_value=True,
                                           strict=False)
        file_name = DEFAULT_CONFIG_PATH / 'viperleed_config.ini'
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
        config.set('measurement_settings', 'measurement_class',
            measurement_cls.__name__)
        with open(file_name, 'w') as configfile:
            config.write(configfile)
        self.measurement = measurement_cls(config)

        self.measurement.new_data_available.connect(self.__on_data_received)

        self._ctrls['abort'].clicked.connect(self.measurement.abort)
        self._ctrls['set_energy'].clicked.connect(self.__on_set_energy)
        self.measurement.error_occurred.connect(self.error_occurred)
        for controller in self.measurement.controllers:
            controller.error_occurred.connect(self.error_occurred)
        for camera in self.measurement.cameras:
            camera.error_occurred.connect(self.error_occurred)
            self.__camera_viewers.append(CameraViewer(camera,
                                                      stop_on_close=False,
                                                      roi_visible=False))
        self.measurement.finished.connect(self.__on_finished)
        self.measurement.finished.connect(self.__on_measurement_finished)

        self._glob['plot'].data_points = self.measurement.data_points
        self._glob['plot'].show()
        self.__delayed_start.start(50)

    def __on_data_received(self):
        """Plot measured data."""
        meas = self.sender()
        if not isinstance(meas, MeasurementABC):
            raise RuntimeError(
                f"Got unexpected sender class {meas.__class__.__name__}"
                "for plotting measurements."
                )

        self._glob['plot'].plot_new_data()

    def __on_set_energy(self):
        """Set energy on primary controller."""
        energy = float(self._ctrls['select'].currentText())
        self.measurement.set_LEED_energy((energy,0))

    def __on_error_occurred(self, error_info):
        """React to an error."""
        sender = self.sender()
        self.__errors.append((sender, *error_info))
        self.__error_report_timer.start(35)
        self.__delayed_start.stop()

    def __report_errors(self):
        if not self.__errors:
            return

        err_text = []
        for sender, error_code, error_message in self.__errors:
            if isinstance(sender, CameraABC):
                source = f"camera {sender.name}"
            elif isinstance(sender, ControllerABC):
                source = f"controller at {sender.name}"
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

    def closeEvent(self, event):
        """Reimplement closeEvent to abort measurements as well."""
        if self.measurement and self.measurement.running:
            self.measurement.abort()
            self.__retry_close.start(50)
            event.ignore()
            return
        if self._glob['plot']:
            self._glob['plot'].close()
        # accept has to be called in order to
        # savely quit the dialog and its threads.
        self._dialogs['bad_px_finder'].accept()
        super().closeEvent(event)

    def __on_change_settings_pressed(self):
        settings = SettingsEditor()
        self._dialogs['change_settings'] = settings
        # Change settings is already open
        if settings:
            settings.show()
            # Move window to front
            settings.raise_()
            # Show as a window, also in case it is minimized
            settings.setWindowState(settings.windowState()
                                  & ~qtc.Qt.WindowMinimized
                                  | qtc.Qt.WindowActive)
            return
        settings.show()

    def __on_read_pressed(self):
        csv_name, _ = qtw.QFileDialog.getOpenFileName(parent=self)
        if not csv_name:
            return
        data = DataPoints()
        data.error_occurred.connect(self.error_occurred)
        try:
            data.read_data(csv_name)
        except RuntimeError as err:                                             # TODO: may want to make it also an error_occurred?
            qtw.QMessageBox.critical(self, "Error", str(err),
                                     qtw.QMessageBox.Ok)

        config_name = csv_name.replace(".csv", ".ini")
        config = configparser.ConfigParser(comment_prefixes='/',
                                           allow_no_value=True,
                                           strict=False)
        try:
            f = open(config_name, 'r')
            f.close()
            config.read(config_name)
        except OSError:
            raise FileNotFoundError(f"Could not open/read file: {config_name}."
                                    " Check if this file exists and is in the "
                                    "correct folder.")
            return
        self._glob['plot'].data_points = data
