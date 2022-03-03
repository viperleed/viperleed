"""Module uimeasurement of viperleed.guilib.measure

===============================================
      ViPErLEED Graphical User Interface
===============================================

Created: 2021-10-12
Author: Michele Riva
Author: Florian Doerr

Defines the Measure class.
"""

from pathlib import Path
import inspect
import time

import PyQt5.QtCore as qtc
import PyQt5.QtWidgets as qtw

# ViPErLEED modules
from viperleed import guilib as gl
from viperleed.guilib.measure.measurement import ALL_MEASUREMENTS
from viperleed.guilib.measure import hardwarebase as base
from viperleed.guilib.basewidgets import QDoubleValidatorNoDot

from viperleed.guilib.measure.camera.abc import CameraABC
from viperleed.guilib.measure.controller.abc import ControllerABC
from viperleed.guilib.measure.measurement.abc import MeasurementABC
from viperleed.guilib.measure.widgets.camerawidgets import CameraViewer
from viperleed.guilib.measure.uimeasurementsettings import SettingsEditor
from viperleed.guilib.measure.datapoints import DataPoints
from viperleed.guilib.measure.widgets.measurement_plot import MeasurementPlot
from viperleed.guilib.measure import dialogs
from viperleed.guilib.measure.classes.settings import (
    ViPErLEEDSettings, MissingSettingsFileError
    )

# temporary solution till we have a system config file
from viperleed.guilib import measure as vpr_measure


DEFAULT_CONFIG_PATH = (Path(inspect.getfile(vpr_measure)).parent                # TODO: Get it from system settings
                       / 'configuration')

TITLE = 'Measurement UI'


class UIErrors(base.ViPErLEEDErrorEnum):
    """Class for errors occurring in the UI."""
    FILE_NOT_FOUND_ERROR = (1000, "Could not find file {}.\n{}")
    FILE_UNSUPPORTED = (1001, "Cannot open {}.\n{}")


# too-many-instance-attributes
class Measure(gl.ViPErLEEDPluginBase):
    """A class that allows to take measurements."""

    error_occurred = qtc.pyqtSignal(tuple)
    _dummy_signal = qtc.pyqtSignal()

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

    def __compose(self):  # too-many-statements
        """Place children widgets and menus."""
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

        self.__compose_menu()

    def __compose_menu(self):
        """Put together menu."""
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
        for cam_name, cam_cls in base.get_devices("camera").items():
            act = cam_devices.addAction(cam_name)
            act.setData(cam_cls)
            act.triggered.connect(self.__on_camera_clicked)
        for ctrl_name, ctrl_cls in base.get_devices("controller").items():
            act = controller_devices.addAction(ctrl_name)
            act.setData(ctrl_cls)
            act.triggered.connect(self.__on_controller_clicked)

        tools_menu = self._ctrls['menus']['tools']
        menu.insertMenu(self.about_action, tools_menu)
        act = tools_menu.addAction("Find bad pixels...")
        act.triggered.connect(self._dialogs['bad_px_finder'].show)
        self._dialogs['bad_px_finder'].setModal(True)

    def __run_measurement(self):
        self.measurement.begin_preparation()
        self.__switch_enabled(False)
        self.statusBar().showMessage('Busy')

    def __on_camera_clicked(self, *_):
        cam_name = self.sender().text()
        cam_cls = self.sender().data()
        print(cam_name, cam_cls)

    def __on_controller_clicked(self, *_):
        action = self.sender()
        ctrl_name = action.text()
        ctrl_cls = action.data()
        ctrl_port = ctrl_name.split("(")[1].replace(")", "")
        ctrl = ctrl_cls(port_name=ctrl_port)                                    # TODO: would be better to get the right config file already

        # TEMP FOR VIPERINO ONLY!                                               # TODO: move to a device info dialog (prob. without QThread)
        thread = qtc.QThread()
        ctrl.moveToThread(thread)
        thread.start()
        self._dummy_signal.connect(ctrl.get_hardware)
        self._dummy_signal.emit()
        ctrl.serial.port.waitForReadyRead(100)

        print(ctrl.hardware)
        ctrl.disconnect_()
        thread.quit()
        time.sleep(0.01)

    def __on_finished(self, *_):
        # After the measurement is done, close the serial ports.
        for controller in self.measurement.controllers:
            controller.disconnect_()
        self.__switch_enabled(True)
        self.statusBar().showMessage('Ready')

    def __print_done(self):
        print("\n#### DONE! ####")

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

    def __on_start_pressed(self):
        text = self._ctrls['select'].currentText()
        if self.measurement:
            for thread in self.measurement.threads:
                thread.quit()
            self.__on_finished()
        for viewer in self.__camera_viewers:                                    # Maybe only those relevant for measurement?
            viewer.close()
        self.__camera_viewers = []
        config = ViPErLEEDSettings()
        file_name = DEFAULT_CONFIG_PATH / 'viperleed_config.ini'                # TODO: will be one "default" per measurement type
        try:
            config.read(file_name)
        except MissingSettingsFileError as err:
            base.emit_error(self, UIErrors.FILE_NOT_FOUND_ERROR,
                            file_name, err)
            return

        measurement_cls = ALL_MEASUREMENTS[text]
        config.set('measurement_settings', 'measurement_class',                 # TODO: should eventually not do this!
                   measurement_cls.__name__)
        config.update_file()

        self.measurement = measurement_cls(config)

        self.measurement.new_data_available.connect(self.__on_data_received)

        self._ctrls['abort'].clicked.connect(self.measurement.abort)
        self._ctrls['set_energy'].clicked.connect(self.__on_set_energy)         # NOT NICE! Keeps live the measurement, and appends data at random.
        self.measurement.error_occurred.connect(self.error_occurred)
        for controller in self.measurement.controllers:
            controller.error_occurred.connect(self.error_occurred)
        for camera in self.measurement.cameras:
            camera.error_occurred.connect(self.error_occurred)
            self.__camera_viewers.append(
                CameraViewer(camera, stop_on_close=False, roi_visible=False)
                )
        self.measurement.finished.connect(self.__on_finished)
        self.measurement.finished.connect(self.__print_done)

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
        self.measurement.set_leed_energy(energy, 0, trigger_meas=False)

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
                source = f"controller {sender.name} at {sender.port_name}"
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
            # TODO: Perhaps would be nicer to ask for confirmation
            # rather than always (silently) aborting the measurement
            self.measurement.abort()
            self.__retry_close.start(50)
            event.ignore()
            return
        if self._glob['plot']:
            self._glob['plot'].close()
        # accept has to be called in order to
        # safely quit the dialog and its threads
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
        except RuntimeError as err:
            base.emit_error(self, UIErrors.FILE_UNSUPPORTED, csv_name, err)
            return

        config_name = csv_name.replace(".csv", ".ini")
        config = ViPErLEEDSettings()
        try:
            config.read(config_name)
        except MissingSettingsFileError as err:
            base.emit_error(self, UIErrors.FILE_NOT_FOUND_ERROR,
                            config_name, err)
            return
        self._glob['plot'].data_points = data
