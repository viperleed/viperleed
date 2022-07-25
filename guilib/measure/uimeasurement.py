"""Module uimeasurement of viperleed.guilib.measure

===============================================
      ViPErLEED Graphical User Interface
===============================================

Created: 2021-10-12
Author: Michele Riva
Author: Florian Doerr

Defines the Measure class.
"""

# TODO: it looks like the largest fraction of the time required by
#       abort_trigger_burst is in fact on the call to _dll_start_live
#       (97% of the 80ms), while .pause() takes only a short time
#       (approx. 3%, i.e., 2.4ms). One could in principle then only
#       .pause() during abort_trigger_burst, and restart the camera
#       at the beginning of the next energy step. It is not that trivial
#       to make it general and thread-safe: one has to decide based on
#       how long the camera_delay (and perhaps i0_settle_time) are whether
#       cameras should be started before setting the energy (if too short)
#       or while setting the energy (if it takes long enough). The additional
#       issue is that one would need to know in advance how long does it
#       take for a camera to restart after being paused, as this determines
#       what 'too short/long' means. It could be potentially an information
#       to be gathered during preparation.
# TODO: progress bar for non-endless
# TODO: quick IV video to find max intensity, and adjust camera
# TODO: auto-scale contrast on camera viewer
# TODO: busy dialog where appropriate
# TODO: fix the documentation in .ini files
# TODO: can we use a decorator on class (or __init__) for reporting init errors
#       with delay? The complication is that one needs to call QObject.__init__
#       before using object as parent (e.g., timers) and before signals can be
#       connected.

from copy import deepcopy
from pathlib import Path
from zipfile import ZipFile
import time

import PyQt5.QtCore as qtc
import PyQt5.QtWidgets as qtw

# ViPErLEED modules
from viperleed.guilib.pluginsbase import ViPErLEEDPluginBase
from viperleed.guilib.widgetslib import move_to_front, AllGUIFonts
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
    ViPErLEEDSettings, MissingSettingsFileError, get_system_config
    )


TITLE = 'Measure LEED-IV'

SYS_CFG = get_system_config()
DEFAULT_CONFIG_PATH = Path(
    SYS_CFG.get("PATHS", 'configuration', fallback='')
    )

_TIME_CRITICAL = qtc.QThread.TimeCriticalPriority


class UIErrors(base.ViPErLEEDErrorEnum):
    """Class for errors occurring in the UI."""
    FILE_NOT_FOUND_ERROR = (1000, "Could not find file {}.\n{}")
    FILE_UNSUPPORTED = (1001, "Cannot open {}.\n{}")
    RUNTIME_ERROR = (1002, "{}")


# too-many-instance-attributes
class Measure(ViPErLEEDPluginBase):
    """A GUI that allows to take measurements."""

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

        # Parent-less dialogs show up in the task bar. They
        # are set to modal where appropriate (in __compose)
        self._dialogs = {
            'change_settings': SettingsEditor(),
            'bad_px_finder': dialogs.BadPixelsFinderDialog(),
            'camera_viewers' : []
            }
        self._glob = {
            'plot': MeasurementPlot(),
            'last_dir': SYS_CFG.get("PATHS", "measurements", fallback=""),
            # Keep track of the last config used for a measurement.
            # Useful if one wants to repeat a measurement.
            'last_cfg': ViPErLEEDSettings(),
            'errors': [],       # Report a bunch at once
            'n_retry_close': 0  # Try at most 50 times, i.e., 2.5 sec
            }
        self._timers = {
            'report_errors': qtc.QTimer(parent=self),
            'retry_close': qtc.QTimer(parent=self),
            'start_measurement': qtc.QTimer(parent=self),
            }

        self.measurement = None
        self.__measurement_thread = qtc.QThread()
        self.__measurement_thread.start(_TIME_CRITICAL)

        for timer, interval in (('report_errors', 35),
                                ('start_measurement', 50)):
            self._timers[timer].setSingleShot(True)
            self._timers[timer].setInterval(interval)
        self._timers['retry_close'].setInterval(50)

        # Set window properties
        self.setWindowTitle(TITLE)
        self.setAcceptDrops(True)

        self.__compose()
        self.__connect()

        self._timestamps = {'start': -1, 'prepared': -1, 'finished': -1}

    def closeEvent(self, event):  # pylint: disable=invalid-name
        """Reimplement closeEvent to abort measurements as well."""
        if self.measurement and self.measurement.running:
            # TODO: Perhaps would be nicer to ask for confirmation
            # rather than always (silently) aborting the measurement
            self._ctrls['abort'].click()  # aborts measurement
            self._timers['retry_close'].start()
            event.ignore()
            return

        retry_later = False
        if self.__measurement_thread.isRunning():
            self.__measurement_thread.quit()
            retry_later = True

        if self._glob['plot']:
            self._glob['plot'].close()
        # accept has to be called in order to
        # safely quit the dialog and its threads
        self._dialogs['bad_px_finder'].accept()

        # Stop all cameras
        for viewer in self._dialogs['camera_viewers']:
            viewer.close()
            camera = viewer.camera
            if camera and camera.is_running:
                print(camera.name, "is running")
                camera.stop()
                retry_later = True

        if retry_later and self._glob['n_retry_close'] <= 50:
            self._glob['n_retry_close'] += 1
            self._timers['retry_close'].start()
            event.ignore()
            return

        super().closeEvent(event)

    def showEvent(self, event):  # pylint: disable=invalid-name
        """Show self."""
        self._glob['n_retry_close'] = 0
        super().showEvent(event)

    def update_device_lists(self):
        """Update entries in "Devices" menu."""
        menu = self._ctrls['menus']['devices']
        cameras, controllers = [a.menu() for a in menu.actions()]
        cameras.clear()
        controllers.clear()

        for cam_name, cam_cls in base.get_devices("camera").items():
            act = cameras.addAction(cam_name)
            act.setData(cam_cls)
            act.triggered.connect(self.__on_camera_clicked)
        for ctrl_name, ctrl_cls in base.get_devices("controller").items():
            act = controllers.addAction(ctrl_name)
            act.setData(ctrl_cls)
            act.triggered.connect(self.__on_controller_clicked)

        # Leave enabled only those containing entries
        cameras.setEnabled(bool(cameras.actions()))
        controllers.setEnabled(bool(controllers.actions()))

    def __can_take_camera_from_viewer(self, cam_name, viewer, cfg_path):
        """Return whether cam_name can be taken from viewer."""
        camera = viewer.camera
        if cam_name != camera.name:
            return False
        viewer.stop_on_close = True
        cfg = deepcopy(camera.settings)
        if cfg_path:                                # TODO: TEMP: re-read config to apply changes. Will not be done when I have a settings editor.
            cfg.read(cfg_path)
        cfg['camera_settings']['mode'] = 'live'     # TODO: this and the next check will not be useful when we have an editor. Can simply check .mode == 'live'
        if cfg != camera.settings:
            if not camera.connected:
                # There's no speed advantage over recreating the
                # camera object with the settings (and its viewer)
                viewer.close()
                self._dialogs['camera_viewers'].remove(viewer)
                return False
            camera.settings = cfg

        # Try starting it. If not possible, we probably lost it
        if not camera.is_running:
            try:
                camera.start()
            except camera.exceptions:
                # Probably lost the device
                viewer.close()
                self._dialogs['camera_viewers'].remove(viewer)
                return False
        if viewer.isVisible():
            move_to_front(viewer)
        return True

    def __compose(self):
        """Place children widgets and menus."""
        self.setStatusBar(qtw.QStatusBar())
        self.setCentralWidget(qtw.QWidget())
        self.centralWidget().setLayout(qtw.QGridLayout())

        self._ctrls['select'].addItems(ALL_MEASUREMENTS.keys())
        self._ctrls['energy_input'].setValidator(QDoubleValidatorNoDot())
        self._ctrls['energy_input'].validator().setLocale(qtc.QLocale.c())
        self._ctrls['set_energy'].setEnabled(False)
        self._ctrls['energy_input'].setEnabled(False)

        font = AllGUIFonts().buttonFont
        for ctrl in ('measure', 'abort', 'select', 'set_energy',
                     'energy_input', 'settings_editor',):
            self._ctrls[ctrl].setFont(font)
            self._ctrls[ctrl].ensurePolished()
        self.__switch_enabled(True)

        layout = self.centralWidget().layout()

        layout.addWidget(self._ctrls['measure'], 1, 1, 1, 1)
        layout.addWidget(self._ctrls['abort'], 1, 2, 1, 1)
        layout.addWidget(self._ctrls['select'], 2, 1, 1, 2)
        layout.addWidget(self._ctrls['set_energy'], 3, 1, 1, 1)
        layout.addWidget(self._ctrls['energy_input'], 3, 2, 1, 1)
        layout.addWidget(self._ctrls['settings_editor'], 4, 1, 1, 2)

        self._glob['plot'].show()

        self.__compose_menu()

    def __compose_menu(self):
        """Put together menu."""
        menu = self.menuBar()

        # File
        file_menu = self._ctrls['menus']['file']
        menu.insertMenu(self.about_action, file_menu)
        act = file_menu.addAction("&Open...")
        act.setShortcut("Ctrl+O")
        act.triggered.connect(self.__on_read_pressed)

        # Devices
        devices_menu = self._ctrls['menus']['devices']                          # TODO: have to update the lists regularly. Use timer to update_device_lists.
        devices_menu.aboutToShow.connect(self.update_device_lists)
        menu.insertMenu(self.about_action, devices_menu)
        devices_menu.addMenu("Cameras")
        devices_menu.addMenu("Controllers")
        self.update_device_lists()

        # Tools
        tools_menu = self._ctrls['menus']['tools']
        menu.insertMenu(self.about_action, tools_menu)
        act = tools_menu.addAction("Find bad pixels...")
        act.triggered.connect(self._dialogs['bad_px_finder'].show)
        self._dialogs['bad_px_finder'].setModal(True)

        act = tools_menu.addAction("Upload/upgrade firmware...")
        act.setEnabled(False)                                                   # TODO: fix when implemented
        # act.triggered.connect(self._dialogs['firmware_upgrade'].show)

    def __connect(self):
        """Connect signals to appropriate slots."""
        # CONTROLS
        self._ctrls['measure'].clicked.connect(self.__on_start_pressed)
        self._ctrls['settings_editor'].clicked.connect(
            self.__on_change_settings_pressed
            )
        self._ctrls['set_energy'].clicked.connect(self.__on_set_energy)

        # OTHERS
        self.error_occurred.connect(self.__on_error_occurred)

        # TIMERS
        slots = (
            ('report_errors', self.__report_errors),
            ('retry_close', self.close),
            ('start_measurement', self.__on_measurement_started),
            )
        for timer, slot in slots:
            self._timers[timer].timeout.connect(slot)

    def __connect_measurement(self):
        measurement = self.measurement
        starter = self._timers['start_measurement'].timeout

        # Errors
        measurement.error_occurred.connect(self.error_occurred)
        for device in measurement.devices:
            device.error_occurred.connect(self.error_occurred)

        # Measurement events and start/stopping
        measurement.new_data_available.connect(self.__on_data_received)
        measurement.prepared.connect(self.__on_measurement_prepared)
        measurement.finished.connect(self.__on_measurement_finished)
        measurement.finished.connect(self.__print_done)
        starter.connect(measurement.begin_preparation)
        self._ctrls['abort'].clicked.connect(measurement.abort)

        for camera in measurement.cameras:
            self._dialogs['camera_viewers'].append(
                CameraViewer(camera, stop_on_close=False, roi_visible=False)
                )

    def __on_camera_clicked(self, *_):                                          # TODO: may want to display a busy dialog with "starting camera <name>..."
        cam_name = self.sender().text()
        cam_cls = self.sender().data()

        cfg_path = base.get_device_config(cam_name,
                                          directory=DEFAULT_CONFIG_PATH,
                                          prompt_if_invalid=False)

        # Decide whether we can take the camera object
        # (and its settings) from the known camera viewers
        viewers = self._dialogs['camera_viewers']
        for viewer in viewers.copy():
            if self.__can_take_camera_from_viewer(cam_name, viewer, cfg_path):
                break
        else:  # Not already available. Make a new camera.
            if not cfg_path:
                print("no config found", cfg_path)                              # TODO: error out here
                return
            cfg = ViPErLEEDSettings.from_settings(cfg_path)
            cfg['camera_settings']['mode'] = 'live'
            camera = cam_cls(settings=cfg)
            camera.error_occurred.connect(self.error_occurred)
            viewer = CameraViewer(camera, stop_on_close=True,
                                  roi_visible=False)
            viewers.append(viewer)
            camera.start()

    def __on_change_settings_pressed(self):
        settings = SettingsEditor()
        self._dialogs['change_settings'] = settings
        # Change settings is already open
        if settings:
            move_to_front(settings)
            return
        settings.show()

    def __on_controller_clicked(self, *_):
        action = self.sender()
        ctrl_name = action.text()
        ctrl_cls = action.data()
        ctrl_port = ctrl_name.split("(")[1].replace(")", "")
        ctrl = ctrl_cls(port_name=ctrl_port)                                    # TODO: would be better to get the right config file already

        # TEMP FOR VIPERINO ONLY!                                               # TODO: move to a device info dialog (prob. without QThread)
        ctrl.get_hardware()
        ctrl.serial.port.waitForReadyRead(100)  # this always times out, irrespective of how long (but the response arrived and was processed correctly!)
        ctrl.disconnect_()

        # TODO: here would be a good place to check if the serial
        # number in the name and the one in ctrl.hardware match.
        # If not, issue an error and update the device list.
        print(ctrl.hardware)

    @qtc.pyqtSlot(dict)
    def __on_data_received(self, new_data):
        """Plot measured data."""
        plot = self._glob['plot']
        plot.data_points.append(new_data)
        plot.data_points.nr_steps_done += 1
        plot.plot_new_data()

    @qtc.pyqtSlot(tuple)
    def __on_error_occurred(self, error_info):
        """React to an error."""
        sender = self.sender()
        self._glob['errors'].append((sender, *error_info))
        self._timers['report_errors'].start()
        self._timers['start_measurement'].stop()

    @qtc.pyqtSlot()
    def __on_measurement_finished(self, *_):
        """Reset all after a measurement is over."""
        self._timestamps['finished'] = time.perf_counter()
        for controller in self.measurement.controllers:
            controller.disconnect_()
        self.__switch_enabled(True)

    @qtc.pyqtSlot()
    def __on_measurement_prepared(self):
        self._timestamps['prepared'] = time.perf_counter()

    @qtc.pyqtSlot()
    def __on_measurement_started(self):
        self._timestamps['start'] = time.perf_counter()
        base.safe_disconnect(self.__measurement_thread.started,
                             self._timers['start_measurement'].start)
        self.__switch_enabled(False)

    def __on_read_pressed(self, *_):
        """Read data from a measurement file."""
        fname, _ = qtw.QFileDialog.getOpenFileName(
            parent=self, directory=self._glob['last_dir']
            )
        if not fname:
            return

        data = DataPoints()
        data.error_occurred.connect(self.error_occurred)
        config = ViPErLEEDSettings()
        fname = Path(fname)
        suffix = fname.suffix

        if suffix not in ('.csv', '.zip'):                                        # TODO: may make sense to rather look for a (single) .csv file in the folder before erroring out
            base.emit_error(self, UIErrors.FILE_UNSUPPORTED, fname, '')
            return
        if (suffix == '.csv' and not self.__read_folder(fname, data, config)):
            return
        if (suffix == '.zip' and not self.__read_archive(fname, data, config)):
            return

        self._glob['last_dir'] = str(fname.parent)
        self._glob['plot'].data_points = data
        self._glob['last_cfg'] = config

        # Select the measurement type that was just loaded.
        # NB: we do not actually instantiate the new measurement
        # because it can take time to handle the communication
        # with the devices. It makes sense to do this only once
        # the user actually wants to start a measurement.
        cls_name = config['measurement_settings']['measurement_class']
        meas_dict = {c.__name__: t for t, c in ALL_MEASUREMENTS.items()}
        self._ctrls['select'].setCurrentText(meas_dict[cls_name])

        # And use the information in the config for
        # correctly updating the DataPoints
        data.time_resolved = (cls_name == "TimeResolved")
        if data.is_time_resolved:
            data.continuous = config.getboolean('measurement_settings',
                                                'is_continuous')

    def __on_set_energy(self):
        """Set energy on primary controller."""
        energy = float(self._ctrls['select'].currentText())
        # self.measurement.set_leed_energy(energy, 0, trigger_meas=False)       # NOT NICE! Keeps live the measurement, and appends data at random.

    def __on_start_pressed(self):
        """Prepare to begin a measurement."""
        print("\n\nSTARTING\n")
        text = self._ctrls['select'].currentText()
        if self.measurement:
            self.__on_measurement_finished()                                    # TODO: necessary?
        for viewer in self._dialogs['camera_viewers']:                          # TODO: Maybe only those relevant for measurement?
            viewer.camera.disconnect_()
            viewer.close()
        self._dialogs['camera_viewers'] = []

        config = ViPErLEEDSettings()                                            # TODO: should decide whether to use the 'last_cfg' or the default here! Probably open dialog.
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
        self.measurement.moveToThread(self.__measurement_thread)

        self.__connect_measurement()

        plot = self._glob['plot']
        plot.data_points = plot_data = DataPoints()
        plot_data.primary_controller = self.measurement.primary_controller
        plot.show()
        self._glob['last_cfg'] = self.measurement.settings

        timer = self._timers['start_measurement']
        if not self.__measurement_thread.isRunning():
            base.safe_connect(self.__measurement_thread.started,
                              timer.start, qct.Qt.UniqueConnection)
            self.__measurement_thread.start(_TIME_CRITICAL)
        else:
            timer.start()

    @qtc.pyqtSlot()
    def __print_done(self):
        print("\n#### DONE! ####")
        start, prep, finish = self._timestamps.values()
        n_steps = self.measurement.data_points.nr_steps_done
        txt = f"Measurement took {finish-start:.2f} s"
        if prep-start > 0:
            txt += f", of which {prep-start:.2f} s for preparation. "
        else:
            txt += ". "
        if n_steps:
            txt += (f"This is on average {1000*(finish-prep)/n_steps:.2f} ms"
                    " per energy step")
        print(txt, "\n")

    def __read_archive(self, fname, datapts, meas_config):
        """Read from a .zip archive containing .csv and .ini."""
        fname = Path(fname)
        with ZipFile(fname, 'r') as arch:
            try:
                csv_lines = arch.read("measurement.csv").decode().split('\n')
            except KeyError as err:
                # not found
                base.emit_error(self, UIErrors.FILE_NOT_FOUND_ERROR,
                                fname / "measurement.csv", err)
                return False

            try:
                cfg_lines = arch.read("measurement.ini").decode()
            except KeyError as err:
                # not found
                base.emit_error(self, UIErrors.FILE_NOT_FOUND_ERROR,
                                fname / "measurement.ini", err)
                return False
        try:
            datapts.read_lines(csv_lines)
        except RuntimeError as err:
            base.emit_error(self, UIErrors.FILE_UNSUPPORTED, fname, err)
            return False

        meas_config.read_string(cfg_lines)
        meas_config.base_dir = fname
        return True

    def __read_folder(self, csv_name, datapts, meas_config):
        """Read data from a folder containing .csv and .ini."""
        try:
            datapts.read_csv(csv_name)
        except RuntimeError as err:
            base.emit_error(self, UIErrors.FILE_UNSUPPORTED, csv_name, err)
            return False

        config_name = csv_name.replace(".csv", ".ini")
        try:
            meas_config.read(config_name)
        except MissingSettingsFileError as err:
            base.emit_error(self, UIErrors.FILE_NOT_FOUND_ERROR,
                            config_name, err)
            return False
        return True

    def __report_errors(self):
        if not self._glob['errors']:
            return

        err_text = []
        for sender, error_code, error_message in self._glob['errors']:
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
        self._glob['errors'] = []

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
        self.statusBar().showMessage('Ready' if idle else 'Busy')
