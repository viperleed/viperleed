"""Module uimeasurement of viperleed.guilib.measure

===============================================
      ViPErLEED Graphical User Interface
===============================================

Created: 2021-10-12
Author: Michele Riva
Author: Florian Doerr

Defines the Measure class, a plug-in for performing LEED(-IV) measurements.
"""
# FIXED? camera error should close viewer --> check that no frames can arrive
#        and reopen the viewer. If this is the case, the .close() can be tied
#        to a timer with a small delay.
# FIXED?: on archive write, .relative_to; fail on network drive
# FIXED?: camera.ini name brackets
# TESTME: make sure multiple controllers work fine with the new ABC.stop

# BUG: bad_px_finder: camera can time out. Reproducible in Prague.
#      Seemed to happen very reproducibly especially while trying to
#      acquire the long-dark movie (no frames at all).
#      May be related to inappropriate (too slow) communication speed.
#      Could perhaps be solved if we try halving the frame rate when
#      a timeout occurs (with limits on the number of retries and/or the
#      slowest sensible frame rate).
# BUG: IS: hardware may be so slow that all frames are lost when estimating
#      frame loss. This means that the camera never starts. Possible solution:
#      start another timeout (see how long it should take at max) that restarts
#      the frame rate opt with half the current rate.
# BUG: list_devices makes TPD COMs stuff go crazy.
# BUG? ROI increments concern only w/h, but should rather be (w/h - min_w/h) % delta!
# BUG: viperleed serial, unknown command error misinterpreted? << not sure what this means
# BUG: timing estimate: does not include hv_settle_time
# BUG: bad pixels: flat goes on forever for wiggly intensities.
#      Should keep track of the adjustments and decide to pick a value
#      at some point
# BUG: bad pixels: Camera too good?
# BUG: update COM and camera name when starting measurement (from known devices)

#   G E N E R I C
# TODO: init errors cause obfuscation of the original object that had problems
# TODO: can we use a decorator on class (or __init__) for reporting init errors
#       with delay? The complication is that one needs to call QObject.__init__
#       before using object as parent (e.g., timers) and before signals can be
#       connected. Perhaps we could use some funkiness on QMetaABC.__init__ and
#       .__call__. QMetaABC.__init__(cls, name, bases, dict_, **clsargs) can be
#       used to add a class-level .error_occurred signal; while executing
#       QMetaABC.__call__(cls, *args, **kwargs), super().__call__(*args, **kwargs)
#       returns an __init__ialized instance that can be modified.
# TODO: use __init_subclass__(subclass) in base classes to register all concrete
#       subclasses (check inspect.is_abstract). Could this allow to avoid
#       explicitly adding new subclasses to the imports? Can also be done in
#       QMetaABC by extending __new__(mcs, name, bases, dict_, *kwargs). Probably
#       neither does solve the problem of explicitly importing subclasses, though!
# TODO: fix the documentation in .ini files
# TODO: proper _ensure_connected instance method decorator for ControllerABC

#   C A M E R A   &  C O.
# TODO: bad pixels finder top progress bar should scale better, with actual
#       duration of tasks
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
# TODO: auto-scale contrast on camera viewer
# TODO: see if possible to not complain about missing bad pixels

#   M E A S U R E M E N T
# TODO: energy ramps are not equivalent for iv == calibration != time_resolved
#       will be solved with the EnergyRamp class
# TODO: Measurement. If primary does not measure, find a better way
#       than sending an empty data_ready for getting the times right
# TODO: Ecal coefficients stored only if user OK with it

#   H A R D W A R E
# TODO: out to I0, measure HV --> not constant??

#   G U I
# TODO: progress bar for non-endless
# TODO: busy dialog where appropriate

#   F E A T U R E S
# TODO: quick IV video to find max intensity, and adjust camera


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
from viperleed.guilib.measure.classes.datapoints import DataPoints
from viperleed.guilib.measure.widgets.measurement_plot import MeasurementPlot
from viperleed.guilib.measure import dialogs
from viperleed.guilib.measure.classes.settings import (
    ViPErLEEDSettings, MissingSettingsFileError, SYSTEM_CONFIG_PATH
    )


TITLE = 'Measure LEED-IV'

_TIME_CRITICAL = qtc.QThread.TimeCriticalPriority
_QMSG = qtw.QMessageBox


def _get_sys_setting_dialog():
    """Return a SettingsDialog for system settings."""
    _sys_settings = ViPErLEEDSettings()
    _sys_settings.read(SYSTEM_CONFIG_PATH)
    _dialog = dialogs.SettingsDialog(settings=_sys_settings,
                                     title="System settings")
    _dialog.setModal(True)

    # Add some information to the entries
    handler = _dialog.handler
    _infos = (
        ('configuration',
         "<nobr>This is the directory that contains all the</nobr> "
         "configuration files for your devices and measurements. "
         "It must be set before you can run any measurement."),
        ('measurements',
         "<nobr>This is the default folder where all your</nobr> "
         "measurements will be automatically saved. IN THE FUTURE "
         "you will be able to decide if you want to be asked each "
         "time a measurement starts."),
        )
    for key, info in _infos:
        handler['PATHS'][key].set_info_text(info)
    return _dialog


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
            'energy_input': qtw.QLineEdit(''),                                  # TODO: QDoubleSpinBox?
            'set_energy': qtw.QPushButton("Set energy"),
            'menus': {
                'file': qtw.QMenu("&File"),
                'devices': qtw.QMenu("&Devices"),
                'tools': qtw.QMenu("&Tools"),
                'sys_settings': qtw.QAction("&Settings"),
                },
            }

        # Parent-less dialogs show up in the task bar. They
        # are set to modal where appropriate (in __compose)
        self._dialogs = {
            'sys_settings': _get_sys_setting_dialog(),
            'bad_px_finder':
                dialogs.badpxfinderdialog.BadPixelsFinderDialog(),
            'camera_viewers': [],
            'error_box': _QMSG(self),                                           # TODO: can look at qtw.QErrorMessage for errors that can be dismissed
            'device_settings': {},     # keys: unique names; No cameras
            }
        self._glob = {
            'plot': MeasurementPlot(),
            'last_dir': self.system_settings.get("PATHS", "measurements",
                                                 fallback=""),
            # Keep track of the last config used for a measurement.
            # Useful if one wants to repeat a measurement.
            'last_cfg': ViPErLEEDSettings(),
            'errors': [],         # Report a bunch at once
            'n_retry_close': 0,   # Try at most 50 times, i.e., 2.5 sec
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
                                ('start_measurement', 50),):
            self._timers[timer].setSingleShot(True)
            self._timers[timer].setInterval(interval)
        self._timers['retry_close'].setInterval(50)

        # Set window properties
        self.setWindowTitle(TITLE)
        self.setAcceptDrops(True)

        self.__compose()
        self.__connect()

        self._timestamps = {'start': -1, 'prepared': -1, 'finished': -1}

    @property
    def system_settings(self):
        """Return a ViPErLEEDSettings with system settings."""
        return self._dialogs['sys_settings'].settings

    def closeEvent(self, event):         # pylint: disable=invalid-name
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

        # Reject all the others (i.e., discard changes)
        for dialog in self._dialogs['device_settings'].values():
            dialog.reject()

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

        self._dialogs['sys_settings'].close()
        super().closeEvent(event)

    def keyPressEvent(self, event):      # pylint: disable=invalid-name
        """Allow copying (Ctrl+C) device name to clipboard when visible."""
        if (event.key() != qtc.Qt.Key_C
                or not event.modifiers() & qtc.Qt.ControlModifier):
            # Not "Ctrl+C"
            super().keyPressEvent(event)
            return

        # See if any of the "Devices" submenus are open
        menus = [a.menu() for a in self._ctrls['menus']['devices'].actions()]
        for menu in menus:
            if menu.isVisible():
                visible = menu.activeAction()
                break
        else:
            # Nothing visible
            visible = None

        if visible and visible.text():
            qtw.qApp.clipboard().setText(visible.text())
        super().keyPressEvent(event)

    def showEvent(self, event):          # pylint: disable=invalid-name
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

    def __can_take_camera_from_viewer(self, cam_name, viewer):
        """Return whether cam_name can be taken from viewer."""
        camera = viewer.camera
        if cam_name != camera.name:
            return False
        viewer.stop_on_close = True
        cfg = deepcopy(camera.settings)
        cfg.read_again()                                                        # TODO: TEMP: re-read config to apply changes. Will not be done when I have a settings editor.
        cfg['camera_settings']['mode'] = 'live'                                 # TODO: this and the next check will not be useful when we have an editor. Can simply check .mode == 'live'
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
                     'energy_input',):
            self._ctrls[ctrl].setFont(font)
            self._ctrls[ctrl].ensurePolished()
        self.__switch_enabled(True)

        layout = self.centralWidget().layout()

        layout.addWidget(self._ctrls['measure'], 1, 1, 1, 1)
        layout.addWidget(self._ctrls['abort'], 1, 2, 1, 1)
        layout.addWidget(self._ctrls['select'], 2, 1, 1, 2)
        layout.addWidget(self._ctrls['set_energy'], 3, 1, 1, 1)
        layout.addWidget(self._ctrls['energy_input'], 3, 2, 1, 1)

        self._glob['plot'].show()

        self.__compose_menu()
        self.__compose_error_box()

    def __compose_error_box(self):
        """Prepare the message box shown when errors happen."""
        err_box = self._dialogs['error_box']
        err_box.setWindowTitle("Error")
        err_box.setTextInteractionFlags(qtc.Qt.TextSelectableByMouse)
        err_box.setIcon(err_box.Critical)

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

        # System settings
        act = self._ctrls['menus']['sys_settings']
        menu.insertAction(self.about_action, act)
        act.triggered.connect(self.__on_sys_settings_triggered)

    def __connect(self):
        """Connect signals to appropriate slots."""
        # CONTROLS
        self._ctrls['measure'].clicked.connect(self.__on_start_pressed)
        self._ctrls['set_energy'].clicked.connect(self.__on_set_energy)

        # OTHERS
        self.error_occurred.connect(self.__on_error_occurred)
        self._dialogs['sys_settings'].settings_changed.connect(
            self.__on_sys_settings_changed
            )

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
        measurement.error_occurred.connect(self.__on_error_occurred)
        for device in measurement.devices:
            device.error_occurred.connect(self.__on_error_occurred)

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

    def __make_ctrl_settings_dialog(self, ctrl_cls, name, port):
        """Make a new settings dialog for a controller."""
        # Find an appropriate settings file, searching in the default
        # configuration folder, and falling back on the base default
        _cfg_dir = self.system_settings['PATHS']['configuration']
        config = base.get_device_config(name, directory=_cfg_dir,
                                        parent_widget=self)
        if config:
            # Found one. Make sure it is in the tree of _cfg_dir,
            # as the user may have selected another folder
            ctrl = ctrl_cls(settings=config, port_name=port)
            try:
                ctrl.settings.last_file.relative_to(_cfg_dir)
            except ValueError:
                # Not in the same tree. Complain, as we expect the
                # configuration tree to contain one default config
                # file per each known device. (Measurements can use
                # device settings from anywhere, though.)
                print("Config not in correct folder tree. "
                      "Consider editing system settings.")                      # TODO
        else:                                                                   # TODO: ask explicitly for making a file for the device
            # We did not find one. Fall back onto _defaults
            ctrl = ctrl_cls(port_name=port)
            # Check validity of default settings loaded
            if not ctrl.settings or not ctrl.serial:                            # TODO: from this on it could be a more general function (also for cameras)
                # Something is wrong with the default configuration file
                print("SOMETHING WRONG WITH DEFAULT CONFIG")                    # TODO
                return

            # Make necessary edits, then save to file
            ctrl.settings['controller']['device_name'] = name
            new_cfg_path = (Path(_cfg_dir)
                            / (name.replace(' ', '_') + '.ini'))
            if new_cfg_path.exists():
                print("Config file name conflict! "
                      "Existing file will be overwritten")                      # TODO: error
            with new_cfg_path.open('w', encoding='utf-8') as fproxy:
                ctrl.settings.write(fproxy)

        dialog = dialogs.SettingsDialog(ctrl, parent=self)                      # TODO: modal?
        ctrl.ready_to_show_settings.connect(dialog.show)
        dialog.finished.connect(ctrl.disconnect_)
        self._dialogs['device_settings'][name] = dialog

    def __on_camera_clicked(self, *_):                                          # TODO: may want to display a busy dialog with "starting camera <name>..."
        cam_name = self.sender().text()
        cam_cls = self.sender().data()

        _cfg_dir = self.system_settings['PATHS']['configuration']
        cfg_path = base.get_device_config(cam_name,
                                          directory=_cfg_dir,
                                          prompt_if_invalid=False)

        # Decide whether we can take the camera object
        # (and its settings) from the known camera viewers
        viewers = self._dialogs['camera_viewers']
        for viewer in viewers.copy():
            if self.__can_take_camera_from_viewer(cam_name, viewer):
                break
        else:  # Not already available. Make a new camera.
            if not cfg_path:
                print("no config found", cfg_path)                              # TODO: error out here
                return
            cfg = ViPErLEEDSettings.from_settings(cfg_path)
            cfg['camera_settings']['mode'] = 'live'
            if cfg['camera_settings']['device_name'] != cam_name:
                cfg['camera_settings']['device_name'] = cam_name
                cfg.update_file()
            camera = cam_cls(settings=cfg)
            camera.error_occurred.connect(self.__on_error_occurred)
            viewer = CameraViewer(camera, stop_on_close=True,
                                  roi_visible=False)
            viewers.append(viewer)
            camera.start()

    def __on_controller_clicked(self, *_):
        """Show settings of the controller selected."""
        action = self.sender()
        ctrl_name = action.text()
        ctrl_cls = action.data()
        ctrl_name, ctrl_port = (s.strip() for s in ctrl_name.split("("))
        ctrl_port = ctrl_port.replace(")", "")

        if ctrl_name not in self._dialogs['device_settings']:
            self.__make_ctrl_settings_dialog(ctrl_cls, ctrl_name, ctrl_port)

        _dialog = self._dialogs['device_settings'].get(ctrl_name, None)
        if not _dialog:
            # Something went wrong with creating the dialog
            return

        if _dialog.isVisible():                                                 # TODO: we should make sure (regularly?) somewhere that the controller is still where it is supposed to be (COM-wise)
            move_to_front(_dialog)
            return

        ctrl = _dialog.handled_object
        if ctrl_port != ctrl.port_name:
            ctrl.port_name = ctrl_port
            ctrl.settings.update_file()
        ctrl.prepare_to_show_settings()

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
        error = (sender, *error_info)
        if error not in set(self._glob['errors']):
            self._glob['errors'].append(error)
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
        self.__switch_enabled(False)
        qtw.qApp.processEvents()

        config = ViPErLEEDSettings()                                            # TODO: should decide whether to use the 'last_cfg' or the default here! Probably open dialog.
        _cfg_dir = Path(self.system_settings['PATHS']['configuration'])
        file_name = _cfg_dir.resolve() / 'measurement.ini'                      # TODO: will be one "default" per measurement type
        try:
            config.read(file_name)
        except MissingSettingsFileError as err:
            base.emit_error(self, UIErrors.FILE_NOT_FOUND_ERROR,
                            file_name, err)
            return

        measurement_cls = ALL_MEASUREMENTS[text]
        config.set('measurement_settings', 'measurement_class',                 # TODO: should eventually not do this!
                   measurement_cls.__name__)

        self.measurement = measurement_cls(config)
        self.measurement.moveToThread(self.__measurement_thread)

        self.__connect_measurement()

        plot = self._glob['plot']
        plot.data_points = deepcopy(self.measurement.data_points)
        plot.show()
        self._glob['last_cfg'] = self.measurement.settings

        timer = self._timers['start_measurement']
        if not self.__measurement_thread.isRunning():
            base.safe_connect(self.__measurement_thread.started,
                              timer.start, qct.Qt.UniqueConnection)
            self.__measurement_thread.start(_TIME_CRITICAL)
        else:
            timer.start()

    def __on_sys_settings_changed(self):
        """Save settings file to disk when system settings change."""
        _dialog = self._dialogs['sys_settings']

        # Check that the paths exist, ask to create them if not
        handler = _dialog.handler
        for option in handler['PATHS'].values():
            path_widg = option.handler_widget
            if path_widg.path and path_widg.path.exists():
                continue
            _reply = _QMSG.question(
                _dialog, "Directory does not exist",
                f"'{option.option_name.title()}' directory "
                "does not exist. Would you like to create it?",
                _QMSG.Yes | _QMSG.No
                )
            if _reply == _QMSG.Yes:
                path_widg.path.mkdir(parents=True)
        self.system_settings.update_file()

        # Since older device dialogs may now be still pointing to               # TODO: do the same for cameras
        # old configuration files, remove them completely. This
        # will create new settings in the new folder when devices
        # are selected (unless the new folder contains settings).
        for dialog in self._dialogs['device_settings'].values():
            dialog.reject()
            dialog.deleteLater()
        self._dialogs['device_settings'] = {}

    def __on_sys_settings_triggered(self):
        """React to a user clicking on 'Settings'."""
        # Update from file, then .show (which updates widgets)
        self.system_settings.read(SYSTEM_CONFIG_PATH)
        self._dialogs['sys_settings'].show()

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

            err_text.append(f"ERROR from {source}\n"
                            f"(Code: {error_code})"
                            f"\n\n{error_message}")

        err_box = self._dialogs['error_box']
        err_box.setText("\n\n".join(err_text))
        err_box.exec_()
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
        self.menuBar().setEnabled(idle)
        self.statusBar().showMessage('Ready' if idle else 'Busy')
