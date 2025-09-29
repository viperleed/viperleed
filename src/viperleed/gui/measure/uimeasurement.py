"""Module uimeasurement of viperleed.gui.measure.

Defines the Measure class, a plug-in for performing LEED(-IV) measurements.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    'Florian Dörr (@FlorianDoerr)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2021-10-12'
__license__ = 'GPLv3+'

# FIXED? camera error should close viewer --> check that no frames can arrive
#        and reopen the viewer. If this is the case, the .close() can be tied
#        to a timer with a small delay.
# FIXED?: measurement camera timeout unless started in live mode before
# FIXED?: measure (timeout), then start viewer --> serial port not open error + wrong camera mode
# FIXED?: QtWarning infinite recursion on err_box.exec_() in _report_errors
#         may be related to the weird errors above?
# TESTME: make sure multiple controllers work fine with the new ABC.stop
# CANT REPRODUCE: badpx file with one line --> not found?? Cannot reproduce with file from J.Mislivicek

# BUG: bad_px_finder: camera can time out. Reproducible in Prague.
#      Did not happen at all to me, nor in Julich (Wuttig). They had
#      fast PC and communication.
#      Seemed to happen very reproducibly especially while trying to
#      acquire the long-dark movie (no frames at all).
#      May be related to inappropriate (too slow) communication speed.
#      Could perhaps be solved if we try halving the frame rate when
#      a timeout occurs (with limits on the number of retries and/or the
#      slowest sensible frame rate).
# BUG: Progress bar "Finding bad pixels" gets reproducibly stuck at 12%
#      in Prague; then timeout error. Cannot reproduce at all here with
#      any of my cameras. Should prepare a debug version with some logging
# BUG: IS: hardware may be so slow that all frames are lost when estimating
#      frame loss. This means that the camera never starts. Possible solution:
#      start another timeout (see how long it should take at max) that restarts
#      the frame rate opt with half the current rate.
# BUG: list_devices makes TPD COMs stuff go crazy.
# ???: viperleed serial, unknown command error misinterpreted? << not sure what this means
# BUG: update COM and camera name when starting measurement (from known devices)
# BUG: measurement start, serial connect failed, attempts to connect three times??
# BUG?: ICCapture open with some settings, settings are retained in viperleed???
# BUG?: IC_SetVideoFormat error -- see 20221123_103752 from Max Buchta -> probably device not open!
# BUG?: (-1/3|1/3) (-1/3|2/3) (-2/3|2/3) not found in beamlist with Max Buchta's correct POSCAR (2022-12-19)

#   G E N E R I C
# TODO: measurement over, restart cameras that were live before?
# TODO: init errors cause obfuscation of the original object that had problems
# TODO: can we use a decorator on class (or __init__) for reporting init errors
#       with delay? The complication is that one needs to call QObject.__init__
#       before using object as parent (e.g., timers) and before signals can be
#       connected. Perhaps we could use some funkiness on QMetaABC.__init__ and
#       .__call__. QMetaABC.__init__(cls, name, bases, dict_, **clsargs) can be
#       used to add a class-level .error_occurred signal; while executing
#       QMetaABC.__call__(cls, *args, **kwargs), super().__call__(*args, **kwargs)
#       returns an __init__ialized instance of cls that can be modified.
# TODO: use __init_subclass__(subclass) in base classes to register all concrete
#       subclasses (check inspect.is_abstract). Could this allow to avoid
#       explicitly adding new subclasses to the imports? Can also be done in
#       QMetaABC by extending __new__(mcs, name, bases, dict_, **clsargs). Probably
#       neither does solve the problem of explicitly importing subclasses, though!
# TODO: fix the documentation in .ini files
# TODO: User event tracking for bug report. Possible with eventFilter on
#       qApp, but see also https://www.qtcentre.org/threads/16182-Logging-Qt-User-Actions
# TODO: error severity in Enum. Non-fatal errors should not stop everything
# TODO: replace all type check complaints with a full_name(object), where
#       full_name does type(object).__module__ (unless None, 'builtins' or
#       AttributeError) + type(object).__qualname__. Possibly falling back
#       to type(object).__name__ in case __qualname__ is missing
# TODO: take inspiration from https://www.the-compiler.org/tmp/qutebrowser_crash.png
#       for a nicer bug-reporting dialog
# TODO: would be very nice to always load the default settings first, and
#       then read in on top of them the ones saved locally. This would allow
#       to seriously reduce the number of lines in many of the "local"
#       configuration files. This also make it possible to have default
#       config files for the same object type, differentiated by some
#       other search criterion.
# TODO: speed up import time by picking more specific module parts
# BUG: PathSelector looks weird on Linux
# BUG: system settings: folder renamed with small/capital not detected as change
# Bug: Low priority bug, about section can be opened multiple times. Use
#       same window in measure and simulation GUI.

#   C A M E R A   &  C O.
# BUG: camera with short exposure consumes an insane amount of memory. Is there a leak?
# BUG: camera viewer & ROI. Something is still not right with the bounds. top-left went
#      to maximum on the 265 when setting a ROI without bad pixels and with properties
#      open.
# BUG: camera lost with viewer & settings open --> cannot open again after power-up
# TODO: Limit frame rate in live view (lag when pixels are saturated)
# TODO: complain if camera lost in live view
# TODO: bad pixels finder top progress bar should scale better, with actual
#       duration of tasks
# TODO: improve progress for preliminary tasks
# TODO: BadPX - Camera without settings file does not prompt to make a new
#       one from defaults
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
# TODO: see if possible to not complain about missing bad pixels
# TODO: moveToThread may be possible if we recreate .driver after we have
#       moved to the new thread (like we do for .serial in ControllerABC).
#       This requires either (i) taking a driver CLASS rather than an
#       instance, and instantiate it (optionally using the extra args/kwargs
#       for it), or (ii) still take an instance, but require its class not
#       take any args/kwargs at __init__.
# TODO: to restore the BadPixelsFinder CameraViewer, see edits in commit
#       7dadce312c775d63487501cb9c5af42b6680df2b. I would need to make
#       CameraABC instances a singleton/Borg (on name?) or find a good
#       way to prevent accessing the same device from different objects.
#       The most difficult part with the singleton/Borg is handling
#       movements to threads.
# TODO?: CameraViewer can still pop up. Try using "with qtc.QSignalBlocker:"
#       while entering the closeEvent? This seems to happen especially when
#       closing the selector with a running camera.
# TODO: Define settings that can be modified without stop-starting the
#       camera. This can speed up some camera startups.
# TODO: find a proper way to CameraABC.moveToThread. This would solve the
#       issue with set_settings holding the UI still.
# TODO: new CameraViewer contains current one + visible controls that
#       allow the same operations as in the context menu.
# TODO: camera returns other data with images. Stuff that comes to mind:
#       max(image); fraction of saturated pixels;
# BUG: exception ignored in ctypes callback: "camera has no exposure"
#      probably a masked AttributeError??
# TODO: bad pixels info not updated on show
# TODO: bad pixels & dark progress bar: text is not always accurate
# TODO: camera extra time should consider that delivery of the
#       first frame in a trigger burst also takes a bit longer
# BUG: Try opening Camera properties with camera lost. OSError? See exceptions_20230119_174155.log
# BUG: imaging source max no. frames should account for the overestimate

#   M E A S U R E M E N T
# TODO: energy ramps are not equivalent for iv == calibration != time_resolved
#       will be solved with the EnergyRamp class
# TODO: Measurement. If primary does not measure, find a better way
#       than sending an empty data_ready for getting the times right
# TODO: Ecal coefficients stored only if user OK with it. Requires plot
#       of calibration residuals.
# TODO: remove default button from Settings dialog
# TODO: filename for measurements: include a progressive number (3 digits?) for each day, ad date_time_NNN

#   H A R D W A R E
# TODO: out to I0, measure HV --> not constant??
# TODO: check that we compare serial number from settings with
#       the one in the hardware when both are available!
# TODO: include datasheet of components
# MUST DO IN 1.0: Add a HWID byte to the info returned by ?
#       for future-proofing to be able to reuse the same protocol
#       and classes for other boxes (e.g., B field compensation)

#   G U I
# TODO: progress bar for non-endless
# TODO: busy dialog where appropriate
# TODO: move_to_front bad-pixels finder dialog when open and click on ui
# TODO: more general than previous one: figure out how to keep a window
#       on top of another set of windows (but NOT with Qt::WindowStaysOnTopHint)
# TODO: consider using a WaitCursor when stuff takes time.
# TODO: QDoubleSpinBox subclass that switches automatically SI unit prefix

#   F E A T U R E S
# TODO: quick IV video to find max intensity, and adjust camera
#       Saturation: oversaturating a bit the very highest maxima gives
#       better S/N on the lower-intensity part of the data, which is most
#       of them. How to choose the overshooting? Perhaps fraction of pixels
#       in saturation?
# TODO: settle time calculation.
#       Look at https://link.springer.com/article/10.1007/s00034-017-0560-3
# TODO: progress report, take inspiration from github.com/verigak/progress/
# QtBUG: should the file dialog issues with network drives appear, try
#        also the QFileDialog::DontUseCustomDirectoryIcons option.
#        Perhaps even QFileDialog::DontResolveSymlinks. The problem was
#        supposed to be solved in Qt 5.15.1, but is still present. It has
#        to do with unreachable network drives, and system calls to populate
#        the list of drives. See bugreports.qt.io/browse/QTBUG-6039


from copy import deepcopy
import functools
from pathlib import Path
import time
from zipfile import ZipFile

import PyQt5.QtCore as qtc
import PyQt5.QtWidgets as qtw

from viperleed.gui.dialogs.errors import DialogDismissedError
from viperleed.gui.measure import hardwarebase as base
from viperleed.gui.measure.camera.abc import CameraABC
from viperleed.gui.measure.classes.decorators import emit_default_faulty
from viperleed.gui.measure.classes.datapoints import DataPoints
from viperleed.gui.measure.classes.settings import DefaultSettingsError
from viperleed.gui.measure.classes.settings import MissingSettingsFileError
from viperleed.gui.measure.classes.settings import NoSettingsError
from viperleed.gui.measure.classes.settings import SystemSettings
from viperleed.gui.measure.classes.settings import ViPErLEEDSettings
from viperleed.gui.measure.controller.abc import ControllerABC
from viperleed.gui.measure.dialogs.badpxfinderdialog import (
    BadPixelsFinderDialog,
    )
from viperleed.gui.measure.dialogs.firmwareupgradedialog import (
    FirmwareUpgradeDialog,
    )
from viperleed.gui.measure.dialogs.new_measurement_dialog import (
    SelectNewMeasurementDialog,
    )
from viperleed.gui.measure.dialogs.settingsdialog import (
    MeasurementSettingsDialog,
    SettingsDialog,
    )
from viperleed.gui.measure.measurement.abc import MeasurementABC
from viperleed.gui.measure.serial.abc import SerialABC
from viperleed.gui.measure.widgets.cameraviewer import CameraViewer
from viperleed.gui.measure.widgets.measurement_plot import MeasurementPlot
from viperleed.gui.pluginsbase import ViPErLEEDPluginBase
from viperleed.gui.widgets.lib import AllGUIFonts
from viperleed.gui.widgets.lib import QDoubleValidatorNoDot
from viperleed.gui.widgets.lib import move_to_front


TITLE = 'Measure LEED-IV'

_QMSG = qtw.QMessageBox
_TIME_CRITICAL = qtc.QThread.TimeCriticalPriority
_UNIQUE = qtc.Qt.UniqueConnection


class UIErrors(base.ViPErLEEDErrorEnum):
    """Class for errors occurring in the UI."""

    FILE_NOT_FOUND_ERROR = (1000, "Could not find file {}.\n{}")
    FILE_UNSUPPORTED = (1001, "Cannot open {}.\n{}")
    RUNTIME_ERROR = (1002, "{}")


# too-many-instance-attributes
class Measure(ViPErLEEDPluginBase):                                             # TODO: Figure out how to inherit error_occurred from QObjectWithError. QObjectMeta hook?
    """A GUI that allows to take measurements."""

    error_occurred = qtc.pyqtSignal(tuple)

    def __init__(self, parent=None):
        """Initialize window."""
        super().__init__(parent, name=TITLE)
        # Keep references to controls, dialogs, and some globals
        self._ctrls = {
            'measure': qtw.QPushButton("New Measurement"),
            'abort': qtw.QPushButton("Abort"),
            'energy_input': qtw.QLineEdit(''),                                  # TODO: QDoubleSpinBox?
            'set_energy': qtw.QPushButton("Set energy"),
            'menus': {
                'file': qtw.QMenu("&File"),
                'devices': qtw.QMenu("&Devices"),
                'tools': qtw.QMenu("&Tools"),
                'views': qtw.QMenu("&View"),
                'sys_settings': qtw.QAction("&Settings"),
                },
            }

        # Parent-less dialogs show up in the task bar. They
        # are set to modal where appropriate (in _compose)
        self._dialogs = {
            'sys_settings':
                SettingsDialog(handled_obj=SystemSettings(),
                               title="System settings"),
            'bad_px_finder': BadPixelsFinderDialog(),
            'camera_viewers': [],
            'error_box': _QMSG(self),                                           # TODO: can look at qtw.QErrorMessage for errors that can be dismissed
            'device_settings': {},     # keys: unique names; No cameras
            'firmware_upgrade': FirmwareUpgradeDialog(self),
            'measurement_selection': SelectNewMeasurementDialog(self),
            # measurement_settings is either None or a
            # MeasurementSettingsDialog created in
            # _create_measurement_settings_dialog.
            'measurement_settings': None,
            }
        self._glob = {
            'plot': MeasurementPlot(),
            'last_dir': self.system_settings.paths["measurements"],
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
            'retry_open_bpx_dialog': qtc.QTimer(parent=self),
            'delay_check_settings': qtc.QTimer(parent=self),
            }

        self.measurement = None
        self._measurement_thread = qtc.QThread()
        self._measurement_thread.start(_TIME_CRITICAL)

        _timer_setup = (
            # key,      interval, single shot
            ('report_errors', 50, True),
            ('retry_close', 50, False),
            ('start_measurement', 50, True),
            ('retry_open_bpx_dialog', 50, True),
            ('delay_check_settings', 5, True),
            )
        for timer, interval, single in _timer_setup:
            self._timers[timer].setSingleShot(single)
            self._timers[timer].setInterval(interval)

        # Set window properties
        self.setWindowTitle(TITLE)
        self.setAcceptDrops(True)

        self._compose()
        self._connect()

        self._timestamps = {'start': -1, 'prepared': -1, 'finished': -1}

    @property
    def system_settings(self):
        """Return a ViPErLEEDSettings with system settings."""
        return self._dialogs['sys_settings'].settings

    def closeEvent(self, event):         # pylint: disable=invalid-name
        """Extend closeEvent to abort measurements as well."""
        if self.measurement and self.measurement.running:
            # TODO: Perhaps would be nicer to ask for confirmation
            # rather than always (silently) aborting the measurement
            self._ctrls['abort'].click()  # aborts measurement
            self._timers['retry_close'].start()
            event.ignore()
            return

        try:
            self.measurement.stop_threads()
        except AttributeError:
            pass

        retry_later = False
        if self._measurement_thread.isRunning():
            self._measurement_thread.quit()
            if not self._measurement_thread.wait(100):
                self._measurement_thread.terminate()
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
        self._dialogs['firmware_upgrade'].close()
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
        self._glob['plot'].show()
        super().showEvent(event)
        if event.spontaneous():  # e.g., Restore after minimize
            return

        # Make sure that the system settings are OK, but use a short
        # timer, as for non-spontaneous events showEvent is called
        # BEFORE the widgets are actually shown.
        self._timers['delay_check_settings'].start()

    def update_device_lists(self):
        """Update entries in "Devices" menu."""
        devices_menu = self._ctrls['menus']['devices']
        cameras, controllers = [a.menu() for a in devices_menu.actions()]
        cameras.clear()
        controllers.clear()

        devices_and_slots = {
            'camera': (cameras, self._on_camera_clicked),
            'controller': (controllers, self._on_controller_clicked),
            }
        for device, (menu, slot) in devices_and_slots.items():
            try:
                detected_devices = self._detect_devices(device)
            except DefaultSettingsError:
                continue
            # The _detect_devices method returns the device name,
            # class and, additional information. The class and
            # additional information are returned as a tuple.
            for device_name, cls_and_info in detected_devices.items():
                act = menu.addAction(device_name)
                act.setData(cls_and_info)
                act.triggered.connect(slot)

        # Leave enabled only those containing entries
        cameras.setEnabled(bool(cameras.actions()))
        controllers.setEnabled(bool(controllers.actions()))

    def _can_take_camera_from_viewer(self, cam_name, viewer):
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
            camera.settings.update_file()

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

    def _check_sys_settings_ok(self):
        """Complain if the system settings are missing entries."""
        if self.system_settings.valid:
            return

        reply = _QMSG.critical(
            self, "Invalid system settings",
            "Missing or invalid system settings. Please fill "
            "in all mandatory (*) fields in the next dialog.",
            _QMSG.Ok | _QMSG.Close
            )
        if reply == _QMSG.Close:
            self.close()
        else:
            self._on_sys_settings_triggered()

    @qtc.pyqtSlot()
    def _cleanup_measurement_settings_dialog(self):
        """Destroy the settings dialog for the current measurement."""
        try:
            self._dialogs['measurement_settings'].deleteLater()
        except AttributeError:
            return  # No dialog yet
        qtw.qApp.processEvents()
        self._dialogs['measurement_settings'] = None

    def _compose(self):
        """Place children widgets and menus."""
        self.setStatusBar(qtw.QStatusBar())
        self.setCentralWidget(qtw.QWidget())
        self.centralWidget().setLayout(qtw.QGridLayout())

        self._ctrls['energy_input'].setValidator(QDoubleValidatorNoDot())
        self._ctrls['energy_input'].validator().setLocale(qtc.QLocale.c())
        self._ctrls['set_energy'].setEnabled(False)
        self._ctrls['energy_input'].setEnabled(False)

        font = AllGUIFonts().buttonFont
        for ctrl in ('measure', 'abort', 'set_energy', 'energy_input'):
            self._ctrls[ctrl].setFont(font)
            self._ctrls[ctrl].ensurePolished()
        self._switch_button_enable(True)

        layout = self.centralWidget().layout()

        layout.addWidget(self._ctrls['measure'], 1, 1, 1, 1)
        layout.addWidget(self._ctrls['abort'], 1, 2, 1, 1)
        layout.addWidget(self._ctrls['set_energy'], 2, 1, 1, 1)
        layout.addWidget(self._ctrls['energy_input'], 2, 2, 1, 1)
        self._compose_menu()
        self._compose_error_box()

        # Take care of dialogs and other windows
        self._dialogs['sys_settings'].setModal(True)
        self._dialogs['measurement_selection'].setModal(True)

    def _compose_error_box(self):
        """Prepare the message box shown when errors happen."""
        err_box = self._dialogs['error_box']
        err_box.setWindowTitle("Error")
        err_box.setTextInteractionFlags(qtc.Qt.TextSelectableByMouse)
        err_box.setIcon(err_box.Critical)

    def _compose_menu(self):
        """Put together menu."""
        menu = self.menuBar()

        # File
        file_menu = self._ctrls['menus']['file']
        menu.insertMenu(self.about_action, file_menu)
        act = file_menu.addAction("&Open...")
        act.setShortcut("Ctrl+O")
        act.triggered.connect(self._on_read_pressed)

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
        act.triggered.connect(self._on_bad_pixels_selected)

        act = tools_menu.addAction("Upload/upgrade firmware...")
        act.setEnabled(True)
        act.triggered.connect(self._dialogs['firmware_upgrade'].open)

        # Views
        views_menu = self._ctrls['menus']['views']
        menu.insertMenu(self.about_action, views_menu)
        act = views_menu.addAction('Show '
                                   + self._glob['plot'].windowTitle().lower())
        act.triggered.connect(self._glob['plot'].show)

        # System settings
        act = self._ctrls['menus']['sys_settings']
        menu.insertAction(self.about_action, act)
        act.triggered.connect(self._on_sys_settings_triggered)

    def _connect(self):
        """Connect signals to appropriate slots."""
        # CONTROLS
        self._ctrls['measure'].clicked.connect(
            self._on_new_measurement_pressed
            )
        self._ctrls['set_energy'].clicked.connect(self._on_set_energy)

        # DIALOGS
        connect_dialogs = (
            ('sys_settings', 'settings_changed',
             self._on_sys_settings_changed),
            ('sys_settings', 'finished', self._check_sys_settings_ok),
            ('bad_px_finder', 'finished',
             functools.partial(self._switch_button_enable, True)),
            ('error_box', 'finished', self._report_errors),
            ('firmware_upgrade', 'error_occurred', self._on_error_occurred),
            ('bad_px_finder', 'error_occurred', self._on_error_occurred),
            ('measurement_selection', 'measurement_selected',
             self._create_measurement),
            ('measurement_selection', 'settings_not_found',
             self._on_measurement_settings_not_found),
            ('measurement_selection', 'rejected',
             self._on_measurement_cancelled),
            )
        for dialog_name, signal_name, slot in connect_dialogs:
            dialog = self._dialogs[dialog_name]
            signal = getattr(dialog, signal_name)
            signal.connect(slot)
        # OTHERS
        self.error_occurred.connect(self._on_error_occurred)
        self._measurement_thread.finished.connect(self._switch_button_enable)

        # TIMERS
        slots = (
            ('report_errors', self._report_errors),
            ('retry_close', self.close),
            ('start_measurement', self._on_measurement_started),
            ('retry_open_bpx_dialog', self._on_bad_pixels_selected),
            ('delay_check_settings', self._check_sys_settings_ok),
            )
        for timer, slot in slots:
            self._timers[timer].timeout.connect(slot)

    def _connect_measurement(self):
        connect = functools.partial(base.safe_connect, type=_UNIQUE)
        measurement = self.measurement
        starter = self._timers['start_measurement'].timeout

        # Errors
        connect(measurement.error_occurred, self._on_error_occurred)
        for device in measurement.devices:
            connect(device.error_occurred, self._on_error_occurred)

        # Measurement events and start/stopping
        connect(measurement.new_data_available, self._on_data_received)
        connect(measurement.prepared, self._on_measurement_prepared)
        connect(measurement.finished, self._on_measurement_finished)
        connect(measurement.finished, self._print_done)
        connect(starter, measurement.start)
        connect(self._ctrls['abort'].clicked, measurement.abort)

        for camera in measurement.cameras:
            self._dialogs['camera_viewers'].append(
                CameraViewer(camera, stop_on_close=False, roi_visible=False,
                             interactions_enabled=False)
                )

    @qtc.pyqtSlot(type, ViPErLEEDSettings)
    def _create_measurement(self, measurement_class, settings):
        """Create measurement object for use in a measurement."""
        self.measurement = measurement_class(settings)
        self.measurement.devices_disconnected.connect(
            self._create_measurement_settings_dialog
            )
        self.measurement.disconnect_devices_and_notify()

    @qtc.pyqtSlot()
    def _create_measurement_settings_dialog(self):
        """Create a new dialog for editing measurement settings."""
        base.safe_disconnect(self.measurement.devices_disconnected,
                             self._create_measurement_settings_dialog)
        dialog = MeasurementSettingsDialog(self.measurement, parent=self)
        dialog.rejected.connect(self._on_measurement_cancelled)
        dialog.accepted.connect(self._on_settings_accepted)
        dialog.error_occurred.connect(self._on_error_occurred)
        dialog.setModal(True)
        self._dialogs['measurement_settings'] = dialog
        dialog.open()

    def _delete_outdated_ctrl_dialog(self, ctrl_name):
        """Remove controller dialogs for ctrl_name."""
        this_ctrl_dialogs = [
            (full_name, dialog)
            for full_name, dialog in self._dialogs['device_settings'].items()
            if full_name.startswith(ctrl_name)
            ]
        for full_name, dialog in this_ctrl_dialogs:
            dialog.reject()  # Also closes the dialog
            dialog.deleteLater()
            del self._dialogs['device_settings'][full_name]

    @emit_default_faulty
    def _detect_devices(self, device_type):
        """Detect and return devices of a certain type."""
        # Notice that self is used by emit_default_faulty.
        return base.get_devices(device_type)

    def _make_ctrl_settings_dialog(self, ctrl_cls, ctrl_info):
        """Make a new settings dialog for a controller."""
        address = ctrl_info.more['address']
        try:
            ctrl = self._make_device(ctrl_cls, ctrl_info, address=address)
        except DefaultSettingsError:
            return
        if not ctrl:
            # Making a controller failed for some reason. We cannot
            # create a SettingsDialog, so we just return. This check is
            # here to prevent the GUI from breaking if the user refuses
            # to make a settings file for a controller without settings
            # when selecting this controller from the devices menu.
            # (See #391)
            return
        dialog = SettingsDialog(ctrl, parent=self)                              # TODO: modal?
        ctrl.ready_to_show_settings.connect(dialog.open)
        dialog.finished.connect(ctrl.disconnect_)
        full_name = ctrl_info.unique_name
        self._dialogs['device_settings'][full_name] = dialog

    @emit_default_faulty
    def _make_device(self, device_cls, settings_info, **other_info):
        """React to the selection of a device."""
        # Find an appropriate settings file, searching in the user
        # configuration folder, and falling back on the base default.
        _cfg_dir = self.system_settings.paths['configuration']
        kwargs = {"directory": _cfg_dir, "parent_widget": self,
                  "third_btn_text": "Create a new settings file"}

        try:
            config = base.get_object_settings(device_cls, settings_info,
                                              **kwargs)
        except NoSettingsError:
            # No settings selected. Will make a new one from defaults.
            config = None
        except DialogDismissedError:
            # Did not find one, and user dismissed the dialog.
            return None

        if config:
            # Found one. Make sure it is in the tree of _cfg_dir,
            # as the user may have selected another folder
            try:
                config.relative_to(_cfg_dir)
            except ValueError:
                # Not in the same tree. Complain, as we expect the
                # configuration tree to contain one default config
                # file per each known device. (Measurements can use
                # device settings from anywhere, though.)
                print("Config not in correct folder tree. "
                      "Consider editing system settings.")                      # TODO
            return device_cls(settings=config, **other_info)

        # Not found, but user wants to make a new one. Use _defaults
        device = device_cls(**other_info)
        if not device.has_valid_settings:
            # Something is wrong with the default configuration file.
            raise DefaultSettingsError

        # Edit the device name in the settings, then save to file
        if issubclass(device_cls, ControllerABC):
            section = "controller"
        elif issubclass(device_cls, CameraABC):
            section = "camera_settings"
        else:
            raise TypeError('Unknown device class detected. Please '
                            'contact the ViPErLEED developers.')

        device_name = (settings_info.more.get('name')
                       or settings_info.unique_name)
        device.settings[section]['device_name'] = device_name
        new_cfg_path = Path(_cfg_dir) / f"{device.name_clean}.ini"
        if new_cfg_path.exists():
            print(f"{section} config file name conflict! Overwriting existing")  # TODO: ask what to do with the (invalid) file
        with new_cfg_path.open('w', encoding='utf-8') as fproxy:
            device.settings.write(fproxy)
        device.uses_default_settings = False
        return device

    def _on_bad_pixels_selected(self):
        """Stop all cameras, then open the dialog."""
        _ok_to_open = True
        for viewer in self._dialogs['camera_viewers']:
            if viewer.camera.is_running:
                _ok_to_open &= viewer.camera.stop()
            # Since we will delete all of the viewers, we can as well
            # make sure they won't show up again after we close them
            viewer.show_auto = False
            viewer.close()
        if not _ok_to_open:
            # Retry soon
            self._timers['retry_open_bpx_dialog'].start()
            return

        for viewer in self._dialogs['camera_viewers']:
            viewer.camera.disconnect_()
        # After disconnecting cameras, it does not really make
        # sense to keep the viewers around: we will not reuse
        # the disconnected cameras anyway
        self._dialogs['camera_viewers'] = []

        self._switch_button_enable(False)
        self._ctrls['abort'].setEnabled(False)
        self._dialogs['bad_px_finder'].open()

    def _on_camera_clicked(self, *_):                                           # TODO: may want to display a busy dialog with "starting camera <name>..."
        cam_name = self.sender().text()

        # Decide whether we can take the camera object
        # (and its settings) from the known camera viewers
        viewers = self._dialogs['camera_viewers']
        for viewer in viewers.copy():
            if self._can_take_camera_from_viewer(cam_name, viewer):
                return

        # Not already available. Make a new camera.
        try:
            camera = self._make_device(*self.sender().data())
        except DefaultSettingsError:
            return

        camera.error_occurred.connect(self._on_camera_error)
        if camera.mode != 'live':
            camera.settings.set('camera_settings', 'mode', 'live')
            camera.settings.update_file()

        viewer = CameraViewer(camera, stop_on_close=True, roi_visible=False)

        # Make sure that the camera is connected, as when we create
        # it from the default file it may not yet have a valid name
        camera.connect_()
        try:
            camera.start()
        except camera.exceptions:
            # Something wrong, but hopefully already
            # reported by the failed camera.start()
            pass
        else:
            viewers.append(viewer)

    @qtc.pyqtSlot(tuple)
    def _on_camera_error(self, error_info):
        """Stop camera and report errors."""
        *_, err_msg = error_info
        if 'bad_pixels' not in err_msg.replace(" ", "_"):
            self.sender().stop()
        self._on_error_occurred(error_info)

    def _on_controller_clicked(self, *_):
        """Show settings of the controller selected."""
        full_name = self.sender().text()
        ctrl_cls, ctrl_info = self.sender().data()

        if full_name not in self._dialogs['device_settings']:
            # Note that ctrl_name may be different from
            # the displayed controller name full_name.
            ctrl_name = ctrl_info.more['name']
            self._delete_outdated_ctrl_dialog(ctrl_name)
            self._make_ctrl_settings_dialog(ctrl_cls, ctrl_info)

        _dialog = self._dialogs['device_settings'].get(full_name, None)
        if not _dialog:
            # Something went wrong with creating the dialog
            return
        if _dialog.isVisible():                                                 # TODO: we should make sure (regularly?) somewhere that the controller is still where it is supposed to be (COM-wise)
            move_to_front(_dialog)
            return
        ctrl = _dialog.handled_object
        ctrl.prepare_to_show_settings()

    @qtc.pyqtSlot(dict)
    def _on_data_received(self, new_data):
        """Plot measured data."""
        plot = self._glob['plot']
        plot.data_points.append(new_data)
        plot.data_points.nr_steps_done += 1
        plot.plot_new_data()

    @qtc.pyqtSlot(tuple)
    def _on_error_occurred(self, error_info):
        """React to an error."""
        sender = self.sender()
        error = (sender, *error_info)
        if error not in set(self._glob['errors']):
            self._glob['errors'].append(error)
        if self._glob['errors']:
            self._timers['report_errors'].start()
        self._timers['start_measurement'].stop()

    @qtc.pyqtSlot()
    def _on_measurement_finished(self, *_):
        """Reset all after a measurement is over."""
        self._timestamps['finished'] = time.perf_counter()
        self.measurement.devices_disconnected.connect(
            self._on_ready_for_next_measurement
            )
        self.measurement.disconnect_devices_and_notify()

    @qtc.pyqtSlot()
    def _on_measurement_prepared(self):
        self._timestamps['prepared'] = time.perf_counter()
        plot_data = self._glob['plot'].data_points
        # After the measurement is prepared, we can expect
        # the total number of steps to have been calculated.
        plot_data.nr_steps_total = self.measurement.data_points.nr_steps_total

    @qtc.pyqtSlot()
    def _on_measurement_started(self):
        self._timestamps['start'] = time.perf_counter()
        base.safe_disconnect(self._measurement_thread.started,
                             self._timers['start_measurement'].start)

    def _on_new_measurement_pressed(self):
        """Prepare to begin a measurement."""
        if self.measurement:
            self._on_measurement_finished()                                     # TODO: necessary?
        for viewer in self._dialogs['camera_viewers']:                          # TODO: Maybe only those relevant for measurement?
            viewer.camera.disconnect_()
            viewer.close()
        self._dialogs['camera_viewers'] = []
        self._switch_button_enable(False)
        qtw.qApp.processEvents()
        self._dialogs['measurement_selection'].cfg_dir = Path(
            self.system_settings['PATHS']['configuration']
            )
        self._dialogs['measurement_selection'].open()

    def _on_read_pressed(self, *_):
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

        if suffix not in ('.csv', '.zip'):                                      # TODO: may make sense to rather look for a (single) .csv file in the folder before erroring out
            base.emit_error(self, UIErrors.FILE_UNSUPPORTED, fname, '')
            return
        if (suffix == '.csv' and not self._read_folder(fname, data, config)):
            return
        if (suffix == '.zip' and not self._read_archive(fname, data, config)):
            return

        self._glob['last_dir'] = str(fname.parent)
        self._glob['plot'].data_points = data
        self._glob['plot'].show()
        self._glob['last_cfg'] = config

    def _on_ready_for_next_measurement(self):
        """Reset controls after a measurement."""
        base.safe_disconnect(self.measurement.devices_disconnected,
                             self._on_ready_for_next_measurement)
        self._switch_button_enable(True)
        for viewer in self._dialogs['camera_viewers']:
            viewer.stop_on_close = True
            viewer.interactions_enabled = True

    def _on_set_energy(self):
        """Set energy on primary controller."""
                                                                                # TODO: implement

    @qtc.pyqtSlot()
    def _on_settings_accepted(self):
        print("\n\nSTARTING\n")
        dialog = self.sender()
        assert dialog is self._dialogs['measurement_settings']
        self._connect_measurement()  # For errors from the measurement.
        settings_ok = self.measurement.set_settings(dialog.settings.last_file)
        self._connect_measurement()  # Again, for the new devices.
        # We will report errors coming from devices separately from the
        # one coming from the measurement. Disconnect the signals from
        # one another here.
        for device in self.measurement.devices:
            base.safe_disconnect(device.error_occurred,
                                 self.measurement.error_occurred)
        # It is necessary to call deleteLater(), otherwise the measurement
        # object will remain alive and the next dialog may attempt to
        # reuse a measurement object for another measurement. This call
        # is done in _cleanup_measurement_settings_dialog().
        self._cleanup_measurement_settings_dialog()
        if not settings_ok:
            return

        self.measurement.moveToThread(self._measurement_thread)
        plot = self._glob['plot']
        plot.data_points = deepcopy(self.measurement.data_points)
        plot.show()
        self._glob['last_cfg'] = self.measurement.settings

        timer = self._timers['start_measurement']
        if not self._measurement_thread.isRunning():
            base.safe_connect(self._measurement_thread.started,
                              timer.start, type=qtc.Qt.UniqueConnection)
            self._measurement_thread.start(_TIME_CRITICAL)
        else:
            timer.start()

        # The next line is approximate: it is only useful in case
        # the measurement errors out before the timer fires
        self._timestamps['start'] = time.perf_counter()

    def _on_sys_settings_changed(self):
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

        # Since older device dialogs may now be still pointing to
        # old configuration files, remove them completely. This
        # will create new settings in the new folder when devices
        # are selected (unless the new folder contains settings).
        for dialog in self._dialogs['device_settings'].values():
            dialog.reject()
            dialog.deleteLater()
        self._dialogs['device_settings'] = {}

        # Do the same for all the camera viewers (and their
        # settings dialogs)
        for viewer in self._dialogs['camera_viewers']:
            viewer.close()
            viewer.deleteLater()
        CameraViewer.clear_cache()
        self._dialogs['camera_viewers'] = []

    def _on_sys_settings_triggered(self):
        """React to a user clicking on 'Settings'."""
        # Update from file, then .open (which updates widgets)
        self.system_settings.read_again()
        self._dialogs['sys_settings'].open()

    @qtc.pyqtSlot()
    def _print_done(self):
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

    def _read_archive(self, fname, datapts, meas_config):
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

        meas_config.read_string(cfg_lines)
        meas_config.base_dir = fname
        cls_name = meas_config['measurement_settings']['measurement_class']
        # Use the information in the config for
        # correctly updating the DataPoints
        datapts.time_resolved = cls_name == "TimeResolved"
        if datapts.is_time_resolved:
            datapts.continuous = meas_config.getboolean('measurement_settings',
                                                        'is_continuous')
        try:
            datapts.read_lines(csv_lines)
        except RuntimeError as err:
            base.emit_error(self, UIErrors.FILE_UNSUPPORTED, fname, err)
            return False
        return True

    def _read_folder(self, csv_name, datapts, meas_config):
        """Read data from a folder containing .csv and .ini."""
        config_name = csv_name.with_suffix(".ini")
        try:
            meas_config.read(config_name)
        except MissingSettingsFileError as err:
            base.emit_error(self, UIErrors.FILE_NOT_FOUND_ERROR,
                            config_name, err)
            return False

        cls_name = meas_config['measurement_settings']['measurement_class']
        datapts.time_resolved = cls_name == 'TimeResolved'
        if datapts.is_time_resolved:
            datapts.continuous = meas_config.getboolean('measurement_settings',
                                                        'is_continuous')
        try:
            datapts.read_csv(csv_name)
        except RuntimeError as err:
            base.emit_error(self, UIErrors.FILE_UNSUPPORTED, csv_name, err)
            return False
        return True

    def _report_errors(self):
        if not self._glob['errors']:
            return

        err_box = self._dialogs['error_box']
        if err_box.isVisible():
            return

        err_text = []
        for sender, error_code, error_message in self._glob['errors']:
            if isinstance(sender, CameraABC):
                source = f'camera {sender.name}'
            elif isinstance(sender, ControllerABC):
                source = f'controller {sender.name} at {sender.address}'
            elif isinstance(sender, MeasurementABC):
                source = f'measurement {type(sender).__name__}'
            elif isinstance(sender, SerialABC):
                # Theoretically we should only receive error
                # messages from controller instances.
                source = f'{type(sender).__name__} at {sender.port_name}'
            elif isinstance(sender, FirmwareUpgradeDialog):
                source = 'firmware upgrade dialog'
            elif isinstance(sender, BadPixelsFinderDialog):
                source = 'bad pixels finder dialog'
            else:
                source = 'system or unknown'

            err_text.append(f'ERROR from {source}\n'
                            f'(Code: {error_code})'
                            f'\n\n{error_message}')

        err_box.setText('\n\n'.join(err_text))
        err_box.open()
        self._glob['errors'] = []

    @qtc.pyqtSlot()
    def _on_measurement_cancelled(self):
        """Handle measurement cancellation."""
        try:
            self.sender().handled_object.stop_threads()
        except AttributeError:
            pass
        self._cleanup_measurement_settings_dialog()
        self._switch_button_enable(True)

    @qtc.pyqtSlot(Path, str)
    def _on_measurement_settings_not_found(self, path, err):
        """Emit an error in response to no settings being detected."""
        base.emit_error(self, UIErrors.FILE_NOT_FOUND_ERROR, path, err)

    def _switch_button_enable(self, idle=False):
        """Switch enabled status of buttons.

        Parameters
        ----------
        idle : bool, optional
            False when a measurement is running. Default is False.
        """
        self._ctrls['measure'].setEnabled(idle)
        self._ctrls['abort'].setEnabled(not idle)
        self.menuBar().setEnabled(idle)
        self.statusBar().showMessage('Ready' if idle else 'Busy')
