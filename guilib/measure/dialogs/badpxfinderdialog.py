"""Module badpxfinderdialog of viperleed.guilib.measure.dialogs.

========================================
   ViPErLEED Graphical User Interface
========================================

Created: 2021-11-22
Author: Michele Riva

Defines the BadPixelsFinderDialog class that handles user interaction
while finding bad pixels for a camera.
"""

import inspect
from pathlib import Path

from PyQt5 import (QtWidgets as qtw,
                   QtCore as qtc,
                   QtGui as qtg)

from viperleed.guilib.measure.hardwarebase import get_devices
from viperleed.guilib.measure.camera import badpixels
from viperleed.guilib.measure.camera import abc as camera_abc
from viperleed.guilib.widgetslib import change_control_text_color

# temporary solution till we have a system config file
from viperleed.guilib import measure as vpr_measure


DEFAULT_CONFIG_PATH = (Path(inspect.getfile(vpr_measure)).parent
                       / 'configuration')


class BadPixelsFinderDialog(qtw.QDialog):
    """Dialog to handle user interaction when finding bad pixels."""

    __start_finder = qtc.pyqtSignal()
    __abort_finder = qtc.pyqtSignal()
    __set_camera_settings = qtc.pyqtSignal(object)

    def __init__(self, parent=None):
        """Initialize dialog."""
        super().__init__(parent=parent)

        self.__ctrls = {
            'camera': qtw.QComboBox(),
            'bad_px_path': qtw.QLabel("None selected"),
            }
        self.__progress = {
            'group' : qtw.QGroupBox("Progress"),
            'total': qtw.QProgressBar(),
            'section_text': qtw.QLabel(),
            'section': qtw.QProgressBar(),
            }
        self.__buttons = {
            'done': qtw.QPushButton("&Done"),
            'abort': qtw.QPushButton("&Abort"),
            'find': qtw.QPushButton("&Find"),
            'set_bad_px_path': qtw.QToolButton(),
            }
        self.__bad_px_info = {
            'date_time': qtw.QLabel("Date/Time: \u2014"),
            'n_bad': qtw.QLabel("No. bad pixels: \u2014"),
            'n_uncorrectable': qtw.QLabel("No. uncorrectable: \u2014")
            }
        self.__timers = {
            # 'update_list' will update the list of camera
            # devices every 2 seconds
            'update_list': (qtc.QTimer(self), 2000),
            # 'process_events' will force processing of currently
            # non-handled GUI events to improve response time
            # of the UI. This is necessary because the camera
            # may take a while to start up, and makes the UI
            # unresponsive.
            'process_events': (qtc.QTimer(self), 10)
            }

        timer, interval = self.__timers['process_events']
        timer.timeout.connect(qtw.qApp.processEvents)
        timer.start(interval)

        self.__available_cameras = {}
        self.active_camera = None
        self.__finder = None
        self.__finder_thread = qtc.QThread()
        self.__finder_thread.start()

        self.setWindowTitle("Find bad pixels")

        self.__compose()
        self.__connect()

    def showEvent(self, event):
        """Extend showEvent to update controls."""
        self.__reset_progress_bars()
        self.__progress['group'].hide()
        timer, interval = self.__timers['process_events']
        timer.timeout.connect(qtw.qApp.processEvents)
        timer.start(interval)
        self.update_available_camera_list()
        self.adjustSize()
        super().showEvent(event)

    def accept(self):
        """Clean up, then accept."""
        self.__clean_up()
        super().accept()

    def reject(self):
        """Clean up, then reject."""
        self.__clean_up()
        super().reject()

    def update_available_camera_list(self, *_):
        """Update the list of available cameras."""
        camera_combo = self.__ctrls['camera']
        old_selection = camera_combo.currentText()
        self.__available_cameras = get_devices('camera')
        old_items = set(camera_combo.itemText(i)
                        for i in range(camera_combo.count()))
        if old_items != set(self.__available_cameras.keys()):
            # Camera list has changed
            camera_combo.clear()
            camera_combo.addItems(self.__available_cameras.keys())
            if old_selection:
                camera_combo.setCurrentText(old_selection)
        self.__enable_controls(True)

    def __abort(self, *_):
        """Abort bad-pixel-finder routine."""
        self.__abort_finder.emit()
        self.__reset_progress_bars()
        self.__progress['group'].hide()
        self.__enable_controls(True)
        self.adjustSize()

    def __clean_up(self):
        """Clean up self before closing."""
        self.__stop_timers()
        self.__ctrls['camera'].clear()
        self.active_camera = None
        if self.__finder_thread.isRunning():
            self.__finder_thread.quit()

    def __compose(self):
        """Place children widgets."""
        for btn in self.__buttons.values():
            try:
                btn.setAutoDefault(False)
            except AttributeError:
                pass

        self.__buttons['set_bad_px_path'].setDefaultAction(
            qtw.QAction(
                self.style().standardIcon(self.style().SP_DialogOpenButton),
                ''
                )
            )

        layout = qtw.QVBoxLayout()
        layout.setSpacing(layout.spacing() + 10)
        layout.addLayout(self.__compose_camera_selector())
        layout.addWidget(self.__compose_bad_pixel_path())
        layout.addWidget(self.__compose_bad_pixel_info())
        layout.addWidget(self.__compose_progress_bars())
        layout.addStretch(1)
        layout.addLayout(self.__compose_bottom_buttons())

        self.setLayout(layout)

    def __compose_bad_pixel_info(self):
        group = qtw.QGroupBox("Bad pixels information")
        layout = qtw.QVBoxLayout()
        layout.setAlignment(qtc.Qt.AlignLeft)
        layout.addWidget(self.__bad_px_info['date_time'])
        layout.addWidget(self.__bad_px_info['n_bad'])
        layout.addWidget(self.__bad_px_info['n_uncorrectable'])
        group.setLayout(layout)
        return group

    def __compose_bad_pixel_path(self):
        group = qtw.QGroupBox("Bad pixel directory")
        layout = qtw.QHBoxLayout()
        layout.addWidget(self.__buttons['set_bad_px_path'])
        layout.addWidget(self.__ctrls['bad_px_path'], stretch=1)

        has_bad_px_path = self.__ctrls['bad_px_path'].text() != 'None selected'
        change_control_text_color(self.__ctrls['bad_px_path'],
                                  'black' if has_bad_px_path else 'red')
        group.setLayout(layout)
        return group

    def __compose_bottom_buttons(self):
        layout = qtw.QHBoxLayout()
        layout.addWidget(self.__buttons['find'])
        layout.addWidget(self.__buttons['abort'])
        layout.addStretch(1)  # flush "Done" to the right
        layout.addWidget(self.__buttons['done'])
        return layout

    def __compose_camera_selector(self):
        layout = qtw.QHBoxLayout()
        layout.addWidget(qtw.QLabel("Camera:"))
        layout.addWidget(self.__ctrls['camera'], stretch=1)
        return layout

    def __compose_progress_bars(self):
        group = self.__progress['group']
        layout = qtw.QVBoxLayout()
        layout.setAlignment(qtc.Qt.AlignLeft)
        layout.addWidget(qtw.QLabel("Finding bad pixels..."))
        layout.addWidget(self.__progress['total'])
        layout.addWidget(self.__progress['section_text'])
        layout.addWidget(self.__progress['section'])
        group.setLayout(layout)
        self.__reset_progress_bars()
        group.hide()
        return group

    def __connect(self):
        """Connect children signals."""
        self.__buttons['done'].clicked.connect(self.accept)
        self.__buttons['find'].clicked.connect(self.__start)
        self.__buttons['abort'].clicked.connect(self.__abort)
        self.__buttons['set_bad_px_path'].triggered.connect(
            self.__on_set_bad_pixel_directory
            )

        self.__ctrls['camera'].currentTextChanged.connect(
            self.__on_camera_selected,
            type=qtc.Qt.QueuedConnection
            )

        timer, _ = self.__timers['update_list']
        timer.setSingleShot(True)
        timer.timeout.connect(self.update_available_camera_list)

    def __enable_controls(self, enabled):
        """Enable or disable controls."""
        for ctrl in (*self.__buttons.values(), *self.__ctrls.values()):
            ctrl.setEnabled(enabled)

        timer, interval = self.__timers['update_list']
        if enabled:
            timer.start(interval)
        else:
            timer.stop()

        # Decide whether controls can be enabled. This
        # depends on whether there is any camera. Cannot
        # find anything if there is no camera...
        has_camera = self.__ctrls['camera'].count() > 0
        has_bad_px_path = self.__ctrls['bad_px_path'].text() != 'None selected'

        bad_pix_enabled = enabled and has_camera
        find_enabled = bad_pix_enabled and has_bad_px_path
        abort_enabled = not enabled and has_camera

        for ctrl in (self.__ctrls['bad_px_path'],
                     self.__buttons['set_bad_px_path']):
            ctrl.setEnabled(bad_pix_enabled)

        change_control_text_color(self.__ctrls['bad_px_path'],
                                  'black' if has_bad_px_path else 'red')

        self.__buttons['abort'].setEnabled(abort_enabled)
        self.__buttons['find'].setEnabled(find_enabled)

    def __get_bad_pixel_info(self):
        """Return bad pixel info from the active camera.

        Returns
        -------
        date_time : str
            Date and time of the most recent bad-pixel file.
            Returns "\u2014" (em-dash) if no bad-pixel file
            was read. Otherwise a string with format
            "YYYY-mm-dd HH:MM:SS".
        n_bad : int
            Number of bad pixels. Returns -1 if no bad-pixels
            file was found.
        bad_fraction : float
            Number of bad pixels as fraction of the sensor
            area in percent. Returns -1 if no bad-pixels file
            was found.
        n_uncorrectable : int
            Number of uncorrectable bad pixels. Returns -1 if
            no bad-pixels file was found.
        uncorrectable_fraction : float
            Number of uncorrectable bad pixels as fraction of
            the sensor area in percent. Returns -1 if no bad-
            pixels file was found.
        """
        cam = self.active_camera
        if not cam.bad_pixels.file_name:
            return "\u2014", -1, -1, -1, -1

        *_, date, time = cam.bad_pixels.file_name.split('_')
        date_time = (f"{date[:4]}-{date[4:6]}-{date[6:]} "
                     f"{time[:2]}:{time[2:4]}:{time[4:]}")

        # We suppose that the camera ROI is currently set
        # to be the whole sensor, as it should whenever
        # this function is called
        width, height, *_ = cam.image_info
        sensor = width*height
        n_bad = cam.bad_pixels.n_bad_pixels_sensor
        n_uncorrectable = cam.bad_pixels.n_uncorrectable_sensor
        bad_fraction = 100*n_bad/sensor
        uncorrectable_fraction = 100*n_uncorrectable/sensor

        return (date_time, n_bad, bad_fraction,
                n_uncorrectable, uncorrectable_fraction)

    def __get_camera_config(self, camera_name):  # TODO: use self to report errors
        """Return the configuration file for a camera with a given name."""
        config_path = DEFAULT_CONFIG_PATH
        config_files = [f for f in config_path.glob('**/*')
                        if f.is_file() and f.suffix == '.ini']
        camera_config_files = []
        for config_path in config_files:
            with open(config_path, 'r') as config_file:
                if camera_name in config_file.read():
                    camera_config_files.append(config_path)

        if not camera_config_files:
            # TODO: report error with a dialog
            raise RuntimeError(
                f"Found no config file for camera {camera_name}."
                )
        if len(camera_config_files) > 1:
            # TODO: report error with a dialog
            raise RuntimeError(
                f"Found multiple ({len(camera_config_files)} configuration"
                f"files for camera {camera_name}."
                )
        return camera_config_files[0]

    def __on_camera_selected(self, camera_name):
        """React to selection of a new camera."""
        if not camera_name:
            # List is empty
            self.active_camera = None
            return
        if self.active_camera and self.active_camera.name == camera_name:
            # Same selection
            return

        # TODO: pop up a small modal window reporting that
        # the camera is being set up. To be closed as soon
        # as the camera object is instantiated correctly

        # New camera selected.
        settings = self.__get_camera_config(camera_name)
        cls = self.__available_cameras[camera_name]
        self.active_camera = cls()
        # self.active_camera.moveToThread(self.__finder_thread)
        self.active_camera.error_occurred.connect(self.__on_error_occurred)
        self.active_camera.started.connect(self.adjustSize)
        self.__set_camera_settings.connect(self.active_camera.set_settings,
                                           type=qtc.Qt.QueuedConnection)
        self.__set_camera_settings.emit(settings)

        self.__update_controls()

    def __on_error_occurred(self, error_info):
        """React to an error situation."""
        error_code, error_msg = error_info
        try:
            error = camera_abc.CameraErrors.from_code(error_code)
        except AttributeError:
            error = None
        if (error is camera_abc.CameraErrors.MISSING_SETTINGS
            or (error is camera_abc.CameraErrors.INVALID_SETTINGS
                and ("could not find a bad pixels file" in error_msg
                     or "No bad_pixel_path found" in error_msg))):
            # We don't want to spam messages if the current
            # folder does not contain a bad pixels file. Also,
            # we will not allow running anything if there is
            # no bad_pixel_path option in the configuration,
            # but we report this 'invalid' situation by having
            # red text.
            return
        qtw.QMessageBox.critical(self, "Error",
                                 f"{error_msg}\n\n(Code: {error_code})")

    def __on_finder_done(self):
        """React to bad pixels being found."""
        bar_total = self.__progress['total']
        bar_total.setValue(bar_total.maximum())

        # Store the old bad pixel information
        keys = ('date_time', 'n_bad', 'bad_fraction',
                'n_uncorrectable', 'uncorrectable_fraction')
        old = {k: v for k, v in zip(keys, self.__get_bad_pixel_info())}
        if old['n_bad'] < 0:
            # No file
            old['date_time'] = "None"
            for key in keys[1:]:
                old['key'] = None

        # Now update the bad pixel info in the camera. This will
        # read the newly created bad-pixels file.
        self.active_camera.update_bad_pixels()
        new = {k: v for k, v in zip(keys, self.__get_bad_pixel_info())}

        # Now prepare strings
        fmt = "{} ({:.2f}% of sensor)"
        fmt_previous = "\t[Previous: {}]"
        date_time_txt = (f"{new['date_time']}"
                         + fmt_previous.format(old['date_time']))

        bad_txt = fmt.format(new['n_bad'], new['bad_fraction'])
        if old['n_bad'] is not None:
            bad_txt += fmt_previous.format(old['n_bad'])

        uncorrectable_txt = fmt.format(new['n_uncorrectable'],
                                       new['uncorrectable_fraction'])
        if old['n_uncorrectable'] is not None:
            uncorrectable_txt += fmt_previous.format(old['n_uncorrectable'])

        self.__bad_px_info['date_time'].setText(f"Date/Time: {date_time_txt}")
        self.__bad_px_info['n_bad'].setText(f"No. bad pixels: {bad_txt}")
        self.__bad_px_info['n_uncorrectable'].setText(
                f"No. uncorrectable: {uncorrectable_txt}"
                )

        self.__enable_controls(True)

    def __on_progress(self, *args):
        """React to a progress report from the finder."""
        sec_txt, sec_no, tot_secs, progress, n_tasks = args
        bar_total = self.__progress['total']
        section = self.__progress['section_text']
        bar_section = self.__progress['section']

        if bar_total.maximum() != tot_secs:
            bar_total.setMaximum(tot_secs)

        section.setText(sec_txt)
        if bar_section.maximum() != n_tasks:
            bar_section.setMaximum(n_tasks)
        bar_section.setValue(progress)

        total_progress = sec_no
        if bar_section.value() < bar_section.maximum():
            total_progress -= 1
        bar_total.setValue(total_progress)

    def __on_set_bad_pixel_directory(self, *_):
        """React to a user request to set the bad pixel directory."""
        new_directory = qtw.QFileDialog.getExistingDirectory(
            parent=self, caption="Set directory for bad pixels files"
            )
        if not new_directory:
            # User exited without selecting
            return
        cam = self.active_camera
        cam.settings.set("camera_settings",  "bad_pixels_path", new_directory)
        cam.update_bad_pixels()
        self.__update_controls()

    def __reset_progress_bars(self):
        for bar in (self.__progress['total'], self.__progress['section']):
            bar.setMinimum(0)
            bar.setValue(0)
        self.__progress['section_text'].setText('')

    def __start(self, *_):
        """Begin finding bad pixels for the selected camera."""
        self.__enable_controls(False)
        self.__reset_progress_bars()
        self.__progress['group'].show()

        # Process events to have the progress
        # bars properly shown immediately.
        qtw.qApp.processEvents()

        self.__finder = badpixels.BadPixelsFinder(self.active_camera)
        self.__finder.progress_occurred.connect(self.__on_progress)
        self.__finder.done.connect(self.__on_finder_done)
        self.__finder.error_occurred.connect(self.__on_error_occurred)
        self.__start_finder.connect(self.__finder.find)
        self.__abort_finder.connect(self.__finder.abort)
        self.__finder.moveToThread(self.__finder_thread)
        self.__start_finder.emit()
        # self.__finder.find()

    def __stop_timers(self):
        """Stop all timers."""
        for timer, _ in self.__timers.values():
            timer.stop()

    def __update_bad_px_info_latest(self):
        """Update the bad pixel info widgets with the latest file only."""
        (date_time, n_bad, bad_fraction,
         n_uncorrectable, uncorrectable_fraction) = self.__get_bad_pixel_info()

        if n_bad < 0:
            # No file found
            date_time_txt = "No file found!"
            bad_txt = uncorrectable_txt = "\u2014"  # em-dash
        else:
            fmt = "{} ({:.2f}% of sensor)"
            date_time_txt = date_time
            bad_txt = fmt.format(n_bad, bad_fraction)
            uncorrectable_txt = fmt.format(n_uncorrectable,
                                           uncorrectable_fraction)
        self.__bad_px_info['date_time'].setText(f"Date/Time: {date_time_txt}")
        self.__bad_px_info['n_bad'].setText(f"No. bad pixels: {bad_txt}")
        self.__bad_px_info['n_uncorrectable'].setText(
                f"No. uncorrectable: {uncorrectable_txt}"
                )

    def __update_controls(self):
        """Update the contents of controls from a camera."""
        cam = self.active_camera
        # See if we have a bad_pixels_path:
        try:
            bad_pixels_path = cam.settings.get("camera_settings",
                                               "bad_pixels_path",
                                               fallback="")
        except AttributeError:
            # settings is None
            bad_pixels_path = ""

        if bad_pixels_path:
            bad_pixels_path = Path(bad_pixels_path)
            if not bad_pixels_path.exists():
                bad_pixels_path = ""
        if not bad_pixels_path:
            # Cannot read bad pixels, nor save them.
            self.__buttons['find'].setEnabled(False)
            self.__ctrls['bad_px_path'].setText("None selected")
            self.__bad_px_info['date_time'].setText("Date/Time: \u2014")
            self.__bad_px_info['n_bad'].setText("No. bad pixels: \u2014")
            self.__bad_px_info['n_uncorrectable'].setText(
                "No. uncorrectable: \u2014"
                )
            return

        self.__ctrls['bad_px_path'].setText(str(bad_pixels_path.resolve()))
        self.__update_bad_px_info_latest()
