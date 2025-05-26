"""Module badpxfinderdialog of viperleed.gui.measure.dialogs.

Defines the BadPixelsFinderDialog class that handles user interaction
while finding bad pixels for a camera.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2021-11-22'
__license__ = 'GPLv3+'

from pathlib import Path

from PyQt5 import (QtWidgets as qtw,
                   QtCore as qtc)

from viperleed.gui.measure import dialogs
from viperleed.gui.measure import hardwarebase as base
from viperleed.gui.measure.camera import badpixels
from viperleed.gui.measure.camera import abc as camera_abc
from viperleed.gui.measure.classes.settings import get_system_config
from viperleed.gui.widgetslib import change_control_text_color


DEFAULT_CONFIG_PATH = Path(
    get_system_config().get("PATHS", 'configuration', fallback='')
    )
NOT_FOUND = "No file found!"
NOT_SET = "\u2014"
NO_BAD_PX_PATH = "None selected"


class BadPixelsFinderDialog(qtw.QDialog):
    """Dialog to handle user interaction when finding bad pixels."""

    __start_finder = qtc.pyqtSignal()
    __abort_finder = qtc.pyqtSignal()

    def __init__(self, parent=None):
        """Initialize dialog."""
        super().__init__(parent=parent)

        self.__children = {
            'ctrls': {
                'camera': qtw.QComboBox(),
                'bad_px_path': qtw.QLabel(NO_BAD_PX_PATH),
                },
            'buttons' : {
                'done': qtw.QPushButton("&Done"),
                'abort': qtw.QPushButton("&Abort"),
                'find': qtw.QPushButton("&Find"),
                'set_bad_px_path': qtw.QToolButton(),
                'save_uncorr_mask': qtw.QPushButton(
                    "Save uncorrectable &mask"
                    ),
                },
            'bad_px_info': {
                'date_time': qtw.QLabel(f"Date/Time: {NOT_SET}"),
                'n_bad': qtw.QLabel(f"No. bad pixels: {NOT_SET}"),
                'n_uncorrectable': qtw.QLabel(f"No. uncorrectable: {NOT_SET}"),
                'date_time_old': qtw.QLabel(""),
                'n_bad_old': qtw.QLabel(""),
                'n_uncorrectable_old': qtw.QLabel(""),
                },
            'progress': {
                'group' : qtw.QGroupBox("Progress"),
                'total': qtw.QProgressBar(),
                'section_text': qtw.QLabel(),
                'section': qtw.QProgressBar(),
                },
            'timers': {
                # 'update_list' will update the list of camera
                # devices every 2 seconds
                'update_list': (qtc.QTimer(self), 2000),
                # 'delay_busy_hide' will delay shortly the hiding
                # of the camera-busy dialog.
                'delay_busy_hide': (qtc.QTimer(self), 100),
                },
            'dialogs': {
                'camera_busy': dialogs.BusyWindow(parent=self,
                                                  text="Preparing camera...",
                                                  max_onscreen_time=10)
                }
            }

        self.__available_cameras = {}
        self.active_camera = None
        self.__finder = None
        self.__finder_thread = qtc.QThread()
        self.__finder_thread.start()

        self.setWindowTitle("Find bad pixels")

        self.__compose()
        self.__connect()

    def showEvent(self, event):  # pylint: disable=invalid-name
        """Extend showEvent to update controls."""
        self.__reset_progress_bars()
        self.__progress['group'].hide()
        self.update_available_camera_list()
        self.__ctrls['camera'].setCurrentIndex(-1)  # TODO: path label not updated correctly
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
        self.__available_cameras = base.get_devices('camera')
        old_items = set(camera_combo.itemText(i)
                        for i in range(camera_combo.count()))
        if old_items != set(self.__available_cameras.keys()):
            # Camera list has changed
            camera_combo.clear()
            camera_combo.addItems(self.__available_cameras.keys())
            if old_selection:
                camera_combo.setCurrentText(old_selection)
        self.__enable_controls(True)

    @property
    def __bad_px_info(self):
        """Return controls for reporting bad-pixel information."""
        return self.__children['bad_px_info']

    @property
    def __camera_busy(self):
        """Return a reference to the modal camera-busy window."""
        return self.__children['dialogs']['camera_busy']

    @property
    def __buttons(self):
        """Return a dictionary of buttons."""
        return self.__children['buttons']

    @property
    def __ctrls(self):
        """Return controls to be enabled/disabled."""
        return self.__children['ctrls']

    @property
    def __progress(self):
        """Return controls for progress reporting."""
        return self.__children['progress']

    @property
    def __timers(self):
        """Return a dictionary of timers."""
        return self.__children['timers']

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
        # Remove the '?' button from the title bar
        self.setWindowFlags(self.windowFlags()
                            & ~qtc.Qt.WindowContextHelpButtonHint)

        # Make all buttons not the 'default' button (i.e.,
        # the one that is automatically triggered when the
        # user presses "Enter".
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
        layout = qtw.QGridLayout()
        layout.addWidget(self.__bad_px_info['date_time'], 0, 0)
        layout.addWidget(self.__bad_px_info['date_time_old'], 0, 1)
        layout.addWidget(self.__bad_px_info['n_bad'], 1, 0)
        layout.addWidget(self.__bad_px_info['n_bad_old'], 1, 1)
        layout.addWidget(self.__bad_px_info['n_uncorrectable'], 2, 0)
        layout.addWidget(self.__bad_px_info['n_uncorrectable_old'], 2, 1)
        layout.addWidget(self.__buttons['save_uncorr_mask'], 3, 0, 1, 2)
        group.setLayout(layout)
        return group

    def __compose_bad_pixel_path(self):
        group = qtw.QGroupBox("Bad pixel directory")
        layout = qtw.QHBoxLayout()
        layout.addWidget(self.__buttons['set_bad_px_path'])
        layout.addWidget(self.__ctrls['bad_px_path'], stretch=1)

        has_bad_px_path = self.__ctrls['bad_px_path'].text() != NO_BAD_PX_PATH
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
        self.__ctrls['camera'].setPlaceholderText("<Select camera>")
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
        self.__buttons['save_uncorr_mask'].clicked.connect(
                self.__on_save_mask_clicked
                )

        self.__ctrls['camera'].textActivated.connect(
            self.__on_camera_selected,
            type=qtc.Qt.QueuedConnection
            )

        timer, _ = self.__timers['update_list']
        timer.setSingleShot(True)
        timer.timeout.connect(self.update_available_camera_list)

        timer, _ = self.__timers['delay_busy_hide']
        timer.setSingleShot(True)
        timer.timeout.connect(self.__camera_busy.hide)

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
        # find anything if there is no camera, or if
        # settings are invalid
        has_camera = (self.__ctrls['camera'].count() > 0
                      and bool(self.active_camera)
                      and bool(self.active_camera.settings))
        has_bad_px_path = self.__ctrls['bad_px_path'].text() != NO_BAD_PX_PATH

        bad_pix_enabled = enabled and has_camera
        find_enabled = bad_pix_enabled and has_bad_px_path
        abort_enabled = not enabled and has_camera
        has_bad_px_info = (
            find_enabled
            and NOT_FOUND not in self.__bad_px_info['date_time'].text()
            and NOT_SET not in self.__bad_px_info['date_time'].text()
            )

        for ctrl in (self.__ctrls['bad_px_path'],
                     self.__buttons['set_bad_px_path']):
            ctrl.setEnabled(bad_pix_enabled)

        change_control_text_color(self.__ctrls['bad_px_path'],
                                  'black' if has_bad_px_path else 'red')

        self.__buttons['abort'].setEnabled(abort_enabled)
        self.__buttons['find'].setEnabled(find_enabled)
        self.__buttons['save_uncorr_mask'].setEnabled(has_bad_px_info)

    def __get_bad_pixel_info(self, previous=False):
        """Return bad pixel info from the active camera.

        Parameters
        ----------
        previous : bool, optional
            If True, return the info from the pre-last bad pixels
            file in the folder. If False, take the info from the
            most recent (i.e., the one loaded in the camera).
            Default is False.

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
        if not previous:
            bad_px = cam.bad_pixels
        else:
            bad_px = badpixels.BadPixels(cam)
            try:
                bad_px.read(self.__ctrls['bad_px_path'].text(),
                            most_recent=False)
            except FileNotFoundError:
                bad_px = None

        if not bad_px or not bad_px.file_name:
            return NOT_SET, -1, -1, -1, -1

        *_, date, time = bad_px.file_name.split('_')
        date_time = (f"{date[:4]}-{date[4:6]}-{date[6:]} "
                     f"{time[:2]}:{time[2:4]}:{time[4:]}")

        width, height = cam.sensor_size
        sensor = width*height
        n_bad = bad_px.n_bad_pixels_sensor
        n_uncorrectable = bad_px.n_uncorrectable_sensor
        bad_fraction = 100*n_bad/sensor
        uncorrectable_fraction = 100*n_uncorrectable/sensor

        return (date_time, n_bad, bad_fraction,
                n_uncorrectable, uncorrectable_fraction)

    def __on_camera_preparing(self, busy):
        """Show a busy dialog while camera prepares to acquire."""
        timer, interval = self.__timers['delay_busy_hide']
        if busy:
            if not self.__camera_busy.isVisible():
                self.__camera_busy.show()
            timer.stop()
        else:
            timer.start(interval)

    def __on_camera_selected(self, camera_name):
        """React to selection of a new camera."""
        if not camera_name:
            # List is empty
            self.active_camera = None
            return
        if self.active_camera and self.active_camera.name == camera_name:
            # Same selection
            return

        # New camera selected.
        settings = base.get_device_config(camera_name,
                                          directory=DEFAULT_CONFIG_PATH,
                                          parent_widget=self)

        # Signal errors by picking an invalid entry
        if not settings:
            self.__ctrls['camera'].setCurrentIndex(-1)
            return

        if not self.__camera_busy.isVisible():
            self.__camera_busy.show()
        cls = self.__available_cameras[camera_name]
        self.active_camera = cls()
        self.active_camera.error_occurred.connect(self.__on_error_occurred)
        self.active_camera.started.connect(self.adjustSize)
        self.active_camera.preparing.connect(self.__on_camera_preparing)
        self.active_camera.settings = settings
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
        old = dict(zip(keys, self.__get_bad_pixel_info(previous=True)))
        if old['n_bad'] < 0:
            # No file
            old['date_time'] = "None"
            for key in keys[1:]:
                old[key] = None

        # Now update the bad pixel info in the camera. This will
        # read the newly created bad-pixels file.
        self.active_camera.update_bad_pixels()
        new = dict(zip(keys, self.__get_bad_pixel_info()))

        # Now prepare strings
        fmt = "{} ({:.2f}% of sensor)"
        fmt_previous = "[Previous: {}]"
        date_time_txt = f"{new['date_time']}"
        date_time_txt_old = fmt_previous.format(old['date_time'])

        bad_txt = fmt.format(new['n_bad'], new['bad_fraction'])
        bad_txt_old = ""
        if old['n_bad'] is not None:
            bad_txt_old = fmt_previous.format(old['n_bad'])

        uncorrectable_txt = fmt.format(new['n_uncorrectable'],
                                       new['uncorrectable_fraction'])
        uncorrectable_txt_old = ""
        if old['n_uncorrectable'] is not None:
            uncorrectable_txt_old = fmt_previous.format(old['n_uncorrectable'])

        self.__bad_px_info['date_time'].setText(f"Date/Time: {date_time_txt}")
        self.__bad_px_info['n_bad'].setText(f"No. bad pixels: {bad_txt}")
        self.__bad_px_info['n_uncorrectable'].setText(
            f"No. uncorrectable: {uncorrectable_txt}"
            )
        self.__bad_px_info['date_time_old'].setText(date_time_txt_old)
        self.__bad_px_info['n_bad_old'].setText(bad_txt_old)
        self.__bad_px_info['n_uncorrectable_old'].setText(
            uncorrectable_txt_old
            )
        self.__report_uncorrectable()
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

    def __on_save_mask_clicked(self, *_):
        """React to a request to save an uncorrectable pixels mask."""
        cam = self.active_camera
        if not cam or not cam.bad_pixels:
            # This should never happen: the button
            # is disabled if this is the case
            return
        mask_path = qtw.QFileDialog.getExistingDirectory(
            parent=self,
            caption="Set directory for uncorrectable pixels mask",
            directory=self.__ctrls['bad_px_path'].text()
            )
        if not mask_path:
            # User exited without selecting
            return
        cam.bad_pixels.save_uncorrectable_mask(mask_path)

    def __on_set_bad_pixel_directory(self, *_):
        """React to a user request to set the bad pixel directory."""
        new_directory = qtw.QFileDialog.getExistingDirectory(
            parent=self, caption="Set directory for bad pixels files"
            )
        if not new_directory:
            # User exited without selecting
            return
        cam = self.active_camera
        cam.settings.set("camera_settings", "bad_pixels_path", new_directory)
        cam.update_bad_pixels()
        cam.settings.update_file()
        self.__update_controls()

    def __reset_progress_bars(self):
        for progress_bar in (self.__progress['total'],
                             self.__progress['section']):
            progress_bar.setMinimum(0)
            progress_bar.setValue(0)
        self.__progress['section_text'].setText('')

    def __report_uncorrectable(self):
        """Report an error if there are many uncorrectable pixels."""
        bad_px = self.active_camera.bad_pixels
        cluster_sizes = bad_px.get_uncorrectable_clusters_sizes()

        # Complain only if there are clusters larger than 1
        cluster_sizes = [s for s in cluster_sizes if s > 1]

        if not cluster_sizes:
            return
        qtw.QMessageBox.critical(
            self,  # parent
            "Too many uncorrectable bad pixels",  # title
            f"Camera {self.active_camera.name} has a significant number "
            "of bad pixels that cannot be corrected. This will have a "
            "serious impact on the quality of your LEED-IV data. Consider "
            f"replacing the camera!\nFound {len(cluster_sizes)} bad-pixel "
            f"clusters with sizes as large as {max(cluster_sizes)} pixels.\n"
            "\nAs a temporary solution you can create an uncorrectable-"
            "pixels mask, process it the same way as your images, and "
            "combine it with the mask for the LEED screen and electron "
            "gun to discard the damaged areas of the sensor."
            )

    def __start(self, *_):
        """Begin finding bad pixels for the selected camera."""
        self.__enable_controls(False)
        self.__reset_progress_bars()
        self.__progress['group'].show()
        self.__camera_busy.show()

        # Process events to have the progress
        # bars properly shown immediately.
        qtw.qApp.processEvents()

        self.__finder = badpixels.BadPixelsFinder(self.active_camera)
        self.__finder.progress_occurred.connect(self.__on_progress)
        self.__finder.done.connect(self.__on_finder_done)
        self.__finder.done.connect(self.__finder.deleteLater)
        self.__finder.aborted.connect(self.__finder.deleteLater)
        self.__finder_thread.finished.connect(self.__finder.deleteLater)
        self.__finder.error_occurred.connect(self.__on_error_occurred)
        self.__start_finder.connect(self.__finder.find)
        self.__abort_finder.connect(self.__finder.abort)
        self.__finder.moveToThread(self.__finder_thread)
        self.__start_finder.emit()

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
            date_time_txt = NOT_FOUND
            bad_txt = uncorrectable_txt = NOT_SET  # em-dash
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
        for key in ('date_time_old', 'n_bad_old', 'n_uncorrectable_old'):
            self.__bad_px_info[key].setText("")

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
            # pylint: disable=redefined-variable-type
            #         Disabled for str->Path as it is easier
            #         than making a new variable name
            bad_pixels_path = Path(bad_pixels_path)
            if not bad_pixels_path.exists():
                bad_pixels_path = ""
        if not bad_pixels_path:
            # Cannot read bad pixels, nor save them.
            self.__buttons['find'].setEnabled(False)
            self.__ctrls['bad_px_path'].setText(NO_BAD_PX_PATH)
            self.__bad_px_info['date_time'].setText(f"Date/Time: {NOT_SET}")
            self.__bad_px_info['n_bad'].setText(f"No. bad pixels: {NOT_SET}")
            self.__bad_px_info['n_uncorrectable'].setText(
                f"No. uncorrectable: {NOT_SET}"
                )
            return

        self.__ctrls['bad_px_path'].setText(str(bad_pixels_path.resolve()))
        self.__update_bad_px_info_latest()
        self.__enable_controls(True)
