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

    def __init__(self, parent=None):
        """Initialize dialog."""
        super().__init__(parent=parent)

        self.__ctrls = {
            'camera': qtw.QComboBox(),
            'bad_px_path': qtw.QLabel("None selected"),
            'total_progress': qtw.QProgressBar(),
            'section_progress': qtw.QProgressBar(),
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
        self.__timers = {'update_list': (qtc.QTimer(self), 2000),}

        self.__available_cameras = {}
        self.active_camera = None
        self.__finder = None

        self.setWindowTitle("Find bad pixels")

        self.__compose()
        self.__connect()

    def showEvent(self, event):
        """Extend showEvent to update list of cameras."""
        self.update_available_camera_list()
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

    def __clean_up(self):
        """Clean up self before closing."""
        self.__stop_timers()
        self.__ctrls['camera'].clear()
        self.active_camera = None

    def __compose(self):
        """Place children widgets."""
        select_cam_layout = qtw.QHBoxLayout()
        select_cam_layout.addWidget(qtw.QLabel("Camera:"))
        select_cam_layout.addWidget(self.__ctrls['camera'], stretch=1)

        for btn in self.__buttons.values():
            try:
                btn.setAutoDefault(False)
            except AttributeError:
                pass

        set_bad_px_path = qtw.QAction(
            self.style().standardIcon(self.style().SP_DialogOpenButton), ''
            )
        self.__buttons['set_bad_px_path'].setDefaultAction(set_bad_px_path)

        bad_px_path_group = qtw.QGroupBox("Bad pixel directory")
        bad_px_path_layout = qtw.QHBoxLayout()
        bad_px_path_layout.addWidget(self.__buttons['set_bad_px_path'])
        bad_px_path_layout.addWidget(self.__ctrls['bad_px_path'],stretch=1)
        bad_px_path_group.setLayout(bad_px_path_layout)

        bad_px_info_group = qtw.QGroupBox("Bad pixels information")
        bad_px_info_layout = qtw.QVBoxLayout()
        bad_px_info_layout.setAlignment(qtc.Qt.AlignLeft)
        bad_px_info_layout.addWidget(self.__bad_px_info['date_time'])
        bad_px_info_layout.addWidget(self.__bad_px_info['n_bad'])
        bad_px_info_layout.addWidget(self.__bad_px_info['n_uncorrectable'])
        bad_px_info_group.setLayout(bad_px_info_layout)

        bottom_buttons_layout = qtw.QHBoxLayout()
        bottom_buttons_layout.addWidget(self.__buttons['find'])
        bottom_buttons_layout.addWidget(self.__buttons['abort'])
        bottom_buttons_layout.addStretch(1)  # flush "Done" to the right
        bottom_buttons_layout.addWidget(self.__buttons['done'])

        layout = qtw.QVBoxLayout()
        layout.setSpacing(layout.spacing() + 10)
        layout.addLayout(select_cam_layout)
        layout.addWidget(bad_px_path_group)
        layout.addWidget(bad_px_info_group)
        layout.addStretch(1)
        layout.addLayout(bottom_buttons_layout)

        self.setLayout(layout)

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
        settings = self.__find_camera_config(camera_name)
        cls = self.__available_cameras[camera_name]
        self.active_camera = cls(settings=settings)
        self.active_camera.error_occurred.connect(self.__on_error_occurred)

        self.__update_controls()

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

    def __start(self, *_):
        """Begin finding bad pixels for the selected camera."""
        self.__enable_controls(False)
        self.__finder = badpixels.BadPixelsFinder(self.active_camera)
        self.__finder.find()

    def __stop_timers(self):
        """Stop all timers."""
        for timer, _ in self.__timers.values():
            timer.stop()

    def __abort(self, *_):
        """Abort bad-pixel-finder routine."""
        self.__finder.abort()
        self.__enable_controls(True)

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

    def __on_error_occurred(self, error_info):
        """React to an error situation."""
        error_code, error_msg = error_info
        try:
            error = camera_abc.CameraErrors.from_code(error_code)
        except AttributeError:
            error = None
        if (error is camera_abc.CameraErrors.INVALID_SETTINGS
            and ("could not find a bad pixels file" in error_msg
                 or "No bad_pixel_path found" in error_msg)):
            # We don't want to spam messages if the current
            # folder does not contain a bad pixels file. Also,
            # we will not allow running anything if there is
            # no bad_pixel_path option in the configuration,
            # but we report this 'invalid' situation by having
            # red text.
            return

        print("\n>>>>   ERROR  <<<<")
        print(error_info)

    def __update_controls(self):
        """Update the contents of controls from a camera."""
        cam = self.active_camera
        # See if we have a bad_pixels_path:
        bad_pixels_path = cam.settings.get("camera_settings",
                                           "bad_pixels_path",
                                           fallback="")
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
        if not cam.bad_pixels.file_name:
            # No file could be read
            self.__bad_px_info['date_time'].setText(
                "Date/Time: No file found!"
                )
            self.__bad_px_info['n_bad'].setText("No. bad pixels: \u2014")
            self.__bad_px_info['n_uncorrectable'].setText(
                "No. uncorrectable: \u2014"
                )
            return
        *_, date, time = cam.bad_pixels.file_name.split('_')
        date_time = (f"{date[:4]}-{date[4:6]}-{date[6:]} "
                     f"{time[:2]}:{time[2:4]}:{time[4:]}")
        width, height, *_ = cam.image_info
        sensor = width*height
        n_bad = cam.bad_pixels.n_bad_pixels_sensor
        n_uncorrectable = cam.bad_pixels.n_uncorrectable_sensor
        self.__bad_px_info['date_time'].setText(f"Date/Time: {date_time}")
        self.__bad_px_info['n_bad'].setText(
            f"No. bad pixels: {n_bad} ({100*(n_bad/sensor):.2f}% of sensor)"
            )
        self.__bad_px_info['n_uncorrectable'].setText(
            f"No. uncorrectable: {n_uncorrectable} "
            f"({100*(n_uncorrectable/sensor):.2}% of sensor)"
            )

    def __find_camera_config(self, camera_name):  # TODO: use self to report errors
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
