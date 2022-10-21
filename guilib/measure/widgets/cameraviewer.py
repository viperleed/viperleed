"""Module cameraviewer of viperleed.guilib.measure.widgets.

========================================
   ViPErLEED Graphical User Interface
========================================

Created: 2021-10-29
Author: Michele Riva

Defines the CameraViewer class, a QScrollArea that allows
visualizing frames from a concrete subclass of CameraABC.
"""

import weakref
from copy import deepcopy
from time import perf_counter as timer                                          # TEMP: till "saturating" is fixed

import numpy as np
from PyQt5 import (QtCore as qtc,
                   QtWidgets as qtw,
                   QtGui as qtg)

from viperleed.guilib.measure import hardwarebase as base
from viperleed.guilib.measure.camera.imageprocess import ImageProcessor
from viperleed.guilib.measure.dialogs.settingsdialog import SettingsDialog
from viperleed.guilib.measure.widgets.imageviewer import ImageViewer
from viperleed.guilib.measure.widgets.roiwidgets import RegionOfInterest
from viperleed.guilib.widgetslib import screen_fraction


# TODO: ImageViewer.optimum_size is not updated when screen is changed
# TODO: ROI show size in image coordinates as it is resized (tooltip?)
# TODO: ROI context menu: precisely set with coordinates
# TODO: context -- properties
# TODO: use ROI position increments
# TODO: If camera is not open at CameraViewer__init__, I never update
#       the limits! Should use camera.started to fetch them if needed


# pylint: disable=too-many-instance-attributes
# Disabled because pylint counts also class attributes, but
# pyqtSignals cannot be grouped into a container like other
# attributes and returned as @property.
class CameraViewer(qtw.QScrollArea):
    """Class for displaying camera frames with scroll bars.

    Zooming can be performed when pressing '+' and '-' or by
    holding down the Control key (Command on Mac) while scrolling.
    A region of interest can be selected and moved/resized.

    Attributes
    ----------
    image_size : QtCore.QSize
        Size in pixel (without scaling) of the image shown.
    roi : RegionOfInterest
        Region of interest object.
    roi_visible : bool
        Whether the ROI is made visible upon user interaction.
    screen : QtGui.QScreen or None
        The screen on which this widget is currently shown.
    show_auto : bool
        Whether the widget is made visible when new frames arrive.
    stop_on_close : bool
        Whether the camera is stopped when the window is closed.

    Public methods
    --------------
    scale_image(by_factor)
        Resize self and the image shown by the given factor.

    Reimplemented methods
    ---------------------
    changeEvent(event)
        Handle window maximization.
    closeEvent(event)
        Stop camera when user closes the window.
    keyPressEvent(event)
        Handle zooming and ROI visibility/scaling/movements.
    mouseDoubleClickEvent(event)
        Prevent propagation to parent if ROI is visible.
    mouseMoveEvent(event)
        Draw a ROI.
    mousePressEvent(event)
        Initiate ROI drawing.
    mouseReleaseEvent(event)
        Prevent propagation to parent if ROI is visible.
    showEvent(event)
        Resize self and image when shown.
    sizeHint()
        Optimal size is bound to image size.
    wheelEvent(event)
        Handle zooming.

    Signals
    -------
    shown : no argument
        Emitted each time this object is made visible.
    image_size_changed : no argument
        Emitted each time the number of pixels acquired by the
        camera changes. Access the new size via self.image_size.
    zoom_changed : float
        Emitted each time the scaling factor used to display
        camera frames changes. Carries the new zoom factor.
    """

    # _viewers is a dict of {camera: viewer} used to prevent
    # creating multiple viewers for the same camera. It holds
    # weak references to the viewer objects, so that we don't
    # risk creating cyclic references.
    _viewers = {}

    shown = qtc.pyqtSignal()
    image_size_changed = qtc.pyqtSignal()
    zoom_changed = qtc.pyqtSignal(float)  # new zoom factor

    def __new__(cls, camera, *args, **kwargs):
        """Return a CameraViewer instance.

        This method is reimplemented such that only one viewer
        can exists for a given camera.

        Parameters
        ----------
        camera : CameraABC
            The camera object whose frames are displayed.
        *args : object
            Other positional arguments not used within__new__.
        **kwargs : object
            Optional arguments not used within __new__.

        Returns
        -------
        viewer : CameraViewer
            The viewer object.
        """
        if camera in cls._viewers:
            viewer = cls._viewers[camera]()
            if viewer is not None:
                # Viewer is alive
                return viewer
            # Viewer has been garbage-collected. Remove it from cache.
            cls._viewers.pop(camera)
        new_viewer = super().__new__(cls, camera, *args, **kwargs)
        new_viewer._initialized = False   # Call __init__ only once
        cls._viewers[camera] = weakref.ref(new_viewer)
        return new_viewer

    def __del__(self):
        """Remove self from the dictionary of viewers."""
        self._viewers.pop(self.camera)
        try:
            super().__del__()
        except AttributeError:
            pass

    def __init__(self, camera, *args, parent=None, stop_on_close=True,
                 show_auto=True, roi_visible=True, interactions_enabled=True,
                 **kwargs):
        """Initialize widget.

        Parameters
        ----------
        camera : CameraABC
            The camera object whose frames are displayed
        *args : object
            Other unused positional arguments, passed to
            QScrollArea.__init__
        parent : QWidget, optional
            Parent widget of self. Default is None.
        stop_on_close : bool, optional
            Stop the camera when this widget is closed.
            Default is True.
        show_auto : bool, optional
            If True this widget is made visible automatically
            whenever a new frame arrives from the camera (if
            it is not visible already). Otherwise, it should
            be made explicitly visible by the caller. This
            behaviour can be changed at any time by setting
            the .show_auto attribute. Default is True.
        roi_visible : bool, optional
            If True, a region of interest can be set and
            manipulated. This behaviour can be changed at any
            time by setting the .roi_visible attribute.
            Default is True.
        interactions_enabled : bool, optional
            Whether actions in the right-click menu are enabled.
            This attribute can be set to False in case the user
            is not supposed to interact with the camera. Default
            is True.
        **kwargs : object
            Other unused optional arguments, passed to
            QScrollArea.__init__

        Raises
        ------
        TypeError
            If camera is not a subclass of CameraABC
        """
        # pylint: disable=access-member-before-definition
        # Attribute is defined in __new__
        if self._initialized:
            return
        # pylint: enable=access-member-before-definition

        self._initialized = True

        super().__init__(*args, parent=parent, **kwargs)

        self.__flags = {"stop_on_close": bool(stop_on_close),
                        "show_auto": bool(show_auto),
                        "roi_visible": bool(roi_visible),
                        "interactions_enabled": bool(interactions_enabled),}
        self.__glob = {"image_size": qtc.QSize(),
                       "camera": camera,
                       "mouse_button": None,
                       "img_array": None,                                       # TODO: may not be necessary
                       "max_intensity": 2**15 - 1,}
        self.__children = {
            "viewer": ImageViewer(),
            "context_menu": qtw.QMenu(parent=self),
            "settings_dialog": SettingsDialog(handled_obj=camera),
            }
        self.__children["roi"] = RegionOfInterest(parent=self.__img_view)

        try:
            self.roi.limits = camera.get_roi_size_limits()
        except camera.exceptions:
            # Most likely camera is not open
            pass
        else:
            self.roi.update_size_limits()

        self.__compose()
        self.__connect()

    @property
    def camera(self):
        """Return the camera whose frames are displayed."""
        return self.__glob["camera"]

    @property
    def interactions_enabled(self):
        """Return whether the user can interact with the camera."""
        return self.__flags["interactions_enabled"]

    @interactions_enabled.setter
    def interactions_enabled(self, enabled):
        """Set whether the user can interact with the camera."""
        self.__flags["interactions_enabled"] = bool(enabled)
        _menu = self.__children["context_menu"]
        for action in _menu.actions():
            action.setEnabled(self.interactions_enabled)

    @property
    def image_size(self):
        """Return the true size in pixel of the image currently shown."""
        return self.__glob["image_size"]

    @image_size.setter
    def image_size(self, new_size):
        """Store the true size in pixel of the image currently shown."""
        if isinstance(new_size, tuple):
            new_size = qtc.QSize(*new_size)
        if (not new_size
                or any(s <=0 for s in (new_size.width(), new_size.height()))):
            return  # Invalid size

        if new_size != self.image_size:
            self.__glob["image_size"] = new_size
            self.__img_view.updateGeometry()
            self.updateGeometry()
            self.__img_view.optimum_size = qtc.QSize()
            self.__img_view.adjustSize()
            self.adjustSize()
            self.image_size_changed.emit()

    @property
    def roi(self):
        """Return the RegionOfInterest object of self."""
        return self.__children["roi"]

    @property
    def roi_visible(self):
        """Return whether the ROI is made visible upon user interaction."""
        return self.__flags["roi_visible"]

    @roi_visible.setter
    def roi_visible(self, visible):
        """Set the ROI to be made visible upon user interaction."""
        self.__flags["roi_visible"] = bool(visible)
        if not self.roi_visible:
            self.roi.hide()

    @property
    def screen(self):
        """Return a the QScreen on which self is shown."""
        try:
            screen = self.window().windowHandle().screen()
        except AttributeError:
            screen = None
        return screen

    @property
    def show_auto(self):
        """Return whether self is automatically shown on new frames."""
        return self.__flags["show_auto"]

    @show_auto.setter
    def show_auto(self, enabled):
        """Set if self should be shown automatically on new frames."""
        self.__flags["show_auto"] = bool(enabled)

    @property
    def stop_on_close(self):
        """Return whether the camera is stopped when closing the widget."""
        return self.__flags["stop_on_close"]

    @stop_on_close.setter
    def stop_on_close(self, stop):
        """Set whether the camera is stopped when closing the widget."""
        self.__flags["stop_on_close"] = bool(stop)

    def changeEvent(self, event):  # pylint: disable=invalid-name
        """Extend changeEvent to react to window maximization."""
        # Keep track of quantities for later adjusting the ROI
        old_scaling = self.__img_view.image_scaling
        old_roi_size = self.roi.size()
        old_roi_pos = self.roi.pos()

        if event.type() == event.WindowStateChange:
            was_maximized = event.oldState() & qtc.Qt.WindowMaximized
            if not was_maximized and self.isMaximized():
                # If window got maximized, adjust image to fit
                self.__img_view.get_scaling_to_fit(self.size()
                                                   - qtc.QSize(2, 2))
                self.__img_view.scale_to_image()

            elif was_maximized and not self.isMaximized():
                # If it is turned from maximized to normal,
                # scale it to optimum (i.e., same as after shown)
                self.__img_view.scale_to_optimum()
                self.adjustSize()

            by_factor = self.__img_view.image_scaling / old_scaling
            if abs(by_factor - 1) > 1e-3:
                # Image has been scaled. Adjust also the ROI.
                new_roi = qtc.QRect(by_factor * old_roi_pos,
                                    by_factor * old_roi_size)
                self.roi.setGeometry(new_roi)

        super().changeEvent(event)

    def closeEvent(self, event):  # pylint: disable=invalid-name
        """Extend to stop camera when window is closed by the user."""
        camera = self.camera
        if event.spontaneous() and camera.is_running and self.stop_on_close:
            # Disconnect signals that may pop up the window
            # again should a frame arrive in the meantime
            if self.show_auto:
                base.safe_disconnect(camera.image_processed, self.__show_image)
                base.safe_disconnect(camera.frame_ready, self.__show_image)
            camera.stop()
            qtw.qApp.processEvents()
        super().closeEvent(event)

    def keyPressEvent(self, event):  # pylint: disable=invalid-name
        """Extend keyPressEvent for zooming and ROI movements.

        Pressing '+' or '-' resizes the image (and self),
        pressing 'Esc' hides the ROI, pressing arrow keys
        while holding down the Control key (Command on Mac)
        affects the ROI: If Alt is also pressed, the ROI is
        resized, otherwise it is moved. If Shift is pressed,
        movements/scalings are 10 times larger.

        Parameters
        ----------
        event : QtGui.QKeyEvent
            The event corresponding to a key press

        Returns
        -------
        None.
        """
        # pylint: disable=no-member
        #    I could not find a solution to the no-member
        #    errors raised by pylint for Key.Key_*
        #    Should be fixed by white-listing PyQt5, but it
        #    is possible that Key_* attributes are added at
        #    runtime depending on the system.
        if event.key() == qtc.Qt.Key.Key_Plus:
            self.__zoom('in')
        elif event.key() == qtc.Qt.Key.Key_Minus:
            self.__zoom('out')
        elif event.key() == qtc.Qt.Key.Key_Escape:
            self.roi.hide()
        elif (self.roi.isVisible()
              and event.modifiers() & qtc.Qt.ControlModifier
              and event.key() in (qtc.Qt.Key_Left, qtc.Qt.Key_Right,
                                  qtc.Qt.Key_Up, qtc.Qt.Key_Down)):
            # Move/resize ROI if visible
            move_resize = qtc.QPoint(0, 0)
            delta = 10 if event.modifiers() & qtc.Qt.ShiftModifier else 1
            if event.key() == qtc.Qt.Key_Left:
                move_resize += qtc.QPoint(-delta, 0)
            if event.key() == qtc.Qt.Key_Right:
                move_resize += qtc.QPoint(delta, 0)
            if event.key() == qtc.Qt.Key_Up:
                move_resize += qtc.QPoint(0, -delta)
            if event.key() == qtc.Qt.Key_Down:
                move_resize += qtc.QPoint(0, delta)
            if event.modifiers() & qtc.Qt.AltModifier:
                self.roi.scale(move_resize)
            else:
                self.roi.translate(move_resize)
            return
        # pylint: enable=no-member
        super().keyPressEvent(event)

    # pylint: disable=invalid-name
    def mouseDoubleClickEvent(self, event):
        """Prevent parent propagation if ROI is visible."""
        if not self.roi.isVisible():
            super().mouseDoubleClickEvent(event)
    # pylint: enable=invalid-name

    def mouseMoveEvent(self, event):     # pylint: disable=invalid-name
        """Extend mouseMoveEvent to draw a ROI (if visible).

        If Control (Command on Mac) is held down, the ROI will
        be square.

        Parameters
        ----------
        event : QtGui.QMouseEvent
            The mouse move event that triggered this call.

        Returns
        -------
        None.
        """
        if (not self.roi.isVisible()
                or self.__mouse_button == qtc.Qt.RightButton):
            super().mouseMoveEvent(event)
            return
        mouse_pos = self.roi.parent().mapFromGlobal(event.globalPos())
        self.roi.setGeometry(
            qtc.QRect(self.roi.origin, mouse_pos)
            )

    def mousePressEvent(self, event):    # pylint: disable=invalid-name
        """Reimplement mousePressEvent to begin drawing a ROI."""
        self.__mouse_button = event.button()
        if not self.roi_visible or event.button() == qtc.Qt.RightButton:
            super().mousePressEvent(event)
            return
        if not self.roi.isVisible():
            self.roi.show()
        self.roi.origin = self.roi.parent().mapFromGlobal(event.globalPos())
        self.roi.setGeometry(qtc.QRect(self.roi.origin, qtc.QSize()))

    def mouseReleaseEvent(self, event):  # pylint: disable=invalid-name
        """Extend mouseReleaseEvent to handle ROI drawing."""
        self.__mouse_button = None
        if not self.roi.isVisible():
            super().mouseReleaseEvent(event)

    def scale_image(self, by_factor):
        """Scale self and image by the given factor."""
        # Keep track of the (scaled) offset and size of the region
        # of interest. NB: this should be done before actually scaling
        # the underlying image, as this also changes ROI size limits.
        scaled_roi_offset = (by_factor - 1) * self.roi.pos()
        scaled_roi_size = by_factor * self.roi.size()

        # Rescale the image viewer
        self.__img_view.image_scaling = (by_factor
                                         * self.__img_view.image_scaling)
        self.__img_view.scale_to_image()

        # Rescale and reposition the region of interest, if visible
        if self.roi.isVisible():
            new_roi = qtc.QRect(self.roi.pos() + scaled_roi_offset,
                                scaled_roi_size)
            self.roi.setGeometry(new_roi)

        # Now decide whether it makes sense to scale also self
        new_self_size = self.__img_view.size() + qtc.QSize(2, 2)
        if self.isMaximized():
            # Do not rescale if maximized
            scr_fraction = -1
        else:
            # See if rescaling would make self go out of the screen
            bottom_right_delta = (
                qtc.QSize(self.pos().x(), self.pos().y())
                + new_self_size
                + qtc.QSize(self.frameWidth(), self.frameWidth())
                - self.screen.availableSize())
            if any(s > 0 for s in (bottom_right_delta.width(),
                                   bottom_right_delta.height())):
                scr_fraction = -1
            else:
                scr_fraction = screen_fraction(self, new_self_size)

        # Use 10% and 80% of the screen area as limits.
        if 0.1 < scr_fraction < 0.8:
            self.resize(new_self_size)

        self.__adjust_scroll_bars(by_factor)

    def showEvent(self, event):          # pylint: disable=invalid-name
        """Extend showEvent to rescale self and the image to optimum."""
        super().showEvent(event)
        # Distinguish whether the showEvent was triggered because
        # the window is shown by the application (.spontaneous()
        # == False), or if the user restored a minimized window.
        # In the latter case, we do not want to rescale the window
        # but merely restore it to what it was last.
        if not event.spontaneous():
            self.__img_view.scale_to_optimum()
            self.adjustSize()
        self.shown.emit()

    def sizeHint(self):                  # pylint: disable=invalid-name
        """Reimplement sizeHint to bind to the underlying image."""
        # Adding (2, 2) ensures that scroll bars are
        # visible only if the image does not fit.
        return self.__img_view.sizeHint() + qtc.QSize(2, 2)

    def wheelEvent(self, event):         # pylint: disable=invalid-name
        """Extend wheelEvent for zooming while Control is pressed."""
        if not event.modifiers() & qtc.Qt.ControlModifier:
            super().wheelEvent(event)
            return
        # Mouse wheel turned with Ctrl down --> zoom
        direction = 'in' if event.angleDelta().y() > 0 else 'out'
        self.__zoom(direction=direction)

    @property
    def __img_view(self):
        """Return the ImageViewer children."""
        return self.__children["viewer"]

    @property
    def __mouse_button(self):
        """Return the mouse button currently pressed, or None."""
        return self.__glob["mouse_button"]

    @__mouse_button.setter
    def __mouse_button(self, button):
        """Set the mouse button currently pressed."""
        self.__glob["mouse_button"] = button

    def __adjust_scroll_bars(self, by_factor):
        """Adjust the position of scrollbars when zooming."""
        if self.__img_view.underMouse():
            # Zoom in at the current mouse position
            mouse_offset = (self.mapFromGlobal(qtg.QCursor.pos())
                            - self.widget().pos()) * (by_factor - 1)
            for scroll_bar, mouse_offs in zip(
                    (self.horizontalScrollBar(), self.verticalScrollBar()),
                    (mouse_offset.x(), mouse_offset.y())
                    ):
                scroll_bar.setValue(scroll_bar.value() + mouse_offs)
            return

        # Zoom centered when the mouse cursor is somewhere else
        for scroll_bar in (self.horizontalScrollBar(),
                           self.verticalScrollBar()):
            scroll_bar.setValue(
                int(by_factor * scroll_bar.value()
                    + ((by_factor - 1) * scroll_bar.pageStep()/2))
                )

    def __apply_roi(self, *_):
        """Set a new ROI in the camera."""
        old_roi_x, old_roi_y, *_ = self.camera.roi
        new_roi_x, new_roi_y, roi_w, roi_h = self.roi.image_coordinates

        # Give the camera new, updated settings. This is the easier way
        # to go, as it will update all the information in the camera.
        settings = deepcopy(self.camera.settings)
        settings.set("camera_settings", "roi",
                     f"({old_roi_x + new_roi_x}, {old_roi_y + new_roi_y}, "
                     f"{roi_w}, {roi_h})")
        was_running = self.camera.is_running
        self.camera.settings = settings
        self.camera.settings.update_file()

        # Remove the rubber-band, and restart the camera.
        self.roi.hide()
        if was_running:
            self.camera.start()

    def __compose(self):
        """Place children widgets."""
        self.setSizePolicy(self.sizePolicy().Preferred,
                           self.sizePolicy().Preferred)
        self.setBackgroundRole(qtg.QPalette.Shadow)
        self.setWidget(self.__img_view)
        self.setAlignment(qtc.Qt.AlignCenter)
        self.adjustSize()
        self.setWindowTitle(self.camera.name)
        self.__compose_context_menu()

    def __update_title(self):
        """Update the window title with camera information."""
        img_w, img_h = self.image_size.width(), self.image_size.height()
        full_w, full_h = self.camera.sensor_size
        bin_factor = self.camera.binning

        img_size = f"{img_w}x{img_h}"
        if bin_factor == 1 and (img_w, img_h) != (full_w, full_h):
            img_size += f" of {full_w}x{full_h}"
        try:
            scale = self.windowTitle().split(' - ')[1]
        except IndexError:
            scale = ""
        title = f"{self.camera.name} ({img_size})"
        if scale:
            title += f" - {scale}"
        self.setWindowTitle(title)

    def __compose_context_menu(self):
        """Set up the context menu."""
        menu = self.__children["context_menu"]
        # Set up to receive context-menu requests (i.e.,
        # right click), and prepare the context menu
        self.setContextMenuPolicy(qtc.Qt.CustomContextMenu)

        menu.addAction("Reset ROI")
        menu.addAction("Snap image")

        # Flags
        menu.addSeparator()
        act = menu.addAction("Allow setting ROI")
        act.setCheckable(True)
        act.setChecked(
            qtc.Qt.Checked if self.roi_visible else qtc.Qt.Unchecked
            )
        act = menu.addAction("Show on new frames")
        act.setCheckable(True)
        act.setChecked(
            qtc.Qt.Checked if self.show_auto else qtc.Qt.Unchecked
            )
        act = menu.addAction("Stop camera when closed")
        act.setCheckable(True)
        act.setChecked(
            qtc.Qt.Checked if self.stop_on_close else qtc.Qt.Unchecked
            )

        # Device settings
        menu.addSeparator()
        menu.addAction("Properties")

        for action in menu.actions():
            action.setEnabled(self.interactions_enabled)

    def __connect(self):
        """Connect signals."""
        self.image_size_changed.connect(self.__update_title)
        self.__img_view.image_scaling_changed.connect(self.__on_image_scaled)
        self.camera.started.connect(self.__on_camera_started)
        self.camera.error_occurred.connect(self.__on_camera_error)
        self.customContextMenuRequested.connect(self.__show_context_menu)
        self.roi.apply_roi_requested.connect(self.__apply_roi)
        self.__children["context_menu"].triggered.connect(
            self.__on_context_menu_triggered
            )
        self.__children["settings_dialog"].settings_changed.connect(
            self.__on_settings_changed
            )
        self.__children["settings_dialog"].settings_saved.connect(
            self.__on_settings_saved
            )

    @qtc.pyqtSlot(tuple)
    def __on_camera_error(self, _):
        """Close viewer when errors occur."""
        self.close()

    @qtc.pyqtSlot()
    def __on_camera_started(self):
        """React to a start of the camera."""
        mode = self.camera.mode
        if mode == 'triggered':
            connect_to = self.camera.image_processed
            disconnect_from = self.camera.frame_ready
        else:
            connect_to = self.camera.frame_ready
            disconnect_from = self.camera.image_processed

        base.safe_disconnect(disconnect_from, self.__show_image)
        base.safe_connect(connect_to, self.__show_image,
                          type=qtc.Qt.UniqueConnection)
        _, self.__glob['max_intensity'] = self.camera.intensity_limits

    def __on_context_menu_triggered(self, action):
        """React to a selection in the context menu."""
        text = action.text().lower()
        if "allow" in text:
            self.roi_visible = action.isChecked()
            return
        if "show" in text:
            self.show_auto = action.isChecked()
            return
        if "reset" in text:  # reset ROI
            self.__reset_roi()
            return
        if "stop" in text:
            self.stop_on_close = action.isChecked()
            return
        if "snap" in text:
            self.__on_snap_image()
            return
        if "properties" in text:
            self.__children['settings_dialog'].open()
            return
        print(action.text(), ": not implemented yet")

    def __on_image_scaled(self):
        """React to a change of zoom factor."""
        new_scaling = self.__img_view.image_scaling
        self.roi.image_scaling = new_scaling
        self.roi.update_size_limits()

        title, *_ = self.windowTitle().split(' - ')
        title += f' - {100*new_scaling:.1f}%'
        self.setWindowTitle(title)
        self.zoom_changed.emit(new_scaling)

    def __on_settings_changed(self):
        """React to a user press of "Apply" in the settings dialog."""
        _dialog = self.__children["settings_dialog"]
        was_running = self.camera.is_running
        self.camera.settings = _dialog.settings
        if was_running:
            self.camera.start()

    def __on_settings_saved(self, saved):
        """Restore original settings in case of unsaved changes."""
        if saved:
            return
        was_running = self.camera.is_running
        old_settings = self.camera.settings
        old_settings.read_again()
        self.camera.settings = old_settings
        if was_running:
            self.camera.start()

    def __on_snap_image(self):
        """Save current frame to file."""
        image = self.__glob['img_array'].copy()
        fname, _ = qtw.QFileDialog.getSaveFileName(                             # TODO: add default filename with date/time and directory?
            parent=self,
            filter="TIFF Image (*.tiff *.tif)"
            )
        if not fname:  # No file selected
            return
        if not (fname.endswith('.tif') or fname.endswith('.tiff')):
            fname += '.tiff'
        process_info = self.camera.process_info.copy()
        process_info.filename = fname
        process_info.n_frames = 1
        processor = ImageProcessor()
        processor.prepare_to_process(process_info, image)
        processor.process_frame(image)

    def __reset_roi(self):
        """Reset camera ROI to full sensor and update settings."""
        was_running = self.camera.is_running
        if was_running:
            self.camera.stop()
        self.camera.set_roi(no_roi=True)
        self.camera.settings.set("camera_settings", "roi",
                                 str(self.camera.get_roi()))
        self.camera.settings.update_file()
        try:
            self.camera.bad_pixels.apply_roi(no_roi=True)
        except RuntimeError:
            # There were no bad pixels to begin with.
            # Notice that we cannot use if .bad_pixels, as
            # this checks only those in the (now old) ROI
            pass
        if was_running:
            self.camera.start()

    def __show_context_menu(self, position):
        """Show a context menu when right-clicking at position."""
        self.__children["context_menu"].popup(self.mapToGlobal(position))

    @qtc.pyqtSlot(np.ndarray)
    def __show_image(self, img_array):
        """Show the gray-scale image in the img_array numpy.ndarray."""
        if img_array.dtype != np.uint16:
            # images coming from triggering are ready for TIFF,
            # and have an improper data type for viewing as a
            # QImage with Grayscale16 format.
            img_array = np.uint16(img_array)

        # If the camera does not support ROI at the hardware
        # level, show only the portion of the image inside
        # the ROI.
        supports_roi = ('roi' in self.camera.hardware_supported_features
                        or self.camera.get_roi())
        if not supports_roi:
            roi_x, roi_y, roi_w, roi_h = self.camera.roi
            img_array = img_array[roi_y:roi_y+roi_h, roi_x:roi_x+roi_w]

        self.__glob['img_array'] = img_array                                    # TODO: used for snapping. Maybe take from ImageViewer?

        # Notice the swapping of width and height!
        height, width = img_array.shape

        max_int = self.__glob['max_intensity']
        if np.sum(img_array > .95*max_int) > 5:
            print("saturating", timer())                                        # TEMP. Will show overlay

        image = qtg.QImage(img_array, width, height,
                           img_array.strides[0],
                           qtg.QImage.Format_Grayscale16)

        self.__img_view.set_image(image)
        self.image_size = image.size()
        if self.show_auto and not self.isVisible():
            self.show()

    def __zoom(self, direction='in'):
        """Make image larger(smaller) by 25%(20%).

        Zooming is limited to:
            direction=='in':
                At most a factor of 8 of the original image size
            direction=='out':
                At most such that the image occupies >=10% of the screen

        Parameters
        ----------
        direction : {'in', 'out'}
            Whether zooming should enlarge the view of
            the image ('in') or make it smaller ('out')

        Returns
        -------
        None.
        """
        if direction == 'in' and self.__img_view.image_scaling < 8:
            self.scale_image(1.25)
        elif direction == 'out':
            # Limit zoom-out to 10% of the screen size
            img_size = 0.8*self.__img_view.scaled_image_size
            if screen_fraction(self, img_size) < 0.1:
                return
            self.scale_image(0.8)
# pylint: enable=too-many-instance-attributes
