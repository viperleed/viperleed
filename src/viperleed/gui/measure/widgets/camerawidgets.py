"""Module camerawidgets of viperleed.gui.measure.widgets.

Contains definition of widgets specific to visualizing frames
from a concrete subclass of CameraABC.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2021-10-29'
__license__ = 'GPLv3+'

import weakref
from copy import deepcopy
from time import perf_counter as timer                                          # TEMP: till "saturating" is fixed

import numpy as np

from PyQt5 import QtCore as qtc
from PyQt5 import QtGui as qtg
from PyQt5 import QtWidgets as qtw

from viperleed.gui.measure import hardwarebase as base
from viperleed.gui.measure.camera import abc
from viperleed.gui.measure.camera.imageprocess import ImageProcessor
from viperleed.gui.widgetslib import screen_fraction


# TODO: ImageViewer.optimum_size is not updated when screen is changed
# TODO: ROI show size in image coordinates as it is resized (tooltip?)
# TODO: ROI context menu: precisely set with coordinates
# TODO: context -- properties
# TODO: use ROI position increments


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
                 show_auto=True, roi_visible=True, **kwargs):
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
            behavior can be changed at any time by setting
            the .show_auto attribute. Default is True.
        roi_visible : bool, optional
            If True, a region of interest can be set and
            manipulated. This behavior can be changed at any
            time by setting the .roi_visible attribute.
            Default is True.
        **kwargs : object
            Other unused optional arguments, passed to
            QScrollArea.__init__

        Raises
        ------
        TypeError
            If camera is not a subclass of CameraABC
        """
        if not isinstance(camera, abc.CameraABC):
            raise TypeError(f"{self.__class__.__name__}: camera argument "
                            "must be a subclass of CameraABC.")

        super().__init__(*args, parent=parent, **kwargs)

        self.__flags = {"stop_on_close": bool(stop_on_close),
                        "show_auto": bool(show_auto),
                        "roi_visible": bool(roi_visible),}
        self.__glob = {"image_size": qtc.QSize(),
                       "camera": camera,
                       "mouse_button": None,
                       "img_array": None,                                       # TODO: may not be necessary
                       "max_intensity": 2**15 - 1,}
        self.__children = {"viewer": ImageViewer(),
                           "context_menu": qtw.QMenu(parent=self)}
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
        if (event.spontaneous() and camera.is_running and self.stop_on_close):
            camera.stop()

            # Disconnect signals that may pop up the window
            # again should a frame arrive in the meantime
            if self.show_auto:
                base.safe_disconnect(camera.image_processed, self.__show_image)
                base.safe_disconnect(camera.frame_ready, self.__show_image)
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

    def mouseMoveEvent(self, event):  # pylint: disable=invalid-name
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

    def mousePressEvent(self, event):  # pylint: disable=invalid-name
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

    def showEvent(self, event):  # pylint: disable=invalid-name
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

    def sizeHint(self):  # pylint: disable=invalid-name
        """Reimplement sizeHint to bind to the underlying image."""
        # Adding (2, 2) ensures that scroll bars are
        # visible only if the image does not fit.
        return self.__img_view.sizeHint() + qtc.QSize(2, 2)

    def wheelEvent(self, event):  # pylint: disable=invalid-name
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

        img_size = f"{img_w}x{img_h}"
        if img_w != full_w or img_h != full_h:
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

    def __connect(self):
        """Connect signals."""
        self.image_size_changed.connect(self.__update_title)
        self.__img_view.image_scaling_changed.connect(self.__on_image_scaled)
        self.camera.started.connect(self.__on_camera_started)
        self.customContextMenuRequested.connect(self.__show_context_menu)
        self.__children["context_menu"].triggered.connect(
            self.__on_context_menu_triggered
            )
        self.roi.apply_roi_requested.connect(self.__apply_roi)

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

    def __on_snap_image(self):
        """Save current frame to file."""
        image = self.__glob['img_array'].copy()
        fname, _ = qtw.QFileDialog.getSaveFileName(                             # TODO: add default filename and directory?
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
        self.camera.bad_pixels.apply_roi(no_roi=True)
        if was_running:
            self.camera.start()

    def __show_context_menu(self, position):
        """Show a context menu when right-clicking at position."""
        self.__children["context_menu"].popup(self.mapToGlobal(position))

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
                At most a factor of 4 of the original image size
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
        if direction == 'in' and self.__img_view.image_scaling < 4:
            self.scale_image(1.25)
        elif direction == 'out':
            # Limit zoom-out to 10% of the screen size
            img_size = 0.8*self.__img_view.scaled_image_size
            if screen_fraction(self, img_size) < 0.1:
                return
            self.scale_image(0.8)
# pylint: enable=too-many-instance-attributes


class ImageViewer(qtw.QLabel):                                                  # TODO: add children (QLabel?) that holds saturation overlay
    """A QLabel that displays images in a CameraViewer.

    The images are scaled with a smooth transform, which makes
    pixels appear blurred at large upscaling factors.

    Attributes
    ----------
    optimum_size : QtCore.QSize
        The size of the image that optimally fits the current
        screen (i.e., occupies roughly 60% of the screen).
    image_scaling : float
        The scaling factor currently used for the image shown.
    scaled_image_size : QtCore.QSize
        Size of the image scaled by self.image_scaling

    Public methods
    --------------
    scale_to_optimum()
        Scale self and pixmap to self.optimum_size.
    scale_to_image()
        Scale self to fully fit the (scaled) pixmap.
    set_image(image)
        Set a QImage to be shown.
    get_scaling_to_fit(size, size_fraction=1)
        Return the scaling factor such that the image
        optimally fits in size*size_fraction.

    Reimplement methods
    -------------------
    sizeHint()
        Return optimal size for self

    Signals
    -------
    image_scaling_changed : no argument
        Emitted every time the scaling factor of the image is changed
    """

    image_scaling_changed = qtc.pyqtSignal()

    def __init__(self, *args, parent=None, **kwargs):
        """Initialize widget."""
        super().__init__(*args, parent=parent, **kwargs)

        self.__image_scaling = 1
        self.optimum_size = qtc.QSize()
        self.__optimum_scaling = -1
        self.__optimal_screen_fraction = 0.6

        self.setBackgroundRole(qtg.QPalette.Base)

        # The QLabel will fill the whole area of the widget containing it,
        # and can scale its contents rather than itself when adjusted.
        self.setSizePolicy(qtw.QSizePolicy.Ignored,
                           qtw.QSizePolicy.Ignored)
        self.setScaledContents(True)

    @property
    def image_scaling(self):
        """Return the currently used pixmap-scaling factor."""
        return self.__image_scaling

    @image_scaling.setter
    def image_scaling(self, new_scaling):
        """Set a new pixmap-scaling factor."""
        if new_scaling != self.image_scaling:
            self.__image_scaling = new_scaling
            self.image_scaling_changed.emit()

    @property
    def scaled_image_size(self):
        """Return the scaled size of self.pixmap()."""
        return self.pixmap().size() * self.image_scaling

    def scale_to_optimum(self):
        """Resize self to the optimal size."""
        if self.optimum_size:
            self.image_scaling = self.__optimum_scaling
            self.resize(self.optimum_size)

    def scale_to_image(self):
        """Resize self to the (scaled) size of the current pixmap."""
        self.resize(self.scaled_image_size)

    def set_image(self, image):
        """Set an image to be shown."""
        self.setPixmap(qtg.QPixmap.fromImage(image))

    def get_scaling_to_fit(self, size, size_fraction=1):
        """Return scaling to best fit pixmap into a fraction of size.

        The scaling is also stored in self.image_scaling.

        Parameters
        ----------
        size : QtCore.QSize
            Reference size into which the pixmap should fit.
        size_fraction : float, optional
            The pixmap is made to fit within size_fraction*size.
            Default is 1.0.

        Returns
        -------
        scaling : float
            Scaling factor as (1.25)**n or (0.8)**m that makes
            the pixmap fit into size_fraction*size.
        """
        scaling = 1
        fraction = max(self.pixmap().width()/size.width(),
                       self.pixmap().height()/size.height())
        while fraction < size_fraction:
            fraction *= 1.25
            scaling *= 1.25
        while fraction > size_fraction:
            fraction *= 0.8
            scaling *= 0.8
        self.image_scaling = scaling
        return scaling

    def sizeHint(self):  # pylint: disable=invalid-name
        """Return optimal size for self."""
        pixmap = self.pixmap()
        if (not pixmap
            or not pixmap.size()
            or not all(s > 0 for s in (pixmap.width(), pixmap.height()))):
            return super().sizeHint()  # invalid pixmap

        if self.optimum_size:
            return self.optimum_size

        # When a pixmap is present, return a scaled version of
        # the pixmap size, trying to keep more or less always
        # the same true size (relative to the screen).
        try:
            screen = self.window().windowHandle().screen()
        except AttributeError:
            # Window does not exist yet
            self.image_scaling = 0.25
            return self.scaled_image_size

        self.get_scaling_to_fit(screen.availableSize(),
                                size_fraction=self.__optimal_screen_fraction)
        self.__optimum_scaling = self.image_scaling
        self.optimum_size = self.scaled_image_size
        return self.optimum_size


class RegionOfInterest(qtw.QWidget):
    """Class for a rectangular rubber-band used for ROI.

    Attributes
    ----------
    image_scaling : float
        Current scaling factor of the image on which this
        rubber-band rectangle is shown. Used to convert
        screen coordinates into image coordinates.
    origin : QtCore.QPoint
        The origin of the ROI while drawing or moving.
    increments : tuple
        Two elements. Minimum change of width and height
        in image coordinates (i.e., pixels).
    limits : tuple
        Three elements, corresponding to self.minimum,
        self.maximum and self.increments.
    maximum : tuple
        Two elements. Maximum width and height in image
        coordinates (i.e., pixels).
    minimum : tuple
        Two elements. Minimum width and height in image
        coordinates (i.e., pixels).

    Reimplement methods
    -------------------
    enterEvent(event)
        Edit mouse cursor when mouse enters the widget.
    leaveEvent(event)
        Reset mouse cursor when mouse exits the widget.
    mouseDoubleClickEvent(event)
        Reimplement to prevent propagation to parent.
    mouseMoveEvent(event)
        Move rubber-band within parent frame.
    mousePressEvent(event)
        Initiate rubber-band movement on mouse click.
    mouseReleaseEvent(event)
        Reimplement to prevent propagation to parent.
    resizeEvent(event)
        Resize also rubber-band.

    Public methods
    --------------
    update_size_limits()
        Update the size limits in current-view coordinates
        using self.limits and self.image_scaling.
    scale(delta_scale)
        Resize self by delta_scale increments in image coordinates.
    translate(delta_pixels)
        Move self by delta_pixels pixels in image coordinates.
    """

    apply_roi_requested = qtc.pyqtSignal()

    def __init__(self, *__args, parent=None, increments=(1, 1), **__kwargs):
        """Initialize RegionOfInterest instance.

        Parameters
        ----------
        parent : QWidget, optional
            The widget on which this region-of-interest rubber-band
            is shown. Typically a CameraViewer. Default is None.
        increments : tuple, optional
            Two elements. Minimum increments of width and height
            in image coordinates (i.e., pixels).

        Returns
        -------
        None.
        """
        super().__init__(parent=parent)
        self.__rubberband = qtw.QRubberBand(qtw.QRubberBand.Rectangle, self)
        self.__limits = {'min': (1, 1),
                         'max': (1000000, 1000000),
                         'increments': increments,
                         'pos_increments': (1, 1)}
        self.__drag_origin = qtc.QPoint(0, 0)
        self.image_scaling = 1
        self.origin = qtc.QPoint(0, 0)
        self.__context_menu = qtw.QMenu(parent=self)

        self.setWindowFlags(qtc.Qt.SubWindow)

        self.__compose()
        self.update_size_limits()

    @property
    def image_coordinates(self):
        """Return position and size in image coordinates.

        Returns
        -------
        top_x : int
            Horizontal position of top-left corner in image
            coordinates (i.e. pixels). Zero is the leftmost
            pixel. This coordinate does not account for an
            already-applied ROI.
        top_y : int
            Vertical position of top-left corner in image
            coordinates (i.e. pixels). Zero is the topmost
            pixel. This coordinate does not account for an
            already-applied ROI.
        width : int
            Width of the selected region of interest in
            image coordinates (i.e. pixels).
        height : int
            Height of the selected region of interest in
            image coordinates (i.e. pixels).
        """
        # Use the bounding box of the rubber-band rather than that of
        # self as the former gets its size normalized on resizeEvent,
        # and will be closer to the correct size when scaled back
        # by 1/self.image_scaling.
        rect = self.__rubberband.rect()

        top_left = self.mapToParent(rect.topLeft())
        top_x, top_y = top_left.x(), top_left.y()
        width, height = rect.width(), rect.height()

        # Scale back to image coordinates (and round as appropriate)
        scale = 1/self.image_scaling
        min_dx, min_dy = self.increments
        min_w, min_h = self.minimum

        top_x = round(top_x * scale)
        top_y = round(top_y * scale)
        width = round((width * scale - min_w) / min_dx) * min_dx + min_w
        height = round((height * scale - min_h) / min_dy) * min_dy + min_h

        return top_x, top_y, width, height

    @property
    def increments(self):
        """Return the smallest change of width/height in image coordinates."""
        return self.__limits['increments']

    @property
    def limits(self):
        """Return .minimum, .maximum., .increments."""
        return self.minimum, self.maximum, self.increments

    @limits.setter
    def limits(self, new_limits):
        """Set .minimum, .maximum., .increments."""
        (self.__limits['min'],
         self.__limits['max'],
         self.__limits['increments'],
         self.__limits['pos_increments']) = new_limits

    @property
    def maximum(self):
        """Return the largest width/height in image coordinates."""
        return self.__limits['max']

    @property
    def minimum(self):
        """Return the smallest width/height in image coordinates."""
        return self.__limits['min']

    def enterEvent(self, event):  # pylint: disable=invalid-name
        """Change mouse cursor when entering the widget."""
        self.setCursor(qtc.Qt.SizeAllCursor)
        super().enterEvent(event)

    def leaveEvent(self, event):  # pylint: disable=invalid-name
        """Reset mouse cursor when exiting the widget."""
        self.unsetCursor()
        super().leaveEvent(event)

    # pylint: disable=invalid-name,no-self-use
    # Disable invalid-name no-self-use as the name
    # and signature must stay unaltered
    def mouseDoubleClickEvent(self, event):
        """Reimplement to prevent propagation to parent."""
        event.accept()
    # pylint: enable=invalid-name,no-self-use

    def mouseMoveEvent(self, event):  # pylint: disable=invalid-name
        """Reimplement mouseMoveEvent to move rubber-band."""
        new_pos = self.origin + event.globalPos() - self.__drag_origin
        # Make sure that self does not go beyond the frame of the parent
        parent = self.parent()
        if parent:
            new_x = max(0, new_pos.x())
            new_x -= max(new_x + self.width() - parent.width(), 0)
            new_y = max(0, new_pos.y())
            new_y -= max(new_y + self.height() - parent.height(), 0)
            new_pos.setX(new_x)
            new_pos.setY(new_y)
        self.move(new_pos)

    def mousePressEvent(self, event):  # pylint: disable=invalid-name
        """Reimplement mousePressEvent to initiate rubber-band move."""
        self.origin = self.pos()
        self.__drag_origin = event.globalPos()

    # pylint: disable=invalid-name,no-self-use
    # Disable invalid-name no-self-use as the name
    # and signature must stay unaltered
    def mouseReleaseEvent(self, event):
        """Reimplement to prevent propagation to parent."""
        event.accept()
    # pylint: enable=invalid-name,no-self-use

    def resizeEvent(self, event):  # pylint: disable=invalid-name
        """Reimplement to resize the rubber-band."""
        self.__rubberband.resize(self.__normalized_size(self.size()))
        super().resizeEvent(event)

    # pylint: disable=invalid-name
    def setGeometry(self, new_rect):
        """Reimplement setGeometry to edit sizes according to increments."""
        # Make sure the new rectangle does not exceed the parent
        # frame. This has to be done on a 'normalized' rectangle
        # with left() < right() and top() < bottom()
        new_size = new_rect.size()
        parent = self.parentWidget()
        if parent:
            norm_rect = new_rect.normalized()
            new_left = max(0, norm_rect.left())
            new_right = min(norm_rect.right(), parent.width())
            new_top = max(0, norm_rect.top())
            new_bottom = min(norm_rect.bottom(), parent.height())
            norm_rect.setLeft(new_left)
            norm_rect.setRight(new_right)
            norm_rect.setTop(new_top)
            norm_rect.setBottom(new_bottom)
            new_size.setWidth(norm_rect.width() * np.sign(new_size.width()))
            new_size.setHeight(norm_rect.height() * np.sign(new_size.height()))

        # Force the size to minimum and maximum constraints
        # as well as to minimum increments in both directions
        new_rect.setSize(self.__normalized_size(new_size))

        # Swap top/left/bottom/right in case the rectangle
        # has width()/height() < 0.
        new_rect = new_rect.normalized()

        # Check again that new_rect fits the parent. This time,
        # however, keep the size constant and rather translate
        # the rectangle.
        if parent:
            new_x = max(0, new_rect.left())
            new_x -= max(new_x + new_rect.width() - parent.width(), 0)
            new_y = max(0, new_rect.top())
            new_y -= max(new_y + new_rect.height() - parent.height(), 0)
            new_rect.translate(new_x - new_rect.x(), new_y - new_rect.y())

        super().setGeometry(new_rect)
    # pylint: enable=invalid-name

    def scale(self, delta_scale):
        """Resize by delta_scale increments, in image coordinates."""
        d_width, d_height = self.increments
        delta_size = qtc.QSize(delta_scale.x()*d_width,
                               delta_scale.y()*d_height) * self.image_scaling
        self.resize(self.__normalized_size(self.size() + delta_size))

    def translate(self, delta_pixels):
        """Translate self by delta_pixels, in image coordinates."""
        delta_pos = delta_pixels*self.image_scaling
        if all(abs(p) < 1e-3 for p in (delta_pos.x(), delta_pos.y())):
            # Image is scaled down so much that no movement
            # would result from this. Rather use delta_pixels
            # as shift in the current-view coordinates
            delta_pos = delta_pixels
        self.move(self.pos() + delta_pos)

    def update_size_limits(self):
        """Update the size limits in viewing pixels.

        Essentially converts self.limits into view pixels by
        scaling with self.image_scaling.

        Returns
        -------
        None.
        """
        self.setMinimumSize(self.image_scaling * qtc.QSize(*self.minimum))
        self.setMaximumSize(self.image_scaling * qtc.QSize(*self.maximum))

    def __compose(self):
        """Place children widgets."""
        # The layout contains only four size grips for resizing
        # (at corners), while the rubber-band is parented directly
        # to self in __init__, and stays on top of all.
        layout = qtw.QGridLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(qtw.QSizeGrip(self), 0, 0,
                         qtc.Qt.AlignLeft | qtc.Qt.AlignTop)
        layout.addWidget(qtw.QSizeGrip(self), 1, 1,
                         qtc.Qt.AlignRight | qtc.Qt.AlignBottom)
        layout.addWidget(qtw.QSizeGrip(self), 0, 1,
                         qtc.Qt.AlignRight | qtc.Qt.AlignTop)
        layout.addWidget(qtw.QSizeGrip(self), 1, 0,
                         qtc.Qt.AlignLeft | qtc.Qt.AlignBottom)
        self.__rubberband.show()
        self.hide()

        # Set up to receive right-clicks, and prepare
        # the corresponding context menu. It only contains
        # a single QAction for applying the ROI
        self.setContextMenuPolicy(qtc.Qt.CustomContextMenu)
        self.customContextMenuRequested.connect(self.__show_context_menu)
        self.__context_menu.triggered.connect(self.__on_context_menu_triggered)
        self.__context_menu.addAction("Apply ROI")
        act = self.__context_menu.addAction("Set ROI coordinates")
        act.setEnabled(False)  # TODO: open up a dialog

    def __normalized_size(self, size):
        """Return a size that fits with increments."""
        # Make sure size is within limits
        signed_width, signed_height = size.width(), size.height()
        width, height = abs(signed_width), abs(signed_height)
        min_w, min_h = (m * self.image_scaling for m in self.minimum)
        max_w, max_h = (m * self.image_scaling for m in self.maximum)

        width = min(max(width, min_w), max_w)
        height = min(max(height, min_h), max_h)

        # Now handle increments
        min_dx, min_dy = (i * self.image_scaling for i in self.increments)
        width = min_w + round((width - min_w)/min_dx) * min_dx
        height = min_h + round((height - min_h)/min_dy) * min_dy

        size.setWidth(round(width) * np.sign(signed_width))
        size.setHeight(round(height) * np.sign(signed_height))
        return size

    def __show_context_menu(self, position):
        """Show a context menu when right-clicking at position."""
        self.__context_menu.popup(self.mapToGlobal(position))

    def __on_context_menu_triggered(self, action):
        """React to a user selection in the context menu."""
        if "apply" in action.text().lower():
            self.apply_roi_requested.emit()
            return
        print("Set requested")
