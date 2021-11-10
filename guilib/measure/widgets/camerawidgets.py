"""Module camerawidgets of viperleed.guilib.measure.widgets.

========================================
   ViPErLEED Graphical User Interface
========================================

Contains definition of widgets specific to visualizing frames
from a concrete subclass of CameraABC.

Created: 2021-10-29
Author: Michele Riva
"""

import numpy as np

from PyQt5 import (QtCore as qtc,
                   QtWidgets as qtw,
                   QtGui as qtg)

from viperleed.guilib.measure.camera.abc import CameraABC
from viperleed.guilib.widgetslib import screen_fraction


# TODO: zooming is still a little weird when done around the mouse
#       cursor, both zoom in and zoom out.
# TODO: ImageViewer.optimum_size is not updated when screen is changed


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
    screen : QtGui.QScreen
        The screen on which this widget is currently shown.

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

    shown = qtc.pyqtSignal()
    image_size_changed = qtc.pyqtSignal()
    zoom_changed = qtc.pyqtSignal(float)  # new zoom factor

    def __init__(self, camera, *args, parent=None, stop_on_close=True,
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
        stop_on_close : bool
            Stop the camera when this widget is closed.
        **kwargs : object
            Other unused optional arguments, passed to
            QScrollArea.__init__
        """
        if not isinstance(camera, CameraABC):
            raise TypeError(f"{self.__class__.__name__}: camera argument "
                            "must be a subclass of CameraABC.")

        super().__init__(*args, parent=parent, **kwargs)

        self.__img_view = ImageViewer()
        self.__img_size = qtc.QSize()
        self.__camera = camera
        self.__stop_on_close = bool(stop_on_close)
        self.__roi = RegionOfInterest(parent=self.__img_view)

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
    def image_size(self):
        """Return the true size in pixel of the image currently shown."""
        return self.__img_size

    @image_size.setter
    def image_size(self, new_size):
        """Store the true size in pixel of the image currently shown."""
        if isinstance(new_size, tuple):
            new_size = qtc.QSize(*new_size)
        if (not new_size
                or any(s <=0 for s in (new_size.width(), new_size.height()))):
            return  # Invalid size

        if new_size != self.image_size:
            self.__img_size = new_size
            self.__img_view.updateGeometry()
            self.updateGeometry()
            self.__img_view.optimum_size = qtc.QSize()
            self.__img_view.adjustSize()
            self.adjustSize()
            self.image_size_changed.emit()

    @property
    def roi(self):
        """Return the RegionOfInterest object of self."""
        return self.__roi

    @property
    def screen(self):
        """Return a the QScreen on which self is shown."""
        try:
            screen = self.window().windowHandle().screen()
        except AttributeError:
            screen = None
        return screen

    def changeEvent(self, event):
        """Extend changeEvent to react to window maximization."""
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
        super().changeEvent(event)

    def closeEvent(self, event):
        """Extend to stop camera when window is closed by the user."""
        if (event.spontaneous()
            and self.__camera.is_running
                and self.__stop_on_close):
            self.__camera.stop()
        super().closeEvent(event)

    def keyPressEvent(self, event):
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
        super().keyPressEvent(event)

    def mouseDoubleClickEvent(self, event):
        """Prevent parent propagation if ROI is visible."""
        if not self.roi.isVisible():
            super().mouseDoubleClickEvent(event)

    def mouseMoveEvent(self, event):
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
        if not self.roi.isVisible():
            super().mouseMoveEvent(event)
            return
        roi_pos = self.roi.parent().mapFromGlobal(event.globalPos())
        new_roi = qtc.QRect(self.roi.origin, roi_pos).normalized()
        if event.modifiers() & qtc.Qt.ShiftModifier:
            # Make square region
            new_side = min(new_roi.width(), new_roi.height())
            min_side = max(*self.roi.minimum)*self.roi.image_scaling
            new_side = round(max(min_side, new_side))
            new_roi.setWidth(new_side)
            new_roi.setHeight(new_side)
        self.roi.setGeometry(new_roi)

    def mousePressEvent(self, event):
        """Reimplement mousePressEvent to begin drawing a ROI."""
        if not self.roi.isVisible():
            self.roi.show()
        self.roi.origin = self.roi.parent().mapFromGlobal(event.globalPos())
        self.roi.setGeometry(qtc.QRect(self.roi.origin, qtc.QSize()))

    def mouseReleaseEvent(self, event):
        """Extend mouseReleaseEvent to handle ROI drawing."""
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

    def showEvent(self, event):
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

    def sizeHint(self):
        """Reimplement sizeHint to bind to the underlying image."""
        # Adding (2, 2) ensures that scroll bars are
        # visible only if the image does not fit.
        return self.__img_view.sizeHint() + qtc.QSize(2, 2)

    def wheelEvent(self, event):
        """Extend wheelEvent for zooming while Control is pressed."""
        if not event.modifiers() & qtc.Qt.ControlModifier:
            super().wheelEvent(event)
            return
        # Mouse wheel turned with Ctrl down --> zoom
        direction = 'in' if event.angleDelta().y() > 0 else 'out'
        self.__zoom(direction=direction)

    def __adjust_scroll_bars(self, by_factor):
        """Adjust the position of scrollbars when zooming."""
        if self.__img_view.underMouse():
            # Zoom in at the current mouse position
            mouse_offset = (self.mapFromGlobal(qtg.QCursor.pos())
                            - self.widget().pos()) * (by_factor - 1)
            for bar, mouse_offs in zip(
                    (self.horizontalScrollBar(), self.verticalScrollBar()),
                    (mouse_offset.x(), mouse_offset.y())
                    ):
                bar.setValue(bar.value() + mouse_offs)
            return

        # Zoom centered when the mouse cursor is somewhere else
        for bar in (self.horizontalScrollBar(), self.verticalScrollBar()):
            bar.setValue(int(by_factor * bar.value()
                             + ((by_factor - 1) * bar.pageStep()/2)))

    def __compose(self):
        """Place children widgets."""
        self.setSizePolicy(self.sizePolicy().Preferred,
                           self.sizePolicy().Preferred)
        self.setBackgroundRole(qtg.QPalette.Shadow)
        self.setWidget(self.__img_view)
        self.setAlignment(qtc.Qt.AlignCenter)
        self.adjustSize()
        self.setWindowTitle(self.__camera.name)

    def __connect(self):
        """Connect signals."""
        self.__img_view.image_scaling_changed.connect(self.__on_image_scaled)
        self.__camera.started.connect(self.__on_camera_started)

    def __on_camera_started(self):
        """React to a start of the camera."""
        mode = self.__camera.mode
        if mode == 'triggered':
            connect_to = self.__camera.image_processed
            disconnect_from = self.__camera.frame_ready
        else:
            connect_to = self.__camera.frame_ready
            disconnect_from = self.__camera.image_processed
        try:
            disconnect_from.disconnect(self.__show_image)
        except TypeError:
            # Not connected
            pass
        try:
            connect_to.connect(self.__show_image, qtc.Qt.UniqueConnection)
        except TypeError:
            # Already connected
            pass

    def __on_image_scaled(self):
        """React to a change of zoom factor."""
        new_scaling = self.__img_view.image_scaling
        self.roi.image_scaling = new_scaling
        self.roi.update_size_limits()

        title, *_ = self.windowTitle().split(' - ')
        title += f' - {100*new_scaling:.1f}%'
        self.setWindowTitle(title)
        self.zoom_changed.emit(new_scaling)

    def __show_image(self, img_array):
        """Show the gray-scale image in the img_array numpy.ndarray."""
        if img_array.dtype != np.uint16:
            # images coming from triggering are ready for TIFF,
            # and have an improper data type for viewing as a
            # QImage with Grayscale16 format.
            img_array = np.uint16(img_array)

        width, height = img_array.shape

        # Notice the swapping of width and height!
        image = qtg.QImage(img_array, height, width,
                           img_array.strides[0],
                           qtg.QImage.Format_Grayscale16)

        self.__img_view.set_image(image)
        self.image_size = image.size()
        if not self.isVisible():
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


class ImageViewer(qtw.QLabel):
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

    def sizeHint(self):
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

    def __init__(self, *args, parent=None, increments=(1, 1), **kwargs):

        super().__init__(parent=parent)
        self.__rubberband = qtw.QRubberBand(qtw.QRubberBand.Rectangle, self)
        self.__min = (1, 1)
        self.__max = (1000000, 1000000)
        self.__incr = increments
        self.__drag_origin = qtc.QPoint(0, 0)
        self.image_scaling = 1
        self.origin = qtc.QPoint(0, 0)

        self.setWindowFlags(qtc.Qt.SubWindow)

        self.__compose()
        self.update_size_limits()

    @property
    def increments(self):
        """Return the smallest change of width/height in image coordinates."""
        return self.__incr

    @property
    def limits(self):
        """Return .minimum, .maximum., .increments."""
        return self.minimum, self.maximum, self.increments

    @limits.setter
    def limits(self, new_limits):
        """Set .minimum, .maximum., .increments."""
        self.__min, self.__max, self.__incr = new_limits

    @property
    def maximum(self):
        """Return the largest width/height in image coordinates."""
        return self.__max

    @property
    def minimum(self):
        """Return the smallest width/height in image coordinates."""
        return self.__min

    def enterEvent(self, event):
        """Change mouse cursor when entering the widget."""
        self.setCursor(qtc.Qt.SizeAllCursor)
        super().enterEvent(event)

    def leaveEvent(self, event):
        """Reset mouse cursor when exiting the widget."""
        self.unsetCursor()
        super().leaveEvent(event)

    def mouseDoubleClickEvent(self, event):
        """Reimplement to prevent propagation to parent."""
        event.accept()

    def mouseMoveEvent(self, event):
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

    def mousePressEvent(self, event):
        """Reimplement mousePressEvent to initiate rubber-band move."""
        self.origin = self.pos()
        self.__drag_origin = event.globalPos()

    def mouseReleaseEvent(self, event):
        """Reimplement to prevent propagation to parent."""
        event.accept()

    def resizeEvent(self, event):
        """Reimplement to resize the rubber-band."""
        self.__rubberband.resize(self.size())

    def scale(self, delta_scale):
        """Resize by delta_scale increments, in image coordinates."""
        d_width, d_height = self.increments
        delta_size = qtc.QSize(delta_scale.x()*d_width,
                               delta_scale.y()*d_height) * self.image_scaling
        self.resize(self.size() + delta_size)

    def translate(self, delta_pixels):
        """Translate self by delta_pixels, in image coordinates."""
        delta_pos = delta_pixels*self.image_scaling
        if all(p == 0 for p in (delta_pos.x(), delta_pos.y())):
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
        self.setSizeIncrement(self.image_scaling * qtc.QSize(*self.increments))

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
