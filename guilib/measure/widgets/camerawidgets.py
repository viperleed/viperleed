"""Module camerawidgets of viperleed.guilib.measure.

========================================
   ViPErLEED Graphical User Interface
========================================

Contains definition of widgets specific to visualizing frames
from a concrete subclass of CameraABC.

Created: 2021-10-29
Author: Michele Riva
"""


from PyQt5 import (QtCore as qtc,
                   QtWidgets as qtw,
                   QtGui as qtg)

from viperleed.guilib.measure.camera.cameraabc import CameraABC
from viperleed.guilib.widgetslib import screen_fraction


# TODO: zooming is still a little weird when done around the mouse
#       cursor, both zoom in and zoom out.
# TODO: ImageViewer.optimum_size is not updated when screen is changed


class CameraViewer(qtw.QScrollArea):
    """Class for displaying camera frames with scroll bars.

    Zooming can be performed when pressing '+' and '-' or by
    holding down the Control key (Command on Mac) while scrolling.

    Attributes
    ----------
    image_size : QtCore.QSize
        Size in pixel (without scaling) of the image shown.
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
        Handle zooming.
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

    def __init__(self, camera, *args, parent=None, **kwargs):
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
                self.__img_view.scale_to_pixmap()

            elif was_maximized and not self.isMaximized():
                # If it is turned from maximized to normal,
                # scale it to optimum (i.e., same as after shown)
                self.__img_view.scale_to_optimum()
                self.adjustSize()
        super().changeEvent(event)

    def closeEvent(self, event):
        """Extend to stop camera when window is closed by the user."""
        if event.spontaneous() and self.__camera.is_running:
            self.__camera.stop()
        super().closeEvent(event)

    def keyPressEvent(self, event):
        """Extend keyPressEvent for zooming."""
        if event.key() == qtc.Qt.Key.Key_Plus:
            self.__zoom('in')
        elif event.key() == qtc.Qt.Key.Key_Minus:
            self.__zoom('out')
        super().keyPressEvent(event)

    def scale_image(self, by_factor):
        """Scale self and image by the given factor."""
        self.__img_view.image_scaling = (by_factor
                                         * self.__img_view.image_scaling)
        self.__img_view.scale_to_pixmap()

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
        """Extend wheelEvent for zooming."""
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
        title, *_ = self.windowTitle().split(' - ')
        title += f' - {100*new_scaling:.1f}%'
        self.setWindowTitle(title)
        self.zoom_changed.emit(new_scaling)

    def __show_image(self, img_array):
        """Show the gray-scale image in the img_array numpy.ndarray."""
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
    scale_to_pixmap()
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
            self.resize(self.optimum_size)
        self.image_scaling = self.__optimum_scaling

    def scale_to_pixmap(self):
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
