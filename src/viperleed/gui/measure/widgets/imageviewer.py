"""Module imageviewer of viperleed.gui.measure.widgets.

Defines the ImageViewer class, a QLabel that can be used to display
images. This class used to be part of the camerawidgets module.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2022-10-21'
__license__ = 'GPLv3+'

from PyQt5 import (QtCore as qtc,
                   QtWidgets as qtw,
                   QtGui as qtg)


class ImageViewer(qtw.QLabel):
    """A QLabel that displays images.

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

    Reimplemented methods
    ---------------------
    sizeHint()
        Return optimal size for self
    paintEvent(event)
        Paint only the visible portion of the image.

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

        # The QLabel will not use sizeHint, but fill
        # the whole area of the widget containing it
        self.setSizePolicy(qtw.QSizePolicy.Ignored,
                           qtw.QSizePolicy.Ignored)

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

    def set_image(self, image, overlay_mask=None):
        """Set an image to be shown."""
        pixmap = qtg.QPixmap.fromImage(image)
        if overlay_mask is not None:
            # Use overlay_mask as clip region, and paint red inside
            # on top of the pixmap generated from image. Assume
            # overlay_mask is a QRegion.
            painter = qtg.QPainter(pixmap)
            painter.setClipRegion(overlay_mask)
            brush = qtg.QBrush(qtc.Qt.red)
            painter.fillRect(pixmap.rect(), brush)
            painter.end()
        self.setPixmap(pixmap)

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

    def paintEvent(self, event):         # pylint: disable=invalid-name
        """Paint only the visible portion of the image."""
        # Pretty much rewrite the relevant parts of qlabel.cpp
        painter = qtg.QPainter(self)
        self.drawFrame(painter)

        if not self.__has_pixmap:
            return

        pix = self.pixmap()
        if not self.isEnabled():
            style = self.style()
            options = qtw.QStyleOption()
            options.initFrom(self)
            pix = style.generateIconPixmap(qtg.QIcon.Disabled, pix, options)

        # Rather than painting via style, as is done in qlabel.cpp,
        # use directly the painter. This has an overload that allows
        # to select only a portion of the image to be painted. The
        # area that requires repaint is event.rect(). It has to be
        # rescaled to image coordinates to pick the correct region.
        image_rect = qtc.QRect(event.rect().topLeft() / self.image_scaling,
                               event.rect().size() / self.image_scaling)
        painter.drawPixmap(event.rect(), pix, image_rect)

    def sizeHint(self):                  # pylint: disable=invalid-name
        """Return optimal size for self."""
        if not self.__has_pixmap:
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

    @property
    def __has_pixmap(self):
        """Return whether self has a valid pixmap."""
        return self.pixmap() and not self.pixmap().isNull()
