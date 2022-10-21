"""Module roiwidgets of viperleed.guilib.measure.widgets.

========================================
   ViPErLEED Graphical User Interface
========================================

Created: 2022-10-21
Author: Michele Riva

Classes defined in this module
------------------------------
RegionOfInterest
    A QWidget that shows up as a rubber-band rectangle, and can be
    used to display a region of interest rectangle. The rubber-band
    sizes can be limited to size limits that the ROI may have, and
    is always contained in the parent widget. Scaling and movement
    can be done in image coordinates (physical pixels) rather than
    in screen coordinates. This class used to be part of camerawidgets.
"""

import numpy as np
from PyQt5 import (QtCore as qtc,
                   QtWidgets as qtw)


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
            is shown. Typically an ImageViewer. Default is None.
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

    # pylint: disable=invalid-name
    # Disable invalid-name no-self-use as the name
    # and signature must stay unaltered
    def mouseDoubleClickEvent(self, event):
        """Reimplement to prevent propagation to parent."""
        event.accept()
    # pylint: enable=invalid-name

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

    # pylint: disable=invalid-name
    # Disable invalid-name no-self-use as the name
    # and signature must stay unaltered
    def mouseReleaseEvent(self, event):
        """Reimplement to prevent propagation to parent."""
        event.accept()
    # pylint: enable=invalid-name

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
