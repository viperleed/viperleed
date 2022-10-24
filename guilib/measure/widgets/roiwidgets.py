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

import ast

import numpy as np
from PyQt5 import (QtCore as qtc,
                   QtWidgets as qtw)


# TODO: resizing with the size grips does not fit increments nicely
# TODO: RegionOfInterest.setGeometry will emit roi_changed twice if
#       both moved and resized

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
    roi_changed = qtc.pyqtSignal()

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
        min_dw, min_dh = self.increments
        min_dx, min_dy = self.position_increments
        min_w, min_h = self.minimum

        top_x = round((top_x * scale)/min_dx) * min_dx
        top_y = round((top_y * scale)/min_dy) * min_dy
        width = round((width * scale - min_w) / min_dw) * min_dw + min_w
        height = round((height * scale - min_h) / min_dh) * min_dh + min_h

        return top_x, top_y, width, height

    @image_coordinates.setter
    def image_coordinates(self, new_roi):
        """Set position and size in image coordinates.

        Parameters
        ----------
        new_roi : Sequence of int
            The form should be (top_x, top_y, width, height).
        """
        screen_roi = (round(v * self.image_scaling) for v in new_roi)
        self.setGeometry(qtc.QRect(*screen_roi))

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

    @property
    def position_increments(self):
        """Return the minimum increments in the position in image coordinates."""
        return self.__limits['pos_increments']

    def enterEvent(self, event):         # pylint: disable=invalid-name
        """Change mouse cursor when entering the widget."""
        self.setCursor(qtc.Qt.SizeAllCursor)
        super().enterEvent(event)

    def leaveEvent(self, event):         # pylint: disable=invalid-name
        """Reset mouse cursor when exiting the widget."""
        self.unsetCursor()
        super().leaveEvent(event)

    # pylint: disable=invalid-name
    def mouseDoubleClickEvent(self, event):
        """Reimplement to prevent propagation to parent."""
        event.accept()
    # pylint: enable=invalid-name

    def mouseMoveEvent(self, event):     # pylint: disable=invalid-name
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

    def mousePressEvent(self, event):    # pylint: disable=invalid-name
        """Reimplement mousePressEvent to initiate rubber-band move."""
        self.origin = self.pos()
        self.__drag_origin = event.globalPos()

    def mouseReleaseEvent(self, event):  # pylint: disable=invalid-name
        """Reimplement to prevent propagation to parent."""
        event.accept()

    def moveEvent(self, event):          # pylint: disable=invalid-name
        """Emit roi_changed when moving."""
        super().moveEvent(event)
        if self.isVisible():
            self.roi_changed.emit()

    def resizeEvent(self, event):        # pylint: disable=invalid-name
        """Reimplement to resize the rubber-band."""
        self.__rubberband.resize(self.__normalized_size(self.size()))
        super().resizeEvent(event)
        if self.isVisible():
            self.roi_changed.emit()

    def setGeometry(self, new_rect):     # pylint: disable=invalid-name
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
        # the rectangle. When translating, update self.origin
        if parent:
            new_x = max(0, new_rect.left())
            new_x -= max(new_x + new_rect.width() - parent.width(), 0)
            new_y = max(0, new_rect.top())
            new_y -= max(new_y + new_rect.height() - parent.height(), 0)
            new_rect.translate(new_x - new_rect.x(), new_y - new_rect.y())
            self.origin = new_rect.topLeft()
        super().setGeometry(new_rect)

    def scale(self, delta_scale):
        """Resize by delta_scale increments, in image coordinates.

        Parameters
        ----------
        delta_scale : QtCore.QPoint
            The number of times the width [.x()] and
            height [.y()] should be increased by
            self.increments[0] and self.increments[1]

        Returns
        -------
        None.
        """
        d_width, d_height = self.increments
        delta_size = qtc.QSize(delta_scale.x()*d_width,
                               delta_scale.y()*d_height) * self.image_scaling
        self.resize(self.__normalized_size(self.size() + delta_size))

    def translate(self, delta_increments):
        """Translate self by delta_increments, in image coordinates.

        Parameters
        ----------
        delta_increments : QtCore.QPoint
            The number of times the left [.x()] and top [.y()] should
            be increased by self.position_increments[0]/[1]. If the
            image is scaled down very much, i.e., so that no movement
            would occur, the ROI is translated by delta_increments
            in screen coordinates instead.

        Returns
        -------
        None.
        """
        d_left, d_top = self.position_increments
        _scale = self.image_scaling

        _delta = qtc.QPoint(delta_increments.x() * d_left,
                            delta_increments.y() * d_top)
        delta_pos = _delta * _scale
        if not any(p for p in (delta_pos.x(), delta_pos.y())):
            # Image is scaled down so much that no movement would
            # result from this. Rather use delta_increments as a
            # shift in the current-view coordinates. NB: this may
            # give an inconsistent ROI! This must be corrected
            # from outside.
            delta_pos = delta_increments
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


class ROIEditor(qtw.QWidget):
    """Collection of four SpinBox widgets for setting a ROI."""

    roi_changed = qtc.pyqtSignal(tuple)               # x, y, w, h
    roi_position_changed = qtc.pyqtSignal(int, int)   # x, y
    roi_size_changed = qtc.pyqtSignal(int, int)       # w, h

    def __init__(self, **kwargs):
        """Initialize instance."""
        super().__init__(kwargs.get('parent', None))

        self.__original_roi = tuple()

        (self.__left,
         self.__top,
         self.__width,
         self.__height) = (qtw.QSpinBox() for _ in range(4))

        self.__compose()
        self.__connect()

    @property
    def increments(self):
        """Return d_left, d_top, d_width, d_height."""
        return tuple(w.singleStep() for w in self.__widgets)

    @property
    def maxima(self):
        """Return min_left, min_top, min_width, min_height."""
        return tuple(w.maximum() for w in self.__widgets)

    @property
    def minima(self):
        """Return min_left, min_top, min_width, min_height."""
        return tuple(w.minimum() for w in self.__widgets)

    @property
    def original_roi(self):
        """Return the reference ROI, before changes were applied."""
        return self.__original_roi

    @original_roi.setter
    def original_roi(self, new_roi):
        """Set a reference ROI that is emitted in roi_changed."""
        self.__original_roi = tuple(new_roi)

    @qtc.pyqtSlot()
    def fix_values(self):
        """Fix values in controls to fit increments."""
        # Since set_roi already takes care of this
        # correctly, it's just a matter of calling it
        self.set_roi(self.roi)

    def get_as_string(self):
        """Return the values in the widgets as a string."""
        return str(self.roi)

    def get_roi(self):
        """Return the region of interest as (left, top, width, height)."""
        return tuple(c.value() for c in self.__widgets)

    def set_increments(self, delta_left, delta_top, delta_width, delta_height):
        """Set minimum steps for position and size.

        Parameters
        ----------
        delta_left : int
            Minimum increment of horizontal position of ROI
        delta_top : int
            Minimum increment of vertical position of ROI
        delta_width : int
            Minimum increment of ROI width
        delta_height : int
            Minimum increment of ROI height
        """
        args = (delta_left, delta_top, delta_width, delta_height)
        for delta, widg in zip(args, self.__widgets):
            widg.setSingleStep(delta)

    def set_ranges(self, size_min, size_max):
        """Set range of all widgets.

        Parameters
        ----------
        size_min : tuple of int
            Minimum width and minimum height of ROI.
        size_max : tuple of int
            Maximum width and maximum height of ROI.

        Returns
        -------
        None.
        """
        width_range, height_range = zip(size_min, size_max)
        self.__width.setRange(*width_range)
        self.__height.setRange(*height_range)
        self.__left.setRange(0, width_range[1])
        self.__top.setRange(0, height_range[1])

    def set_from_string(self, string_roi):
        """Set widgets from a string value."""
        if string_roi in (None, 'None', 'none'):
            new_roi = (0, 0, self.__width.maximum(), self.__height.maximum())
        else:
            new_roi = ast.literal_eval(string_roi)
        self.roi = new_roi

    @qtc.pyqtSlot(tuple)
    def set_roi(self, new_roi):
        """Set the values in the widgets from new_roi.

        Parameters
        ----------
        new_roi : tuple of int
            Form should be (left, top, width, height).

        Returns
        -------
        None.
        """
        # Make sure values fit the increments. Min and
        # max should already be coerced correctly
        new_roi = tuple(
            (round((v - _min)/d)) * d + _min
            for v, _min, d in zip(new_roi, self.minima, self.increments)
            )
        old_roi = self.roi

        # Block signals and emit them selectively
        # later. This prevents repeated emissions
        with qtc.QSignalBlocker(self):
            for widg, value in zip(self.__widgets, new_roi):
                widg.setValue(value)

        # Now emit the right signals
        pos_changed = old_roi[:2] != new_roi[:2]
        size_changed = old_roi[2:] != new_roi[2:]
        if pos_changed:
            self.__on_pos_changed()
        if size_changed:
            self.__on_size_changed()
        if pos_changed or size_changed:
            self.__on_roi_changed()

    roi = property(get_roi, set_roi)

    @property
    def __widgets(self):
        """Return left, top, width, and height widgets."""
        return self.__left, self.__top, self.__width, self.__height

    def __compose(self):
        """Place children widgets."""
        layout = qtw.QGridLayout()
        layout.setContentsMargins(0, 0, 0, 0)
        _labels = ('Left', 'Top', 'Width', 'Height')
        for col, (label, widg) in enumerate(zip(_labels, self.__widgets)):
            layout.setColumnStretch(col, 1)  # Same width for all
            layout.addWidget(qtw.QLabel(label),
                             0, col, 1, 1, qtc.Qt.AlignHCenter)
            layout.addWidget(widg, 1, col)
            widg.setAccelerated(True)
        self.setLayout(layout)

    def __connect(self):
        """Connect signals."""
        for widg in self.__widgets:
            widg.valueChanged.connect(self.__on_roi_changed)
            widg.editingFinished.connect(self.fix_values)
        for widg in (self.__left, self.__top):
            widg.valueChanged.connect(self.__on_pos_changed)
        for widg in (self.__width, self.__height):
            widg.valueChanged.connect(self.__on_size_changed)

    @qtc.pyqtSlot()
    @qtc.pyqtSlot(int)
    def __on_pos_changed(self, *_):
        """Emit roi_position_changed."""
        if self.__value_fits_increment(self.sender()):
            self.roi_position_changed.emit(*self.roi[:2])

    @qtc.pyqtSlot()
    @qtc.pyqtSlot(int)
    def __on_size_changed(self, *_):
        """Emit roi_size_changed."""
        if self.__value_fits_increment(self.sender()):
            self.roi_size_changed.emit(*self.roi[2:])

    @qtc.pyqtSlot()
    @qtc.pyqtSlot(int)
    def __on_roi_changed(self, *_):
        """Check that new values fit with increments."""
        if self.__value_fits_increment(self.sender()):
            self.roi_changed.emit(self.roi)

    def __value_fits_increment(self, ctrl):
        """Check that the value in ctrl fits its increment.

        Parameters
        ----------
        ctrl : QSpinBox or None
            The control to check. If the value does not fit,
            the text color is changed to red. Notice that we
            do not force the value to fit the increments
            because this can be very annoying to the user.
            However, users can call the fix_values() method
            to force this on all controls.

        Return
        ------
        value_ok : bool
            Whether the value does fit the increment.
        """
        try:
            idx = self.__widgets.index(ctrl)
        except ValueError:
            return True

        _min, _delta = self.minima[idx], self.increments[idx]
        value = ctrl.value()
        corrected = round((value - _min)/_delta)*_delta + _min

        value_fits = (value == corrected)
        color = qtc.Qt.black if value_fits else qtc.Qt.red
        palette = ctrl.palette()
        palette.setColor(palette.Text, color)
        ctrl.setPalette(palette)
        return value_fits
