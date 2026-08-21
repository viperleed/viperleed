"""Module roieditor of viperleed.gui.measure.widgets.

This module defines the ROIEditor class, a QWidget that allows
graphically editing position and size of a region of interest
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2022-10-21'
__license__ = 'GPLv3+'

import ast

from PyQt5 import QtCore as qtc
from PyQt5 import QtWidgets as qtw


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
    def current_roi_valid(self):
        """Return whether the current ROI is acceptable."""
        return all(getattr(w, "valid_value", False) for w in self.__widgets)

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
        """Set a reference ROI."""
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
        if not self.__value_fits_limits(self.sender()):
            return
        self.roi_position_changed.emit(*self.roi[:2])

        # Changing the position may make the size exceed its maximum
        for widg in self.__widgets[2:]:
            self.__value_fits_limits(widg, check_max=True)

    @qtc.pyqtSlot()
    @qtc.pyqtSlot(int)
    def __on_size_changed(self, *_):
        """Emit roi_size_changed."""
        if not self.__value_fits_limits(self.sender()):
            return
        self.roi_size_changed.emit(*self.roi[2:])

        # Changing the size may make the position exceed its maximum
        for widg in self.__widgets[:2]:
            self.__value_fits_limits(widg, check_max=True)

    @qtc.pyqtSlot()
    @qtc.pyqtSlot(int)
    def __on_roi_changed(self, *_):
        """Check that new values fit with increments."""
        if self.__value_fits_limits(self.sender()):
            self.roi_changed.emit(self.roi)

    def __value_fits_limits(self, ctrl, check_max=False):
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
        corrected = round((value - _min) / _delta) * _delta + _min

        if check_max:
            # The maximum allowed value depends on the control
            _buddies = {self.__top: self.__height, self.__height: self.__top,
                        self.__left: self.__width, self.__width: self.__left}
            _max = self.maxima[idx]
            corrected = min(corrected, _max - _buddies[ctrl].value())

        value_fits = (value == corrected)
        color = qtc.Qt.black if value_fits else qtc.Qt.red
        palette = ctrl.palette()
        palette.setColor(palette.Text, color)
        ctrl.setPalette(palette)
        ctrl.valid_value = value_fits
        return value_fits
