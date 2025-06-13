"""Module roi of viperleed.gui.measure.widgets.

Classes defined in this module
------------------------------
Corner
    An enumeration of QtCore.Qt.Corner with some convenience functions.
RegionOfInterest
    A QWidget that shows up as a rubber-band rectangle, and can be
    used to display a region of interest rectangle. The rubber-band
    sizes can be limited to size limits that the ROI may have, and
    is always contained in the parent widget. Scaling and movement
    can be done in image coordinates (physical pixels) rather than
    in screen coordinates. This class used to be part of camerawidgets.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2022-10-21'
__license__ = 'GPLv3+'

from enum import Enum
import sys

from PyQt5 import QtCore as qtc
from PyQt5 import QtGui as qtg
from PyQt5 import QtWidgets as qtw


def _fixed_corner_rect(size, reference_rect, fixed_corner, normalize=False):
    """Return a new QRect with a given size and a fixed corner.

    Parameters
    ----------
    size : QtCore.QSize
        The size of new_rect.
    reference_rect : QtCore.QRect
        The reference rectangle from which one should
        take the position of the corner to be kept fixed.
    fixed_corner : Corner
        The corner that should be kept fixed. It should not
        be Corner.NO_CORNER
    normalize : bool, optional
        Whether the rectangle returned should be normalized, i.e.,
        should have left < right and top < bottom. Notice that the
        normalization procedure is different from the one that Qt
        implements, as that one CHANGES THE SIZE (by 2 pixels) for
        invalid rectangles. The implementation used here will keep
        fixed_corner fixed in place. Default is False.

    Returns
    -------
    new_rect : QtCore.QRect
        The rectangle with size size and the fixed_corner
        fixed to the corresponding one of reference_rect

    Raises
    ------
    TypeError
        If fixed_corner is not an instance of Corner
    ValueError
        If fixed_corner is Corner.NO_CORNER
    """
    if not isinstance(fixed_corner, Corner):
        raise TypeError("fixed_corner should be a 'Corner', not "
                        f"{type(fixed_corner).__name__!r}")
    if fixed_corner is Corner.NO_CORNER:
        raise ValueError("Cannot handle Corner.NO_CORNER")

    new_rect = qtc.QRect(qtc.QPoint(), size)
    if fixed_corner is Corner.TOP_LEFT:
        new_rect.moveTopLeft(reference_rect.topLeft())
    elif fixed_corner is Corner.TOP_RIGHT:
        new_rect.moveTopRight(reference_rect.topRight())
    elif fixed_corner is Corner.BOTTOM_LEFT:
        new_rect.moveBottomLeft(reference_rect.bottomLeft())
    else:   # corner == Corner.BOTTOM_RIGHT
        new_rect.moveBottomRight(reference_rect.bottomRight())

    if not normalize or new_rect.isValid():
        return new_rect
    return _normalized_rect_fixed_size(new_rect, fixed_corner)


def _normalized_rect_fixed_size(rect, fixed_corner):
    """Return a normalized version of rect with size fixed."""
    # Decide what needs to be swapped, and what to fix
    _left, _right = rect.left(), rect.right()
    _top, _bot = rect.top(), rect.bottom()
    fix_left, fix_top = fixed_corner.at_left, fixed_corner.at_top

    if _left > _right and fix_left:
        _left, _right = _right, _left - 2
    elif _left > _right and not fix_left:  # fix right edge
        _left, _right = _right + 2, _left

    if _top > _bot and fix_top:
        _top, _bot = _bot, _top - 2
    elif _top > _bot and not fix_top:  # fix bottom edge
        _top, _bot = _bot + 2, _top

    return qtc.QRect(qtc.QPoint(_left, _top), qtc.QPoint(_right, _bot))


def _get_max_edge_positions(widget, available_geo=None):
    """Return the extrema for edge positions given an available_geo.

    Parameters
    ----------
    widget : QWidget
        The widget for which the edge positions should be returned.
    available_geo : QtCore.QRect, optional
        The rectangle outlining the region of space available for
        the expansion of widget. If not given or None, automatically
        find the available geometry for widget. Default is None.

    Returns
    -------
    min_left : int
        Minimum horizontal position of the left edge
    max_right : int
        Maximum horizontal position of the right edge
    min_top : int
        Minimum vertical position of the top edge
    max_bot : int
        Maximum vertical position of the bottom edge
    """
    if available_geo is None:
        available_geo, *_ = _get_resize_limits(widget)

    # Detect sizes of frame elements
    frame_geo = widget.frameGeometry()
    geo = widget.geometry()
    title_bar_h = max(0, geo.top() - frame_geo.top())
    bottom_h = max(0, frame_geo.bottom() - geo.bottom())
    left_right_width = max(0, (frame_geo.width() - geo.width()) / 2)

    # Determine the extreme positions for left/right/top/bottom
    # edges, in screen coordinates relative to the "parent"
    min_left = available_geo.left() + left_right_width
    max_right = available_geo.right() - left_right_width
    min_top = available_geo.top() + title_bar_h
    max_bot = available_geo.bottom() - bottom_h

    return min_left, max_right, min_top, max_bot


_NOBAR = qtc.Qt.ScrollBarAlwaysOff


def _get_resize_limits(widget):
    """Return the available geometry and whether size is constrained.

    Parameters
    ----------
    widget : QWidget
        The widget for which available geometry and size
        constraints should be returned

    Returns
    -------
    available_geo : QtCore.QRect
        The geometry available to widget for resizing.
    has_size_constraint : tuple
        Two bool elements, representing whether the width and height
        of widget have a geometry-related constraint. The only case
        in which one of the elements is False is if widget is part
        of a QAbstractScrollArea with potentially active scroll bars.
    """
    has_size_constraint = (True, True)
    if widget.isWindow():
        available_geo = qtw.QDesktopWidget().availableGeometry(widget)
    else:
        parent = widget.parentWidget()
        available_geo = parent.contentsRect()
        if isinstance(parent, qtw.QAbstractScrollArea):
            has_size_constraint = (
                parent.horizontalScrollBarPolicy() == _NOBAR,
                parent.verticalScrollBarPolicy() == _NOBAR
                )
    return available_geo, has_size_constraint


def _normalized_size(size, limits):
    """Return a modified size that fits limits.

    Parameters
    ----------
    size : QtCore.QSize
        The size to be normalized.
    limits : dict
        Should have (at least) the following key/values:
        {'min': (min_width, min_height),
         'max': (max_width, max_height),
         'increments': (delta_width, delta_height)}

    Return
    ------
    size : QtCore.QSize
        The modified size
    """
    ((min_w, min_h),
     (max_w, max_h),
     (min_dw, min_dh)) = (limits[k] for k in ('min', 'max', 'increments'))

    w_sign = 1 if size.width() >=0 else -1
    h_sign = 1 if size.height() >=0 else -1

    # Make sure size does not exceed the min and max
    _, width, _ = sorted((min_w, abs(size.width()), max_w))
    _, height, _ = sorted((min_h, abs(size.height()), max_h))

    # And handle increments
    width = min_w + round((width - min_w) / min_dw) * min_dw
    height = min_h + round((height - min_h) / min_dh) * min_dh

    return qtc.QSize(round(width) * w_sign, round(height) * h_sign)


def _corner_pos(rect, corner):
    """Return the position of Corner corner for the given rect."""
    if corner is Corner.TOP_LEFT:
        return rect.topLeft()
    if corner is Corner.TOP_RIGHT:
        return rect.topRight()
    if corner is Corner.BOTTOM_LEFT:
        return rect.bottomLeft()
    if corner is Corner.BOTTOM_RIGHT:
        return rect.bottomRight()
    raise ValueError(f"Invalid corner {corner}")


def _rect_fits(widget, rect):
    """Return whether rect fits the region available to widget."""
    min_left, max_right, min_top, max_bot = _get_max_edge_positions(widget)
    return all((rect.left() >= min_left, rect.right() <= max_right,
                rect.top() >= min_top, rect.bottom() <= max_bot))


class Corner(Enum):
    """Enumeration of QtCore.Qt.Corner."""

    BOTTOM_LEFT = qtc.Qt.BottomLeftCorner
    BOTTOM_RIGHT = qtc.Qt.BottomRightCorner
    TOP_LEFT = qtc.Qt.TopLeftCorner
    TOP_RIGHT = qtc.Qt.TopRightCorner
    NO_CORNER = None

    @staticmethod
    def from_position(at_left, at_bottom):
        """Return a Corner from knowledge of its position."""
        if at_left:
            return Corner.BOTTOM_LEFT if at_bottom else Corner.TOP_LEFT
        return Corner.BOTTOM_RIGHT if at_bottom else Corner.TOP_RIGHT

    @property
    def opposite(self):
        """Return the corner opposite to self."""
        if self is Corner.BOTTOM_LEFT:
            return Corner.TOP_RIGHT
        if self is Corner.BOTTOM_RIGHT:
            return Corner.TOP_LEFT
        if self is Corner.TOP_LEFT:
            return Corner.BOTTOM_RIGHT
        if self is Corner.TOP_RIGHT:
            return Corner.BOTTOM_LEFT
        return Corner.NO_CORNER

    @property
    def at_left(self):
        """Return whether this is a left-side corner."""
        return self in (Corner.BOTTOM_LEFT, Corner.TOP_LEFT)

    @property
    def at_top(self):
        """Return whether this is a bottom-side corner."""
        return self in (Corner.TOP_LEFT, Corner.TOP_RIGHT)

# pylint: disable=too-many-public-methods
# Most methods (11) are reimplementations of base-class ones.
# Also, it looks like pylint is counting @property decorated
# methods as public methods (rather than attributes)
class RegionOfInterest(qtw.QWidget):
    """Class for a rectangular rubber-band used for ROI."""

    roi_changed = qtc.pyqtSignal()

    def __init__(self, *, parent=None, increments=(1, 1)):
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
        self.__limits = {'min': (0, 0),
                         'max': (1000000, 1000000),
                         'increments': increments,
                         'pos_increments': (1, 1)}
        self.__image_scaling = 1
        self.__rect = qtc.QRect()   # The actual ROI (in screen coordinates)

        # delta_max dict:
        #   "offset": origin shift when resizing left/right/top/bottom
        #   "expansion": max. resizing towards left/right/top/bottom
        self.__drawing_info = {
            "origin": qtc.QPoint(),             # Where drawing starts
            "moving_corner": Corner.NO_CORNER,  # Corner of grip moved
            "delta_max": {},   # Resizing "offset" and max "expansion"
            "drag_origin": qtc.QPoint(),    # Initial pos when dragged
            }
        self.__edit_flags = {"moved": False, "resized": False}

        self.is_being_edited = False

        self.setWindowFlags(qtc.Qt.SubWindow)
        self.__compose()
        self.__update_size_limits()

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
        top_x, top_y = self.__rect.x(), self.__rect.y()
        width, height = self.__rect.width(), self.__rect.height()

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
    def image_scaling(self):
        """Return the current scaling factor of the image."""
        return self.__image_scaling

    @image_scaling.setter
    def image_scaling(self, new_image_scaling):
        """Set a new scaling factor for the image."""
        if new_image_scaling == self.image_scaling:
            return
        self.__image_scaling = new_image_scaling
        self.__update_size_limits()
        scaled_lims = self.__scaled_limits
        for idx in range(self.layout().count()):
            grip = self.layout().itemAt(idx).widget()
            grip.set_limits(scaled_lims)

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

        self.__update_size_limits()
        scaled_lims = self.__scaled_limits
        for idx in range(self.layout().count()):
            grip = self.layout().itemAt(idx).widget()
            grip.set_limits(scaled_lims)

    @property
    def maximum(self):
        """Return the largest width/height in image coordinates."""
        return self.__limits['max']

    @property
    def minimum(self):
        """Return the smallest width/height in image coordinates."""
        return self.__limits['min']

    @property
    def origin(self):
        """Return the origin in parent coordinates."""
        return self.__drawing_info["origin"]

    @origin.setter
    def origin(self, new_origin):
        """Set a new origin ensuring it fits increments."""
        # Notice that we cannot really check whether new_origin
        # is too close to one of the parent's edges as this
        # depends on which direction we will draw the ROI
        min_dx, min_dy = self.__scaled_limits['pos_increments']
        if isinstance(new_origin, tuple):
            new_origin = qtc.QPoint(*new_origin)
        new_x = round(new_origin.x() / min_dx) * min_dx
        new_y = round(new_origin.y() / min_dy) * min_dy
        self.__drawing_info["origin"] = qtc.QPoint(round(new_x), round(new_y))

    @property
    def position_increments(self):
        """Return the minimum increments in the position in image coordinates."""
        return self.__limits['pos_increments']

    def enterEvent(self, event):         # pylint: disable=invalid-name
        """Change mouse cursor when entering the widget."""
        self.setCursor(qtc.Qt.SizeAllCursor)
        super().enterEvent(event)

    def eventFilter(self, obj, event):   # pylint: disable=invalid-name
        """Return whether event of obj should be filtered.

        Currently, no filtering occurs. This method is used to
        initiate/finish_draw when a size grip is clicked/released.

        Parameters
        ----------
        obj : QObject
            The object whose event is checked
        event : QEvent
            The event being filtered

        Returns
        -------
        filtered : bool
            Whether event of obj should not be handled.
        """
        if isinstance(obj, _IncrementalSizeGrip):
            if (event.type() == event.MouseButtonPress
                    and event.button() == qtc.Qt.LeftButton):
                self.__drawing_info["moving_corner"] = moving = obj.corner
                fixed = moving.opposite
                self.initiate_draw(_corner_pos(self.geometry(), fixed))
            elif event.type() == event.MouseButtonRelease:
                self.finish_draw()
        return super().eventFilter(obj, event)

    def finish_draw(self):
        """Signal that drawing the rubber-band rectangle is over."""
        self.is_being_edited = False
        self.__drawing_info["moving_corner"] = Corner.NO_CORNER
        self.__drawing_info["delta_max"] = {}

    def initiate_draw(self, origin):
        """Begin drawing the rubber-band rectangle at origin."""
        self.origin = origin
        self.is_being_edited = True
        self.__calculate_max_resize()

    def leaveEvent(self, event):         # pylint: disable=invalid-name
        """Reset mouse cursor when exiting the widget."""
        self.unsetCursor()
        super().leaveEvent(event)

    # pylint: disable=invalid-name
    def mouseDoubleClickEvent(self, event):
        """Override to prevent propagation to parent."""
        event.accept()
    # pylint: enable=invalid-name

    def mouseMoveEvent(self, event):     # pylint: disable=invalid-name
        """Override to move rubber-band consistently with increments."""
        new_pos = (self.origin
                   + event.globalPos() - self.__drawing_info["drag_origin"])
        # Make sure that self does not go beyond the frame of the parent
        parent = self.parentWidget()
        if parent:
            new_x = max(0, new_pos.x())
            new_x -= max(new_x + self.width() - parent.width(), 0)
            new_y = max(0, new_pos.y())
            new_y -= max(new_y + self.height() - parent.height(), 0)
            new_pos.setX(new_x)
            new_pos.setY(new_y)

        # Now make sure only allowed position increments are used
        incr_left, incr_top = self.__scaled_limits['pos_increments']
        delta_pos = new_pos - self.origin
        d_left = round(delta_pos.x() / incr_left) * incr_left
        d_top = round(delta_pos.y() / incr_top) * incr_top

        delta_pos.setX(round(d_left))
        delta_pos.setY(round(d_top))

        self.move(self.origin + delta_pos)

    def mousePressEvent(self, event):    # pylint: disable=invalid-name
        """Override to initiate rubber-band move."""
        self.origin = self.pos()
        self.__drawing_info["drag_origin"] = event.globalPos()
        self.is_being_edited = True

    def mouseReleaseEvent(self, __event):  # pylint: disable=invalid-name
        """Override to conclude rubber-band move."""
        self.is_being_edited = False

    def move(self, new_pos):
        """Move self to new_pos."""
        self.__edit_flags['moved'] = True
        self.__rect.moveTo(new_pos)
        super().move(new_pos)

    def moveEvent(self, event):          # pylint: disable=invalid-name
        """Emit roi_changed when moving and not being resized."""
        self.__edit_flags['moved'] = False
        super().moveEvent(event)
        # Since resize always happens after move, testing the attribute
        # ensures only one roi_changed is emitted for each setGeometry
        if self.isVisible() and not self.__edit_flags['resized']:
            self.roi_changed.emit()

    def resize(self, new_size):
        """Resize self and rubber-band to new_size after normalizing it."""
        _norm_size = self.__normalized_size(new_size)
        self.__edit_flags['resized'] = True
        super().resize(_norm_size)
        self.__rubberband.resize(_norm_size)
        self.__rect.setSize(_norm_size)

    def resizeEvent(self, event):        # pylint: disable=invalid-name
        """Emit roi_changed when resized."""
        self.__edit_flags['resized'] = False
        super().resizeEvent(event)
        if self.isVisible():
            self.roi_changed.emit()

    def scale(self, delta_scale):
        """Resize by delta_scale increments, in image coordinates.

        Parameters
        ----------
        delta_scale : QtCore.QPoint
            The number of times the width [.x()] and height
            [.y()] should be increased by self.increments[0]
            and self.increments[1]. If delta_scale used as
            above would result in no change at all (usually
            when the image is scaled down a lot), delta_scale
            is used as an increment in screen coordinates.

        Returns
        -------
        None.
        """
        if delta_scale.isNull():  # No scaling
            return
        d_width, d_height = self.increments
        _scale = self.image_scaling
        delta_size = qtc.QSize(delta_scale.x()*d_width,
                               delta_scale.y()*d_height) * _scale

        new_size = self.__normalized_size(self.size() + delta_size)
        if new_size == self.size():
            # Size is not changed, probably scaled down too much. Use
            # delta_scale as an increment in screen coordinates rather
            # than in image coordinates. If the reason size did not
            # change because it's already at its min/max, this will
            # not change it either.
            delta_size = qtc.QSize(delta_scale.x(), delta_scale.y()) / _scale
            new_size = self.__normalized_size(self.size() + delta_size)
        if not _rect_fits(self, qtc.QRect(self.pos(), new_size)):
            # New size would exceed the available geometry
            return
        self.resize(new_size)

    def setGeometry(self, new_rect):     # pylint: disable=invalid-name
        """Edit position and size consistently with increments."""
        norm_rect = self.__normalized_rect(new_rect)

        moved = norm_rect.topLeft() != self.__rect.topLeft()
        self.__edit_flags['moved'] = moved
        self.__edit_flags['resized'] = norm_rect.size() != self.__rect.size()

        self.__rect = norm_rect
        super().setGeometry(norm_rect)
        self.__rubberband.resize(norm_rect.size())

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
        if delta_pos.isNull():
            # Image is scaled down so much that no movement would
            # result from this. Rather use delta_increments as a
            # shift in the current-view coordinates. Make it such
            # that it still fits the position increments, though:
            delta_pos = qtc.QPoint(
                round(delta_increments.x() / (_scale * d_left)) * d_left,
                round(delta_increments.y() / (_scale * d_top)) * d_top
                )
        new_pos = self.pos() + delta_pos
        if not _rect_fits(self, qtc.QRect(new_pos, self.__rect.size())):
            # New position would move outside the available geometry
            return
        self.move(new_pos)

    @property
    def __scaled_limits(self):
        """Return a version of self.__limits in screen coordinates."""
        return {k: tuple(vi * self.image_scaling for vi in v)
                for k, v in self.__limits.items()}

    # pylint: disable=too-many-locals
    # Disabled because it seems clearer to treat each edge explicitly
    def __calculate_max_resize(self):
        """Calculate and store resizing constraints in all directions."""
        ((min_w, min_h),
         (max_w, max_h),
         (d_w, d_h),
         (d_x, d_y)) = self.__scaled_limits.values()
        min_left, max_right, min_top, max_bot = _get_max_edge_positions(self)

        # Determine the maximum expansions in all directions,
        # starting at the minimum allowed values
        expand_left = expand_right = min_w
        expand_top = expand_bot = min_h

        # Shift the origin if it is too close to the edges, otherwise
        # determine how much we can expand as multiples of d_w/d_h.
        # Notice that the "expand" will be used for sizes, while
        # edge positions are computed in QRect as x2 = x1 + w - 1.
        # If we have to shift, it will be done when drawing happens
        # as only then we can decide which shift to pick, depending
        # on whether we go left, right, top or bottom.
        delta_x_left = delta_x_right = delta_y_top = delta_y_bot = 0

        # Too far left when moving left (w == -min_w)
        edge_left = self.origin.x() - min_w - 1
        while round(edge_left + delta_x_left) < min_left:
            delta_x_left += d_x

        # Too far right when moving right (w == min_w)
        edge_right = self.origin.x() + min_w - 1
        while round(edge_right + delta_x_right) > max_right:
            delta_x_right -= d_x

        # Too far top when moving up (h == -min_h)
        edge_top = self.origin.y() - min_h - 1
        while round(edge_top + delta_y_top) < min_top:
            delta_y_top += d_y

        # Too far bottom when moving down (h == min_h)
        edge_bot = self.origin.y() + min_h - 1
        while round(edge_bot + delta_y_bot) > max_bot:
            delta_y_bot -= d_y

        # In range
        w_max_right = min(max_right - self.origin.x() + 1, round(max_w))
        while round(expand_right + d_w) <= w_max_right:
            expand_right += d_w

        w_max_left = min(self.origin.x() - min_left + 1, round(max_w))
        while round(expand_left + d_w) <= w_max_left:
            expand_left += d_w

        h_max_top = min(self.origin.y() - min_top + 1, round(max_h))
        while round(expand_top + d_h) <= h_max_top:
            expand_top += d_h

        h_max_bot = min(max_bot - self.origin.y() + 1, round(max_h))
        while round(expand_bot + d_h) <= h_max_bot:
            expand_bot += d_h

        self.__drawing_info["delta_max"] = {
            "offset": (delta_x_left, delta_x_right, delta_y_top, delta_y_bot),
            "expansion": (expand_left, expand_right, expand_top, expand_bot)
            }
    # pylint: enable=too-many-locals

    def __compose(self):
        """Place children widgets."""
        # The layout contains only four size grips for resizing
        # (at corners), while the rubber-band is parented directly
        # to self in __init__, and stays on top of all.
        layout = qtw.QGridLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(_IncrementalSizeGrip(self), 0, 0,
                         qtc.Qt.AlignLeft | qtc.Qt.AlignTop)
        layout.addWidget(_IncrementalSizeGrip(self), 1, 1,
                         qtc.Qt.AlignRight | qtc.Qt.AlignBottom)
        layout.addWidget(_IncrementalSizeGrip(self), 0, 1,
                         qtc.Qt.AlignRight | qtc.Qt.AlignTop)
        layout.addWidget(_IncrementalSizeGrip(self), 1, 0,
                         qtc.Qt.AlignLeft | qtc.Qt.AlignBottom)
        for idx in range(layout.count()):
            grip = layout.itemAt(idx).widget()
            grip.installEventFilter(self)
        self.__rubberband.show()
        self.hide()

    def __normalized_rect(self, rect):
        """Return a QRect that fits all size/position constraints.

        Returns
        -------
        new_rect : QtCore.QRect
            The rectangle that complies with all constraints: does
            not exceed the parent (if assigned), fits size limits
            and increments, fits position increments
        """
        # Make sure we know which corner is moving; work always with a
        # "valid" rect (left <= right and top <= bot). If not valid,
        # fix the corner opposite to the one that moves. This ensures
        # that we always obtain consistent results between size-grips
        # resize and other drawing operations.
        moving = self.__drawing_info["moving_corner"]
        if moving is Corner.NO_CORNER:
            # Not being resized by size grips.
            # Figure out what is moving
            at_left = rect.width() < 0
            at_bottom = rect.height() >= 0
            moving = Corner.from_position(at_left, at_bottom)
            rect = _fixed_corner_rect(rect.size(), rect, moving.opposite,
                                      normalize=True)

        # Clip the new_rect size to the maximum resizing distances,
        # if available, and offset in the right direction, if needed.
        # This ensures the rect is always inside the parent.
        size = rect.size()
        offset = qtc.QPoint()
        if self.__drawing_info["delta_max"]:
            idx_x = 0 if moving.at_left else 1
            idx_y = 2 if moving.at_top else 3
            _offs, _maxs = self.__drawing_info["delta_max"].values()
            size = qtc.QSize(min(size.width(), round(_maxs[idx_x])),
                             min(size.height(), round(_maxs[idx_y])))
            offset = qtc.QPoint(round(_offs[idx_x]), round(_offs[idx_y]))

        # Force the size to minimum and maximum constraints
        # as well as to minimum increments in both directions
        norm_size = self.__normalized_size(size)

        # Translate in the right direction if necessary
        rect.translate(offset)

        # Apply the size change fixing the non-moving corner
        return _fixed_corner_rect(norm_size, rect, moving.opposite)

    def __normalized_size(self, size):
        """Return a size that fits with increments."""
        return _normalized_size(size, self.__scaled_limits)

    def __update_size_limits(self):
        """Update the size limits in viewing pixels.

        Essentially converts self.limits into view pixels by
        scaling with self.image_scaling. This has any effect
        only when self is inserted in a layout.

        Returns
        -------
        None.
        """
        self.setMinimumSize(self.image_scaling * qtc.QSize(*self.minimum))
        self.setMaximumSize(self.image_scaling * qtc.QSize(*self.maximum))


# Most methods here are brute rewriting of the ones in qsizegrip.cpp
class _IncrementalSizeGrip(qtw.QSizeGrip):
    """A QSizeGrip that handles finite position/size increments."""

    def __init__(self, parent):
        """Initialize instance."""
        super().__init__(parent)
        self.__got_mouse_press = False
        self.__origin = qtc.QPoint()
        self.__original_rect = qtc.QRect()
        self.__delta_max = (sys.maxsize, sys.maxsize)
        self.__limits = {'min': (0, 0),             # width, height
                         'max': (sys.maxsize, sys.maxsize),
                         'increments': (1, 1),
                         'pos_increments': (1, 1)}  # left, top

    @property
    def corner(self):
        """Return the Corner of this size grip."""
        # The very same implementation as in qsizegrip.cpp but
        # there the property is available only for the private
        top_level = self.__top_level_widget(self)
        pos = self.mapTo(top_level, qtc.QPoint())
        at_left = pos.x() <= top_level.width() / 2
        at_bottom = pos.y() >= top_level.height() / 2
        return Corner.from_position(at_left, at_bottom)

    def mousePressEvent(self, event):    # pylint: disable=invalid-name
        """Initiate a resize of the top-level parent widget on left button."""
        if event.button() != qtc.Qt.LeftButton:
            super().mousePressEvent(event)
            return
        self.__origin = event.globalPos()
        top_level = self.__top_level_widget(self)
        self.__original_rect = top_level.geometry()
        self.__got_mouse_press = True
        self.__calculate_max_resize()

    def mouseMoveEvent(self, event):     # pylint: disable=invalid-name
        """Resize the top-level parent widget."""
        if not event.buttons() & qtc.Qt.LeftButton:
            super().mouseMoveEvent(event)
            return

        top_level = self.__top_level_widget(self)
        if (not self.__got_mouse_press
                or top_level.testAttribute(qtc.Qt.WA_WState_ConfigPending)):
            return

        delta = event.globalPos() - self.__origin
        dx_max, dy_max = self.__delta_max
        new_height = self.__original_rect.height()
        new_width = self.__original_rect.width()
        _at_top, _at_left = self.corner.at_top, self.corner.at_left
        if _at_top:
            new_height -= max(delta.y(), dy_max)
        else:
            new_height += min(delta.y(), dy_max)

        if _at_left:
            new_width -= max(delta.x(), dx_max)
        else:
            new_width += min(delta.x(), dx_max)

        new_size = qtc.QSize(new_width, new_height)
        new_size = qtw.QLayout.closestAcceptableSize(top_level, new_size)

        # Now apply size limits and constraints
        new_size = _normalized_size(new_size, self.__limits)

        # Move new_rect to always keep the opposite corner fixed
        new_rect = _fixed_corner_rect(new_size, self.__original_rect,
                                      self.corner.opposite)
        top_level.setGeometry(new_rect)

    def mouseReleaseEvent(self, event):  # pylint: disable=invalid-name
        """Emit released if the event is accepted."""
        if event.button() != qtc.Qt.LeftButton:
            super().mouseReleaseEvent(event)
            return
        self.__got_mouse_press = False
        self.__origin = qtc.QPoint()

    def moveEvent(self, event):          # pylint: disable=invalid-name
        """Pick cursor shape according to location."""
        if not self.__origin.isNull():   # Not being dragged
            super().moveEvent(event)
            return
        if self.corner in (Corner.TOP_LEFT, Corner.BOTTOM_RIGHT):
            cursor = qtc.Qt.SizeFDiagCursor
        else:
            cursor = qtc.Qt.SizeBDiagCursor
        self.setCursor(cursor)

    def paintEvent(self, __event):       # pylint: disable=invalid-name
        """Draw size grip widget."""
        # This reimplementation is needed to draw the grip with
        # the correct orientation depending on self.corner
        painter = qtg.QPainter(self)
        opt = qtw.QStyleOptionSizeGrip()
        opt.initFrom(self)
        opt.corner = self.corner.value
        self.style().drawControl(self.style().CE_SizeGrip, opt, painter, self)

    def set_limits(self, new_limits):
        """Set size and position limits for this size grip."""
        self.__limits = new_limits

    @staticmethod
    def __top_level_widget(widg):
        """Return the widget that this grip is controlling."""
        # This is the same implementation as in qsizegrip.cpp
        # for the static qt_sizegrip_topLevelWidget function
        while (widg and not widg.isWindow()
               and widg.windowType() != qtc.Qt.SubWindow):
            widg = widg.parentWidget()
        return widg

    def __calculate_max_resize(self):
        """Store maximum resizing constraints in self.__delta_max."""
        top_level = self.__top_level_widget(self)
        (available_geo,
         (w_constrained, h_constrained)) = _get_resize_limits(top_level)

        (min_left, max_right,
         min_top, max_bot) = _get_max_edge_positions(top_level, available_geo)

        _at_top, _at_left = self.corner.at_top, self.corner.at_left
        sign_y = -1 if _at_top else 1
        dy_max = sys.maxsize
        if not _at_top and h_constrained:
            dy_max = max_bot - self.__original_rect.bottom()
        elif _at_top and h_constrained:
            dy_max = self.__original_rect.top() - min_top

        dx_max = sys.maxsize
        sign_x = -1 if _at_left else 1
        if _at_left and w_constrained:
            dx_max = self.__original_rect.left() - min_left
        elif not _at_left and w_constrained:
            dx_max = max_right - self.__original_rect.right()

        # Here we could also trim dx_max/dy_max by increments,
        # but since we're using the grips only for the ROI, we
        # can have the ROI take care of this.
        self.__delta_max = (dx_max * sign_x, dy_max * sign_y)
