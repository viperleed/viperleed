"""Module hoverannot of viperleed.gui.leedsim.widgets.

Defines the HoverAnnot widget used to mark selected LEED spots.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2020-01-12'
__license__ = 'GPLv3+'

import numpy as np
import PyQt5.QtCore as qtc
import PyQt5.QtGui as qtg
import PyQt5.QtWidgets as qtw

from viperleed.gui.base import check_type
from viperleed.gui.widgetslib import AllGUIFonts


class HoverAnnot(qtw.QWidget):
    # Arrow properties:
    headW = 6  # total width of arrowhead
    headL = 1.2*headW  # total length of arrowhead
    arrowL = 40  # total length of arrow segment

    # quantities used for painting the frame around the text
    framePad = (5, 3)  # amount of space added around the text

    # Painting styles and colors:
    bgColor = qtg.QColor(qtc.Qt.white)
    bgColor.setAlphaF(0.8)
    lineColor = qtg.QColor.fromHsvF(0, 0, 0.15)  # somewhat gray
    lineWidth = 1.5
    gap = 3 + lineWidth

    paintPen = qtg.QPen()
    paintPen.setWidthF(lineWidth)
    paintPen.setBrush(lineColor)
    paintPen.setCapStyle(qtc.Qt.RoundCap)
    paintPen.setJoinStyle(qtc.Qt.RoundJoin)

    def __init__(self, parent=None, text='Test some longer text', center=(500, 500)):
        super().__init__()

        self.setAttribute(qtc.Qt.WA_ShowWithoutActivating)
        self.setAttribute(qtc.Qt.WA_TranslucentBackground)
        # With the next attribute, even if the widget is too large mouse
        # events will be passed on to the widgets directly below it (spatially)
        self.setAttribute(qtc.Qt.WA_TransparentForMouseEvents)
        self.setParent(parent)

        self.setWindowFlags(qtc.Qt.Tool
                            | qtc.Qt.CustomizeWindowHint
                            | qtc.Qt.FramelessWindowHint
                            | qtc.Qt.WindowDoesNotAcceptFocus)

        self.setFont(AllGUIFonts().labelFont)

        # Global position of the point that, together with the arrow head_pos
        # defines the direction of the arrow (from head_pos to center)
        self.__center = np.asarray(center)

        self.__head_pos = None  # Global position of the arrowhead
        self.__tip_to_center = None  # distance from tip to center of text box

        self.__text = ''
        self.frame = qtg.QPainterPath()  # Frame that will contain the text,
        self.framePos = qtc.QPoint()     # its local position (top-left),
        self.frameRadius = 0             # the curvature radius of its corners
        self.textBbox = qtc.QRectF()     # and the boundingRect of text

        self.text = text

    @property
    def center(self):
        """
        Global position of the point towards which the arrow points
        """
        return self.__center

    @center.setter
    def center(self, ctr):
        """
        Set the global position of the point towards which the arrow points to
        ctr
        """
        if (not check_type(ctr, 'arraylike')
            and not isinstance(ctr, (qtc.QPoint, qtc.QPointF))):
            raise TypeError("head_pos must be array-like, QPoint or QPointF."
                            f"Found {type(ctr)} instead.")
        if isinstance(ctr, (qtc.QPoint, qtc.QPointF)):
            ctr = (ctr.x(), ctr.y())

        try:  # this accounts for the case in which head_pos is not defined yet
            delta_head = self.head_pos - self.center
        except TypeError:
            delta_head = None
        self.__center = np.asarray(ctr)

        if delta_head is not None:
            self.head_pos = self.center + delta_head

    @property
    def text(self):
        return self.__text

    @text.setter
    def text(self, text):
        if not isinstance(text, str):
            raise
        self.__text = text
        self.updateFrameOutline()
        self.updateFramePosition()

    @property
    def head_pos(self):
        """
        Global position of the arrowhead
        """
        return self.__head_pos

    @head_pos.setter
    def head_pos(self, pos):
        """
        Set the global position of the arrowhead to pos
        """
        if (not check_type(pos, 'arraylike')
            and not isinstance(pos, (qtc.QPoint, qtc.QPointF))):
            raise TypeError("head_pos must be array-like, QPoint or QPointF."
                            f"Found {type(pos)} instead.")
        if isinstance(pos, (qtc.QPoint, qtc.QPointF)):
            pos = (pos.x(), pos.y())

        if not len(pos) == 2:
            raise ValueError("Invalid position. An array-like with length == "
                             f"2 is needed. Found {len(pos)} instead.")

        self.__head_pos = np.asarray(pos)

        self.updateFramePosition()

    @property
    def tip_pos(self):
        """
        Returns the correct local position of the tip of the arrow such that
        the whole annotation can be painted within the area of the widget.
        Painting of the arrow will start at this position.
        """
        theta = np.radians(self.angle)
        sin_theta = np.sin(theta)
        cos_theta = np.cos(theta)
        frame_h = self.frame.boundingRect().height() + self.lineWidth
        frame_w = self.frame.boundingRect().width() + self.lineWidth

        tip_to_center = self.__tip_to_center + self.lineWidth

        delta_h = frame_h/2 + tip_to_center*sin_theta
        delta_w = frame_w/2 + tip_to_center*cos_theta

        # set the minimum delta to 2 to take care of the arrowhead triangle
        # +1 px for anti-alias, and one more pixel on h to avoid issues when the
        # arrow is in the bottom-left quadrant and exiting from the short side
        delta_w = max(2, delta_w) + 1
        delta_h = max(2, delta_h) + 2

        # notice that I return rounded values so that, later on, the widget can
        # be moved correctly to this position: only int values are acceptable
        # for widget movements
        return qtc.QPointF(np.round(delta_w), np.round(delta_h))

    @property
    def direction(self):
        try:
            dir_vec = self.center - self.head_pos
        except TypeError:
            # one of the two above is undefined
            return None
        return dir_vec/np.linalg.norm(dir_vec)

    @property
    def angle(self):
        """
        Returns the angle in degrees
        """
        dir_vec = self.direction
        if dir_vec is None:
            return None
        return np.degrees(np.arctan2(dir_vec[1], dir_vec[0]))

    def sizeHint(self):
        self.ensurePolished()
        rect_size = self.frame.boundingRect().size()
        delta = self.__tip_to_center + self.lineWidth
        w = np.round(rect_size.width())
        h = np.round(rect_size.height())
        # notice that the return value is slightly over-sized even in the worst
        # case scenario, as one would strictly need to add only w/2 and h/2
        return qtc.QSize(w + delta, h + delta)

    def minimumSizeHint(self):
        return self.sizeHint()

    def updateFrameOutline(self):
        self.frame = qtg.QPainterPath()
        self.textBbox = qtg.QFontMetricsF(self.font()).boundingRect(self.text)

        bbox = qtc.QRectF(self.textBbox)
        bbox.adjust(-self.framePad[0], -self.framePad[1],
                     self.framePad[0], self.framePad[1])
        bbox.moveTo(self.framePad[0]//2, self.framePad[1]//2)
        self.frameRadius = round(0.2*bbox.height())
        self.frame.addRoundedRect(bbox, self.frameRadius, self.frameRadius)

        self.textBbox.adjust(-0.5, -0.5, 0.5, 0.5)
        self.textBbox.moveTo(self.framePad[0] + 1, self.framePad[1] + 1)

    def updateFramePosition(self):
        # Get the right position for self.frame such that the center of the
        # text box is coincident with the prolongation of the arrow

        if self.angle is not None:
            bbox = self.frame.boundingRect()
            w = bbox.width()
            h = bbox.height()

            if np.abs(np.sin(np.radians(self.angle))) > np.radians(1):
                # away from the horizontal direction by more than 1 degree
                dy = -h/2 * np.sign(self.direction[1])
            else:
                dy = 0

            if dy == 0:
                dx = -w/2 * np.sign(self.direction[0])
            elif np.abs(np.cos(np.radians(self.angle))) > np.radians(1):
                # away from the vertical direction by more than 1 degree
                dx = dy / np.tan(np.radians(self.angle))
            else:
                dx = 0

            if np.abs(dx) > w/2:
                # Shallow horizontal angle, the arrow should come out from the
                # vertical edge rather than from the horizontal one
                dx = -w/2 * np.sign(self.direction[0])
                dy = dx * np.tan(np.radians(self.angle))

            # Set up the gap between the arrow and the outer edge of the box,
            # including a correction for the rounded corner
            gap = self.gap  # takes into account line width
            tt = np.abs(np.tan(np.radians(self.angle)))
            ct = np.abs(np.cos(np.radians(self.angle)))
            st = np.abs(np.sin(np.radians(self.angle)))
            r = self.frameRadius
            if tt > max((h - 2*r)/w, h/(w - 2*r)):
                (ddx, ddy) = (0.5*gap*ct, 0.5*gap*(1 + st))
            elif tt < min((h - 2*r)/w, h/(w - 2*r)):
                (ddx, ddy) = (0.5*gap*(1 + ct), 0.5*gap*st)
            else:  # arrow intersects at rounded corner -> different gap
                (ddx, ddy) = (0.8*gap*ct, 0.8*gap*st)

            dx -= ddx*np.sign(self.direction[0])
            dy -= ddy*np.sign(self.direction[1])

            dc = -(self.arrowL)*self.direction + np.array([dx, dy])
            dc -= np.array([2, 1.5])
            self.__tip_to_center = np.linalg.norm(dc)

            deltaCenter = self.tip_pos + qtc.QPointF(*dc)

            self.framePos = deltaCenter - qtc.QPointF(w/2, h/2)

    def paintEvent(self, event):
        painter = qtg.QPainter()

        painter.begin(self)

        painter.setRenderHint(qtg.QPainter.Antialiasing)
        painter.setFont(self.font())
        painter.setPen(self.paintPen)

        self.paintArrow(painter)
        self.paintTextBox(painter)

        painter.end()

    def paintArrow(self, painter):

        painter.save()
        painter.translate(self.tip_pos)
        painter.rotate(self.angle)

        painter.drawPolyline(qtc.QPointF(-self.headL, -self.headW/2),
                             qtc.QPointF(0, 0),
                             qtc.QPointF(-self.headL, self.headW/2))
        painter.drawLine(qtc.QPointF(0, 0), qtc.QPointF(-self.arrowL, 0))
        painter.restore()

    def paintTextBox(self, painter):
        painter.translate(self.framePos)

        painter.fillPath(self.frame, self.bgColor)
        painter.drawPath(self.frame)
        painter.drawText(self.textBbox, self.text)

    def update_position(self):
        """
        Moves the annotation such that self.tip_pos is at the same global
        position as self.head_pos. This function should be called before showing
        the annotation to make sure its position on the screen is correct
        """
        tip = self.tip_pos  # local position
        tip_global = self.mapToGlobal(qtc.QPoint(tip.x(), tip.y()))
        new_pos = self.pos()
        delta = qtc.QPoint(*self.head_pos.round()) - tip_global
        new_pos += delta
        self.move(new_pos)
