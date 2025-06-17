"""Module basewidgets of viperleed.gui.widgets.

Contains custom QWidget subclasses that do not fit other modules
of the widgets package.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2020-01-12'
__license__ = 'GPLv3+'

import numpy as np
from PyQt5 import QtCore as qtc
from PyQt5 import QtGui as qtg
from PyQt5 import QtWidgets as qtw

from viperleed.gui.widgets.lib import AllGUIFonts
from viperleed.gui.widgets.lib import drawText


class PainterMatrix(qtw.QWidget):  ## --> use it in a special QPushButton with a tick-box? or checkable?
    def __init__(self, matrix, fs=None, color=(0, 0, 0, 1), parent=None):
        if not isinstance(matrix, (np.ndarray, list)):
            raise

        if not len(np.shape(matrix)) == 2:
            # check that it only has 2 dimensions
            raise

        super().__init__(parent)
        self.setSizePolicy(qtw.QSizePolicy.Expanding, qtw.QSizePolicy.Expanding)
        self.setParent(parent)

        self.color = color

        #----- set up the font -----#
        font = AllGUIFonts().mathFont
        if fs is None:
            fs = AllGUIFonts().mathFont.pointSize()
        font.setPointSize(fs)
        self.setFont(font)

        # Some details about the CMU Serif font (fs = fontSize in points): #
        # lineSpacing() = 1.58*fs
        # leading() = 0
        # descent() = fs/3
        # ascent() = 3.74/3 * fs
        # parenthesis(L,R) tight height = 4*descent() = 4/3 * fs
        # parenthesis(L,R) tight width = .3085 * fs ~ 25/81 * fs
        # parenthesis(L) bbox width = .4405 * fs
        # parenthesis(R) bbox width = 3.47/9 * fs
        # parenthesis(L,R) baseline is 3/4 h from top of tight box,
        #                              or exactly fs
        # parenthesis(L,R) bottom point is at baseline+descent
        # parenthesis(L) horizontalAdvance() = .5174*fs = 6.2*h/16
        #                from left edge of bbox

        # Process the numeric matrix into an array of strings with the same
        # shape also replacing '-' with the right character '\u2212'
        # everywhere
        self.matrix = self.textMatrix(matrix)

        #----- Initialize the QTransforms for painting the matrix -----#
        self.initTransforms()

    def textMatrix(self, matrix):
        """
        Returns a string version of the matrix in which occurrences of '-'
        (minus) are replaced by the Unicode minus sign (u+2212)
        """
        strMatr = np.array(matrix).astype(str)
        fixMinus = np.array([el.replace('-','\u2212')
                             for el in strMatr.ravel()])

        return fixMinus.reshape(strMatr.shape)

    def colSep(self):
        """
        Returns the correct distance (in pixels) between:
        * the right-hand side of the bounding box of elements of column j
        AND
        * the left-hand side of the bounding box of elements of column j+1
        """
        fs = self.font().pointSize()
        return 13.6 + ((1.4)**5 + (1.6*fs/27)**5)**(1/5)

    def colWidths(self):
        """
        Returns the largest width (in pixels) of each column that fits all
        entries in the column
        """
        return [max(self.textWidth(el) for el in col)
                 for col in self.matrix.transpose()]

    def interLine(self):
        """
        Returns the interline spacing, i.e., the distance in pixels between the
        baselines of two adjacent lines
        """
        fs = self.font().pointSize()
        return 1.73*fs

    def parenthesisVGeo(self):
        """
        Returns the geometrical parameters (in pixels) that allow correct
        positioning of the parentheses. These are returned as a tuple
        (h, deltaTop), where
        - h = total height of the parenthesis
        - deltaTop = Distance between top of parenthesis tightRect and baseline
          of top row
        """

        # Analysis for 2 lines:
        # - htight = (3.57418+-0.01904)*fs
        # - interline = 1.73 * fs
        # --> htight - interline = (1.84418+-0.01904)*fs
        # - dTop = (1.27362+-0.00951)*fs ~ 3.82/3 * fs
        # --> dBottom = htight - interline - dTop = 0.57 * fs

        nrows = len(self.matrix)
        fs = self.font().pointSize()

        # - Distance between top of parenthesis tightRect and baseline
        #   of top row
        deltaTop = 3.82*fs/3

        # - Contribution of the number of rows
        deltaRows = self.interLine()*(nrows-1)

        # - Distance between baseline of bottom row and bottom of
        #   parenthesis tightRect
        deltaBottom = 0.57*fs

        h = deltaTop + deltaRows + deltaBottom

        return (h, deltaTop)

    def initTransforms(self):
        if not hasattr(self, 'lParTransform'):   # left parenthesis
            self.lParTransform = qtg.QTransform()
        if not hasattr(self, 'elemsTransform'):  # matrix elements
            self.elemsTransform = qtg.QTransform()
        if not hasattr(self, 'rParTransform'):   # right parenthesis
            self.rParTransform = qtg.QTransform()

        (h, dTop) = self.parenthesisVGeo()   # h and dTop are in the
                                             # unscaled coordinate system
        (sx, sy) = self.scaleFactor()

        # 1) Get the correct starting position for the painter to get the
        #    left parenthesis drawn at the right location
        dxPar = -.1*h * (sx/sy)
        # NB:
        # .1*h comes from
        #    width - tightWidth = .4405*fs - .3085*fs
        #                       = .132*fs = .132*3/4 * h
        #                       = approx. .1*h
        # it is then scaled by (sx/sy) because the parenthesis is narrower
        # than normal (h/sy is the height of a normal parenthesis)
        dyPar = 3*h/4
        self.lParTransform.reset()
        self.lParTransform.translate(dxPar, dyPar)  # First translate to the
                                                    # baseline of the big
                                                    # parenthesis.
        self.lParTransform.scale(sx, sy)            # Then scale

        # 2) Set painter position to start drawing the first element
        #    of the first line of the matrix, relative to the one found
        #    above for the left parenthesis
        # - horizontal: shift by the scaled horizontal advance of the
        #               parenthesis + dxPar
        # - vertical: shift by the vertical distance between the
        #             tightBoundingRect of the large bracket and the
        #             baseline of the first line
        self.elemsTransform.reset()
        self.elemsTransform.translate((6.2*h/16) * (sx/sy) + dxPar, dTop)

        # 3) Set painter position for right parenthesis
        # - horizontal: start at left side of top line, and add the total
        #   width of the matrix, considering both the columns with text and
        #   the inter-column separation.
        # - vertical: same as left parenthesis
        colWs = self.colWidths()
        ncols = len(colWs)
        dxPar = self.elemsTransform.dx()
        dxPar += sum(colWs)
        dxPar += self.colSep()*(ncols-1)

        self.rParTransform.reset()
        # first translate to the baseline:
        self.rParTransform.translate(dxPar, dyPar)
        # then scale
        self.rParTransform.scale(sx, sy)

    def scaleFactor(self):
        """
        Returns the correct scaling factors that make a standard parenthesis
        assume the correct appearance
        """
        nrows = len(self.matrix)
        fs = self.font().pointSize()

        h = self.parenthesisVGeo()[0]  # correct tight height
        h0 = 4/3 * fs
        sy = h/h0

        # w comes from a fit of the dependence of the width of the
        # parenthesis on its height for a 2-line matrix: ~ linear at
        # beginning, saturates to 19 px (from which 1px of anti-alias should
        # be removed)
        w = (19-1) / (1 + (100/h)**10)**(1/12)
        w0 = 25/81 * fs  # tight width from fit
        sx = w/w0

        return (sx, sy)

    def textWidth(self, text):
        fm = qtg.QFontMetricsF(self.font())
        return fm.boundingRect(text).width()

    def minimumSizeHint(self):
        return self.sizeHint()

    def sizeHint(self):
        h = self.parenthesisVGeo()[0]
        (sx, sy) = self.scaleFactor()
        w = self.rParTransform.dx() + (3.47*h/12) * (sx/sy)

        # Note that I'm using ceil() for the return value because sizeHint
        # should return a QSize(int, int), but w and h are floats and would
        # be typecast without rounding otherwise
        return qtc.QSize(int(np.ceil(w)), int(np.ceil(h)))

    def paintEvent(self, event):
        self.adjustSize()
        painter = qtg.QPainter()

        pen = qtg.QPen(qtg.QColor.fromRgbF(*self.color))

        painter.begin(self)
        painter.setFont(self.font())
        painter.setPen(pen)
        painter.setRenderHint(qtg.QPainter.Antialiasing)
        painter.setRenderHint(qtg.QPainter.TextAntialiasing)

        #-- Draw first the entries of the matrix
        self.paintElements(painter)

        #-- Then draw the parentheses (with correct translation/scaling)
        drawText(painter, '(', self.lParTransform, combine = True)
        drawText(painter, ')', self.rParTransform, combine = True)

        painter.end()

    def paintElements(self, painter):
        widths = self.colWidths()

        painter.save()
        painter.setWorldTransform(self.elemsTransform, True)
        for row in self.matrix:  # run through rows
            painter.save()   # save at each run
            for (el, w) in zip(row, widths):
                painter.save()   # save each element
                dx = w - self.textWidth(el)
                painter.translate(dx, 0)
                drawText(painter, el)
                painter.restore()  # restore to go back to the left side
                painter.translate(self.colSep() + w, 0)  # and go to the
                                                         # next column
            painter.restore()  # go back to the left-hand side of the line
            painter.translate(0, self.interLine())  # and switch to the new
                                                    # line
        painter.restore()  # and back to the initial state
