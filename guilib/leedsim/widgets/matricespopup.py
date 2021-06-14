"""
=======================================
   ViPErLEED Graphical User Interface
=======================================
 *** module guilib.leedsim.widgets ***

Created: 2020-01-12
Author: Michele Riva

Blah blah TODO
"""

import PyQt5.QtCore as qtc
import PyQt5.QtWidgets as qtw

from viperleed import guilib as gl


class MatricesPopup(qtw.QWidget):
    def __init__(self, matrices, colors,
                 fs=gl.AllGUIFonts().mathFont.pointSize(),
                 parent=None):
        self.parent = parent
        super().__init__()

        if colors is not None:  # there's more than one domain
            self.matrices = [gl.PainterMatrix(matrix, color=color, fs=fs)
                             for (matrix, color) in zip(matrices, colors)]
            # add one matrix to be displayed in black when the user toggles
            # the visibility of the domains
            self.firstDomMatrix = gl.PainterMatrix(matrices[0], fs=fs)
        elif len(matrices) == 1:
            self.matrices = [gl.PainterMatrix(matrices[0], fs=fs)]
        else:  # some input error: colors are missing
            raise

        self.shown = False  # used for choosing whether is should pop up at
                            # the standard position or not
        self.initPopup()

    def initPopup(self):
        # The next line makes the widget not steal focus from the main
        # window when it is shown
        self.setAttribute(qtc.Qt.WA_ShowWithoutActivating)
        self.setParent(self.parent)  # This should come before the next
                                     # lines to have the window stay on top
                                     # of the main, but not on top of other
                                     # applications

        #-------- Edit the appearance of the window --------#
        # Qt.Tool: thinner frame
        # Qt.CustomizeWindowHint: Remove titlebar, icon, and
        #                         maximize/minimize/close buttons
        # Qt.WindowStaysOnTopHint: always keep it on top - do not use.
        #                          Rather parent it to the MainWindow so
        #                          that it stays only on top of it (and not
        #                          on top of other apps)
        # Qt.WindowDoesNotAcceptFocus: prevent switching focus to the popup
        self.setWindowFlags(qtc.Qt.Tool
                            | qtc.Qt.CustomizeWindowHint
                            | qtc.Qt.WindowDoesNotAcceptFocus)
        gl.editStyleSheet(self, 'background-color: white;')

        # The next attribute is used in the mousePressEvent and
        # mouseMoveEvent for dragging the window around
        self.dragPosition = self.frameGeometry().topLeft()

        #--------- Prepare layout for the matrices ---------#
        nMatrices = len(self.matrices)  # this is 1, 2, 3, 4, 6, 8, or 12
                                        # TRUE?
        if nMatrices in range(1, 4):  # 1 -- 3
            nRows = 1
        elif nMatrices in range(4, 10):  # 4 -- 8
            nRows = 2
        else:
            nRows = 3
        nCols = int(nMatrices/nRows)

        lay = qtw.QGridLayout()
        [lay.addWidget(matrix, index//nCols, index % nCols, qtc.Qt.AlignCenter)
         for (index, matrix) in enumerate(self.matrices)]

        if hasattr(self, 'firstDomMatrix'):
            # More than one domain. Add anyway the matrix of the first domain
            # for changing its color when domains are toggled. Initially
            # hide (toggling the button will show it).
            lay.addWidget(self.firstDomMatrix, 0, 0, qtc.Qt.AlignCenter)
            self.firstDomMatrix.hide()

        sizes = [(matrix.sizeHint().width(), matrix.sizeHint().height())
                 for matrix in self.matrices]
        (widths, heights) = zip(*sizes)

        [lay.setColumnMinimumWidth(col, max(widths))
         for col in range(lay.columnCount())]
        [lay.setRowMinimumHeight(row, max(heights))
         for row in range(lay.rowCount())]
        lay.setSizeConstraint(qtw.QLayout.SetFixedSize)  # This makes the window
                                                         # not re-sizable
        self.setLayout(lay)
        self.layout().activate()

    def mousePressEvent(self, event):
        if event.button() == qtc.Qt.LeftButton:
            self.dragPosition = (event.globalPos()
                                - self.frameGeometry().topLeft())

    def mouseReleaseEvent(self, event):
        event.accept()

    def mouseDoubleClickEvent(self, event):
        event.accept()

    def mouseMoveEvent(self, event):
        if event.buttons() == qtc.Qt.LeftButton:
            self.move(event.globalPos() - self.dragPosition)

    def updateMatrices(self, domsHidden):
        if domsHidden:
            toHide = self.matrices[0]
            toShow = self.firstDomMatrix
        else:
            toHide = self.firstDomMatrix
            toShow = self.matrices[0]
        toHide.hide()
        toShow.show()
