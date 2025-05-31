"""Module domainsblock of viperleed.gui.leedsim.widgets.

Defines the DomsBlock class for displaying information about
symmetry-induced domains.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2020-01-12'
__license__ = 'GPLv3+'

import PyQt5.QtCore as qtc
import PyQt5.QtWidgets as qtw

from viperleed.gui.leedsim.widgets.matricespopup import MatricesPopup
from viperleed.gui.leedsim.widgets.togglebutton import ToggleButton
from viperleed.gui.widgets.lib import AllGUIFonts


class DomsBlock(qtw.QWidget):
    toggleText = ['Hide &Domains', 'Show &Domains']
    # These are the two possible strings that appear on the toggle button.
    # Pressing alt+d fires the toggle button

    def __init__(self, parent=None):
        super().__init__(parent)
        self.hide = False

        # The next line is useful only if I ever want to place this in a
        # layout
        self.setSizePolicy(qtw.QSizePolicy.Minimum, qtw.QSizePolicy.Minimum)
        self.text = qtw.QLabel('\u2014 inequivalent domain(s)')
        self.toggle = ToggleButton(self.toggleText[self.hide])

        self.subWidgs = (self.text, self.toggle)

        self.setTips(text='Open or drag-drop a LEED pattern file to'
                          'determine how many domains!',
                     but='Toggle visibility of symmetry-equivalent domains')
        self.compose()

    def setTips(self, **kwargs):
        if 'text' in kwargs.keys():
            self.text.setStatusTip(kwargs['text'])
        if 'but' in kwargs.keys():
            tip = kwargs['but']
            self.toggle.setToolTip(tip)
            self.toggle.setStatusTip(tip)

    def compose(self):
        #set fonts
        self.text.setFont(AllGUIFonts().labelFont)
        self.toggle.setFont(AllGUIFonts().buttonFont)

        #set sizes
        self.toggle.setSizePolicy(qtw.QSizePolicy.Minimum,
                                  qtw.QSizePolicy.Minimum)
        self.toggle.adjustSize()

        self.text.setAlignment(qtc.Qt.AlignRight | qtc.Qt.AlignBaseline)
        self.text.setSizePolicy(qtw.QSizePolicy.Expanding,
                                qtw.QSizePolicy.Minimum)
        self.text.adjustSize()
        self.text.setMinimumSize(self.text.size())
        self.toggle.setMinimumSize(self.toggle.size())

        #----- Set up the layout -----#
        lay = qtw.QVBoxLayout()
        lay.setSpacing(2)
        lay.setContentsMargins(0, 0, 0, 0)

        lay.addStretch(1)  # the stretches keep things together
        lay.addWidget(self.text)
        lay.addWidget(self.toggle)
        lay.addStretch(1)
        lay.setAlignment(self.text, self.text.alignment())
        lay.setAlignment(self.toggle, qtc.Qt.AlignRight | qtc.Qt.AlignVCenter)

        self.setLayout(lay)

        self.adjustSize()

    def width(self):
        self.ensurePolished()
        return max(sub.width() for sub in self.subWidgs)

    def height(self):
        self.ensurePolished()
        return (sum(sub.height() for sub in self.subWidgs)
                + self.layout().spacing())

    def minimumSizeHint(self):
        return self.sizeHint()

    def sizeHint(self):
        return qtc.QSize(self.width(), self.height())

    def togglePressed(self):
        self.hide = not self.hide
        self.toggle.setText(self.toggleText[self.hide])
        self.matricesPopup.updateMatrices(self.hide)

    def updateText(self, txt):
        if not isinstance(txt, str):
            raise

        self.text.setText(txt)
        self.text.setMinimumSize(self.text.sizeHint())

        trOld = self.geometry().topRight()
        self.adjustSize()
        newPos = (self.geometry().topLeft()
                  - self.geometry().topRight()
                  + trOld)
        self.move(newPos)

    def initPopup(self):
        leed = self.window().leed
        self.matricesPopup = MatricesPopup(leed.superlattices,
                                           leed.domColors,
                                           parent=self.window())

    #----------------------- Re-implement mouseEvents -----------------------#
    def mouseReleaseEvent(self, event):
        if self.underMouse():
            child = self.childAt(event.pos())
            if child == self.text:
                # Position and show matricesPopup
                if not self.matricesPopup.shown:
                    # not shown before --> place it at the standard position
                    tR = self.mapToGlobal(self.text.geometry().topRight())
                    sizePopup = self.matricesPopup.frameSize()
                    newPos = qtc.QPoint(round(tR.x() - 0.7*sizePopup.width()),
                                        tR.y() - sizePopup.height() - 25)
                    self.matricesPopup.move(newPos)
                    self.matricesPopup.dragPosition = newPos
                    self.matricesPopup.shown = True
                if self.matricesPopup.isVisible():
                    self.matricesPopup.hide()
                else:
                    self.matricesPopup.show()
        super().mouseReleaseEvent(event)
