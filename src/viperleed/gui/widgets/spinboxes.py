"""Module spinboxes of viperleed.gui.widgets."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2020-01-12'
__license__ = 'GPLv3+'

import os

import numpy as np
from PyQt5 import QtCore as qtc
from PyQt5 import QtGui as qtg
from PyQt5 import QtWidgets as qtw

from viperleed.gui.widgets.decorators import receive_mouse_broadcast
from viperleed.gui.widgets.lib import AllGUIFonts


@receive_mouse_broadcast
class TextBox(qtw.QLineEdit):

    # This is the signal emitted when the user edits the text and
    # then changes focus, clicks away or Enter is pressed after the edit
    textModified = qtc.pyqtSignal(float, float)  # (before, after)

    # upDownPressed is a signal that is emitted whenever the up or down
    # keyboard keys are pressed
    upDownPressed = qtc.pyqtSignal(qtc.QEvent)

    def __init__(self, *args, **kwargs):
        # args = contents, parent in any order
        contens = ''
        parent = None
        try:
            if isinstance(args[0], str):
                contents = args[0]
        except:
            pass
        try:
            if isinstance(args[0], qtw.QWidget):
                parent = args[0]
            elif isinstance(args[1], qtw.QWidget):
                parent = args[1]
        except:
            pass

        super().__init__(contents, parent)
        self.returnPressed.connect(self.checkText)
        self._before = contents
        self.precision = 1
        if 'precision' in kwargs.keys():
            self.precision = kwargs['precision']

    def setStep(self, baseStep, stepType):
        if stepType not in ['add', 'scale']:
            raise

        self.baseStep = baseStep
        self.stepType = stepType
        if self.stepType == 'add':
            self.smallStep = self.baseStep/4
        elif self.stepType == 'scale':
            self.smallStep = 1 + (self.baseStep-1)/4

    def setLimits(self, minimum, maximum, limTyp='coerce'):
        self.lims = (minimum, maximum)
        self.limTyp = limTyp

    def textToFloat(self, s):
        try:
            ret = float(s)
        except ValueError:
            return np.nan
        else:
            return ret

    def getFloatText(self):
        return self.textToFloat(self.text())

    def roundFloat(self, value):
        return np.round(value, decimals=self.precision)

    def processValue(self, value):
        if not isinstance(value, float):
            raise

        value = self.roundFloat(value)
        if hasattr(self, 'lims'):
            mi = min(self.lims)
            ma = max(self.lims)
            if value < mi or value > ma:
                if self.limTyp == 'cyclic':
                    # this processes the value with periodic boundaries.
                    # NB: since python modulo returns always the same sign as
                    # the denominator the following works
                    value = self.roundFloat((value-mi) % (ma-mi) + mi)
                elif self.limTyp == 'coerce':  # this coerces the value
                    # when v < mi sorted() returns (v, mi, ma)
                    # when v > ma sorted() returns (mi, ma, v)
                    # when mi < v < ma sorted() returns (mi, v, ma)
                    value = sorted((mi, ma, value))[1]
        return value

    def updateText(self, new):
        self.setText(str(new))
        self.checkText()

    def checkText(self):  # also rounds value to self.precision
        if self._before != self.text():
            new = self.getFloatText()
            if not np.isnan(new):  # input is float
                new = self.processValue(new)
                self.setText(str(new))
                self.textModified.emit(self.textToFloat(self._before), new)
                self._before = self.text()
            else:
                self.setText(self._before)

    # The next three are used to check whether the textbox is edited by the
    # user or if the user simply clicked away without modifying the text. In
    # the latter case no signal is emitted.
    def focusInEvent(self, event):
        if event.reason() != qtc.Qt.PopupFocusReason:
            self._before = self.text()
        super().focusInEvent(event)

    def focusOutEvent(self, event):
        if event.reason() != qtc.Qt.PopupFocusReason:
            self.checkText()
        super().focusOutEvent(event)

    def mousePressEvent(self, event):
        # always react, on mousePress as we want to acknowledge the text
        # edit even if someone clicks away, without changing the focus
        self.checkText()
        if self.underMouse():
            super().mousePressEvent(event)  # call normal implementation

    def keyPressEvent(self, event):
        if event.key() in [qtc.Qt.Key_Up, qtc.Qt.Key_Down]:
            self.upDownPressed.emit(event)
        super().keyPressEvent(event)

    def event(self, event):
        """Extend event implementation to catch Ctrl+E key press in Linux."""
        if (isinstance(event, qtg.QKeyEvent)
            and os.name == 'posix'
            and event.key() == qtc.Qt.Key_E
                and event.modifiers() == qtc.Qt.ControlModifier):
            # Returning False causes the event to be propagated to
            # the parent, de facto disabling the 'go to end of line'
            # behavior it would normally have on Linux.
            return False
        return super().event(event)

    def wheelEvent(self, event):
        if self.hasFocus():
            # event.angleDelta() is the amount the 'wheel' is moved along
            # the two directions (works also for trackpad with two fingers).
            # For true mouse wheel, angleDelta() returns the amount in
            # multiples of 1/8 of a degree; for trackpad, the number
            # returned corresponds to the number of pixels of scrolling
            movement = event.angleDelta().y()/8
            if event.source() == qtc.Qt.MouseEventSynthesizedBySystem:
                # track-pad
                movement /= 1  # fractional movement
            else:  # true mouse wheel
                movement /= 20

            step = self.baseStep
            # make steps smaller if Ctrl is held down during wheelEvent
            # --> one could use bitwise & to do the following in case any
            # combination of modifiers is pressed that includes the Ctrl button
            if event.modifiers() == qtc.Qt.ControlModifier:
                step = self.smallStep

            new = self.getFloatText()
            if self.stepType == 'add':
                new += step*movement
            elif self.stepType == 'scale':
                new *= 1 + (step-1)*movement
            new = self.processValue(new)
            self.updateText(new)
        super().wheelEvent(event)


class TextBoxWithButtons(qtw.QWidget):
    def __init__(self, **kwargs):
        parent = kwargs.get('parent', None)

        super().__init__(parent=parent)

        self.label = qtw.QLabel(kwargs.get('labelText', ''), parent)
        self.text = TextBox(kwargs.get('textBoxText', ''), parent)
        self.topBut = qtw.QPushButton(kwargs.get('topButText', ''), parent)
        self.botBut = qtw.QPushButton(kwargs.get('botButText', ''), parent)
        self.bottomWidget = qtw.QWidget(parent)
        # self.bottomWidget is just a container at the bottom that can
        # be filled with different stuff in subclasses

        self.subWidgs = [self.label, self.text, self.topBut,
                         self.botBut, self.bottomWidget]
        self.setFonts()

        self.textBoxWidth = kwargs.get('textBoxWidth',
                                       self.text.sizeHint().width())

        self.smallButtonDims = int(self.textBoxWidth/5)                         # TODO: was .label.width()/6
        # print(f"{self.textBoxWidth=}, {self.smallButtonDims=}")

        # setBuddy allows to use shortcuts if there's a single & prepended
        # to a letter in the text of the QLabel
        self.label.setBuddy(self.text)

        self.compose()

        # connect the upDownPressed signal of TextBox: arrows up/down will
        # activate buttons (without animation)
        self.text.upDownPressed.connect(self.on_upDown)

        # connect the textModified signal of TextBox: update state of buttons
        # if limits are present
        self.text.textModified.connect(self.on_textModified)

    def setFonts(self):
        self.label.setFont(AllGUIFonts().labelFont)
        self.text.setFont(AllGUIFonts().labelFont)
        for but in [self.topBut,self.botBut]:
            but.setFont(AllGUIFonts().buttonFont)
        self.ensurePolished()

    def setTips(self, **kwargs):
        if 'labelTip' in kwargs.keys():
            self.label.setToolTip(kwargs['labelTip'])
            self.label.setStatusTip(kwargs['labelTip'])
        if 'textBoxTip' in kwargs.keys():
            self.text.setToolTip(kwargs['textBoxTip'])
            self.text.setStatusTip(kwargs['textBoxTip'])
        if 'topButTip' in kwargs.keys():
            self.topBut.setToolTip(kwargs['topButTip'])
            self.topBut.setStatusTip(kwargs['topButTip'])
        if 'botButTip' in kwargs.keys():
            self.botBut.setToolTip(kwargs['botButTip'])
            self.botBut.setStatusTip(kwargs['botButTip'])

    def compose(self):
        #----- Set the sizes and policies before setting up the layout ----#
        # Allow the text box to expand vertically
        self.text.setSizePolicy(qtw.QSizePolicy.Fixed,
                                qtw.QSizePolicy.Expanding)
        self.bottomWidget.setSizePolicy(qtw.QSizePolicy.Fixed,
                                        qtw.QSizePolicy.Expanding)
        for but in [self.topBut, self.botBut]:
            but.setSizePolicy(qtw.QSizePolicy.Fixed, qtw.QSizePolicy.Fixed)

        self.text.setMinimumHeight(2*self.smallButtonDims)
        self.text.setMaximumWidth(self.textBoxWidth)
        for but in [self.topBut, self.botBut]:
            but.setFixedSize(self.smallButtonDims, self.smallButtonDims)
        [widg.adjustSize() for widg in self.subWidgs]

        lay = qtw.QGridLayout()
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(0)

        '''
               0         1        2
          |---------|--------|---------|
        0 | label   | label  | label   |
          |---------|--------|---------|
        1 |         |        |         |
          |---------|--------|---------|
        2 | text    |        | topBut  |
          |---------|--------|---------|
        3 | text    |        | botBut  |
          |---------|--------|---------|
        4 |         |        |         |
          |---------|--------|---------|
        5 | bottom  | bottom | bottom  |
          |---------|--------|---------|

        #rows 1 and 4 are a spacer
        #col 1 is a spacer
        '''
        lay.addWidget(self.label, 0, 0, 1, 3)
        lay.addWidget(self.text, 2, 0, 2, 1)
        lay.addWidget(self.topBut, 2, 2, 1, 1)
        lay.addWidget(self.botBut, 3, 2, 1, 1)
        lay.addWidget(self.bottomWidget, 5, 0, 1, 3)

        #----- set the columnWidths and rowHeights -----#
        lay.setColumnMinimumWidth(0, self.textBoxWidth)
        lay.setColumnMinimumWidth(1, 1)
        lay.setColumnMinimumWidth(2, self.smallButtonDims)

        lay.setRowMinimumHeight(0, self.label.height())
        lay.setRowMinimumHeight(1, 3)
        lay.setRowMinimumHeight(2, self.smallButtonDims)
        lay.setRowMinimumHeight(3, self.smallButtonDims)
        lay.setRowMinimumHeight(4, 3)

        lay.setSizeConstraint(qtw.QLayout.SetFixedSize)

        self.setLayout(lay)
        self.adjustSize()

        self.bottomWidget.setMaximumWidth(self.width())

    def makeBottomWidget(self):
        # dummy function to be reimplemented in child classes. The
        # reimplementation should set the row height!
        pass

    def on_upDown(self, event):
        if event.key() == qtc.Qt.Key_Up:
            self.topBut.click()
        else:
            self.botBut.click()

    def on_textModified(self, old, new):
        # handle enabling state of increment/decrement buttons
        # This should have an effect only if there are limits for the values
        # and if the textbox limits are to be coerced
        if not hasattr(self.text, 'lims'):
            return None
        if self.text.limTyp == 'coerce':
            if new == min(self.text.lims):
                enable = [True, False]
            elif new == max(self.text.lims):
                enable = [False, True]
            else:
                enable = [True, True]
            self.topBut.setEnabled(enable[0])
            self.botBut.setEnabled(enable[1])
        return None

    def on_buttonPressed(self, pressed):
        if pressed in [self.topBut, self.botBut]:
            step = self.text.baseStep
            if qtw.qApp.keyboardModifiers() == qtc.Qt.ControlModifier:
                step = self.text.smallStep
            if self.text.stepType == 'scale':
                step -= 1
            if pressed == self.botBut:
                step *= -1

            new = self.text.getFloatText()
            if self.text.stepType == 'add':
                new += step
            else:
                new *= 1 + step
            new = self.text.processValue(new)
            self.text.updateText(new)
        elif any(isinstance(ch, qtw.QPushButton)
                 for ch in self.bottomWidget.children()):
            # There are buttons in the bottomWidget -> defer action to the
            # right handler
            self.on_botWidgPressed(pressed)

    def on_botWidgPressed(self, pressed):
        # dummy function to be reimplemented in subclasses that have
        # buttons in the bottomWidget
        pass
