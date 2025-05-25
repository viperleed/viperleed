"""
======================================
  ViPErLEED Graphical User Interface
======================================
  *** module guilib.basewidgets ***


Created: 2020-01-12
Author: Michele Riva

Contains QWidget subclasses that re-implement/enhance/add some features.
Typically they would be further subclassed to create more complicated widgets
"""

# print("this is guilib.basewidgets")

import os

import numpy as np
import PyQt5.QtCore as qtc
import PyQt5.QtGui as qtg
import PyQt5.QtWidgets as qtw

from viperleed import guilib as gl

if gl.BACKEND == 'mplcairo':
    try:
        import mplcairo  # requires import before matplotlib
    except ImportError:
        gl.BACKEND = 'agg'  # fallback on agg if mplcairo is not there
else:
    gl.BACKEND = 'agg'

import matplotlib as mpl  # this should be the first import ever of matplotlib

if gl.BACKEND == 'mplcairo':
    mpl.use("module://mplcairo.qt")
    from mplcairo.qt import FigureCanvasQTCairo as FigureCanvas
else:  # 'agg'
    mpl.use('qt5agg')
    from matplotlib.backends.backend_qt5agg import FigureCanvas

from matplotlib.figure import Figure

from viperleed.guilib.widgetdecorators import receive_mouse_broadcast
from viperleed.guilib.widgetslib import AllGUIFonts
from viperleed.guilib.widgetslib import drawText
from viperleed.guilib.widgetslib import editStyleSheet


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


class MPLFigureCanvas(FigureCanvas):
    """
    Subclass of the matplotlib (MPL) FigureCanvas class
    """

    def __init__(self, figure=None, **kwargs):
        if figure is None:
            figure = Figure()
        super().__init__(figure)

        self.title = None
        self.__ax = self.figure.subplots()  # matplotlib axes
        self.__wheel_buddy = None  # buddy widget for wheel events

        self.set_transparent_background(True)
        self.show_ticks(False)
        self.fill_canvas()
        self.windowFraction = 0.8
        self.mplLayout = qtw.QVBoxLayout()  # This is the layout containing the
                                            # figure and its title. It is not
                                            # set to any widget, but should be
                                            # added to the layout of the
                                            # centralWidget
        # Change the size policy such that it hasHeightForWidth. This keeps
        # MPLFigureCanvas square when inserted in a layout
        size_policy = self.sizePolicy()
        size_policy.setHeightForWidth(True)
        self.setSizePolicy(size_policy)

        self.setParent(kwargs.get('parent', None))

        if 'titleFont' in kwargs.keys() or 'fontSize' in kwargs.keys():
            self.titleFont = qtg.QFont()
            if 'titleFont' in kwargs.keys():
                self.titleFont.setFamily(kwargs['titleFont'])
            if 'titleFontSize' in kwargs.keys():
                self.titleFont.setPointSize(kwargs['titleFontSize'])
        else:
            self.titleFont = AllGUIFonts().plotTitleFont

        self.setTitle(kwargs.get('title', ''))

        if 'plotCanvasSize' in kwargs.keys():
            self.setSize(kwargs['plotCanvasSize'])
        else:
            self.adjustSize()

    @property
    def ax(self):
        """
        Returns a reference to the matplotlib axes on which stuff is plotted
        """
        return self.__ax

    @property
    def wheel_buddy(self):
        """
        This property holds a reference to a QWidget whose wheelEvent is called
        whenever a wheelEvent occurs on the canvas
        """
        return self.__wheel_buddy

    @wheel_buddy.setter
    def wheel_buddy(self, buddy):
        if not isinstance(buddy, qtw.QWidget):
            raise TypeError("Argument of wheel_buddy must be a QWidget. "
                            f"Found {type(buddy)} instead")
        self.__wheel_buddy = buddy

    # def sizeHint(self):
        # s = round(
            # min(self.window().centralWidget().height()*self.windowFraction,
            # self.window().centralWidget().width()*0.45)
            # )
        # sH = qtc.QSize()
        # sH.setWidth(s)
        # sH.setHeight(s)
        # return sH

    def minimumSizeHint(self):
        min_width = 350
        return qtc.QSize(min_width, self.heightForWidth(min_width))

    # def minimumSize(self):
        # return self.minimumSizeHint()

    def rightSize(self):
        sWindow = self.sizeHint()
        if sWindow.width() > self.minimumSizeHint().width():
            s = sWindow
        else:
            s = self.minimumSizeHint()
        l = min(s.width(), s.height())
        return qtc.QSize(l, l)

    def heightForWidth(self, width):
        """
        Re-implementation of QWidget heightForWidth that allows to keep the
        widget square when inserted in a QLayout
        """
        # print("hfw", end=' ', flush=True)
        return width

    # TODO: this causes the fixed sizes issue. Perhaps would be better to
    # actually set mplLayout as self.layout (may be impossible), or have
    # a widget around the whole thing? Maybe LEEDPattern and RealSpace could
    # just be QWidget and contain MPLFigureCanvas. Removing this resizeEvent
    # makes heightForWidth kinda work, but only when reducing the width of
    # the container widget. Making it larger or changing height makes the
    # plots become rectangular.
    # def resizeEvent(self, event):
        # event = qtg.QResizeEvent(self.rightSize(), event.oldSize())
        # super().resizeEvent(event)  # call the original one

        # # the next updates notify the layout of the change in size,
        # # and keep the alignment correct
        # self.updateGeometry()
        # self.title.updateGeometry()
    # def resizeEvent(self, event):
        # new_len = min(event.size().width(), event.size().height())
        # new_adjusted_size = qtc.QSize(new_len, self.heightForWidth(new_len))
        # resize_event = qtg.QResizeEvent(new_adjusted_size, event.oldSize())
        # super().resizeEvent(resize_event)

        # # the next updates notify the layout of the change in size,
        # # and keep the alignment correct
        # self.updateGeometry()
        # self.title.updateGeometry()


    def getStretchFactor(self):
        return qtc.QSize(
            int(self.rightSize().width()
                / self.window().centralWidget().size().width()),
            int(self.rightSize().height()
                / self.window().centralWidget().size().height())
            )

    def set_transparent_background(self, transparent):
        if transparent:
            editStyleSheet(self, "background-color:transparent;")
            self.figure.patch.set_alpha(0)
        else:
            editStyleSheet(self, "background-color:white;")
            self.figure.patch.set_alpha(1)

    def show_ticks(self, show):
        ax = self.ax
        [x.set_visible(show) for x in (ax.get_xaxis(), ax.get_yaxis())]

    def fill_canvas(self):
        eps = 5e-3
        self.ax.set_position((eps, eps, 1-2*eps ,1-2*eps))

    def setSpinesOn(self, on):
        if on == True:
            self.ax.axis('on')
        elif on == False:
            self.ax.axis('off')
        else:
            raise ValueError("setSpines requires a bool input")

    def setAxLimits(self, limit):
        axlims = [-limit*1.01, limit*1.01]
        self.ax.set_xlim(axlims)
        self.ax.set_ylim(axlims)

    def setTitle(self, title):
        initRun = False
        if self.title is None:
            self.title = qtw.QLabel('', self.parentWidget())
            # I might want to set a SizePolicy here later
            initRun = True
        self.title.setText(title)
        self.title.setFont(self.titleFont)
        if initRun:
            self.mplLayout.setSpacing(3)
            self.mplLayout.setContentsMargins(0, 0, 0, 0)

            # This stretch and the last one are used to keep the figure and
            # title together
            self.mplLayout.addStretch(1)
            self.mplLayout.addWidget(self.title)#,# 0,
                                     # qtc.Qt.AlignHCenter | qtc.Qt.AlignBaseline)
            # self.mplLayout.addWidget(self, 20, qtc.Qt.AlignHCenter)
            self.mplLayout.addWidget(self)#, qtc.Qt.AlignHCenter)
            self.mplLayout.addStretch(2)

            # NB: for unknown reasons setting the alignment directly at the time
            # of adding the widgets to the layout messes up the spacing between
            # items
            # Right now, the title and self do not have a stretch factor defined
            # This might be an issue for resizing
            self.mplLayout.setAlignment(
                self.title,
                qtc.Qt.AlignHCenter | qtc.Qt.AlignBaseline
                )
            self.mplLayout.setAlignment(self, qtc.Qt.AlignHCenter)

            # self.mplLayout.setStretchFactor(self.title, 1)
            # self.mplLayout.setStretchFactor(self, 20)

    def setSize(self, size=None):
        self.setFixedSize(size, size)

    def getSpinesThickness(self):
        return self.ax.spines['top'].get_linewidth()

    def setEnabled(self, enable):
        # reimplementation of QWidget.setEnabled that
        # enables/disables also the title if present
        super().setEnabled(enable)
        if self.title is not None:
            self.title.setEnabled(enable)

    def wheelEvent(self, event):
        """
        Re-implement the standard wheelEvent to give keyboard focus and send a
        wheelEvent to the wheel_buddy.
        """
        if self.wheel_buddy is not None:
            if self.underMouse():
                if not self.wheel_buddy.hasFocus():
                    # Note: if the buddy does not have focus, sending the
                    # wheelEvent might have no effect (this is the case, e.g.,
                    # if TextBox is the buddy)
                    self.wheel_buddy.setFocus()
                self.wheel_buddy.wheelEvent(event)
        super().wheelEvent(event)

    def mpl_to_qt_coordinates(self, *args):
        """
        This is the inverse of the mouseEventCoords function implemented in
        matplotlib.backend_qt5 that returns QWidget (x, y) coordinates given the
        corresponding matplotlib ones. Qt works in logical coordinates with
        origin at top-left corner and y increasing from top to bottom.
        matplotlib uses physical coordinates, has origin at bottom left, and
        y increases from bottom to top.

        Returns a QPoint
        """
        if len(args) != 2:
            raise ValueError('mpl_to_qt_coordinates requires exactly '
                             'two arguments')

        x_mpl, y_mpl = args

        dpi_ratio = self._dpi_ratio
        x_qt = x_mpl/dpi_ratio
        y_qt = (self.figure.bbox.height - y_mpl)/dpi_ratio
        return qtc.QPoint(np.round(x_qt), np.round(y_qt))


class MeasurementFigureCanvas(FigureCanvas):
    """
    Subclass of the matplotlib FigureCanvas class
    used to display measurement datay.
    """

    def __init__(self, figure=None, **kwargs):
        if figure is None:
            figure = Figure()
        super().__init__(figure)

        self.__ax = self.figure.subplots()  # matplotlib axes

        self.setParent(kwargs.get('parent', None))

        self.adjustSize()
    
    @property
    def ax(self):
        """
        Returns a reference to the matplotlib axes on which stuff is plotted
        """
        return self.__ax


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

class QDoubleValidatorNoDot(qtg.QDoubleValidator):
    def validate(self, text:str, cursor_pos:int):
        text = text.replace(',', '.')
        return super().validate(text, cursor_pos)

