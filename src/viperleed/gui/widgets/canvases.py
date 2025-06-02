"""Module canvases of viperleed.gui.widgets.

Defines subclasses of matplotlib FigureCanvas for plotting.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2020-01-12'
__license__ = 'GPLv3+'

from matplotlib.figure import Figure
import numpy as np
import PyQt5.QtCore as qtc
import PyQt5.QtGui as qtg
import PyQt5.QtWidgets as qtw

from viperleed.gui.mpl_graphics import import_figure_canvas
from viperleed.gui.widgets.lib import AllGUIFonts
from viperleed.gui.widgets.lib import editStyleSheet

FigureCanvas = import_figure_canvas()


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

