"""Module busywindow of viperleed.gui.measure.dialogs.

Defines the BusyWindow class, a modal (i.e., blocking
user events) window that can be shown while an operation
(of unknown duration) is running.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2021-11-27'
__license__ = 'GPLv3+'

from PyQt5 import QtCore as qtc
from PyQt5 import QtGui as qtg
from PyQt5 import QtWidgets as qtw


class BusyWindow(qtw.QWidget):
    """A window to signal that an operation is running."""

    timed_out = qtc.pyqtSignal()

    def __init__(self, parent=None, text='', max_onscreen_time=-1):
        """Initialize window.

        Parameters
        ----------
        parent : QObect, optional
            The parent window of this modal window. If given, this
            window will block user input to the parent (and any
            other ancestors) while visible. Default is None.
        text : str, optional
            The text that will be visualized. Default is an
            empty string.
        max_onscreen_time : float, optional
            Maximum time (in seconds) that the modal window will
            stay visible on screen, if a positive value is given.
            If no value (or a negative one) is given, the window
            will stay forever visible. Since the window has no
            way to be closed by the user, when max_onscreen_time
            is not given the caller should explicitly close() this
            window! Default is -1.

        Returns
        -------
        None.
        """
        super().__init__(parent=parent)

        on_screen_timer = qtc.QTimer()
        on_screen_timer.setSingleShot(True)
        on_screen_timer.timeout.connect(self.close)
        on_screen_timer.timeout.connect(self.timed_out)

        self.__ctrls = {'text': qtw.QLabel(text),
                        'circle': _RotatingCircle()}
        self.__timers = {'on_screen': (on_screen_timer,
                                       round(max_onscreen_time * 1000))}

        self.setWindowFlags(qtc.Qt.Tool | qtc.Qt.CustomizeWindowHint)
        self.setWindowModality(qtc.Qt.WindowModal)

        self.__compose()

    def __compose(self):
        """Place children widgets."""
        layout = qtw.QVBoxLayout()
        layout.setSpacing(layout.spacing() + 15)
        layout.addWidget(self.__ctrls['text'])
        layout.addWidget(self.__ctrls['circle'],
                         qtc.Qt.AlignCenter, qtc.Qt.AlignCenter)
        layout.setSizeConstraint(layout.SetFixedSize)

        self.setLayout(layout)
        self.adjustSize()

    def showEvent(self, event):  # pylint: disable=invalid-name
        """Reimplement to have a max visible time."""
        timer, interval = self.__timers['on_screen']
        if interval >= 0:
            timer.start(interval)
        # Center onto parent widget, if any
        parent = self.parent()
        if parent:
            delta_size = (parent.geometry().size()
                          - self.geometry().size()) / 2
            delta_pos = qtc.QPoint(delta_size.width(), delta_size.height())
            self.move(parent.pos() + delta_pos)
        super().showEvent(event)


class _RotatingCircle(qtw.QWidget):
    """A blue circle that continuously rotates."""

    def __init__(self, size=60, parent=None):
        """Initialize widget."""
        super().__init__(parent=parent)
        self.__size = size
        stroke = int((0.1 * size + 1) // 2) * 2
        self.__stroke_width = max(stroke, 2)
        self.__angle = 0

        self.__rotation = qtc.QPropertyAnimation(self, b'angle', parent=parent)
        self.__rotation.setDuration(1200)
        self.__rotation.setStartValue(0)
        self.__rotation.setEndValue(360)
        self.__rotation.finished.connect(self.__rotation.start)

        self.__pen = qtg.QPen()
        self.__pen.setColor(qtg.QColor(0, 120, 204, 200))  # viperleed blue
        self.__pen.setWidth(self.__stroke_width)
        self.__pen.setCapStyle(qtc.Qt.RoundCap)

        self.setMinimumSize(self.sizeHint())
        self.resize(self.sizeHint())

    @qtc.pyqtProperty(float)
    def angle(self):
        """Return the starting angle of self (in degrees)."""
        return self.__angle

    @angle.setter
    def angle(self, new_angle):
        """Set the starting angle of self (in degrees)."""
        self.__angle = new_angle
        self.repaint()

    def hideEvent(self, event):  # pylint: disable=invalid-name
        """Stop rotating when hidden."""
        super().hideEvent(event)
        self.__rotation.stop()

    def minimumSizeHint(self):  # pylint: disable=invalid-name
        """Return a reasonable minimum size for self."""
        return self.sizeHint()

    def paintEvent(self, *_):  # pylint: disable=invalid-name
        """Draw self."""
        circle_size = self.__size - self.__stroke_width
        circle_margin = round(self.__stroke_width / 2)

        painter = qtg.QPainter()
        painter.begin(self)
        painter.setRenderHint(painter.Antialiasing)

        rect = qtc.QRect(0, 0, self.__size, self.__size)
        # 'Draw' a blank, border-less rectangle
        # just to keep the widget size constant
        painter.setPen(qtc.Qt.NoPen)
        painter.drawRect(rect)

        painter.setPen(self.__pen)
        # The * 16 is because angles are expressed in 16th of degree
        painter.drawArc(circle_margin, circle_margin,
                        circle_size, circle_size,
                        round(self.angle * 16), 240 * 16)
        painter.end()

    def showEvent(self, event):  # pylint: disable=invalid-name
        """Start rotating when shown."""
        self.__rotation.start()
        super().showEvent(event)

    def sizeHint(self):  # pylint: disable=invalid-name
        """Return a reasonable size for self."""
        return qtc.QSize(self.__size, self.__size)
