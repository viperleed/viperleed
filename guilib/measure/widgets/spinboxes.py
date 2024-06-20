"""Module spinboxes of viperleed.guilib.measure.widgets.

========================================
   ViPErLEED Graphical User Interface
========================================

Created: 2021-12-06
Author: Michele Riva

Defines convenience subclasses of QSpinBox and QDoubleSpinBox.
"""

from PyQt5 import (QtCore as qtc,
                   QtWidgets as qtw,
                   QtGui as qtg)


class TolerantCommaSpinBox(qtw.QDoubleSpinBox):
    """A QDoubleSpinBox that tolerates both '.' and ',' as separators.

    Each time the user presses the "comma" key, it is interpreted
    as a if a "dot" was pressed instead. This is useful to keep
    intact the user experience of pressing the "comma" key as a
    decimal separator while enforcing "dot" as a decimal separator
    everywhere.

    Reimplemented methods
    ---------------------
    keyPressEvent(event)
        Replace comma key-presses with dot key presses.
    keyReleaseEvent(event)
        Replace comma key-releases with dot key releases.
    """

    def keyPressEvent(self, event):      # pylint: disable=invalid-name
        """Replace commas with dots."""
        super().keyPressEvent(self.__make_dot_key_event(event))

    def keyReleaseEvent(self, event):    # pylint: disable=invalid-name
        """Replace commas with dots."""
        super().keyReleaseEvent(self.__make_dot_key_event(event))

    @staticmethod
    def __make_dot_key_event(event):
        """Return a KeyEvent with comma replaced by dot."""
        if event.key() != qtc.Qt.Key_Comma:
            return event
        return qtg.QKeyEvent(
            event.type(), qtc.Qt.Key_Period, event.modifiers(),
            text=event.text().replace(',', '.'),
            autorep=event.isAutoRepeat(), count=event.count()
            )


class InfIntSpinBox(qtw.QDoubleSpinBox):
    """A spin-box that allows inputting an infinite integer number."""

    def __init__(self, parent=None):
        """Initialize instance."""
        super().__init__(parent)
        self.setRange(0, float('inf'))
        self.setDecimals(0)


class CoercingDoubleSpinBox(TolerantCommaSpinBox):
    """Coercing QDoubleSpinBox that sets limits after edit is done."""

    def __init__(self, decimals=1, range_=tuple(), step=1, suffix='', **kwargs):
        """Initialise widget.

        Parameters
        ----------
        decimals : int, optional
            The amount of decimals. Default is 1.
        range_ : tuple, optional
            The soft minimum and maximum.
        step : int or float, optional
            The increment of the SpinBox value. Default is 1.
        suffix : str, optional
            The suffix of the SpinBox. Default is '', no suffix.

        Returns
        -------
        None
        """
        super().__init__(**kwargs)
        self.setRange(-float('inf'), float('inf'))
        self._soft_min = -float('inf')
        self._soft_max = float('inf')
        self.setDecimals(decimals)
        self.setSingleStep(step)
        if range_:
            self.range_ = range_
        if suffix:
            self.setSuffix(suffix)
        self.editingFinished.connect(self._adjust_value)

    @property
    def soft_minimum(self):
        """Return soft minimum."""
        return self._soft_min

    @soft_minimum.setter
    def soft_minimum(self, new_minimum):
        """Set soft minimum."""
        if new_minimum > self.soft_maximum:
            raise ValueError('The minimum cannot be larger than the maximum.')
        self._soft_min = new_minimum

    @property
    def soft_maximum(self):
        """Return soft maximum."""
        return self._soft_max

    @soft_maximum.setter
    def soft_maximum(self, new_maximum):
        """Set soft maximum."""
        if new_maximum < self.soft_minimum:
            raise ValueError('The maximum cannot be lower than the minimum.')
        self._soft_max = new_maximum

    @property
    def range_(self):
        """Return soft limits."""
        return self.soft_minimum, self.soft_maximum

    @range_.setter
    def range_(self, values):
        """Set soft limits."""
        new_minimum, new_maximum = values
        if new_minimum > new_maximum:
            new_maximum, new_minimum = new_minimum, new_maximum
        self._soft_max = new_maximum
        self._soft_min = new_minimum

    @qtc.pyqtSlot()
    def _adjust_value(self):
        """Check if value is whithin the limits and adjust it if necessary."""
        value = self.value()
        if value > self.soft_maximum:
            self.setValue(self.soft_maximum)
        if value < self.soft_minimum:
            self.setValue(self.soft_minimum)

    @qtc.pyqtSlot(int)
    def stepBy(self, steps):
        """Adjust set vaöue through steps according to soft limits."""
        _, value, _ = sorted((self.soft_minimum, self.soft_maximum,
                             self.value() + steps*self.singleStep()))
        self.setValue(value)


class CoercingSpinBox(CoercingDoubleSpinBox):
    """Coercing QSpinBox that sets limits after edit is done."""

    def __init__(self, step=1, range_=tuple(), suffix='', **kwargs):
        """Initialise widget."""
        super().__init__(decimals=0, range_=range_, step=step,
                         suffix=suffix, **kwargs)

