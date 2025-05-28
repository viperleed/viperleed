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


class CoercingDoubleSpinBox(TolerantCommaSpinBox):
    """Coercing QDoubleSpinBox that sets limits after edit is done."""

    def __init__(self, decimals=None, soft_range=tuple(), step=1, suffix='',
                 **kwargs):
        """Initialise widget.

        Parameters
        ----------
        decimals : int, optional
            The number of decimal places. If not given or None,
            use the default decimal places. Default is None.
        soft_range : tuple, optional
            The soft minimum and maximum. While it is not impossible
            to set values below and above these respectively, the input
            will be set to either the minimum and maximum value after
            editing is finished.
        step : int or float, optional
            The increment of the SpinBox value. Default is 1.
        suffix : str, optional
            The suffix of the SpinBox. Default is '', no suffix.

        Returns
        -------
        None.
        """
        super().__init__(**kwargs)
        self.setRange(-float('inf'), float('inf'))
        self._soft_min = -float('inf')
        self._soft_max = float('inf')
        if decimals is not None:
            self.setDecimals(decimals)
        self.setSingleStep(step)
        if soft_range:
            self.soft_range = soft_range
        if suffix:
            self.setSuffix(suffix)
        self.editingFinished.connect(self._coerce_value)

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
            raise ValueError('The maximum cannot be smaller than the minimum.')
        self._soft_max = new_maximum

    @property
    def soft_range(self):
        """Return soft limits."""
        return self.soft_minimum, self.soft_maximum

    @soft_range.setter
    def soft_range(self, values):
        """Set soft limits."""
        new_minimum, new_maximum = values
        if new_minimum > new_maximum:
            new_maximum, new_minimum = new_minimum, new_maximum
        self._soft_max = new_maximum
        self._soft_min = new_minimum

    @qtc.pyqtSlot(int)
    def stepBy(self, steps):    # pylint: disable=invalid-name
        """Adjust set value through steps according to soft limits."""
        _, value, _ = sorted((*self.soft_range,
                              self.value() + steps*self.singleStep()))
        self.setValue(value)

    @qtc.pyqtSlot()
    def _coerce_value(self):
        """Check if value is whithin the limits and adjust it if necessary."""
        _, value, _ = sorted((self.value(), *self.soft_range))
        self.setValue(value)


class CoercingSpinBox(CoercingDoubleSpinBox):
    """Coercing QSpinBox that sets limits after edit is done."""

    def __init__(self, step=1, soft_range=tuple(), suffix='', **kwargs):
        """Initialise widget."""
        super().__init__(soft_range=soft_range, step=step,
                         suffix=suffix, **kwargs)
        super().setDecimals(0)

    def setDecimals(self, _):   # pylint: disable=invalid-name
        """Disable setting decimals."""
        raise AttributeError('CoercingSpinBox cannot setDecimals.')
