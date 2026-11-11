"""Module checkboxes of viperleed.gui.widgets.

Defines custom check boxes and radio buttons.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    'Florian Dörr (@FlorianDoerr)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2024-06-13'
__license__ = 'GPLv3+'

from PyQt5 import QtCore as qtc
from PyQt5 import QtWidgets as qtw


class QCheckBoxInvertedSignal(qtw.QCheckBox):
    """A QCheckBox with an extra unchecked signal.

    Signals
    -------
    unchecked : bool
        Emitted right after the stateChanged one. It carries
        a True value if the check box was just unchecked,
        False if it was just checked.
    """

    unchecked = qtc.pyqtSignal(bool)

    def __init__(self, **kwargs):
        """Initialise widget."""
        super().__init__(**kwargs)
        self.stateChanged.connect(self._emit_inverted_signal)

    @qtc.pyqtSlot(int)
    def _emit_inverted_signal(self, value):
        """Emit unchecked signal."""
        self.unchecked.emit(not bool(value))
