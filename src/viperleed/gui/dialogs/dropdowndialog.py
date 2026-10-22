"""Module dropdowndialog of viperleed.gui.measure.dialogs.

Defines the DropdownDialog class, a modal dialog that can be used
when execution requires a user choice from a list of options.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2021-12-08'
__license__ = 'GPLv3+'

from collections.abc import Sequence

from PyQt5 import QtWidgets as qtw


_ICON_ALIAS = {'q': 'Question',
               'i': 'Information',
               'w': 'Warning',
               'c': 'Critical'}


class DropdownDialog(qtw.QMessageBox):
    """A modal dialog that presents options to choose from.

    Typical use of the dialog:
    if dialog.exec_() == dialog.Apply:
        selection = dialog.selection
        <do stuff with selection>
    else:
        # Closed, or "Abort" was pressed
        <abort operations>
    """

    def __init__(self, title, text, options, **kwargs):
        """Initialize dialog.

        Parameters
        ----------
        title : str
            The title to be given to the dialog.
        text : str
            Informative text that will be shown to the user.
        options : list of str
            The options the user has to choose from.
        severity : {'question', 'information',
                    'warning', 'critical'}, optional
            Defines which icon (and system sound) will be used
            when the dialog executes. Only the first character is
            checked, i.e., 'q' is equivalent to 'question'.
            Default is 'q'.
        parent : QWidget, optional
            The parent window of this dialog. The parent and all
            its ancestors will be blocked when .exec_() is called
            until the user successfully terminates the dialog.

        Returns
        -------
        None.

        Raises
        ------
        TypeError
            If severity is not a string.
        TypeError
            If options is not a sequence.
        ValueError
            If severity is an unknown severity level.
        """
        severity = kwargs.get('severity', 'q')
        if not isinstance(severity, str):
            raise TypeError(
                f"{self.__class__.__name__}: invalid type "
                f"{type(severity).__name__!r} for severity argument. "
                "Expected 'str'"
                )
        if not isinstance(options, Sequence):
            raise TypeError(
                f"{self.__class__.__name__}: invalid type "
                f"{type(options).__name__!r} for severity argument. "
                "Expected a sequence."
                )
        if severity[0].lower() not in _ICON_ALIAS:
            raise ValueError(
                f"{self.__class__.__name__}: Unkown severity {severity!r}"
                )

        icon = getattr(self, _ICON_ALIAS[severity[0]])
        options = [str(o) for o in options]
        super().__init__(icon, title, text, buttons=self.Apply|self.Cancel,
                         parent=kwargs.get('parent', None))

        self.options = qtw.QComboBox()
        self.options.addItems(options)

        layout = self.layout()
        # The layout (grid) is as follows:
        # col 0: icon
        # col 1: spacer
        # col 2: text
        # The last row is filled with with dialog buttons (flushed right).
        # I will insert the dropdown at the pre-last row, last column,
        # but have to move down the last row by one index
        for idx in reversed(range(layout.count())):
            row, col, row_span, col_span = layout.getItemPosition(idx)
            if row < 2:
                # Anything except the buttons
                break
            row = row + 1
            item = layout.takeAt(idx)
            layout.addItem(item, row, col, row_span, col_span)
        layout.addWidget(self.options, layout.rowCount() - 2, 2, 1, -1)

    @property
    def selection(self):
        """Return the selected option."""
        return self.options.currentText()
