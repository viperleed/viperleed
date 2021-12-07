"""
"""

import traceback

from PyQt5 import QtWidgets as qtw


class ErrorBox(qtw.QMessageBox):
    """A pop-up dialog reporting errors.

    This class should be used to catch Exceptions and report
    them to the user in a graphical manner. It would be typically
    used when overriding sys.excepthook.
    """

    def __init__(self, error_while="", text="", silent=False, parent=None):
        """Initialize dialog.

        Parameters
        ----------
        error_while : str
            Portion of title for the dialog. The full title
            will be "Error while <error_while>".
        text : str, optional
            Descriptive text that will be used every time
            the dialog is shown.
        silent : bool, optional
            If True, no dialog is shown when an error occurs.
            Can be accessed via the .silent attribute. Default
            is False.
        parent : QWidget, optional
            The parent widget. Default is None.

        Returns
        -------
        None

        Raises
        ------
        TypeError
            If the mandatory argument error_while is missing
        """
        if not error_while:
            raise TypeError(f"{self.__class__.__name__} missing "
                            "mandatory argument error_while")
        title = f"Error while {error_while}"

        self.base_text = text
        self.silent = silent

        super().__init__(qtw.QMessageBox.Critical, title, text,
                         buttons=qtw.QMessageBox.Ok, parent=parent)

    def exec_(self, extra_text=""):
        """Open the dialog as modal if there is an exception.

        Parameters
        ----------
        extra_text : str
            Text to be appended to the one given at instantiation.
        """
        if self.silent:
            return self.Ok

        trace = traceback.format_exc()
        if trace.startswith('NoneType'):
            # Do not execute if there is no exception
            return self.Ok

        self.setText(self.base_text + extra_text)
        self.setDetailedText(f"Error details:\n{trace}")
        return super().exec_()
