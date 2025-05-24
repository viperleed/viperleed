"""Module editablematrix of viperleed.guilib.leedsim.widgets.

======================================
  ViPErLEED Graphical User Interface
======================================

Defines the EditableMatrix widget, an interactively
editable matrix of integers or floats.

Created: 2021-06-01
Author: Michele Riva
"""

import numpy as np
import PyQt5.QtCore as qtc
import PyQt5.QtGui as qtg
import PyQt5.QtWidgets as qtw

from viperleed.guilib.widgetslib import AllGUIFonts

from viperleed.guilib import decorators as dev_


class EditableMatrix(qtw.QWidget):
    """An interactively editable matrix of integers or floats."""

    matrix_edited = qtc.pyqtSignal(np.ndarray)

    def __init__(self, parent=None, shape=(2, 2), dtype=int, label=True):
        """Initialize widget.

        Parameters
        ----------
        parent : PyQt5.QtWidgets.QWidget, default=None
            Parent widget that 'contains' this instance.
        shape : tuple, default=(2, 2)
            Shape of the matrix. Should have at most two
            entries.
        label : bool, default=True
            Whether an "M = " label should appear
            on the left of the matrix
        dtype : {int, float}
            The data type that will be used for
            validating the interactive input

        Returns
        -------
        None.

        Raises
        ------
        ValueError
            If shape is not a sequence of length 1 or 2
        ValueError
            If dtype is not int or float
        """
        super().__init__(parent)

        if not hasattr(shape, '__len__') or not 0 < len(shape) <= 2:
            raise ValueError("EditableMatrix: invalid shape. Expected "
                             "tuple with one or two elements, found "
                             f"{shape}")
        if dtype not in (int, float):
            raise ValueError("EditableMatrix: invalid dtype. Expected "
                             f"'int' or 'float' found {dtype}")

        self._dtype = dtype

        # Prepare the controls for the matrix elements
        # NB: one cannot use np.full to create the
        # array because otherwise all the elements
        # refer to the same object in memory
        _ctrls = [qtw.QLineEdit() for _ in range(np.prod(shape))]
        self._ctrls = np.asarray(_ctrls).reshape(shape)

        self._compose(label)
        self._connect()

    @property
    def matrix(self):
        """Return an array of the values in the controls.

        Returns
        -------
        input_values : numpy.ndarray
        """
        matrix = np.array([self._dtype(m_ij.text())
                           for m_ij in self._ctrls.ravel()])
        return matrix.reshape(self._ctrls.shape)

    @matrix.setter
    def matrix(self, matrix):
        """Set control values from a given matrix.

        Parameters
        ----------
        matrix : array-like
            Matrix to be set. It should have the same
            shape and dtype given at instantiation

        Raises
        ------
        ValueError
            If the shape of matrix is inconsistent
            with the one given at instantiation
        TypeError
            If the data in the matrix have a different
            type than the one given at instantiation

        Emits
        -----
        matrix_edited
            After setting the controls
        """
        self.set_matrix(matrix)
        self.matrix_edited.emit(matrix)

    def set_matrix(self, matrix):
        """Set control values from a given matrix.

        This does the same as assigning to the .matrix
        attribute, but does not cause emission of the
        matrix_edited signal.

        Parameters
        ----------
        matrix : array-like
            Matrix to be set. It should have the same
            shape and dtype given at instantiation

        Raises
        ------
        ValueError
            If the shape of matrix is inconsistent
            with the one given at instantiation
        TypeError
            If the data in the matrix have a different
            type than the one given at instantiation
        """
        matrix = np.asarray(matrix)
        if matrix.shape != self._ctrls.shape:
            raise ValueError(f"EditableMatrix: shape {matrix.shape} is"
                             f"incompatible. Expected {self._ctrls.shape}")
        if self._dtype == int != matrix.dtype:
            raise TypeError("EditableMatrix: Invalid elements with dtype="
                            f"{matrix.dtype}. Expected dtype=int")

        for ctrl_ij, m_ij in zip(self._ctrls.ravel(), matrix.ravel()):
            ctrl_ij.setText(str(self._dtype(m_ij)))

    def set_text_color(self, color):
        """Set the color of the text in the controls.

        The setting is done via the internal QPalette, and
        may not necessarily be honored in all systems.

        Parameters
        ----------
        color : QColor
        """
        palette = self._ctrls.ravel()[0].palette()
        palette.setColor(palette.Text, color)
        for ctrl in self._ctrls.ravel():
            ctrl.setPalette(palette)

    def _compose(self, label):
        """Place children widgets.

        Parameters
        ----------
        shape : tuple, default=(2, 2)
            Shape of the matrix. Should have at most two
            entries.
        label: bool, default=True
            Whether an "M = " label should appear
            on the left of the matrix

        Returns
        -------
        None.
        """
        # Fonts
        label_font = AllGUIFonts().labelFont
        bracket_font = qtg.QFont(label_font)
        bracket_font.setPointSize(30)
        # The bracket_font size is not great. It is currently OK
        # for a 2x2 matrix, but would probably look weird otherwise

        labels = []

        # (1) Label, if needed
        if label:
            labels.append(qtw.QLabel('M = '))
            labels[0].setFont(label_font)

        # (1) Brackets
        labels.extend([qtw.QLabel('('), qtw.QLabel(')')])
        for bracket in labels[-2:]:
            bracket.setFont(bracket_font)

        # (3) matrix elements, with their validator
        if self._dtype == int:
            validator = qtg.QIntValidator()
        else:
            validator = qtg.QDoubleValidator()
        for m_ij in self._ctrls.ravel():
            m_ij.setFont(label_font)
            m_ij.setMaximumWidth(35 if self._dtype == int else 50)
            m_ij.setSizePolicy(qtw.QSizePolicy.Fixed,
                               qtw.QSizePolicy.Preferred)
            m_ij.setValidator(validator)
            m_ij.ensurePolished()

        for qlabel in labels:
            qlabel.ensurePolished()

        # Lay out the widgets
        # (i) Matrix elements
        matrix_layout = qtw.QGridLayout()
        matrix_layout.setSpacing(2)
        matrix_layout.setContentsMargins(0, 0, 0, 0)

        if len(self._ctrls.shape) == 2:
            matrix = self._ctrls
        else:
            matrix = [self._ctrls]
        for i, row in enumerate(matrix):
            for j, m_ij in enumerate(row):
                matrix_layout.addWidget(m_ij, i, j)

        # (ii) put it all together
        layout = qtw.QHBoxLayout()
        for qlabel in labels:
            layout.addWidget(qlabel)
            layout.setAlignment(qlabel, qtc.Qt.AlignCenter)
        layout.insertLayout(len(labels) - 1, matrix_layout)
        layout.addStretch(1)

        self.setLayout(layout)

    def _connect(self):
        """Connect relevant signals to handlers."""
        # NB: textEdited is emitted only for user edits,
        # not for programmatic ones [like setText(...)]
        for m_ij in self._ctrls.ravel():
            m_ij.textEdited.connect(self._on_element_edited)

    @dev_.print_call
    def _on_element_edited(self, new_element):
        """React on a change of one of the matrix elements.

        This is the slot associated with the 'textEdited'
        signal of the matrix element QLineEdit(s).

        Returns
        -------
        None.

        Emits
        -----
        matrix_edited(new_matrix) if the edited entry is valid.
        """
        validate = self._ctrls[0, 0].validator().validate
        # The '0' below is the cursor position, unused by both
        # QIntValidator and QDoubleValidator, but mandatory from
        # the signature. Both validators return a tuple:
        # (state, text, cursor_pos), where state is what we need
        if not validate(new_element, 0)[0] == qtg.QValidator.Acceptable:
            return
        self.matrix_edited.emit(self.matrix)
