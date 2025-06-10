"""Module latticeinput of viperleed.gui.leedsim.widgets.

Defines the LatticeInput class, a widget for the interactive input of
lattice parameters and plane group that define a lattice.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2021-06-01'
__license__ = 'GPLv3+'

import re

import numpy as np
import PyQt5.QtCore as qtc
import PyQt5.QtWidgets as qtw

from viperleed.gui.classes.lattice2d import Lattice2D
from viperleed.gui.classes.planegroup import PlaneGroup
from viperleed.gui.widgets.lib import AllGUIFonts
from viperleed.gui.widgets.lib import change_control_text_color

from viperleed.gui import decorators as dev_

ANGSTROM = ' \u212b'
DEGREES = '\u00b0'

# TODO: replace this with a math parser
MATCH_LATTICE_PARAM = re.compile(r'^\s*(\d+(\.\d+)?)\s*[\u212b\u00b0]?\s*$')


class LatticeInput(qtw.QWidget):
    """Generic input for a lattice.

    The underlying Lattice2D instance, used upon
    construction, is updated according to the
    interactive user input.

    The controls are not initialized with up-to-date values.
    Call update_controls_from_lattice() to update from the data.
    """

    shape_changed = qtc.pyqtSignal(str)  # Selected shape
    group_changed = qtc.pyqtSignal(str)  # Selected group
    lattice_parameters_changed = qtc.pyqtSignal(np.ndarray)  # Basis
    need_high_sym_reduction = qtc.pyqtSignal()
    high_sym_pressed = qtc.pyqtSignal()

    def __init__(self, **kwargs):
        """Initialize widget.

        Parameters
        ----------
        parent : PyQt5.QtWidgets.QWidget or None, optional
            Parent widget that 'contains' this instance.
            Default is None.
        lattice : Lattice2D or None, optional
            The lattice instance that will be modified when the
            controls in this instance are edited. If None or not
            given, it should be set via the .lattice property
            before any functionality is available. Default is None.
        title : str, optional
            Title to be used for this widget. It is placed
            centered horizontally, on top. Default is ''.
        with_labels : bool, optional
            If True, labels for the controls are placed on the left.
            Default is False.
        edit_enabled : bool, optional
            If False, all the controls, except for the one
            of the plane group, are not editable, and serve
            only for display purposes. Default is False.

        Returns
        -------
        None.
        """
        super().__init__(kwargs.get('parent', None))

        self.title = kwargs.get('title', '')
        self._lattice = kwargs.get('lattice', None)
        self.edit_enabled = kwargs.get('edit_enabled', False)

        # 'hi_sym' is a button to transform the lattice to the highest
        # possible symmetry, i.e., with shortest vectors, as close to
        # orthogonal as possible.
        self._ctrls = {'shape': qtw.QComboBox(),
                        'a': qtw.QLineEdit(''),
                        'b': qtw.QLineEdit(''),
                        'alpha': qtw.QLineEdit(''),
                        'group': qtw.QComboBox(),
                        'hi_sym': qtw.QPushButton()}

        # Keep track if the user input is OK or not. Used
        # to signal the user that something is wrong. For
        # now this is always True.
        self.__valid_input = True

        # Keep track of the last acceptable
        # lattice basis given by the user
        self.__last_acceptable_basis = np.zeros((2,2))

        self._compose(kwargs.get('with_labels', False))
        self._connect()

    @property
    def basis(self):
        """Return an array of real-space basis vectors.

        Values are taken from the text in the controls.

        Returns
        -------
        basis : numpy.ndarray or None
            When a tuple, a == basis[0], b == basis[1].
            Returns None if the text input is invalid.
        """
        params = [MATCH_LATTICE_PARAM.match(self._ctrls[k].text())
                  for k in ('a', 'b', 'alpha')]
        if not all(params):
            return None

        try:
            a_len, b_len, alpha = [float(match.group(1)) for match in params]
        except ValueError:
            # Some parameter can't be converted to float
            return None
        if a_len > 0 and b_len > 0:
            alpha = np.radians(alpha)
            return np.array([[a_len, 0],
                             [b_len*np.cos(alpha), b_len*np.sin(alpha)]])
        return None

    @property
    def group(self):
        """Return the combo box used as plane group input."""
        return self._ctrls['group']

    @property
    def lattice(self):
        """Return the underlying Lattice2D."""
        if self._lattice is None:
            raise RuntimeError("LatticeInput: underlying lattice never set.")
        return self._lattice

    @lattice.setter
    def lattice(self, new_lattice):
        """Set the underlying lattice.

        Parameters
        ---------
        new_lattice : viperleed.Lattice

        Raises
        ------
        TypeError
            If new_lattice is not a viperleed.Lattice
        """
        if not isinstance(new_lattice, Lattice2D):
            raise TypeError("LatticeInput: invalid lattice argument "
                            f"type {type(new_lattice).__name__}. "
                            "Expected Lattice2D")
        self._lattice = new_lattice
        self.__last_acceptable_basis = new_lattice.basis
        self.update_controls_from_lattice()

    @property
    def is_high_symmetry(self):
        """Return whether the lattice has the highest symmetry."""
        return self.valid_input

    @property
    def top_left_global_handle(self):
        """Return the global position of the top-left alignment handle.
        
        The alignment handle is the global position of the top-left
        corner of the 'shape' control. Can be used for aligning
        with respect to the 'bulk'.

        Returns
        -------
        top_left : QPoint
        """
        shape = self._ctrls['shape']
        return shape.parent().mapToGlobal(shape.pos())

    @property
    def valid_input(self):
        """Return whether the user input is OK."""
        return self.__valid_input

    @dev_.print_call
    def check_high_symmetry(self):
        """Check if lattice is highest symmetry and update controls.

        This only makes the 'high symmetry reduction' button visible
        or not, depending on whether the underlying lattice is at
        the highest possible symmetry.

        Returns
        -------
        None.
        """
        need_high_sym = not self.lattice.is_high_symmetry()
        self._ctrls['hi_sym'].setVisible(need_high_sym)
        if need_high_sym:
            self.__valid_input = False
            # self.need_high_sym_reduction.emit()
        else:
            self.__valid_input = True

    def compatible_groups(self):
        """Return a tuple of groups compatible with lattice."""
        shape = self.lattice.cell_shape
        return PlaneGroup.groups_compatible_with(shape)

    def flash_high_symm_button(self):
        """Flash on/off the high-symmetry button.
        
        This is useful to get the sizeHint of the control
        adjust to the largest it can ever get. Can be called
        if the dialog has never been shown before with the
        purpose of alignment.
        """
        self._ctrls['hi_sym'].show()
        self._ctrls['hi_sym'].hide()

    @dev_.print_call
    def group_from_lattice_and_update_options(self):
        """Update the group combo box with appropriate entries.

        Emits
        -----
        group_changed
            If the selected group has actually changed.
        """
        lattice_group = self.lattice.group.group
        compatible_groups = self.compatible_groups()
        if lattice_group not in compatible_groups:
            self.lattice.group = 'p1'

        group_combo = self._ctrls['group']
        old_group = group_combo.currentText()
        group_combo.clear()
        group_combo.addItems(compatible_groups)
        group_combo.setCurrentText(lattice_group)
        group_combo.setEnabled(group_combo.count() > 1)

        if lattice_group != old_group:
            print("###     o--->", self.__class__.__name__,
                  "-- about to emit group_changed")
            self.group_changed.emit(lattice_group)

    @dev_.print_call
    def update_controls_from_lattice(self):
        """Use the underlying lattice to update the controls.

        Emits
        -----
        [need_high_sym_reduction()]
        [group_changed(new_group)]
        [shape_changed(new_shape)]
        """
        self._update_lattice_params_ctrls_from_lattice()
        self.check_high_symmetry()
        shape = self.lattice.cell_shape
        _ctrl_shape = self._ctrls['shape'].currentText()
        if _ctrl_shape == '_None':
            # It's the first time we run after setting a lattice.
            # Remove the special '_None' from the shape combo
            self._ctrls['shape'].removeItem(
                self._ctrls['shape'].currentIndex()
                )
        if shape != _ctrl_shape:
            self._ctrls['shape'].setCurrentText(self.lattice.cell_shape)
            self.update_lattice_restrictions()
            self.group_from_lattice_and_update_options()
            print("###     o--->", self.__class__.__name__,
                  "-- about to emit shape_changed")
            self.shape_changed.emit(shape)

    @dev_.print_call
    def update_lattice_basis(self, basis=None):
        """Update the basis of the underlying lattice.

        Parameters
        ----------
        basis : numpy.ndarray or None, optional
            The basis to be stored in the underlying lattice.
            Shape == (2, 2), with a == basis[0] and b == basis[1].
            If None or not given, the basis is fetched from the
            last acceptable input. If given and acceptable, the
            last acceptable basis is updated accordingly.

        Returns
        -------
        None.
        """
        if basis is None:
            basis = self.__last_acceptable_basis
        self.lattice.basis = basis
        self.__last_acceptable_basis = self.lattice.basis

    @dev_.print_call
    def update_lattice_restrictions(self):
        """Update restrictions on the lattice, given the current shape.

        Enable/disable controls that would be editable given
        the currently selected cell shape. The value of
        self.edit_enabled determines whether the controls
        that could allow user input are actually enabled.
        """
        # self._on_shape_changed()
        shape = self._ctrls['shape'].currentText()

        # Enable controls that make sense for the
        # current shape, and set their values.
        # (1) a and b
        if shape in ('Square', 'Rhombic', 'Hexagonal'):
            self._ctrls['b'].setEnabled(False)
            self._ctrls['b'].setText(self._ctrls['a'].text())
        else:
            self._ctrls['b'].setEnabled(self.edit_enabled)

        # (2) alpha
        if shape in ('Square', 'Rectangular', 'Hexagonal'):
            self._ctrls['alpha'].setEnabled(False)
            if shape in ('Square', 'Rectangular'):
                self._ctrls['alpha'].setText(f'90{DEGREES}')
            else:
                self._ctrls['alpha'].setText(f'120{DEGREES}')
        else:
            self._ctrls['alpha'].setEnabled(self.edit_enabled)

    def _compose(self, with_labels):
        """Place children widgets.

        Parameters
        ----------
        with_labels : bool
            If True, labels for the controls are placed on the left.

        Returns
        -------
        None.
        """
        # Fonts
        label_font = AllGUIFonts().labelFont

        # Title
        title = qtw.QLabel(self.title)
        title.setFont(label_font)
        title.ensurePolished()

        # Labels of controls, if needed
        labels = []
        if with_labels:
            labels = [qtw.QLabel('Shape'), qtw.QLabel('a = '),
                      qtw.QLabel('b = '), qtw.QLabel('\u03b1 = '),
                      qtw.QLabel('Group')]
        for label in labels:
            label.setFont(label_font)
            label.ensurePolished()
            label.adjustSize()    # Not sure if this is actually needed

        if with_labels:
            label_width = max(label.width() for label in labels)
            for label in labels:
                label.setMaximumWidth(label_width)

        # Actual controls
        self._compose_ctrl_shape()
        self._compose_ctrl_a_b_alpha()
        self._compose_ctrl_group()
        self._compose_ctrl_high_sym_button()

        # Lay controls out:
        # (i) Set up layout
        layout = qtw.QGridLayout()
        layout.setSpacing(5)
        layout.setContentsMargins(5, 5, 5, 5)  # L, T, R, B
        layout.setSizeConstraint(qtw.QLayout.SetFixedSize)

        # (ii) Title and labels
        layout.addWidget(title, 0, 1, 1, 2)  # Row, col, n_rows, n_cols
        layout.setAlignment(title, qtc.Qt.AlignCenter)
        for i, label in enumerate(labels):
            layout.addWidget(label, i + 1, 0, 1, 1)
            layout.setAlignment(label, qtc.Qt.AlignRight | qtc.Qt.AlignVCenter)

        # (iii) Controls
        ctrl_column = 1 if with_labels else 0
        for i, ctrl in enumerate(self._ctrls.values()):
            # Leave one empty row to host the 'extra bulk symmetry'
            # button right below the 'group' combo box, i.e., above
            # the 'hi_sym' button
            if ctrl == self._ctrls['hi_sym']:
                i += 1
            layout.addWidget(ctrl, i + 1, ctrl_column, 1, 1)

        self.setLayout(layout)

    def _compose_ctrl_a_b_alpha(self):
        """Set up controls of lattice parameters."""
        for key in ('a', 'b', 'alpha'):
            self._ctrls[key].setFont(AllGUIFonts().labelFont)
            self._ctrls[key].setSizePolicy(qtw.QSizePolicy.Fixed,
                                           qtw.QSizePolicy.Fixed)
            self._ctrls[key].setMaximumWidth(self._ctrls['shape'].width())
            self._ctrls[key].ensurePolished()
            self._ctrls[key].adjustSize()
            self._ctrls[key].setEnabled(self.edit_enabled)

    def _compose_ctrl_group(self):
        """Set up control for plane group selection."""
        ctrl = self._ctrls['group']

        ctrl.setFont(AllGUIFonts().smallTextFont)
        ctrl.setSizePolicy(qtw.QSizePolicy.Fixed, qtw.QSizePolicy.Preferred)
        ctrl.setMinimumWidth(self._ctrls['a'].width())
        ctrl.ensurePolished()

    def _compose_ctrl_high_sym_button(self):
        """Set up control to reduce to highest symmetry."""
        ctrl = self._ctrls['hi_sym']

        ctrl.setText('Make high\nsymmetry')
        ctrl.setFont(AllGUIFonts().buttonFont)
        ctrl.setSizePolicy(qtw.QSizePolicy.Fixed, qtw.QSizePolicy.Preferred)
        ctrl.setMinimumWidth(self._ctrls['shape'].width())
        ctrl.ensurePolished()
        ctrl.setAutoDefault(False)
        ctrl.hide()
        change_control_text_color(ctrl, qtc.Qt.red)

    def _compose_ctrl_shape(self):
        """Set up 'shape' control."""
        ctrl = self._ctrls['shape']

        ctrl.setFont(AllGUIFonts().smallTextFont)
        ctrl.addItems(['Oblique', 'Rectangular', 'Square',
                       'Rhombic', 'Hexagonal'])

        # The special item '_None' is used as an initial value
        # only if a lattice was not given at instantiation, as
        # we need to emit a shape_changed once the lattice
        # is set the first time. See update_controls_from_lattice
        if self._lattice is None:
            ctrl.addItem('_None')
            ctrl.setCurrentText('_None')
        ctrl.setCurrentText('')
        ctrl.setSizePolicy(qtw.QSizePolicy.Fixed, qtw.QSizePolicy.Preferred)
        ctrl.ensurePolished()
        ctrl.adjustSize()
        ctrl.setEnabled(self.edit_enabled)

    def _connect(self):
        """Connect relevant signals to handlers."""
        self._ctrls['shape'].activated.connect(self._on_shape_changed)
        self._ctrls['group'].activated.connect(self._on_group_changed)

        for param in ('a', 'b', 'alpha'):
            param = self._ctrls[param]
            # NB: Using textEdited that is not
            # emitted on programmatic text changes
            param.textEdited.connect(self._on_lattice_param_edited)

            # editingFinished is emitted after Enter/Return
            # or a loss of focus. Hence, if there was an
            # edit, it is emitted after textEdited.
            param.editingFinished.connect(self._on_lattice_param_edit_finished)
        self.group_changed.connect(self._update_lattice_group)
        self._ctrls['hi_sym'].clicked.connect(self._on_high_sym_pressed)

    @dev_.print_call
    def _on_group_changed(self, new_group_index):
        """React on change of plane group.

        Triggers an update of the plane group of the lattice.

        This is the slot associated with the 'activated'
        signal of the group combo box.

        Parameters
        ----------
        new_group_index : int or None, default=None
            Index of the selected plane group. No action
            occurs if new_group_index is None.

        Returns
        -------
        None.

        Emits
        -----
        group_changed(new_group_name)
        """
        # NB: the 'activated' signal is emitted
        # even if the user choice did not change.
        new_group = self._ctrls['group'].itemText(new_group_index)
        if new_group == self.lattice.group.group:
            return
        print("###     o--->", self.__class__.__name__,
              "-- about to emit group_changed")
        self.group_changed.emit(new_group)

    @dev_.print_call
    def _on_high_sym_pressed(self, __checked):
        """React on a user press of the 'make high symmetry button.

        This is the slot connect with the 'clicked' signal of
        the button.

        Parameters
        ----------
        __checked : int
            Checked state of button. Unused but necessary for
            the signature of the clicked signal of QPushButton.

        Returns
        -------
        None

        Emits
        -----
        high_sym_pressed()
            The signal is emitted before executing the function,
            so that it is possible to manipulate the values of
            controls before reducing to high symmetry, if needed
        [need_high_sym_reduction()]
            Should never happen
        [shape_changed(new_shape)]
            If shape has actually changed
        [group_changed(new_group)]
            If group has actually changed
        lattice_parameters_changed
            After reduction to high symmetry is done
        """
        print("###     o--->", self.__class__.__name__,
              "-- about to emit high_sym_pressed")
        self.high_sym_pressed.emit()
        self.lattice.make_high_symmetry()
        self.update_controls_from_lattice()
        print("###     o--->", self.__class__.__name__,
              "-- about to emit lattice_parameters_changed")
        self.lattice_parameters_changed.emit(self.lattice.basis)

    @dev_.print_call
    def _on_lattice_param_edited(self, new_param):
        """React on change of a lattice parameter input.

        If the text input is valid, it causes update of
        the lattice shape control, as well as of the
        underlying lattices.

        This is the slot associated with the 'textEdited'
        signal of the lattice parameter QLineEdit(s). It
        can also be called without arguments to cause a
        similar behavior.

        Parameters
        ----------
        new_param : str
            String representation of the text in the input widget

        Returns
        -------
        None.
        """
        if MATCH_LATTICE_PARAM.match(new_param):
            # This call allows to enforce parameters to
            # update in a way that conserves the shape
            # currently selected
            self.update_lattice_restrictions()

        basis = self.basis
        if basis is not None:
            # Text entered is acceptable. The underlying lattice
            # will be updated only after the editing is finished
            self.__last_acceptable_basis = basis

    @dev_.print_call
    def _on_lattice_param_edit_finished(self):
        """React on a completed change of a lattice parameter input.

        Causes a reformatting of the text to also include units. The
        values in the controls are fetched from the underlying lattice.

        This is the slot associated with the 'editingFinished' signal
        of the lattice parameter QLineEdit(s).

        Returns
        -------
        None.

        Emits
        -----
        group_changed(new_group)
            If group actually changed
        shape_changed(new_shape)
            If shape actually changed
        lattice_parameters_changed(new_basis)
            Always
        need_high_sym_reduction()
            If the lattice does not have the highest possible symmetry
        """
        if np.allclose(self.lattice.basis, self.__last_acceptable_basis):
            return
        self.update_lattice_basis()
        self.update_controls_from_lattice()
        print("###     o--->", self.__class__.__name__,
              "-- about to emit lattice_parameters_changed")
        self.lattice_parameters_changed.emit(self.lattice.basis)
        if not self.is_high_symmetry:
            print("###     o--->", self.__class__.__name__,
                  "-- about to emit need_high_sym_reduction")
            self.need_high_sym_reduction.emit()

    @dev_.print_call
    def _on_shape_changed(self, shape_index):
        """Set lattice parameters editable or fixed depending on shape.

        Also updates the underlying lattice, including its lattice
        parameters and its plane group, and picks a new list of
        plane groups that are compatible with the selected shape.

        This is the slot associated with the 'activated' signal
        of the shape combo box. It can also be called without
        arguments to cause a similar behavior.

        Parameters
        ----------
        shape_index : int
            Index of the combo-box entry selected

        Returns
        -------
        None.

        Emits
        -----
        [need_high_sym_reduction()]
        [group_changed(new_group)]
        shape_changed(new_shape).
        """
        # The 'activated' signal is emitted even if
        # no change in the actual selection occurred
        new_shape = self._ctrls['shape'].itemText(shape_index)
        if new_shape == self.lattice.cell_shape:
            return

        # Enable/disable controls that should be editable
        # given the current shape. Also restrict values.
        self.update_lattice_restrictions()

        # Update the underlying lattice
        self.update_lattice_basis(self.basis)

        # Check whether it is at highest symmetry
        self.check_high_symmetry()  # May emit need_high_sym_reduction

        # And update the group control
        self.group_from_lattice_and_update_options()  # May emit group_changed

        print("###     o--->", self.__class__.__name__,
              "-- about to emit shape_changed")
        self.shape_changed.emit(new_shape)

        if not self.is_high_symmetry:
            print("###     o--->", self.__class__.__name__,
                  "-- about to emit need_high_sym_reduction")
            self.need_high_sym_reduction.emit()

    @dev_.print_call
    def _update_lattice_group(self, new_group):
        """Update the PlaneGroup of the lattice, if needed.

        This method can be used as the slot of group_changed.

        Parameters
        ----------
        new_group : str
            Hermann-Maugin name of the new plane group
        """
        if new_group != self.lattice.group.group:
            self.lattice.group = new_group

    @dev_.print_call
    def _update_lattice_params_ctrls_from_lattice(self):
        """Update controls of lattice parameters from lattice.

        Also causes a reformatting of the text to include units.

        Returns
        -------
        None.
        """
        a_len, b_len, alpha = self.lattice.lattice_parameters
        self._ctrls['a'].setText(f"{a_len:.4f}{ANGSTROM}")
        self._ctrls['b'].setText(f"{b_len:.4f}{ANGSTROM}")
        self._ctrls['alpha'].setText(f"{alpha:.4f}{DEGREES}")
