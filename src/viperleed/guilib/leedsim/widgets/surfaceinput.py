"""Module surfaceinput of viperleed.guilib.leedsim.widgets.

======================================
  ViPErLEED Graphical User Interface
======================================

Defines the SurfaceStructureInput class, a widget for the interactive
input of lattice parameters, plane group, superlattice matrix, and
Wood's notation super-periodicity of a reconstructed surface.

Created: 2021-06-01
Author: Michele Riva
"""
import numpy as np
import PyQt5.QtCore as qtc
import PyQt5.QtWidgets as qtw

from viperleed.guilib.classes.lattice2d import Lattice2D
from viperleed.guilib.classes.planegroup import PlaneGroup
from viperleed.guilib.leedsim.classes.woods import Woods
from viperleed.guilib.leedsim.classes.woods import WoodsSyntaxError
from viperleed.guilib.leedsim.classes.woods import MatrixIncommensurateError
from viperleed.guilib.leedsim.classes.woods import WoodsInvalidForBasisError
from viperleed.guilib.leedsim.classes.woods import WoodsNotRepresentableError
from viperleed.guilib.leedsim.widgets.editablematrix import EditableMatrix
from viperleed.guilib.leedsim.widgets.latticeinput import LatticeInput
from viperleed.guilib.widgetslib import AllGUIFonts
from viperleed.guilib.widgetslib import change_control_text_color

from viperleed.guilib import decorators as dev_


# TODO: handle equivalent Woods
# TODO: unformatted input string stays in Woods combo! They should be
#       removed whenever a item is selected. Invalid ones should also
#       be not re-selectable


def _wrap_compatible_groups(lattice_input, surface_input):
    """Wrap reimplementation of compatible_groups.

    This wrapper allows reimplementing the method
    LatticeInput.compatible_groups, allowing each
    call to the method -- internally done in
    LatticeInput instances -- to access also
    information contained in surface_input.

    This is necessary to handle changes of the bulk
    group that affect the list of surface groups.

    Parameters
    ----------
    lattice_input : LatticeInput
        The LatticeInput instance for which the method
        should be reimplemented
    surface_input : SurfaceStructureInput
        The SurfaceStructureInput that contains
        lattice_input, from which the information
        on the bulk group is drawn.

    Returns
    -------
    _reimplemented : callable
        Reimplemented method

    Raises
    ------
    TypeError
        If lattice_input is not a LatticeInput instance
    TypeError
        If surface_input is not a SurfaceStructureInput instance
    """
    if not isinstance(lattice_input, LatticeInput):
        raise TypeError("Reimplementation of compatible_groups was not called "
                        "with a LatticeInput instance as first argument.")
    if not isinstance(surface_input, SurfaceStructureInput):
        raise TypeError("Reimplementation of compatible_groups was "
                        "not called with a SurfaceStructureInput "
                        "instance as second argument.")

    def _reimplemented():
        """Reimplement compatible_groups.

        The reimplementation picks the list of groups taking
        into account also the bulk.group of the surface_input
        instance.

        Returns
        -------
        compatible_groups : tuple
            tuple of strings of the Hermann-Maugin names of
            groups compatible with both the lattice_input
            LatticeInput instance and the operations of the
            bulk group of the surface_input instance.

        Emits
        -----
        lattice_input.group_changed
            If the group of the lattice underlying
            lattice_input has actually changed
        """
        print(f"{lattice_input.__class__.__name__}",
              "compatible_groups. Reimplemented by SurfaceStructureInput")

        # (1) Project the bulk operations to the
        #     surface via the superlattice matrix
        #     Rather recalculate the matrix than
        #     taking it from the controls, that
        #     may not be up to date
        transform = np.dot(surface_input.lattice.basis,
                           np.linalg.inv(surface_input.bulk.basis))
        #     Exclude 3D operations, as they only
        #     contribute to which domains exist
        projected_bulk_ops = surface_input.bulk.group.transform(
            transform,
            include_3d=False
            )

        # Now projected_bulk_ops is a tuple of (2, 2)-arrays of floats.
        # Round off those elements that have only 'integers', and
        # throw away all the others (surface operations only have ints)
        projected_bulk_ops = tuple(op.round().astype(int)
                                   for op in projected_bulk_ops
                                   if np.all(np.abs(op - op.round()) < 1e-3)
                                   )

        # (2) Pick the list of compatible groups, also
        #     accounting for the operations just found
        shape = lattice_input.lattice.cell_shape
        return PlaneGroup.groups_compatible_with(shape, projected_bulk_ops)

    return _reimplemented


class SurfaceStructureInput(qtw.QWidget):
    """A widget for interactive input of a surface structure.

    The controls are not initialized with up-to-date values.
    Call update_controls() to cause an update from the data.
    """

    user_gave_invalid_input = qtc.pyqtSignal(str, int)  # msg, timeout
    surface_changed = qtc.pyqtSignal()
    need_high_sym_reduction = qtc.pyqtSignal()

    def __init__(self, parent=None, bulk_lattice=None, surface_lattice=None):
        """Initialize widget.

        Parameters
        ----------
        parent : PyQt5.QtWidgets.QWidget or None, default=None
            Parent widget that 'contains' this instance.
        bulk_lattice : viperleed.Lattice
            The real-space lattice of the bulk. Can
            be read/written via the .bulk property
        surface_lattice : Lattice2D or None, default=None
            The real-space lattice of the surface structure.
            Can be read/written via the .lattice property.
            If not given or None, it should be set using
            .lattice before any functionality is available.

        Returns
        -------
        None.

        Raises
        ------
        ValueError
            If bulk_lattice is not given.
        TypeError
            If bulk_lattice is not a viperleed.Lattice.
        """
        if bulk_lattice is None:
            raise TypeError("SurfaceStructureInput: missing required "
                            "(keyword) argument bulk_lattice")

        super().__init__(parent)
        self.__bulk = None         # Set by next line, for extra checks
        self.bulk = bulk_lattice
        self.__woods = Woods(bulk_basis=self.bulk.basis, style='unicode')

        self._ctrls = {'woods': qtw.QComboBox(),
                       'superlattice': EditableMatrix(),
                       'lattice': LatticeInput(lattice=surface_lattice)}
        if surface_lattice:
            self.update_woods_list_and_selection()

        # Replace .compatible_groups() of the underlying
        # LatticeInput to handle changes of the bulk group
        self._ctrls['lattice'].compatible_groups = (
            _wrap_compatible_groups(self._ctrls['lattice'], self)
            )

        # Keep track if the user input is OK or not. Used
        # to signal the user that something is wrong
        self.__valid_input = True

        self._compose()
        self._connect()


    @property
    def bulk(self):
        """Return the Lattice2D of the bulk."""
        return self.__bulk

    @bulk.setter
    def bulk(self, new_lattice):
        """Set the Lattice2D of the bulk to new_lattice."""
        if isinstance(new_lattice, Lattice2D):
            self.__bulk = new_lattice
            return
        raise TypeError("SurfaceStructureInput: invalid bulk lattice "
                        f"type {type(new_lattice).__name__}. Expected "
                        "Lattice2D or None")

    @property
    def lattice(self):
        """Return the underlying Lattice2D."""
        return self._ctrls['lattice'].lattice

    @lattice.setter
    def lattice(self, new_lattice):
        """Set the underlying lattice.

        Parameters
        ---------
        new_lattice : viperleed.Lattice

        Raises
        ------
        TypeError
            If new_lattice is not a viperleed.Lattice.
        ValueError
            If the basis of new_lattice is incommensurate
            with respect to the one of the bulk.
        """
        if not isinstance(new_lattice, Lattice2D):
            raise TypeError("SurfaceStructureInput: invalid lattice "
                            f"argument type {type(new_lattice).__name__}. "
                            "Expected Lattice2D")
        self._update_superlattice_from_surf_basis(new_lattice.basis)
        self._ctrls['lattice'].lattice = new_lattice
        self.update_woods_list_and_selection()

    @property
    def superlattice(self):
        """Return the superlattice matrix as a numpy.ndarray."""
        return self._ctrls['superlattice'].matrix

    @superlattice.setter
    def superlattice(self, new_matrix):
        """Set superlattice to a new matrix.

        Does not emit EditableMatrix.matrix_edited!
        """
        self._ctrls['superlattice'].set_matrix(new_matrix)

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
        return self._ctrls['lattice'].top_left_global_handle

    @property
    def valid_input(self):
        """Return whether the user input is OK."""
        return (self.__valid_input
                and self._ctrls['lattice'].valid_input)

    @property
    def woods(self):
        """Return the Woods object used to handle Wood's notation."""
        return self.__woods

    @dev_.print_call
    def on_bulk_group_changed(self, *__args):
        """Pick the appropriate list of groups.

        'Appropriate' means all those groups that are compatible
        with the cell shape and the new bulk group, in the sense
        of PlaneGroup.groups_compatible_with().

        This slot should be connected to the group_changed signal
        of the BulkInput widget that holds the same bulk lattice
        as in self.bulk.

        Parameters
        ----------
        *__args : tuple
            Unused but necessary to comply with the signature
            of signals connected to this method. The new bulk
            group is always taken from the underlying bulk
            lattice, which should be up to date.

        Emits
        -----
        lattice.group_changed
            If the group of the underlying lattice actually changed
        """
        self._ctrls['lattice'].group_from_lattice_and_update_options()

    @dev_.print_call
    def on_bulk_basis_changed(self, *__args):
        """React to a change of the bulk basis.

        The surface basis is changed accordingly,
        then the controls are updated.

        This slot should be connected to the lattice_parameters_changed
        signal of the BulkInput widget that holds the same bulk lattice
        as in self.bulk. It can also be connected to the shape_changed
        signal of the same widget.

        Parameters
        ----------
        *__args : tuple
            Unused but necessary to comply with the signature
            of signals connected to this method. The new bulk
            basis is always taken from the underlying bulk
            lattice, which should be up to date.

        Emits
        -----
        lattice.need_high_sym_reduction()
            If the underlying lattice is not highest symmetry
        need_high_sym_reduction()
            If the underlying lattice is not highest symmetry.
            This signal is just a convenience to prevent the
            need to directly look into the underlying lattice
            control
        lattice.group_changed(new_group)
            If group of underlying lattice actually changed
        lattice.shape_changed(new_shape)
            If shape of underlying lattice actually changed
        """
        lattice = self._ctrls['lattice']
        superlattice = self._ctrls['superlattice'].matrix
        lattice.update_lattice_basis(np.dot(superlattice, self.bulk.basis))
        self.__woods.bulk_basis = self.bulk.basis
        lattice.update_controls_from_lattice()

    @dev_.print_call
    def pick_right_woods(self):
        """Choose the appropriate Wood's representation.

        Causes a change of the text in the Wood's combo
        box. The new text is a valid Wood's notation if
        the superlattice matrix is Wood's-representable,
        and 'None' otherwise.

        Returns
        -------
        None.
        """
        woods_combo = self._ctrls['woods']

        # Set color to black whenever this function is called,
        # as it always will pick the correct Wood's text. This
        # undoes 'unprocessed' color changes that may have
        # happened as a result of an improper text input (set
        # in _on_woods_selected).
        palette = woods_combo.lineEdit().palette()
        palette.setColor(palette.Text, qtc.Qt.black)
        woods_combo.lineEdit().setPalette(palette)

        # See if the superlattice matrix is Wood's-representable
        try:
            new_woods = self.__woods.from_matrix(
                self._ctrls['superlattice'].matrix
                )
        except WoodsNotRepresentableError:
            new_woods = 'None'
        else:
            # Representable. See if we already have it in the
            # list (.findText returns index >= 0 if found)
            found = woods_combo.findText(new_woods,
                                         flags=qtc.Qt.MatchExactly) >= 0
            if not found:
                woods_combo.addItem(new_woods)
                self.__woods.add_example(new_woods, self.bulk.cell_shape)

        woods_combo.setCurrentText(new_woods)

    @dev_.print_call
    def update_controls(self):
        """Set the values in the controls from stored attributes.

        Emits
        -----
        lattice.group_changed
            If lattice group is in fact changed
        lattice.shape_changed
            If lattice cell shape is in fact changed
        """
        # Update the superlattice control from self.lattice.basis.
        # NB: _update_superlattice_from_surf_basis() uses set_matrix()
        # that does not emit matrix_edited. Emission of matrix_edited
        # would also cause the lattice constants stored in self.lattice
        # to be changed from the (not necessarily up to date) controls
        self._update_superlattice_from_surf_basis()

        # Update lattice parameters and group. May emit
        # lattice.group_changed and/or lattice.shape_changed
        self._ctrls['lattice'].update_controls_from_lattice()

    @dev_.print_call
    def update_woods_list_and_selection(self, bulk_shape=''):
        """Fetch examples of Wood's given a bulk-lattice shape.

        Also update the current selection, should the current
        superlattice be representable in Wood's notation.

        This slot should be connected to the shape_changed signal
        of the BulkInput widget that holds the same bulk lattice
        as in self.bulk.

        Parameters
        ----------
        bulk_shape : {'Oblique', 'Rectangular', 'Square',
                      'Rhombic', 'Hexagonal', ''}
            If empty, bulk_shape is taken from the underlying
            bulk lattice. Default is an empty string.
        """
        if not bulk_shape:
            bulk_shape = self.bulk.cell_shape

        woods_combo = self._ctrls['woods']
        woods_combo.clear()
        woods_combo.addItems(sorted(
            ex.string for ex in self.__woods.get_examples(bulk_shape))
            )

        self.pick_right_woods()

    def _compose(self):
        """Place children widgets."""
        # Fonts
        small_txt_font = AllGUIFonts().smallTextFont
        label_font = AllGUIFonts().labelFont

        # (1) Labels
        labels = {'title': qtw.QLabel('Reconstruction periodicity'),
                  'woods': qtw.QLabel("Wood's:")}
        labels['title'].setFont(label_font)
        labels['woods'].setFont(small_txt_font)
        labels['woods'].setSizePolicy(qtw.QSizePolicy.Fixed,
                                      qtw.QSizePolicy.Preferred)

        # (2) Woods combo box
        self._ctrls['woods'].setEditable(True)
        self._ctrls['woods'].setInsertPolicy(
            qtw.QComboBox.InsertAlphabetically
            )
        self._ctrls['woods'].setFont(small_txt_font)
        self._ctrls['woods'].setSizePolicy(qtw.QSizePolicy.Fixed,
                                            qtw.QSizePolicy.Preferred)
        self._ctrls['woods'].setSizeAdjustPolicy(
            qtw.QComboBox.AdjustToContents
            )
        self._ctrls['woods'].setMinimumContentsLength(12)

        # (3) Lay out widgets
        for widg in (*labels.values(), *self._ctrls.values()):
            widg.ensurePolished()

        # (3.i) woods input
        woods_layout = qtw.QHBoxLayout()
        woods_layout.addWidget(labels['woods'])
        woods_layout.addWidget(self._ctrls['woods'])
        woods_layout.addStretch(1)  # Flush all to left

        # (3.ii) superstructure input (= title, then woods, and matrix below)
        superstructure_layout = qtw.QVBoxLayout()
        superstructure_layout.addWidget(labels['title'])
        superstructure_layout.addLayout(woods_layout)
        superstructure_layout.addWidget(self._ctrls['superlattice'])

        # (3.iii) put it together: lattice on the left,
        # superstructure right, both aligned to the top.
        layout = qtw.QHBoxLayout()
        layout.addWidget(self._ctrls['lattice'])
        layout.addLayout(superstructure_layout)
        layout.setAlignment(self._ctrls['lattice'], qtc.Qt.AlignTop)
        layout.setAlignment(superstructure_layout, qtc.Qt.AlignTop)

        self.setLayout(layout)

    def _connect(self):
        """Connect relevant signals to handlers."""
        woods_combo = self._ctrls['woods']
        woods_combo.activated.connect(self._on_woods_selected)
        woods_combo.lineEdit().editingFinished.connect(self._on_woods_selected)

        superlattice = self._ctrls['superlattice']
        superlattice.matrix_edited.connect(self._on_superlattice_changed)

        # For the surface lattice, lattice_parameters_changed is
        # emitted only when the 'reduce to high symmetry' button
        # is pressed (after reduction is already done), as the
        # controls are not editable (would normally be emitted
        # on editingFinished).
        self._ctrls['lattice'].lattice_parameters_changed.connect(
            self._on_high_sym_pressed
            )

        self.need_high_sym_reduction.connect(
            self._ctrls['lattice'].need_high_sym_reduction
            )

    @dev_.print_call
    def _on_high_sym_pressed(self, *__args):
        """React on a request to make the lattice high symmetry.

        Most of the correct updating is already done intrinsically
        in the underlying LatticeInput controls. Here we only
        update the superlattice matrix and pick the appropriate
        Woods.

        Returns
        -------
        None

        Emits
        -----
        surface_changed
        """
        self._update_superlattice_from_surf_basis()
        self.pick_right_woods()
        print("###     o--->", self.__class__.__name__,
              "-- about to emit surface_changed")
        self.surface_changed.emit()

    @dev_.print_call
    def _on_superlattice_changed(self, new_superlattice):
        """React to a change of the superlattice matrix.

        This is the slot associated with the matrix_edited
        signal of the superlattice EditableMatrix.

        Parameters
        ----------
        new_superlattice : numpy.ndarray
            The new superlattice matrix. Shape == (2, 2).

        Returns
        -------
        None.

        Emits
        -----
        user_gave_invalid_input :
            If superlattice matrix is singular
        lattice.need_high_sym_reduction :
            If the the new superlattice is acceptable, but the
            underlying lattice is not highest symmetry
        need_high_sym_reduction :
            If the the new superlattice is acceptable, but the
            underlying lattice is not highest symmetry.
            This is a convenience signal to prevent the need
            to access the underlying LatticeInput widget
        lattice.group_changed(new_group) :
            If the the new superlattice is acceptable, and the
            group of the underlying lattice actually changed
        lattice.shape_changed(new_shape) :
             If the the new superlattice is acceptable, and the
             shape of the underlying lattice actually changed
        surface_changed :
            If the new superlattice is acceptable
        """
        if abs(np.linalg.det(new_superlattice)) < 1e-3:
            self._ctrls['superlattice'].set_text_color(qtc.Qt.red)
            self.__valid_input = False
            self.user_gave_invalid_input.emit('Matrix is singular!', 2000)
            return

        self.__valid_input = True
        self._ctrls['superlattice'].set_text_color(qtc.Qt.black)
        # Update the basis of the underlying lattice
        # and the related controls
        lattice = self._ctrls['lattice']
        lattice.update_lattice_basis(np.dot(new_superlattice, self.bulk.basis))
        lattice.update_controls_from_lattice()

        # And see if the matrix is woods representable
        self.pick_right_woods()

        print("###     o--->", self.__class__.__name__,
              "-- about to emit surface_changed")
        self.surface_changed.emit()

    @dev_.print_call
    def _on_woods_selected(self, *args):
        """React on the selection of a new entry in the combo.

        This is the slot associated with the 'activated' signal
        of the woods combo box and with the 'editingFinished'
        signal of the editable entry of the combo box itself
        (which is in fact a QLineEdit). The signals below are
        emitted only in case the text in the combo box is not
        'None'.

        Emits
        -----
        user_gave_invalid_input
            If woods syntax is wrong, or if selected woods would
            give an incommensurate lattice. Notice that self.valid
            is not changed even if this input is invalid, as the
            lattice is not changed if this input is invalid (taken
            from the unchanged superlattice matrix)
        superlattice.matrix_edited
            Emitted only if the new Wood's is acceptable.
        [lattice.group_changed]
            Emitted only if the surface-lattice group has changed
        [lattice.shape_changed]
            Emitted only if the surface-lattice shape has changed
        surface_changed
            If new superlattice is acceptable (indirectly,
            via self._on_superlattice_changed)
        """
        woods_combo = self._ctrls['woods']
        if args:
            # The call came from a user click on a QComboBox item
            # NB: This happens also if the element chosen did not
            #     actually change!
            woods_txt = woods_combo.itemText(args[0])
        else:
            # The call to the function came from
            # the editingFinished of the QLineEdit
            woods_txt = woods_combo.currentText()

        if woods_txt == 'None':
            return

        # Reformat input
        try:
            self.__woods.string = woods_txt
        except (ValueError, WoodsSyntaxError):
            # Text is not an acceptable Wood's notation
            self.user_gave_invalid_input.emit("Invalid Wood's syntax.", 1000)
            change_control_text_color(woods_combo.lineEdit(), qtc.Qt.red)
            return
        except WoodsInvalidForBasisError:
            self.user_gave_invalid_input.emit(
                f"{woods_txt} incommensurate for this bulk basis.",
                7000
                )
            change_control_text_color(woods_combo.lineEdit(), qtc.Qt.red)
            return

        woods_combo.setCurrentText(self.__woods.string)
        self._ctrls['superlattice'].matrix = self.__woods.matrix
        change_control_text_color(woods_combo.lineEdit(), qtc.Qt.black)

    @dev_.print_call
    def _update_superlattice_from_surf_basis(self, surf_basis=None):
        """Update superlattice control from the given surface basis.

        Parameters
        ----------
        surf_basis : numpy.ndarray or None
            Shape == (2, 2). If None, surf_basis is
            taken from self.lattice. Default is None.

        Returns
        -------
        None.

        Raises
        ------
        MatrixIncommensurateError
            If surf_basis gives an incommensurate superlattice
        """
        if surf_basis is None:
            surf_basis = self._ctrls['lattice'].lattice.basis
        superlattice = np.dot(surf_basis, np.linalg.inv(self.bulk.basis))
        if np.any(np.abs(superlattice - superlattice.round()) > 1e-3):
            raise MatrixIncommensurateError(superlattice,
                                            message="SurfaceStructureInput: "
                                            "bulk and surface bases are "
                                            "incommensurate")
        superlattice = superlattice.round().astype(int)

        # Use .set_matrix() to not emit matrix_edited
        self._ctrls['superlattice'].set_matrix(superlattice)
