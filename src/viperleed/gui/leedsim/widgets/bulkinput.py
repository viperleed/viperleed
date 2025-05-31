"""Module bulkinput of viperleed.gui.leedsim.widgets.

Defines the BulkInput widget used within the NewFileDialog class.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2021-06-01'
__license__ = 'GPLv3+'

import PyQt5.QtCore as qtc
import PyQt5.QtWidgets as qtw

from viperleed.gui.leedsim.dialogs.dialogbulk3dsym import Bulk3DSymDialog
from viperleed.gui.leedsim.widgets.latticeinput import LatticeInput
from viperleed.gui.widgets.lib import AllGUIFonts

from viperleed.gui import decorators as dev_


class BulkInput(LatticeInput):
    """Input of bulk lattice parameters and symmetry.

    The controls are not initialized with up-to-date values.
    Call update_controls_from_lattice() to update from the data.
    """

    bulk_3d_operations_changed = qtc.pyqtSignal()

    def __init__(self, parent=None, bulk_lattice=None):
        """Initialize widget.

        Parameters
        ---------
        parent : PyQt5.QtWidgets.QWidget, default=None
            Parent widget that 'contains' this instance.
        bulk_lattice : Lattice2D
            The lattice instance that will be modified
            when the controls in this instance are edited
        """
        kwargs = {'parent': parent,
                  'title': 'Bulk Lattice',
                  'lattice': bulk_lattice,
                  'with_labels': True,
                  'edit_enabled': True}
        super().__init__(**kwargs)

        self.bulk_3d_sym_dialog = Bulk3DSymDialog(self.window())

    def _compose(self, with_labels):
        """Extend _compose to handle the 3D symmetry input.

        Parameters
        ----------
        with_labels : bool
            If True, labels for the controls are placed on the left.

        Returns
        -------
        None.
        """
        super()._compose(True)

        # Extra control for the bulk 3D symmetry
        self._ctrls['3d_sym'] = qtw.QPushButton('Add bulk sym.\noperations')
        self._ctrls['3d_sym'].setFont(AllGUIFonts().buttonFont)
        self._ctrls['3d_sym'].setSizePolicy(qtw.QSizePolicy.Fixed,
                                             qtw.QSizePolicy.Preferred)
        self._ctrls['3d_sym'].setMinimumWidth(self._ctrls['a'].width())
        self._ctrls['3d_sym'].ensurePolished()
        # Prevent the button to be fired automatically when
        # pressing 'Enter' or 'Return', unless it has focus
        self._ctrls['3d_sym'].setAutoDefault(False)

        layout = self.layout()
        layout.addWidget(self._ctrls['3d_sym'], 6, 1, 1, 1)

    def _connect(self):
        """Extend _connect to handle the 3D symmetry input."""
        super()._connect()
        self._ctrls['3d_sym'].clicked.connect(self._on_bulk_3d_pressed)

    @dev_.print_call
    def _on_bulk_3d_pressed(self, __checked):
        """Open the dialog for editing screws/glides.

        This is the slot connected to the
        clicked signal of the '3d_sym' button.

        Parameters
        ----------
        __checked : bool
            Checked state of the button. Unused, but necessary
            from the signature of the clicked signal.

        Returns
        -------
        None

        Emits
        -----
        bulk_3d_operations_changed
            If the user selected a different set of operations
        """
        self.bulk_3d_sym_dialog.update_operations(self.lattice)

        if self.bulk_3d_sym_dialog.exec() == qtw.QDialog.Accepted:
            old_3d = self.lattice.group.screws_glides
            operations = self.bulk_3d_sym_dialog.extra_operations()
            self.lattice.group.set_screws_glides(operations,
                                                 self.lattice.cell_shape)
            new_3d = self.lattice.group.screws_glides
            if set(new_3d) != set(old_3d):  # Order does not matter
                print("###     o--->", self.__class__.__name__,
                      "-- about to emit bulk_3d_operations_changed")
                self.bulk_3d_operations_changed.emit()

    @dev_.print_call
    def _update_lattice_group(self, new_group):
        """Update the PlaneGroup of the lattice, if needed.

        This is a re-implementation of the base-class method,
        needed to correctly handle the 3D bulk operations.

        Parameters
        ----------
        new_group : str
            Hermann-Maugin name of the new plane group
        """
        print("    (reimplemented)")
        if new_group == self.lattice.group.group:
            return

        # Store the value of screws_glides to set it back
        bulk_3d = self.lattice.group.screws_glides

        # Edit the group
        self.lattice.group = new_group
        self.lattice.group.set_screws_glides(bulk_3d, self.lattice.cell_shape)

        # Update the options for extra bulk operations
        self.bulk_3d_sym_dialog.update_operations(self.lattice)

        # And deactivate the button in case there is nothing to add
        self._ctrls['3d_sym'].setEnabled(self.bulk_3d_sym_dialog.n_extra_ops)
