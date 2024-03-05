"""Module newfiledialog of viperleed.guilib.leedsim.dialogs.

===============================================
      ViPErLEED Graphical User Interface
===============================================

Defines the NewFileDialog class, which allows user input of bulk
and surface periodicities and space groups

Created: 2020-01-11
Author: Michele Riva
"""
import os

import numpy as np
import PyQt5.QtCore as qtc
import PyQt5.QtGui as qtg
import PyQt5.QtWidgets as qtw

from viperleed import guilib as gl
from viperleed.gui import resources_path
from viperleed.guilib.classes.lattice2d import Lattice2D
from viperleed.guilib.leedsim import LEEDParametersList
from viperleed.guilib.widgetslib import change_control_text_color
# from viperleed.guilib.leedsim.widgets import BulkInput

from viperleed.guilib import decorators as dev_


DEFAULT_STARTUP = (
    {'SUPERLATTICE': ((1, 0), (0, 1)),
     'surfBasis': ((2, 0), (3*np.cos(np.radians(100)),
                            3*np.sin(np.radians(100)))),
     'surfGroup': 'p1', 'bulkGroup': 'p1', 'eMax': 700}
    )


class NewFileDialog(qtw.QDialog):
    """Dialog to create or edit the input for a LEED pattern."""

    editing_finished = qtc.pyqtSignal(LEEDParametersList)
    leed_parameters_changed = qtc.pyqtSignal(LEEDParametersList)

    def __init__(self, parent=None, old_params=None):
        """Initialize dialog.

        Parameters
        ----------
        parent : qtw.QWidget, default=None
            Parent widget
        old_params : viperleed.LEEDParametersList or None, optional
            LEED parameters used for initialization of the controls.
            They will also be used to restore the old pattern if
            the user presses the 'cancel' button. If not given or
            None, some default parameters will be used.

        Returns
        -------
        None.
        """
        super().__init__(parent)

        # Set the appearance of the window
        self.setWindowModality(qtc.Qt.WindowModal)
        flags = self.windowFlags()
        flags &= ~qtc.Qt.WindowCloseButtonHint  # Disable close button
        self.setWindowFlags(flags)
        self.setWindowTitle('Create/edit LEED input')

        old_params = LEEDParametersList(old_params)
        if old_params:
            self.__startup_parameters = old_params
        else:
            self.__startup_parameters = LEEDParametersList(DEFAULT_STARTUP)
        self.__lattices = []  # List of Lattice2D
        self._ctrls = {       # Control widgets
            'bulk': gl.BulkInput(),
            'surfaces': qtw.QTabWidget(),
            'e_max': qtw.QLineEdit(),
            'status_bar': qtw.QStatusBar(),
            'done': qtw.QPushButton('&Done'),
            'cancel': qtw.QPushButton('&Cancel'),
            'edit_surf_name': qtw.QLineEdit(self),
            'add_surface_tab': qtw.QToolButton(),
            'live_view': qtw.QCheckBox('Live update')
            }
        self._glob = {      # Miscellaneous attributes
            # Keep track of whether the widget has ever been shown
            '__never_shown': True,
            # 'surface_being_edited' is the index of the tab
            # whose name is currently being edited (or -1 if
            # no tab is being edited). Used in self.edit_surf_name
            # and self._on_surface_name_changed
            'surface_being_edited' : -1,
            # 'e_max_valid' contains information on whether
            # the input in the e_max control is a valid
            # maximum LEED energy or not.
            'e_max_valid': True,
            # Keep track of the last valid LEED parameters
            # to decide whether to fire leed_parameters_changed
            'last_leed_parameters': self.__startup_parameters,
            }

        self._compose()
        self._connect()

    @property
    def bulk_lattice(self):
        """Return the bulk lattice."""
        if self.__lattices:
            return self.__lattices[0]
        return None

    @property
    def is_live_view_on(self):
        """Return whether live-view mode is active."""
        return self._ctrls['live_view'].isChecked()

    def accept(self):
        """Reimplement QDialog.reject for some clean-up.

        Emits
        -----
        editing_finished
        """
        params = self.__pack_leed_parameters()
        self.editing_finished.emit(params)
        self._cleanup_before_close()
        super().accept()

    def keyPressEvent(self, event):  # pylint: disable=invalid-name
        """Reimplement keyPressEvent.

        The reimplementation fires Done or Cancel when they have
        keyboard focus and the user presses 'Enter' or 'Return'.
        Also, it catches a press to 'Esc' while editing one of
        the surface structure names (that would otherwise close
        the dialog).

        Parameters
        ----------
        event : QKeyPressEvent

        Returns
        -------
        None.
        """
        focus_widget = self.focusWidget()
        if event.key() in (qtc.Qt.Key_Enter, qtc.Qt.Key_Return):
            if focus_widget in (self._ctrls['cancel'], self._ctrls['done']):
                focus_widget.click()
                return
        if (event.key() == qtc.Qt.Key_Escape
                and focus_widget == self._ctrls['edit_surf_name']):
            self._glob['surface_being_edited'] = -1
            focus_widget.hide()
            return
        super().keyPressEvent(event)

    def open(self, params=None):
        """Open dialog.

        Parameters
        ----------
        params : LEEDParametersList or None, optional
            LEED parameters when the dialog is opened. Will be used for
            initializing controls, and to undo changes when 'Cancel' is
            pressed. If not given or None, the dialog is opened with
            default startup parameters. Default is None.

        Returns
        -------
        None.
        """
        if params:
            self.__startup_parameters = LEEDParametersList(params)
            self._glob['last_leed_parameters'] = self.__startup_parameters
        else:
            # Always have live-view mode on when no parameters
            # are given. This way we can force a plot of the
            # default lattices and LEED
            self._ctrls['live_view'].setChecked(qtc.Qt.Checked)

        self._init_lattices()
        self._init_controls()
        self.__fix_bulk_surf_spacing()

        # When no parameters are given, we should force an update.
        if not params:
            print("Open with no parameters")
            print("###     o--->", self.__class__.__name__,
                  "-- about to emit leed_parameters_changed")
            self.leed_parameters_changed.emit(self.__startup_parameters)
        super().open()

    def reject(self):
        """Reimplement QDialog.reject for some clean-up.

        Emits
        -----
        editing_finished
        """
        self.editing_finished.emit(self.__startup_parameters)
        self._cleanup_before_close()
        super().reject()

    @dev_.print_call
    def _add_surface_structure_tab(self, lattice=None, name=''):
        """Add one tab for a new lattice.

        Parameters
        ----------
        lattice : Lattice2D or None, optional
            The lattice to use to initialize the controls.
            If not given or None, a new 1x1 lattice is
            created and added to the list of lattices
        name : str, optional
            The string to be used as a name. If not given
            or '', the name is chosen from the number of
            tabs currently present. Default is ''.

        Returns
        -------
        None.
        """
        surf_tabs = self._ctrls['surfaces']

        if lattice is None:
            lattice = Lattice2D(basis=self.bulk_lattice.basis)
            self.__lattices.append(lattice)
        if not name:
            name = f"S{surf_tabs.count() + 1}"

        surf = gl.SurfaceStructureInput(bulk_lattice=self.bulk_lattice)
        surf.lattice = lattice
        surf.update_controls()

        # Disconnect some signals that need to be reconnected
        # after the ones of the surface are connected. This is
        # because, upon emission, the slots are called in the
        # order they were originally connected.
        bulk = self._ctrls['bulk']
        bulk.lattice_parameters_changed.disconnect(self._update_done_enabled)
        bulk.group_changed.disconnect(self._on_structures_changed)             # Live-view
        bulk.bulk_3d_operations_changed.disconnect(self._on_structures_changed)# Live-view
        bulk.shape_changed.disconnect(self._on_structures_changed)             # Live-view
        bulk.lattice_parameters_changed.disconnect(self._on_structures_changed)# Live-view

        # Connect signals to get the surface updated
        # when the bulk changes
        bulk.lattice_parameters_changed.connect(surf.on_bulk_basis_changed)
        bulk.group_changed.connect(surf.on_bulk_group_changed)
        bulk.shape_changed.connect(surf.update_woods_list_and_selection)

        # And those that signal user input errors or not
        surf.user_gave_invalid_input.connect(
            self._ctrls['status_bar'].showMessage
            )
        surf.surface_changed.connect(               # Valid edit
            self._ctrls['status_bar'].clearMessage
            )
        surf.user_gave_invalid_input.connect(self._update_done_enabled)
        surf.need_high_sym_reduction.connect(self._update_done_enabled)
        surf.surface_changed.connect(self._update_done_enabled)
        surf.surface_changed.connect(self._on_structures_changed)              # Live-view

        # Now reconnect the signals disconnected a few lines before
        bulk.lattice_parameters_changed.connect(self._update_done_enabled)
        bulk.group_changed.connect(self._on_structures_changed)                # Live-view
        bulk.bulk_3d_operations_changed.connect(self._on_structures_changed)   # Live-view
        bulk.shape_changed.connect(self._on_structures_changed)                # Live-view
        bulk.lattice_parameters_changed.connect(self._on_structures_changed)   # Live-view

        surf_tabs.insertTab(surf_tabs.count(), surf, name)

        # Select the new added tab.
        surf_tabs.setCurrentIndex(surf_tabs.count() - 1)

        # And update some controls and the live view
        self._update_done_enabled()
        self._on_structures_changed()                                          # Live-view

    def _cleanup_before_close(self):
        """Clean up the dialog before closing.

        In essence, removes all the tabs. Otherwise at the
        next call to .open() we would still have the previous
        tabs present.

        Returns
        -------
        None.
        """
        while self._ctrls['surfaces'].count():
            self._ctrls['surfaces'].removeTab(0)
        self._glob['__never_shown'] = True

    def _compose(self):
        """Place children widgets."""
        # (1) Initialize the widgets
        # (1.1) surfaces tab widget
        self._ctrls['surfaces'].setFont(gl.AllGUIFonts().smallTextFont)
        self._ctrls['surfaces'].setTabsClosable(True)

        # Change default behavior, such that we never
        # select the 'add' tab. This would otherwise
        # happen when removing the 'last' structure
        self._ctrls['surfaces'].tabBar().setSelectionBehaviorOnRemove(
            self._ctrls['surfaces'].tabBar().SelectLeftTab
            )

        self._ctrls['add_surface_tab'].setDefaultAction(
            qtw.QAction(
                qtg.QIcon(
                    os.path.join(resources_path('guilib/icons'),
                                 'plus_button_48x48.png')
                    ),
                "Add structural domain"
                )
            )

        self._ctrls['surfaces'].setCornerWidget(
            self._ctrls['add_surface_tab'],
            qtc.Qt.TopLeftCorner
            )

        # (1.2) Maximum LEED energy
        emax_label = qtw.QLabel('Max. Energy:')
        for ctrl in [emax_label, self._ctrls['e_max']]:
            ctrl.setSizePolicy(qtw.QSizePolicy.Fixed,
                               qtw.QSizePolicy.Preferred)
            ctrl.setFont(gl.AllGUIFonts().labelFont)
        self._ctrls['e_max'].setMaximumWidth(70)
        to_be_polished = [emax_label, self._ctrls['e_max']]

        # (1.3) Live-view mode check-box
        self._ctrls['live_view'].setChecked(qtc.Qt.Checked)
        self._ctrls['live_view'].setFont(gl.AllGUIFonts().labelFont)
        to_be_polished.append(self._ctrls['live_view'])

        # (1.4) 'Done' and 'Cancel' buttons
        for key in ('done', 'cancel'):
            self._ctrls[key].setFont(gl.AllGUIFonts().buttonFont)
            self._ctrls[key].setSizePolicy(qtw.QSizePolicy.Fixed,
                                           qtw.QSizePolicy.Fixed)
            to_be_polished.append(self._ctrls[key])

        # (1.5) Status bar
        self._ctrls['status_bar'].setSizeGripEnabled(False)

        # (2) Make sure all styles are up to date
        for widg in to_be_polished:
            widg.ensurePolished()

        # (3) Place controls into layouts
        emax_layout = qtw.QHBoxLayout()
        emax_layout.setSpacing(4)
        emax_layout.setContentsMargins(5, 5, 5, 5)
        emax_layout.addWidget(emax_label)
        emax_layout.addWidget(self._ctrls['e_max'])
        emax_layout.setAlignment(emax_label,
                                 qtc.Qt.AlignRight | qtc.Qt.AlignVCenter)

        status_and_buttons_layout = qtw.QHBoxLayout()
        for key in ('status_bar', 'done', 'cancel'):
            status_and_buttons_layout.addWidget(self._ctrls[key])

        # (3.1) Assemble the layout of the QDialog
        dialog_layout = qtw.QGridLayout()
        # Leave free one row above and one below 'bulk', such that
        # they can be used in _init_controls to align appropriately
        # the lattice parameter controls for bulk and surface.
        dialog_layout.addWidget(self._ctrls['bulk'], 1, 0, 1, 1)
        dialog_layout.addWidget(self._ctrls['surfaces'], 0, 1, 3, 1)
        dialog_layout.addWidget(self._ctrls['live_view'], 3, 0, 1, 1)
        dialog_layout.addLayout(emax_layout, 3, 1, 1, 1)
        dialog_layout.addLayout(status_and_buttons_layout, 4, 0, 1, 2)
        dialog_layout.setAlignment(emax_layout,
                                   qtc.Qt.AlignRight | qtc.Qt.AlignVCenter)
        self.setLayout(dialog_layout)

        # (4) Take care of widgets that do not belong to a layout
        self._ctrls['edit_surf_name'].setFont(gl.AllGUIFonts().smallTextFont)
        self._ctrls['edit_surf_name'].setAlignment(qtc.Qt.AlignCenter)
        self._ctrls['edit_surf_name'].hide()

    def _connect(self):
        """Connect relevant signals to handlers."""
        for key in ('done', 'cancel'):
            button = self._ctrls[key]
            # The next two lines prevent the dialog to be closed
            # automatically by firing one of the buttons if the
            # user presses enter/return while focus is on widgets
            # other than the buttons.  The reimplementation of
            # keyPressEvent takes care of the behavior when the
            # user presses 'Enter'/'Return' while they have focus
            button.setAutoDefault(False)
            button.setDefault(False)

        self._ctrls['done'].clicked.connect(self.accept)
        self._ctrls['cancel'].clicked.connect(self.reject)

        self._ctrls['e_max'].editingFinished.connect(self._on_energy_changed)

        # Handle the surface-structure tabs
        self._ctrls['surfaces'].tabCloseRequested.connect(
            self._on_surface_tab_closed
            )
        self._ctrls['surfaces'].tabBarDoubleClicked.connect(
            self._edit_surface_name
            )
        self._ctrls['add_surface_tab'].clicked.connect(
            self._on_add_structure_pressed
            )
        self._ctrls['edit_surf_name'].editingFinished.connect(
            self._on_surface_name_changed
            )

        # Connect signals that cause updates of the surfaces
        # when the bulk lattice is changed
        bulk = self._ctrls['bulk']
        bulk.shape_changed.connect(self._on_bulk_shape_changed)
        bulk.high_sym_pressed.connect(self._on_bulk_high_sym_pressed)

        # And signals that give feedback to the user about invalid input
        bulk.need_high_sym_reduction.connect(self._update_done_enabled)
        bulk.lattice_parameters_changed.connect(self._update_done_enabled)

        # Finally stuff for live-view
        bulk.group_changed.connect(self._on_structures_changed)
        bulk.bulk_3d_operations_changed.connect(self._on_structures_changed)
        bulk.shape_changed.connect(self._on_structures_changed)
        bulk.lattice_parameters_changed.connect(self._on_structures_changed)
        self._ctrls['live_view'].stateChanged.connect(
            self._on_structures_changed
            )

    @dev_.print_call
    def _edit_surface_name(self, index):
        """Show input for editing name of surface structure.

        This is the slot connected to the tabBarDoubleClicked
        signal of the tab widget for the surface structure.

        Pop up a QLineEdit where the new name can be entered.

        Parameters
        ----------
        index : int
            The index of the tab that was double-clicked

        Returns
        -------
        None.
        """
        edit = self._ctrls['edit_surf_name']
        tabs = self._ctrls['surfaces']

        edit.show()
        edit.setText(tabs.tabText(index))
        edit.resize(tabs.tabBar().tabRect(index).size())
        edit.move(tabs.x() + tabs.tabBar().tabRect(index).x(), tabs.y())
        edit.setFocus()
        edit.selectAll()
        edit.raise_()  # Move on top of all widgets

        self._glob['surface_being_edited'] = index

    def __fix_bulk_surf_spacing(self):
        """Fix vertical offset between bulk and surface controls.

        The vertical offset between the lattice-parameter
        controls of bulk and surface is there because the
        latter is inside the tab-bar widget. We add some
        space above the bulk to have the controls line-up.

        The fix is done only once.
        """
        if self._glob['__never_shown']:
            # To get the correct positions of widgets there
            # is no way around showing and hiding the dialog
            # This should be almost invisible though.
            self.show()
            self._ctrls['bulk'].flash_high_symm_button()
            self.hide()
            self._glob['__never_shown'] = False

        # Now positions are correct. Calculate the current
        # vertical offset between the bulk and surface controls
        surf = self._ctrls['surfaces'].currentWidget()
        offset = surf.top_left_global_handle
        offset -= self._ctrls['bulk'].top_left_global_handle
        offset = offset.y()

        # Finally set the space as the height of the 1st row
        top_row_height = self.layout().rowMinimumHeight(0)
        top_row_height += offset
        top_row_height = max(0, top_row_height)
        self.layout().setRowMinimumHeight(0, top_row_height)

    def _init_controls(self):
        """Set the values in the controls from startup parameters."""
        # Maximum LEED energy
        e_max = self.__startup_parameters[0]['eMax']
        self._ctrls['e_max'].setText(f"{e_max} eV")

        # Bulk lattice
        self._ctrls['bulk'].lattice = self.bulk_lattice

        # Surface lattices
        for i, lattice in enumerate(self.__lattices[1:]):
            # Add one tab per surface lattice
            name = self.__startup_parameters[i].get('name', '')
            if not name:
                name = f"S{i + 1}"
            self._add_surface_structure_tab(lattice, name)

        self._ctrls['surfaces'].setCurrentIndex(0)

    def _init_lattices(self):
        """Make lattices from the parameters.

        The lattices are very minimal and
        should not be used for plotting.
        """
        superlattice = self.__startup_parameters[0]['SUPERLATTICE']
        bulk_basis = np.dot(np.linalg.inv(superlattice),
                            self.__startup_parameters[0]['surfBasis'])
        self.__lattices = [Lattice2D(
            bulk_basis,
            group=self.__startup_parameters[0]['bulkGroup']
            )]
        for param in self.__startup_parameters:
            self.__lattices.append(Lattice2D(
                param['surfBasis'],
                group=param['surfGroup']
                ))

    @dev_.print_call
    def _on_add_structure_pressed(self, __action=None):
        """React to a user requesting a new structure to be added.

        This is the slot connected to the '+' clicked signal
        of the button in the 'add surface structure' tab.

        Parameters
        ----------
        __action : QAction
            The action that has been triggered. Unused.
        """
        self._add_surface_structure_tab()

    @dev_.print_call
    def _on_bulk_high_sym_pressed(self):
        """React to a request to reduce the bulk to high symmetry.

        Opens an extra pop-up to choose whether the surface
        lattices should be correspondingly changed by either
        keeping the superlattice matrix constant (default)
        or keeping the lattice parameters constant.

        Returns
        -------
        None.
        """
        diag = BulkHighSymReductionDialog()
        error_txt = ("Keeping constant lattice parameters for structure {} "
                     "makes the lattice incommensurate. Can only transform "
                     "with constant superlattice matrix.")
        error_msg = qtw.QMessageBox(qtw.QMessageBox.Critical,
                                    "Cannot keep constant lattice parameters",
                                    error_txt)
        surf_tabs = self._ctrls['surfaces']

        if diag.exec_() == diag.keep_constant_parameters:
            # Since the default behavior is to keep the superlattice
            # matrix constant, we rather try to set the superlattice
            # matrix to the correct values. Since Bs = M Bb, when Bb
            # is transformed to Bb' = T Bb, to keep the same Bs we
            # need Bs = M Bb = M * T^-1 Bb' = M' Bb'
            bulk_transform = self.bulk_lattice.high_symm_transform()
            bulk_transform_inv = np.linalg.inv(bulk_transform)

            for i in range(surf_tabs.count()):
                surf_i = surf_tabs.widget(i)
                superlattice = np.dot(surf_i.superlattice, bulk_transform_inv)  # TODO: use ensure_integer_matrix
                if np.any(np.abs(superlattice - superlattice.round()) > 1e-3):
                    # Lattice would be incommensurate. Signal the problem
                    error_msg.setText(error_txt.format(surf_tabs.tabText(i)))
                    error_msg.exec_()
                    continue
                # Matrix is ok. Set it and update the Wood's notation
                surf_i.superlattice = superlattice.round().astype(int)
                surf_i.woods.bulk_basis = np.dot(bulk_transform,
                                                 self.bulk_lattice.basis)
                surf_i.pick_right_woods()

    @dev_.print_call
    def _on_bulk_shape_changed(self, __new_shape=None):
        """React to a change of the shape of the bulk lattice.

        This is the slot connected to the shape_changed
        signal of the bulk control.

        Parameters
        ----------
        __new_shape : str
            Unused, but necessary due to the signature
            of the shape_changed signal.

        Return
        ------
        None.
        """
        # This is not ideal, as shape_changed may be already
        # followed by a lattice_parameters_changed signal.
        # This should not hurt on the final outcome, but it
        # may process twice the same information.
        print("###     o--->", self.__class__.__name__,
              "-- about to emit bulk.lattice_parameters_changed")
        self._ctrls['bulk'].lattice_parameters_changed.emit(
            self.bulk_lattice.basis
            )

    @dev_.print_call
    def _on_energy_changed(self):
        """React to a change of energy.

        This is the handler of the editingFinished signal
        of the 'maximum LEED energy' control.

        For now, it checks that the text entered is a float,
        formats it with one decimal digit, and adds an 'eV'
        unit. Also resizes the widget to fit possibly longer
        text.

        Parameters
        ----------
        new_energy : str
            The new value of energy.

        Returns
        -------
        None.
        """
        energy = self._ctrls['e_max'].text()
        valid = True
        if 'eV' in energy:
            energy = energy.replace('eV', '')
        try:
            energy = float(energy)
        except ValueError:
            # Not a valid float
            self._ctrls['status_bar'].showMessage("Maximum LEED energy should "
                                                  "be a number", 2000)
            valid = False

        if valid and energy <= 0:
            self._ctrls['status_bar'].showMessage("Maximum LEED energy should "
                                                  "be positive", 2000)
            valid = False

        self._glob['e_max_valid'] = valid
        change_control_text_color(self._ctrls['e_max'],
                                  qtc.Qt.black if valid else qtc.Qt.red)
        self._update_done_enabled()

        if not valid:
            return

        # Prepare the new text, and get its size
        energy = f"{energy:.1f} eV"
        font_metrics = self._ctrls['e_max'].fontMetrics()
        # Energies up to 4 digits (999.9 eV) fit into width = 57,
        # but it looks nicer if there is a little space to the right
        width = max(font_metrics.width(energy) + 13, 70)

        # Set the new text, and resize the widget
        self._ctrls['e_max'].setText(energy)
        self._ctrls['e_max'].setMaximumWidth(width)
        self._ctrls['e_max'].adjustSize()

    @dev_.print_call
    def _on_structures_changed(self, *__args, **__kwargs):
        """React on any user change of the structures.

        The leed_parameters_changed signal emitted can be
        used for live updates.

        Returns
        -------
        None

        Emits
        -----
        leed_parameters_changed
            Emitted only if the current input is valid, the
            parameters have actually changed, and live-view
            mode is on.
        """
        # TODO: evaluate whether changing the emax makes things
        # too slow to recalculate. In that case, I may choose to
        # fix the energy to some low value here to limit the number
        # of beams that need to be calculated.
        if not self._ctrls['done'].isEnabled():
            return

        new_params = self.__pack_leed_parameters()
        if not new_params or new_params == self._glob['last_leed_parameters']:
            print("Parameters unchanged")
            return

        print("old:", self._glob['last_leed_parameters'])
        print("new:", new_params)
        self._glob['last_leed_parameters'] = new_params

        if not self.is_live_view_on:
            return
        print("###     o--->", self.__class__.__name__,
              "-- about to emit leed_parameters_changed")
        self.leed_parameters_changed.emit(new_params)

    @dev_.print_call
    def _on_surface_name_changed(self):
        """React on completed edits of a surface-structure name.

        This is the slot to which the editingFinished signal of
        the QLineEdit for editing a surface name is connected.
        """
        new_name = self._ctrls['edit_surf_name'].text()
        tab_names = [self._ctrls['surfaces'].tabText(i)
                     for i in range(self._ctrls['surfaces'].count())]
        if self._glob['surface_being_edited'] != -1:
            # A name is being edited. Make sure the new name
            # does not conflict with others, except the one
            # of the edited tab.
            if any(new_name == tab_name
                   for i, tab_name in enumerate(tab_names)
                   if i != self._glob['surface_being_edited']):
                self._ctrls['status_bar'].showMessage(
                    f"Surface structure name {new_name!r} already used.",
                    3000
                    )
                return
            self._ctrls['surfaces'].setTabText(
                self._glob['surface_being_edited'],
                new_name
                )
            self._glob['surface_being_edited'] = -1
        self._ctrls['edit_surf_name'].hide()

    @dev_.print_call
    def _on_surface_tab_closed(self, index):
        """React on a request to remove a surface structure.

        This is the slot connected to the tabCloseRequested
        signal of the QTabWidget of the surface structures.

        Parameters
        ----------
        index : int
            The index of the tab that the user wants to be closed

        Returns
        -------
        None.
        """
        tabs = self._ctrls['surfaces']

        # There should always be at least one surface structure tab.
        # The new tab that will appear will just be a 1x1 surface
        # structure, i.e., a standard lattice.
        n_structures = tabs.count()
        if n_structures == 1:
            lattice = Lattice2D(self.bulk_lattice.basis)
            self.__lattices[1] = lattice
            tabs.removeTab(index)
            self._add_surface_structure_tab(lattice, "S1")
            return

        # When there are multiple structures, we can remove the
        # tab right away, but we should also rename any of the
        # tabs with the generic name f"S{i + 1}" following it,
        # since their position will be changed
        for i in range(index + 1, n_structures):
            name = tabs.tabText(i)
            if name == f"S{i + 1}":
                tabs.setTabText(i, f"S{i}")
        self.__lattices.pop(index + 1)     # 0 is the bulk
        tabs.removeTab(index)

    @dev_.print_call
    def __pack_leed_parameters(self):
        """Return a LEEDParametersList for the current structures.

        This function should be called only when the 'Done' button
        is enabled.

        Returns
        -------
        parameters : LEEDParametersList
            The LEED parameters corresponding
            to the current user input.
        """
        # Prepare the stuff that is common to all structures
        e_max = float(self._ctrls['e_max'].text().replace('eV', ''))
        bulk_group = self.bulk_lattice.group
        bulk_3d_sym = bulk_group.screws_glides

        # And get the parameters that you cannot change with
        # this dialog from those we got at the beginning
        beam_incidence = self.__startup_parameters[0]['beamIncidence']
        screen_aperture = self.__startup_parameters[0]['screenAperture']

        parameters = LEEDParametersList()
        for i in range(self._ctrls['surfaces'].count()):
            surf = self._ctrls['surfaces'].widget(i)
            parameters.append(
                {'eMax': e_max,
                 'surfBasis': surf.lattice.basis,
                 'SUPERLATTICE': surf.superlattice,
                 'surfGroup': surf.lattice.group,
                 'bulkGroup': bulk_group,
                 'bulk3Dsym': bulk_3d_sym,
                 'beamIncidence': beam_incidence,
                 'screenAperture': screen_aperture,
                 'name': self._ctrls['surfaces'].tabText(i)}
                )

        return parameters

    @dev_.print_call
    def _update_done_enabled(self, *__args):
        """Update enabled/disabled state of Done button.

        The button will be disabled if the user gave any
        invalid input. Also, the color of tabs with valid
        or invalid input is updated.

        Returns
        -------
        None.
        """
        enabled = self._ctrls['bulk'].valid_input
        enabled &= self._glob['e_max_valid']

        for i in range(self._ctrls['surfaces'].count()):
            surf = self._ctrls['surfaces'].widget(i)
            valid = surf.valid_input
            self._ctrls['surfaces'].tabBar().setTabTextColor(
                i,
                qtc.Qt.black if valid else qtc.Qt.red
                )
            enabled &= valid
        self._ctrls['done'].setEnabled(enabled)


class BulkHighSymReductionDialog(qtw.QDialog):
    """Pick how surfaces are transformed when symmetrizing the bulk."""

    keep_constant_parameters = qtw.QDialog.Accepted
    keep_constant_superlattice = qtw.QDialog.Accepted

    def __init__(self, parent=None):
        """Initialize dialog."""
        super().__init__(parent)

        # Set the appearance of the window
        self.setWindowModality(qtc.Qt.WindowModal)
        flags = self.windowFlags()
        flags &= ~qtc.Qt.WindowCloseButtonHint  # Disable close button
        self.setWindowFlags(flags)
        self.setWindowTitle('Make bulk lattice high symmetry')

        self.__compose_and_connect()

    def __compose_and_connect(self):
        """Place children widgets."""
        buttons = [qtw.QPushButton("Constant\nlattice\nparameters"),
                   qtw.QPushButton("Constant\nsuperlattice\nmatrix")]
        btn_layout = qtw.QHBoxLayout()
        for button in buttons:
            button.setFont(gl.AllGUIFonts().buttonFont)
            button.ensurePolished()
            btn_layout.addWidget(button)
        buttons[0].clicked.connect(self.accept)
        buttons[1].clicked.connect(self.reject)

        txt = qtw.QLabel("Pick how to transform surface structures")
        txt.setFont(gl.AllGUIFonts().labelFont)

        layout = qtw.QVBoxLayout()
        layout.addWidget(txt)
        layout.addLayout(btn_layout)

        self.setLayout(layout)
