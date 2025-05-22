"""Module mainwindow of viperleed.guilib.leedsim

======================================
  ViPErLEED Graphical User Interface
======================================
    *** module guilib.leedsim ***

Module: mainwindow.py

Created: 2020-01-13
Author: Michele Riva

Defines the LEEDPatternSimulator class, i.e., the main window of
the ViPErLEED pattern simulator plug-in. Allows simulating LEED
patterns from multiple structures, accounting for symmetry.

2021-07-04: major restyling started, including
            support for multiple structures
"""

import re
import copy
import warnings
from pathlib import Path

import numpy as np
import PyQt5.QtCore as qtc
import PyQt5.QtWidgets as qtw

from viperleed import guilib as gl
from viperleed.guilib.widgetdecorators import broadcast_mouse


DEFAULT_INPUT_FILE_EXTENSION = '*.tlm'
DEFAULT_OPEN_FILE = ('./guilib/leedsim/input examples/', 'PatternInfo.tlm')
DEFAULT_EXPORT_FILE = ('./guilib/leedsim/exported/', 'LEEDBeams.csv')
TITLE = 'LEED Pattern Simulator'
SCREEN_FRACTION = 0.63


def default_input_file_extensions():
    """Return a string that can be used in Open/Save dialogs for filtering."""
    txt = ("LEED Pattern simulation input files ("
           f"{DEFAULT_INPUT_FILE_EXTENSION}, "
           f"{DEFAULT_INPUT_FILE_EXTENSION.upper()})")
    return f"{txt};;All files (*)"


# TODO: this may easily go into a settings.ini file
def default_file_menu():
    """Return settings for the 'File' menu.

    Returns
    -------
    settings : dict
        keys : str
            Shorthand name of menu entry. Typical names are
            'new', 'open', 'save', 'save_as', 'exit'
        values : tuple
            Settings for the menu entry, in the following order:
            icon : QIcon or None
                Icon to be used in the menu and tool bar. If
                None, the corresponding entry will not be added
                to the tool bar.
            menu_text : str
                Text appearing in the menu entry
            shortcut : str
                Keyboard shortcut that activates the action
            tool_tip : str
                Text that appears in a tool tip box when
                hovering with the mouse. Shortcut will be
                added in parentheses.
            status_tip : str
                Text that appears in the status bar when
                hovering with the mouse. Shortcut will be
                added in parentheses.
            handler : str
                Attribute of the class that will be invoked
                when the menu entry/tool-bar icon is activated
            to_be_enabled : bool
                True if the entry should be enabled only when
                a valid input file is loaded. If False, the
                entry is always enabled.
    """
    style = qtw.QWidget().style()
    settings = {
        'new': (style.standardIcon(qtw.QStyle.SP_FileIcon),
                '&New/edit...',
                'Ctrl+N',
                'New/Edit LEED pattern',
                'New/Edit LEED pattern input file',
                '_on_file_new_or_edit_pressed',
                False),
        'open': (style.standardIcon(qtw.QStyle.SP_DialogOpenButton),
                 '&Open...',
                 'Ctrl+O',
                 'Open LEED pattern input',
                 'Open LEED pattern input file',
                 '_on_file_open_pressed',
                 False),
        'save': (style.standardIcon(qtw.QStyle.SP_DialogSaveButton),
                 '&Save',
                 'Ctrl+S',
                 'Save LEED pattern input',
                 'Save LEED pattern input file',
                 '_on_file_save_pressed',
                 True),
        'save_as': (None,
                    '&Save As...',
                    'Ctrl+Shift+S',
                    'Save LEED pattern input as',
                    'Save LEED pattern input file as',
                    '_on_file_save_as_pressed',
                    True),
        'exit': (None,
                 '&Exit',
                 'Ctrl+W',
                 'Exit LEED Pattern simulator',
                 'Exit LEED Pattern simulator plug-in',
                 'close',
                 False)
    }
    return settings


def show_use_betatest_version_popup():
    """Show a pop-up dialog hinting at using the betatest version."""
    _betatest = 'https://github.com/viperleed/viperleed-betatest'
    txt = (
        f'The ViPErLEED graphical user interface for v{gl.GLOBALS["version"]} '
        'is currently under development.<p>'
        'Please use the pre-packed version available from the <a href='
        f'{_betatest}/releases/latest>viperleed-betatest'
        '</a> GitHub repository (file gui.zip).</p>'
        )
    msgBox = qtw.QMessageBox(qtw.QMessageBox.Information,
                             'Use betatest version',
                             txt)
    msgBox.exec_()


# TODO: remember last directories, also across sessions

@broadcast_mouse
class LEEDPatternSimulator(gl.ViPErLEEDPluginBase):
    """A class that allows simulating LEED Patterns."""

    def __init__(self, parent=None):
        """Initialize window."""
        super().__init__(parent, name=TITLE)

        # Keep references to controls, dialogs, and some globals
        self._ctrls = {
            # We have to keep a reference to the top-level
            # menus in the menu bar to have them show up.
            'file_menu': qtw.QMenu("&File"),
            # Real-space and LEED pattern plotting canvases
            'lattices_canvas': gl.RealCanvas(title='Real-Space Lattices'),      # TODO: was self.realSpace
            'leed_canvas': gl.LEEDCanvas(title='LEED Pattern'),                 # TODO: was self.recSpace
            # Rotation, just as a view, i.e., not the azimuthal
            # angle of incidence of the primary beam
            'rotation': gl.RotationBlock(),                                     # TODO: was self.rotWidg
            # Current primary electron energy
            'energy': gl.EnergyBlock(),                                         # TODO: was self.enWidg
            # Domains selector
            'domains': qtw.QPushButton('Select domains')
            }
        # __enabled_on_valid_input is a list of actions
        # and controls that are enabled only when the user
        # gave a valid input for a surface structure.
        self.__enabled_on_valid_input = [
            v for k, v in self._ctrls.items()
            if k not in ('file_menu',)
            ]

        self._dialogs = {'file_new': gl.NewFileDialog(self),}
        self._glob = {
            # parameters holds the LEED parameters currently
            # used for plotting the lattices and LEED pattern
            'parameters': gl.LEEDParametersList(),
            # is_saved keeps track of whether the structure
            # input is saved to disk. Used to keep track of
            # whether the saved file is up to date
            'is_saved': False,
            # filename keeps track of the name of the
            # input file that is currently open, if any.
            # Used to save edits and for exporting.
            'filename': None,
            }

        # Set window properties
        self.setWindowTitle(TITLE)
        self.setAcceptDrops(True)

        self.__compose()
        self.__connect()
        self.adjustSize()
        self.center_on_screen()

    @property
    def view_rotation(self):
        """Return the current viewing rotation as a float.

        Returns
        -------
        angle : float
            The current angle in degrees between the horizontal
            direction and the first basis vector of the bulk,
            modulo 360.                                                         # TODO: true??
        """
        return float(self._ctrls['rotation'].text.text())

    @property
    def energy(self):
        """Return the current energy as a float."""
        return float(self._ctrls['energy'].text.text())

    @property
    def filename(self):
        """Return the name of the open file, if any, None otherwise."""
        return self._glob['filename']

    @filename.setter
    def filename(self, new_filename):
        """Set the name of the currently open file."""
        self._glob['filename'] = new_filename

    def __get_parameters(self):
        """Return the current LEED parameters used.

        This is the getter for the .parameters @property.

        Returns
        -------
        parameters : LEEDParametersList
            The current LEED parameters used for plotting
        """
        return self._glob['parameters']

    def __set_parameters(self, new_parameters):
        """Set new LEED parameters.

        This is the setter for the .parameters @property.

        Parameters
        ----------
        new_parameters : Sequence
            The new LEED parameters used for plotting. Each element
            should be a dict, a ConfigParser, a LEEDParameters, or
            a LEEDParametersList. Also a single LEEDParametersList
            is an acceptable Sequence.
        """
        if new_parameters == self.parameters:
            # No need to recalculate if things are unchanged
            print("same")
            return
        self._glob['parameters'] = gl.LEEDParametersList(new_parameters)
        self.enable_ctrls(True)
        # self.initRealAndLEED(True)                                            # TODO: call appropriate method after reimplementation

        #                                                                       # TODO: was also destroying the exportDialog, but would probably
        #                                                                       # simply be better to reprepare the dialog and its contents each
        #                                                                       # time it is opened

    parameters = property(__get_parameters, __set_parameters)

    @property
    def saved(self):
        """Return whether the current input is saved."""
        return self._glob['is_saved']

    @saved.setter
    def saved(self, is_saved):
        """Set whether the current input is saved."""
        #                                                                       # TODO: also change active state of save button and menu actions
        self._glob['is_saved'] = bool(is_saved)

    def closeEvent(self, event):
        """Extend closeEvent to ask for saving changes."""
        reply = qtw.QMessageBox.No

        # self.filename is None only if no file was ever loaded
        # and if no edit was ever done when starting from scratch
        if self.filename is not None and not self.saved:
            reply = qtw.QMessageBox.question(
                self, 'Unsaved changes',
                f"Save changes to file {self.filename} before closing?",
                qtw.QMessageBox.Yes
                | qtw.QMessageBox.No
                | qtw.QMessageBox.Cancel,
                qtw.QMessageBox.Yes
                )
        if reply == qtw.QMessageBox.Cancel:
            event.ignore()
            return
        if reply == qtw.QMessageBox.Yes:
            self._on_file_save_pressed()
        super().closeEvent(event)

    def dragEnterEvent(self, event):
        """Extend dragEnterEvent to accept drag-drop of input files."""
        mime_data = event.mimeData()

        # Use 'urls' as this returns an iterable of QUrl objects
        # that point to the file paths of the files that are dragged
        if not mime_data.hasUrls():
            # Not dropping a file
            super().dragEnterEvent(event)
            return

        # Determine whether the files dropped are acceptable inputs
        fnames = [url.toLocalFile() for url in mime_data.urls()]
        params = self.__read_files_and_report_errors(fnames, silent=True)
        if not params:
            # Invalid input
            super().dragEnterEvent(event)
            return

        # Files are OK. Enable receiving QDropEvent
        event.acceptProposedAction()

    def dropEvent(self, event):
        """Extend dropEvent to accept drag-drop of input files.

        The reimplementation of dragEnterEvent guarantees that
        event.mimeData().hasUrls(), and that the files being
        dropped are acceptable.

        Parameters
        ----------
        event : QDropEvent
            The drop event to be processed.
        """
        fnames = [url.toLocalFile() for url in event.mimeData().urls()]
        params = self.__read_files_and_report_errors(fnames)
        if not params:
            # This should never happen, but better safe than sorry
            return

        event.acceptProposedAction()
        if params:
            self.parameters = params

        if len(fnames) > 1:
            self.filename = ''
            self.saved = False
        else:
            self.filename = fnames[0]
            self.saved = True

    def enable_ctrls(self, enabled):
        """Enable or disable controls and actions.

        This method should typically be called when a valid
        structural input is loaded.

        Parameters
        ----------
        enabled : bool
            True to enable controls. In fact, any object that
            can evaluate to a bool is accepted. Its truth value
            will be used.
        """
        enabled = bool(enabled)
        for obj in self.__enabled_on_valid_input:
            obj.setEnabled(enabled)

    def mouseReleaseEvent(self, event):                                         # TODO: fixme
        """Extend mouseReleaseEvent.

        Hide the pop-up with matrices when clicking anywhere
        on the window.
        TODO: this will rather close the domain selector window

        Parameters
        ----------
        event : QMouseEvent
            The mouse-release event to be handled.

        Returns
        -------
        None
        """
        if False:
            doms = self._ctrls['domains']
            if (doms.isEnabled()
                and doms.matricesPopup.isVisible()
                and not doms.text.underMouse()):
                # the reason for not acting when doms.text is underMouse() is
                # that the event handling in this case is done within DomsBlock
                doms.matricesPopup.hide()
            # TESTING                                                           # TODO: fix
            # for annot in self.recSpace.annots:
                # annot.show()
            # END TESTING
        super().mouseReleaseEvent(event)

    def moveEvent(self, event):
        """Extend moveEvent to have non-children follow the window."""
        super().moveEvent(event)
        delta = event.pos() - event.oldPos()

        for annot in self._ctrls['leed_canvas'].annots:
            # NB: one cannot use += for in-place assignment since
            # .center is a property, which would then become a
            # standard attribute.
            annot.center = qtc.QPoint(*annot.center) + delta
            if annot.isVisible():
                # Move only visible annotations, which will move
                # together with the window. If one would have
                # this call even for invisible annotations, the move
                # would be issued before the next show. This messes
                # the position of the annotation.
                annot.move(annot.pos() + delta)

    def save_input(self):
        """Save the current input to file."""
        gl.LEEDParser.write_structures(self.parameters, self.filename)

    def sizeHint(self):                                                         # TODO: consider if it makes sense to: (1) also reimplement minimumSizeHint(); (2) move this to base
        """Reimplement sizeHint to scale with the current screen."""
        self.ensurePolished()
        return SCREEN_FRACTION*self.screen().size()                             # TODO: for MPLCanvas it may make sense to look at QSize.scale(..., aspect_ratio)

    def _on_export_beams_pressed(self):                                         # TODO: fix + doc
        # Create a dialog to select which domains should be exported
        # and to process the data to export
        if not hasattr(self, 'exportDialog'):
            self.exportDialog = gl.ExportCSVDialog(self.leed,
                                                   parent=self)
            self.exportDialog.exportSelected.connect(self.exportCSV)
        else:
            self.exportDialog.open()

    def _on_file_new_or_edit_pressed(self):
        """React to a request to produce a new input file or edit one."""
        self._dialogs['file_new'].open(self.parameters)

    def _on_file_open_pressed(self):
        """React to a request to open existing input files."""
        fnames = qtw.QFileDialog.getOpenFileNames(
            self, 'Open LEED pattern input',
            DEFAULT_OPEN_FILE[0],
            default_input_file_extensions()
            )
        if not fnames[0]:
            # No files selected
            return

        params = self.__read_files_and_report_errors(fnames[0], silent=False)
        if not params:
            # Files selected were invalid
            return

        # When multiple files are loaded, we force the
        # user to save the combined input to a new file
        if len(fnames[0]) > 1:
            self.filename = ''
            self.saved = False
        else:
            self.saved = True
            self.filename = fnames[0][0]

        self.parameters = params

    def _on_file_save_as_pressed(self):
        """React to a request to save the input with a new name."""
        fname = qtw.QFileDialog.getSaveFileName(
            self, 'Save Pattern simulator input',
            '', default_input_file_extensions()
            )
        if fname and fname[0]:
            # File selected correctly
            self.filename = fname[0]
            self.save_input()
            self.saved = True

    def _on_file_save_pressed(self):
        """React to a request to save the input file."""
        if not self.filename:  # No file is currently open
            self._on_file_save_as_pressed()
            return
        #                                                                       # TODO: add here a message box asking if you want to overwrite
        self.save_input()
        self.saved = True

    def _on_screen_changed(self, old_screen, new_screen):
        """React to a change of screen.

        This is the slot connected to the screen_changed
        signal, inherited from ViPErLEEDPluginBase.

        Parameters
        ----------
        old_screen : QScreen
            The old screen of the window
        new_screen : QScreen
            The new screen of the window
        """
        print(old_screen, new_screen)
        self.adjustSize()

    def _on_structure_edit_finished(self, leed_parameters):
        """React to a change of the structure input.

        This is the slot connected to the editing_finished of
        NewFileDialog. The sole purpose of this is to have the
        GUI appropriately ask for saving changes to file when
        the first structure is created from scratch.

        Parameters
        ----------
        leed_parameters : LEEDParametersList
            The new LEEDParametersList to be used
        """
        if leed_parameters and self.filename is None:
            print("edited")
            self.filename = ''
        self.parameters = leed_parameters
        self.saved = False

    def __compose(self):
        """Prepare menu, bars, and children widgets."""
        self.__compose_menu_and_toolbar()

        self.setStatusBar(qtw.QStatusBar())
        self.statusBar().showMessage('Ready')

        self.setCentralWidget(qtw.QWidget())
        self.centralWidget().setLayout(qtw.QGridLayout())

        self.__compose_ctrls()
        self.enable_ctrls(False)

    def __compose_ctrls(self):
        """Set up and place children widgets."""
        # (1) Real-space view of the lattices:
        #     Scrolling while the mouse cursor is on the canvas
        #     causes visual rotation of the lattices (an the LEED)
        lattices = self._ctrls['lattices_canvas']
        lattices.wheel_buddy = self._ctrls['rotation'].text

        # (2) LEED pattern
        # * Draw the LEED screen as a circle. It will also act as a
        #   background, and will be used as a clip path for the pattern
        leed = self._ctrls['leed_canvas']
        leed.initLEEDScreen(lineThick=lattices.getSpinesThickness())            # TODO: snake_case. Probably make getSpinesThickness a spine_thickness @property

        # * Scrolling while the mouse cursor is on the LEED plot
        #   will change the energy
        leed.wheel_buddy = self._ctrls['energy'].text

        # (3) Button that opens the domain selector
        domains = self._ctrls['domains']
        domains.setFont(gl.AllGUIFonts().labelFont)
        domains.ensurePolished()
        domains.setSizePolicy(qtw.QSizePolicy.Fixed, qtw.QSizePolicy.Fixed)

        # (4) Place controls in the layout
        layout = self.centralWidget().layout()
        layout.setHorizontalSpacing(15)
        layout.setContentsMargins(20, 10, 20, 10)

        layout.addLayout(lattices.mplLayout,                                    # TODO: snake
                         0, 1,
                         alignment=qtc.Qt.AlignCenter)
        layout.addLayout(leed.mplLayout,                                        # TODO: snake
                         0, 2,
                         alignment=qtc.Qt.AlignCenter)
        layout.addWidget(self._ctrls['rotation'],
                         0, 2,
                         alignment=qtc.Qt.AlignTop | qtc.Qt.AlignLeft)
        layout.addWidget(self._ctrls['energy'],
                         0, 2,
                         alignment=qtc.Qt.AlignTop | qtc.Qt.AlignRight)
        layout.addWidget(domains,
                         0, 2,
                         alignment=qtc.Qt.AlignBottom | qtc.Qt.AlignRight)
        # layout.setColumnStretch(0, 2)
        # layout.setColumnStretch(3, 2)

    def __compose_menu_and_toolbar(self):                                       # TODO: Common stuff that can go to base?: New, Open, Save, Save As, Exit
        """Set up the menu and a matching tool bar."""
        menu_bar = self.menuBar()
        tool_bar = self.addToolBar('')
        tool_bar.setFloatable(False)
        tool_bar.setIconSize(0.7*tool_bar.iconSize())

        file_menu = self._ctrls['file_menu']
        menu_bar.insertMenu(self.about_action,  # before 'About'
                            file_menu)          # insert the rest

        # Need to keep 'exit' to insert 'export' and a separator.
        # All other actions are connected to their handlers.
        exit = None
        for name, (icon, text, shortcut,
                   tooltip, statustip,
                   handler, to_enable) in default_file_menu().items():
            if icon:
                action = qtw.QAction(icon, text, self)
            else:
                action = qtw.QAction(text, self)
            action.setShortcut(shortcut)
            shortcut = action.shortcut().toString()
            action.setToolTip(tooltip + f" ({shortcut})")
            action.setStatusTip(statustip + f" ({shortcut})")
            action.triggered.connect(getattr(self, handler))
            if name == 'exit':
                exit = action
            file_menu.addAction(action)
            if icon:
                tool_bar.addAction(action)
            if to_enable:
                 self.__enabled_on_valid_input.append(action)

        # File -> Export sub-menu...
        export_menu = qtw.QMenu('Export')
        self.__enabled_on_valid_input.append(export_menu)

        # ...and its contents:
        # * Export list of beams to *.csv
        export_csv = qtw.QAction(
            self.style().standardIcon(qtw.QStyle.SP_DriveFDIcon),
            '&Export Beam List...', self)
        export_csv.setShortcut('Ctrl+E')
        shortcut = export_csv.shortcut().toString()
        export_csv.setToolTip('Export list of LEED beams' + f" ({shortcut})")
        export_csv.setStatusTip("Export list of LEED beams "
                                + f"to *.csv ({shortcut})")
        export_csv.triggered.connect(self._on_export_beams_pressed)
        self.__enabled_on_valid_input.append(export_csv)

        export_menu.addAction(export_csv)

        file_menu.insertMenu(exit, export_menu)
        file_menu.insertSeparator(exit)

    def __connect(self):
        """Connect relevant controls and signals."""
        self.screen_changed.connect(self._on_screen_changed)
        self._dialogs['file_new'].leed_parameters_changed.connect(
            self.__set_parameters
            )
        self._dialogs['file_new'].editing_finished.connect(
            self._on_structure_edit_finished
            )

    def exportCSV(self, params):                                                # TODO: rename, fix & doc
        # Ask if the user wants to save the input if it isn't
        if not self.saved:
            self.unsavedPopup()

        # Then open a normal file dialog to save the data to file
        fname = qtw.QFileDialog.getSaveFileName(self, 'Export list of beams',
                                                self.default_export[0],
                                                'Comma-separated file (*.csv)')
        if not fname or not fname[0]:
            return

        # set up the other parameters needed for export_pattern_csv
        if self.filename:
            params['source'] = self.filename
        params['version'] = gl.GLOBALS['version']

        gl.export_pattern_csv(fname[0], (self.leed,), **params)

    def unsavedPopup(self):                                                     # TODO: fix & doc
        reply = qtw.QMessageBox.question(self, 'Edits unsaved',
                                        'Would you like to save changes?',
                                        qtw.QMessageBox.Yes
                                        | qtw.QMessageBox.No,
                                        qtw.QMessageBox.Yes)
        if reply == qtw.QMessageBox.Yes:
            self._on_file_save_pressed()

    def initRealAndLEED(self, active):                                          # TODO: reimplement, rename & doc
        if active:
            # An acceptable LEED input has been loaded.
            # Update all controls with the correct stuff
            self.real = gl.RealSpace(self.leedParams)
            self.leed = gl.LEEDPattern(self.leedParams)
            for ctrl in self.allCtrls:
                self.updateCtrlAfterFileLoad(ctrl)
            self._ctrls['domains'].initPopup()
            # initially set focus to a widget that does not respond to
            # wheelEvent.
            self._ctrls['domains'].toggle.setFocus()

        self.connectControlEvents(active)

    def connectControlEvents(self, active):
        self.buts = [self._ctrls['energy'].enUp, self._ctrls['energy'].enDown,
                     self._ctrls['rotation'].cw, self._ctrls['rotation'].ccw,
                     self._ctrls['rotation'].h10, self._ctrls['rotation'].v10,
                     self._ctrls['rotation'].h01, self._ctrls['rotation'].v01,
                     self._ctrls['domains'].toggle]
        # First disconnect all controls, as this prevents multiple connections  # LOOK AT Qt::UniqueConnection. Should not require to disconnectAll()
        # from being established every time a new file is opened
        self.disconnectAll()

        if active:
            self._ctrls['energy'].text.textModified.connect(self.energyChanged)
            self._ctrls['rotation'].text.textModified.connect(self.rotationChanged)
            [but.clicked.connect(self.buttonPressed) for but in self.buts]

    def disconnectAll(self):
        for (widg, slt) in zip([self._ctrls['energy'].text,
                                self._ctrls['rotation'].text],
                               [self.energyChanged, self.rotationChanged]):
            try:
                widg.textModified.disconnect(slt)
            except TypeError:
                pass
        for but in self.buts:
            try:
                but.clicked.disconnect()
            except TypeError:
                pass

    def updateCtrlAfterFileLoad(self, ctrl):                                    # TODO: make it a list of things to do rather than 1kg of if...elif
        if ctrl not in self.allCtrls:
            raise ValueError('Unknown control')

        leed = self.leed
        real = self.real

        if ctrl == self._ctrls['rotation'].text:
            ctrl.setText('0.0')
        elif ctrl == self._ctrls['energy'].text:
            ctrl.setText(f'{leed.max_energy/2:.1f}')
            ctrl.setLimits(10, leed.max_energy)
        elif ctrl == self._ctrls['energy'].limits:
            ctrl.setText('Min = 10 eV\n'
                         f'Max = {leed.max_energy:.1f} eV')
        elif ctrl == self._ctrls['domains']:
            ctrl.updateText(f'{leed.n_domains} inequivalent domain(s)')
            ctrl.setTips(text='Click to see all the superlattice '
                              'matrices that generate the domains.')
            ctrl.hide = False
            ctrl.toggle.setEnabled(leed.n_domains != 1)
        elif False and ctrl == self._ctrls['cell_shapes']:
            self.insertText(ctrl,
                            f"Slab: {real.surf.cell_shape}, {real.surf.group}. "
                            f"Bulk: {real.bulk.cell_shape}, {real.bulk.group}",
                            'center')
            ctrl.setStatusTip('Cell shapes and plane groups '
                              'for the whole slab and its bulk.')
        elif ctrl == self._ctrls['lattices_canvas']:
            ctrl.plotLattices()
        elif ctrl == self._ctrls['leed_canvas']:
            ctrl.plotLEED(self._ctrls['domains'].hide)
        else:
            pass
            # Nothing should be done for the others
            # (except enabling and connecting events, which is done elsewhere)

        self.statusBar().showMessage('Ready')

    # The next event handlers take care of user-induced changes
    # on the controls
    def energyChanged(self, eOld, eNew):
        self._ctrls['leed_canvas'].plotLEED(self._ctrls['domains'].hide)

    def rotationChanged(self, rotOld, rotNew):
        self._ctrls['lattices_canvas'].plotLattices()
        self._ctrls['leed_canvas'].plotLEED(self._ctrls['domains'].hide)

    def buttonPressed(self, event):
        pressed = self.sender()
        if pressed in self._ctrls['energy'].subWidgs:
            self._ctrls['energy'].on_buttonPressed(pressed)
        elif pressed in self._ctrls['rotation'].subWidgs:
            self._ctrls['rotation'].on_buttonPressed(pressed)
        elif pressed == self._ctrls['domains'].toggle:
            self._ctrls['domains'].togglePressed()
            self._ctrls['leed_canvas'].plotLEED(self._ctrls['domains'].hide)

    def __read_files_and_report_errors(self, fnames, silent=False):
        """Read input files and report errors and warnings.

        Errors and warnings will be reported in a dedicated
        message box if not silenced.

        Parameters
        ----------
        fnames : Sequence
            An iterable of the file names of the files to be read
        silent : bool, optional
            If True, errors are swallowed, and are not reported
            with a message box. Errors are always reported in the
            status bar. Default is False.

        Returns
        -------
        leed_parameters : LEEDParametersList
            The LEED parameters ready for being loaded in the
            controls. leed_parameters is an empty list if fnames
            are an invalid input.
        """
        error_msg = gl.ErrorBox(error_while="opening LEED pattern input",
                                parent=self, silent=silent)
        parser = gl.LEEDParser()

        # Use a context manager to catch DeprecationWarning, which
        # is 'raised' for old-style files (without section headers)
        with warnings.catch_warnings(record=True) as warns:
            warnings.simplefilter("always", category=DeprecationWarning)
            parser.read(fnames)
            if warns and not silent:
                self.statusBar().showMessage(
                    "Old-style input file(s): Will be "
                    "reformatted and overwritten.",
                    5000
                    )

        try:
            params = gl.LEEDParametersList(parser)
        except NameError as err:
            # Missing parameters
            self.statusBar().showMessage(
                f"Invalid LEED input file(s): {err.args[0]}.", 5000
                )
            error_msg.exec_(self.statusBar().currentMessage())
            return gl.LEEDParametersList()
        except (ValueError, RuntimeError) as err:
            # Invalid data
            self.statusBar().showMessage(
                f"Invalid LEED input file(s): {err.args[0]}", 5000
                )
            error_msg.exec_(self.statusBar().currentMessage())
            return gl.LEEDParametersList()

        if not params:
            try:
                raise EOFError(f"The selected file(s) contain no data.")
            except EOFError as err:
                self.statusBar().showMessage("File(s) are empty", 5000)
                error_msg.exec_(f"Invalid LEED input file: {err.args[0]}")
                return gl.LEEDParametersList()
        return params

    def __zz_test_annots_compose(self):                                        # TODO: remove
        # self.ensurePolished()
        # tmpscreencenter = self.recSpace.geometry().center()
        # ctmp = winWidg.mapToGlobal(self.recSpace.geometry().center())
        # print('\n\n', self.pos(), tmpscreencenter,
              # ctmp, self.recSpace.parentWidget()==winWidg, '\n')
        tmpC = np.array([500, 500])#np.array([530.5, 500.5])
        # print("\ncreating")
        self.tmps = [gl.HoverAnnot(winWidg) for i in range(12)]
        deltas = [np.array([np.cos(x), np.sin(x)]) for x in np.linspace(0, 2*np.pi, num=len(self.tmps), endpoint=False)]
        # print("assign deltas")
        for (tmp, delta) in zip(self.tmps, deltas):
            tmp.head_pos = delta*5 + tmp.center
            # tmp.show()

    def __zz_annots_showEvent(self, event):
        # print("begin showEvent")
        super().showEvent(event)

        # print("super().showEvent finished")
        tmpscreencenter = self.recSpace.geometry().center()
        tmpLayCtr = self.recSpace.mplLayout.geometry().center()
        ctmp = self.window().centralWidget().mapToGlobal(
        # ctmp = self.window().mapToGlobal(
               self.recSpace.geometry().center())
        # print('\n\n', tmpscreencenter, tmpLayCtr, ctmp, self.recSpace.parentWidget().mapToGlobal(
               # self.recSpace.geometry().center()), self.recSpace.parentWidget().mapToGlobal(qtc.QPoint(0,0)))
        # print('\n\n', self.pos(), tmpscreencenter,
              # ctmp, self.recSpace.parentWidget()==winWidg, '\n')
        for tmp in self.tmps:
            # print(f"before: {tmp.center}", end=' ')
            tmp.center = ctmp
            tmp.show()
            # print(f"after: {tmp.center}")
