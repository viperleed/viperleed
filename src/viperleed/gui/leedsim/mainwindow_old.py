"""Module mainwindow_old of viperleed.gui.leedsim.

This is the main window of the ViPErLEED pattern-simulation plug-in.
This is the legacy version that only supports a single structural
domain.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2020-01-13'
__license__ = 'GPLv3+'

#TESTING
from time import perf_counter
#TESTING

import re
import copy

import numpy as np
import PyQt5.QtCore as qtc
import PyQt5.QtWidgets as qtw

from viperleed import __version__
from viperleed.gui.classes.planegroup import PlaneGroup
from viperleed.gui.helpers import string_matrix_to_numpy
from viperleed.gui.leedsim.classes.oldleedpatterns import LEEDPattern
from viperleed.gui.leedsim.classes.realspace import RealSpace
from viperleed.gui.leedsim.dialogs.exportcsvdialog import ExportCSVDialog
from viperleed.gui.leedsim.dialogs.newfiledialog_old import NewFileDialog
from viperleed.gui.leedsim.exportcsv import export_pattern_csv
from viperleed.gui.leedsim.widgets.domainsblock import DomsBlock
from viperleed.gui.leedsim.widgets.energyblock import EnergyBlock
from viperleed.gui.leedsim.widgets.hoverannot import HoverAnnot
from viperleed.gui.leedsim.widgets.leedcanvas import LEEDCanvas
from viperleed.gui.leedsim.widgets.realcanvas import RealCanvas
from viperleed.gui.leedsim.widgets.rotationblock import RotationBlock
from viperleed.gui.pluginsbase import ViPErLEEDPluginBase
from viperleed.gui.widgetdecorators import broadcast_mouse
from viperleed.gui.widgetslib import AllGUIFonts


@broadcast_mouse
class LEEDPatternSimulator(ViPErLEEDPluginBase):

    extension = '*.tlm'
    version = __version__

    extStr = 'LEED input files ({})'.format(' '.join([extension,
                                                      extension.upper()]))
    extStr = ';;'.join([extStr, 'All files (*)'])
    default_open = ('./gui/leedsim/input examples/', 'PatternInfo.tlm')
    default_export = ('./gui/leedsim/exported/', 'LEEDSpots.csv')

    # these are the parameters to read from the LEED file
    inputParams = ('eMax', 'surfBasis', 'SUPERLATTICE', 'surfGroup',
                   'bulkGroup')
    optionalParams = ('bulk3Dsym', 'screenAperture')

    def __init__(self):
        super().__init__()
        self.fonts = AllGUIFonts()
        self.leedParams = dict()
        self.openedFile = None  # Keep track of the name of the currently
                                # open file. This is used for exporting.
        self.initUI()
        self.centerOnScreen()

    def centerOnScreen(self):
        qr = self.frameGeometry()
        cp = qtw.QDesktopWidget().availableGeometry(self).center()
        qr.moveCenter(cp)
        self.move(qr.topLeft())

    def setMenu(self):
        menubar = self.menuBar()

        self.actionsToEnable = []

        # "File"
        self.fileMenu = qtw.QMenu('&File')
        menubar.insertMenu(self.about_action, self.fileMenu)

        # File -> New
        newF = qtw.QAction('&New...', self)
        newF.setShortcut('Ctrl+N')
        newF.setStatusTip('New LEED pattern file')
        newF.triggered.connect(self.fileNewDialog)

        # File -> Open
        openF = qtw.QAction('&Open...', self)
        openF.setShortcut('Ctrl+O')
        openF.setStatusTip('Open LEED pattern file')
        openF.triggered.connect(self.fileOpenDialog)

        # File -> Save
        saveF = qtw.QAction('&Save', self)
        saveF.setShortcut('Ctrl+S')
        saveF.setStatusTip('Save LEED pattern file')
        saveF.triggered.connect(self.fileSavePressed)

        self.actionsToEnable.append(saveF)

        # File -> Save As...
        svAs = qtw.QAction('&Save As...', self)
        svAs.setShortcut('Ctrl+Shift+S')
        svAs.setStatusTip('Save LEED pattern file as...')
        svAs.triggered.connect(self.fileSaveAsPressed)
        self.actionsToEnable.append(svAs)

        # File -> Export sub-menu...
        exportMenu = qtw.QMenu('Export')
        self.actionsToEnable.append(exportMenu)

        # ...and its contents:
        # * Export list of beams to *.csv
        exportCSV = qtw.QAction('&Export Beam List...', self)
        exportCSV.setShortcut('Ctrl+E')
        exportCSV.setStatusTip('Export list of LEED beams to *.csv')
        exportCSV.triggered.connect(self.exportBeamsPressed)
        self.actionsToEnable.append(exportCSV)

        exportMenu.addAction(exportCSV)

        # File -> Exit
        exit = qtw.QAction('&Exit', self)
        exit.setShortcut('Ctrl+Q')
        exit.setStatusTip('Exit application')
        exit.triggered.connect(qtw.qApp.quit)
        #exitAct.triggered.connect(self.close)
        # Replacing line 109 with 110 makes the popup appear also when
        # closing with CTRL+Q

        self.fileMenu.addAction(newF)
        self.fileMenu.addAction(openF)
        self.fileMenu.addAction(saveF)
        self.fileMenu.addAction(svAs)
        self.fileMenu.addMenu(exportMenu)
        self.fileMenu.addSeparator()
        self.fileMenu.addAction(exit)

    def setToolBar(self):
        self.tBar = self.addToolBar('')
        self.tBar.setFloatable(False)
        self.tBar.setIconSize(0.7*self.tBar.iconSize())

        newF = qtw.QAction(self.style().standardIcon(qtw.QStyle.SP_FileIcon),
                       'New LEED pattern', self)
        newF.setStatusTip('New LEED pattern file (Ctrl+N)')
        newF.triggered.connect(self.fileNewDialog)

        openF = qtw.QAction(
            self.style().standardIcon(qtw.QStyle.SP_DialogOpenButton),
            'Open LEED pattern', self)
        openF.setStatusTip('Open LEED pattern file (Ctrl+O)')
        openF.triggered.connect(self.fileOpenDialog)

        saveF = qtw.QAction(
            self.style().standardIcon(qtw.QStyle.SP_DialogSaveButton),
            'Save LEED pattern', self)
        saveF.setStatusTip('Save LEED pattern file (Ctrl+S)')
        saveF.triggered.connect(self.fileSavePressed)

        exportAct = qtw.QAction(
            self.style().standardIcon(qtw.QStyle.SP_DriveFDIcon),
            'Export LEED pattern', self)
        exportAct.setStatusTip('Export list of LEED beams to '
                               'comma-separated file (Ctrl+E)')
        exportAct.triggered.connect(self.exportBeamsPressed)

        self.tBar.addAction(newF)
        self.tBar.addAction(openF)
        self.tBar.addAction(saveF)
        self.tBar.addAction(exportAct)

        self.actionsToEnable.append(saveF)
        self.actionsToEnable.append(exportAct)

    def fileOpenDialog(self):
        fname = qtw.QFileDialog.getOpenFileName(self, 'Open LEED pattern',
                                                self.default_open[0],
                                                self.extStr)
        loadedParams = self.parseLEEDFile(fname)
        if loadedParams is not None:
            self.openedFile = fname
            self.leedParams = loadedParams
            self.loadLEEDInput(loadedParams, True)
            self.isSaved = True
            self.openFile = fname[0]
        else:
            self.fileReadError(fname=fname[0])

    def fileNewDialog(self):
        if not hasattr(self, 'newDialog'):
            self.newDialog = NewFileDialog(self, self.leedParams)
            self.newDialog.dialogEdited.connect(self.loadLEEDInput)
            self.isSaved = False
        else:
            self.newDialog.open(self.leedParams)

    def fileSavePressed(self):
        if (not hasattr(self, 'openFile')) or (self.openFile is None):
            # no file is currently open
            self.fileSaveAsPressed()
        else:
            # TODO: ADD HERE A MESSAGE BOX ASKING IF YOU WANT TO OVERWRITE?
            self.saveToFile(self.openFile)
            self.isSaved = True

    def fileSaveAsPressed(self):
        fname = qtw.QFileDialog.getSaveFileName(self, 'Save TLM input', '',
                                                self.extStr)
        if fname is not None and fname[0]:
            # file selected correctly
            self.saveToFile(fname[0])
            self.openFile = fname[0]
            self.isSaved = True

    def saveToFile(self, fname):
        # Sort out the parameters, converting the values to strings
        params = copy.deepcopy(self.leedParams)
        for key, val in params.copy().items():  # BUG: dictionary changes size upon saving
            if key == 'SUPERLATTICE' or key == 'surfBasis':
                params[key] = np.array2string(
                                        val,
                                        separator=',',
                                       # precision=5,
                                        suppress_small=True
                                        ).replace('\n', '').replace(' ', '')
            elif key in ('eMax', 'screenAperture'):
                params[key] = str(val)
            elif key in ('bulkGroup', 'surfGroup'):
                params[key] = val.group
            elif key == 'bulk3Dsym' and val is None:
                del params[key]
        # now params contains correctly formatted strings
        params = list(params.items())
        paramsTxt = '\n'.join('{}'.format('='.join(param))
                              for param in params)
        with open(fname, 'w+') as f:
            f.write(paramsTxt)

    def exportBeamsPressed(self):
        # Create a dialog to select which domains should be exported
        # and to process the data to export
        if not hasattr(self, 'exportDialog'):
            self.exportDialog = ExportCSVDialog(self.leed, parent=self)
            self.exportDialog.exportSelected.connect(self.exportCSV)
        else:
            self.exportDialog.open()

    def parseLEEDFile(self, fname, checkOnly=False):
        # returns: None if file is unreadable/no file was selected; empty
        # dict if file does not contain all parameters/incorrect ones; dict
        # with parameters otherwise
        self.fileErr = ''   # keeps track of what happened

        if fname[0]:
            with open(fname[0], 'r') as f:
                try:
                    flines = f.readlines()
                except:
                    ext = qtc.QFileInfo(fname).suffix().lower()
                    self.fileErr = 'ERROR: File type *.%s not supported.'%(ext)
                    return None
        else:
            self.fileErr = 'No file selected.'
            return None

        if checkOnly:
            return dict()

        if flines is None:
            return dict()

        # Something has been successfully read
        # get rid of white spaces (i.e, ' ', '\t', '\n', '\r')
        flines = [re.sub(r'\s+', '', line) for line in flines]

        #get rid of comments, that begin with \# \! \%
        commentChars = ['#', '!', '%']
        for ch in commentChars:
            flines = [line.partition(ch)[0].lower() for line in flines]
            # i.e., take whatever block is on the left side of any of the
            # comment characters

        # now keep only lines that contain exactly one '=' sign
        flines = [line for line in flines if line.count('=') == 1]

        if len(flines) == 0:
            self.fileErr = 'File does not contain any LEED parameter.'
            return dict()

        # and split them in two:
        # _param contains what's on the left of '=',
        # _values what's on the right
        flines_param, flines_values = zip(
                         *[np.array(line.partition('='))[[0, 2]].tolist()
                           for line in flines
                           if len(line) > 1]
                         )
        # make a dictionary out of these. This also removes duplicate
        # parameters, taking the last occurrence as the value
        par_dict = dict(zip(flines_param, flines_values))

        if len(par_dict) < len(flines_param):
            self.fileErr = ('Warning: Repeated entries found;'
                            ' using last occurrence.')

        # remove all entries that are not acceptable parameters
        acceptable = [p.lower()
                      for p in (*self.inputParams, *self.optionalParams)]
        par_dict = {key.lower(): val
                    for key, val in par_dict.items()
                    if key.lower() in acceptable}

        # replace the dictionary keys in par_dict with their standard
        # capitalization, found in self.inputParams and in self.optionalParams.
        # Will make a difference between keys that are mandatory parameters
        # and keys that are optional ones.
        for par in (*self.inputParams, *self.optionalParams):
            try:
                par_dict[par.lower()]
            except KeyError:
                if par in self.inputParams:
                    self.fileErr = (f'ERROR: Mandatory parameter {par} '
                                    'not found. '
                                    + self.fileErr)
                    return dict()
                # otherwise the parameter is optional, and will be set to a
                # default of 'None'
                par_dict[par] = 'None'
            else:
                par_dict[par] = par_dict.pop(par.lower())

        par_dict['eMax'] = float(par_dict['eMax'])
        par_dict['surfBasis'] = string_matrix_to_numpy(
            par_dict['surfBasis'], dtype=float, needs_shape=(2, 2)
            )
        par_dict['SUPERLATTICE'] = string_matrix_to_numpy(
            par_dict['SUPERLATTICE'], int, needs_shape=(2, 2)
            )

        try:
            par_dict['screenAperture'] = float(par_dict['screenAperture'])
        except ValueError:
            self.fileErr += ("Warning: screenAperture not found, or value "
                             " not acceptable. Will use default 110 deg.")
            par_dict['screenAperture'] = 110.0

        # Now handle the exit conditions
        try:
            b_group = PlaneGroup(par_dict['bulkGroup'])
        except ValueError:
            self.fileErr = ('ERROR: Unknown bulk group type '
                           + par_dict['bulkGroup']
                           + self.fileErr)
            return dict()
        else:
            par_dict['bulkGroup'] = b_group

        try:
            s_group = PlaneGroup(par_dict['surfGroup'])
        except ValueError:
            self.fileErr = ('ERROR: Unknown surface group type '
                            + par_dict['surfGroup']
                            + self.fileErr)
            return dict()
        else:
            par_dict['surfGroup'] = s_group

        if ((par_dict['surfBasis'] is None)
            or (par_dict['SUPERLATTICE'] is None)):
            self.fileErr  = ('ERROR: Too many entries in '
                            + 'surface unit vectors or SUPERLATTICE. '
                            + self.fileErr)
            return dict()
        return par_dict

    def exportCSV(self, params):
        # Ask if the user wants to save the input if it isn't
        if not self.isSaved:
            self.unsavedPopup()

        # Then open a normal file dialog to save the data to file
        fname = qtw.QFileDialog.getSaveFileName(self, 'Export list of beams',
                                                self.default_export[0],
                                                'Comma-separated file (*.csv)')
        if not fname[0]:
            return

        # set up the other parameters needed for export_pattern_csv
        if hasattr(self, 'openFile') and self.openFile:
            params['source'] = self.openFile

        export_pattern_csv(fname[0], (self.leed,), **params)

    def unsavedPopup(self):
        reply = qtw.QMessageBox.question(self, 'Edits unsaved',
                                        'Would you like to save changes?',
                                        qtw.QMessageBox.Yes
                                        | qtw.QMessageBox.No,
                                        qtw.QMessageBox.Yes)
        if reply == qtw.QMessageBox.Yes:
            self.fileSavePressed()

    def loadLEEDInput(self, par_dict, readFile=False):
        if par_dict:
            self.leedParams = par_dict
            self.set_ControlsActive(True)
            self.initRealAndLEED(True)
            self.enableActions(True)
        if readFile:
            self.fileReadError()
        if hasattr(self, 'exportDialog'):
            # delete exportDialog as the number of its widgets (that depends
            # on the number of domains) is most likely not accurate anymore.
            self.exportDialog.destroy()
            del self.exportDialog

    def insertText(self, label, text, align='left'):
        # sets, reshapes, and realigns horizontally a QWidget after changing
        # its text

        if align not in ['left', 'center', 'right']:
            raise ValueError("Alignment of QWidget can only be 'left',"
                             "'center' or 'right'")

        pos0 = label.geometry().topLeft()
        if align == 'left':
            delta0 = 0
        elif align == 'right':
            delta0 = label.width()
        else:
            delta0 = label.width()//2

        label.setText(text)
        label.adjustSize()

        if align == 'left':
            delta1 = 0
        elif align == 'right':
            delta1 = label.width()
        else:
            delta1 = label.width()//2

        label.move(pos0.x() - delta1 + delta0, pos0.y())

    def initControls(self):
        # Initialize all controls and labels, and set their dimensions
        # w is the centralWidget

        w = self.window().centralWidget()

        #### REAL SPACE AND LEED PATTERN PLOTTING CANVASES
        self.realSpace = RealCanvas(parent=w, title='Real Space Lattice')
        self.recSpace = LEEDCanvas(parent=w, title='LEED Pattern')

        # draw the LEED screen as a circle. It will also act as a background,
        # and will be used as a clip path for the pattern
        self.recSpace.initLEEDScreen(
                               lineThick=self.realSpace.getSpinesThickness())

        #### ROTATION ####
        self.rotWidg = RotationBlock(w)
        self.realSpace.wheel_buddy = self.rotWidg.text

        #### ENERGY ####
        self.enWidg = EnergyBlock(w)
        self.recSpace.wheel_buddy = self.enWidg.text

        #### DOMAINS ####
        self.doms = DomsBlock(w)

        #### Unit cells shape and group ####
        self.cellShapes = qtw.QLabel("Slab: \u2014, \u2014. "
                                     "Bulk: \u2014, \u2014",
                                     w)
        self.cellShapes.setFont(AllGUIFonts().largeTextFont)
        self.cellShapes.setStatusTip('Open or drag-drop a LEED pattern file'
                                     'to print the shapes and symmetry'
                                     'groups of the system!')
        self.cellShapes.setAlignment(qtc.Qt.AlignHCenter | qtc.Qt.AlignBaseline)
        self.cellShapes.setSizePolicy(qtw.QSizePolicy.Expanding,
                                      qtw.QSizePolicy.Preferred)
        self.cellShapes.adjustSize()

        # Prepare a list of controls
        # Notice that the plot controls are at the end on purpose,
        # so all the others are updated first after loading a file
        self.allCtrls = [self.rotWidg, *(self.rotWidg.subWidgs),
                         self.enWidg, *(self.enWidg.subWidgs),
                         self.doms, self.cellShapes,
                         self.realSpace, self.recSpace]

        # Disable all controls at initialization.
        # They will be enabled when a file is loaded
        self.set_ControlsActive(False)

    def set_ControlsActive(self, active):
        if isinstance(active, bool):
            [ctrl.setEnabled(active) for ctrl in self.allCtrls]
        else:
            raise ValueError('Invalid datatype for set_ControlsActive')

    def composeUI(self):
        # Position all controls correctly

        winWidg = self.window().centralWidget()

        g0 = winWidg.geometry()   # Origin and size of centralWidget
                                  # of the window, i.e. window size minus
                                  # size of menu and status bars
        w0 = g0.width()
        h0 = g0.height()
        self.area = w0*h0

        topLayout = winWidg.layout()
        '''
        Layout:
            Legend: sh/sv(n)= horizontal/vertical stretch;
                    RS= real space (VBox)
                    LD= leed (Vbox)
                    Ro= rotation block (excluding bulk buttons)
                    bu= bulk alignment buttons
                    En= energy block (excluding text below)
                    Et= energy limits text
                    Db= domains toggle button
                    Dt= domains text
                    Sh= slab shape (sh+Hbox+sh)

                    if more symbols are present in a cell, the left ones correspond to widgets/items below the right ones

              0      1        2        3        4        5        6        7
            -------------------------------------------------------------
        0    | sv |    RS    |  sv    |  LD    |  sv    |        |        |        |
            -------------------------------------------------------------
        1    | sv |        |  sv    |        |  sv    |        |        |        |
            -------------------------------------------------------------
        2    | sv |        |  sv    |        |  sv    |        |        |        |
            -------------------------------------------------------------
        3    | sv |        |  sv    |        |  sv    |        |        |        |
            -------------------------------------------------------------
        4    | sv |        |  sv    |        |  sv    |        |        |        |
            -------------------------------------------------------------
        5    | sv |        |  sv    |        |  sv    |        |        |        |
            -------------------------------------------------------------
        6    | sv |        |  sv    |        |  sv    |        |        |        |
            -------------------------------------------------------------
        7    | sv |        |  sv    |        |  sv    |        |        |        |
            -------------------------------------------------------------
        '''

        ## REAL SPACE AND LEED PATTERN
        if self.testConfig == 'layout':
            topLayout.setSpacing(3)
            topLayout.setContentsMargins(0, 0, 0, 0)
            basicStretch = 1
            plotStretch = (topLayout.columnStretch(0)
                           * self.realSpace.getStretchFactor().width())

            topLayout.addLayout(self.realSpace.mplLayout, 0, 1, 1, 1)
            topLayout.setColumnStretch(0, basicStretch)
            topLayout.setColumnStretch(1, plotStretch)
            topLayout.setColumnStretch(2, basicStretch)

            topLayout.addLayout(self.recSpace.mplLayout, 0, 3, 1, 1)
            topLayout.setColumnStretch(3, plotStretch)
            topLayout.setColumnStretch(4, basicStretch)

            #[topLayout.setColumnStretch(col,1) for col in [3,4,10,11]]
            #topLayout.setRowStretch(3,1)
            # topLayout.setRowStretch(0,1)
            # topLayout.setRowStretch(1,2)
            #topLayout.setRowStretch(5,1)
            #topLayout.setRowStretch(6,1)
        else:
            plotsHBox = qtw.QHBoxLayout()
            plotsHBox.addStretch(1)
            plotsHBox.addLayout(self.realSpace.mplLayout)
            plotsHBox.addStretch(1)
            plotsHBox.addLayout(self.recSpace.mplLayout)
            plotsHBox.addStretch(1)

            winWidg.layout().addLayout(plotsHBox)

        labelToTextOffset = -qtc.QPoint(0,3)

        #### ROTATION ####
        if self.testConfig == 'layout':
            pass
        else:
            rotx = round((w0 - self.rotWidg.width())/2 + 0.032*w0)
            self.rotWidg.move(qtc.QPoint(rotx, 8))

        #### ENERGY ####
        if self.testConfig == 'layout':
            pass
        else:
            enx = round(w0*.98 - self.enWidg.width())
            self.enWidg.move(qtc.QPoint(enx,8))

        #### TOGGLE DOMAINS ###
        dx = self.enWidg.geometry().topRight().x() - self.doms.width()
        dy = round(h0*.97 - self.doms.height())
        self.doms.move(qtc.QPoint(dx,dy))

        #### CELLS SHAPES/GROUPS ####
        if self.testConfig == 'layout':
            self.cellShapes.setAlignment(qtc.Qt.AlignCenter)
                        #(qtc.Qt.AlignHCenter | qtc.Qt.AlignTop)
            shapesLay = qtw.QHBoxLayout()
            shapesLay.addStretch(1)
            shapesLay.addWidget(self.cellShapes,
                                qtc.Qt.AlignHCenter | qtc.Qt.AlignTop)
            shapesLay.addStretch(1)
            topLayout.addLayout(shapesLay, 1, 0, 1, topLayout.columnCount())#,qtc.Qt.AlignCenter)
            # topLayout.setRowStretch(0,100)
            # topLayout.setRowStretch(1,5)
            # vSpacer=QSpacerItem(1, 10, hPolicy = qtw.QSizePolicy.Expanding, vPolicy = qtw.QSizePolicy.Minimum)
            # topLayout.addItem(vSpacer,2,0,1,topLayout.columnCount())
            #topLayout.addWidget(self.cellShapes,1,0,1,1,qtc.Qt.AlignHCenter|qtc.Qt.AlignTop)
            #pass
        else:
            self.cellShapes.move(
                round((w0-self.cellShapes.geometry().width())/2),
                self.doms.geometry().center().y()
                )

        ##TESTING
        # self.ensurePolished()
        # tmpscreencenter = self.recSpace.geometry().center()
        # ctmp = winWidg.mapToGlobal(self.recSpace.geometry().center())
        # print('\n\n', self.pos(), tmpscreencenter,
              # ctmp, self.recSpace.parentWidget()==winWidg, '\n')
        tmpC = np.array([500, 500])#np.array([530.5, 500.5])
        # print("\ncreating")
        self.tmps = [HoverAnnot(winWidg) for i in range(12)]
        deltas = [np.array([np.cos(x), np.sin(x)]) for x in np.linspace(0, 2*np.pi, num=len(self.tmps), endpoint=False)]
        # print("assign deltas")
        for (tmp, delta) in zip(self.tmps, deltas):
            tmp.head_pos = delta*5 + tmp.center
            # tmp.show()

        ## END TESTING

    def enableActions(self, enable):
        if not isinstance(enable, bool):
            raise
        [act.setEnabled(enable) for act in self.actionsToEnable]

    def initUI(self):
        # print("Initializing UI...")
        # print("... Menus")
        self.setMenu()

        # print("... Tool bar")
        self.setToolBar()
        self.enableActions(False)

        self.setAcceptDrops(True)

        # print("... Title and icon")
        self.setWindowTitle('LEED Pattern Indexing Helper')

        # print("... Status bar")
        self.sBar = qtw.QStatusBar()
        self.setStatusBar(self.sBar)
        self.sBar.showMessage('Ready')

        self.testConfig = 'fix'#'layout'#
        cw = qtw.QWidget()
        if self.testConfig == 'layout':
            topLayout = qtw.QGridLayout(cw)
        else:
            topLayout = qtw.QHBoxLayout(cw)
        # print(f"... Central widget ({cw})")
        self.setCentralWidget(cw)

        self.adjustWindowSize()

        # print("... Initialize controls")
        self.initControls()
        # print("... Position controls")
        self.composeUI()

    def sizeHint(self):
        self.ensurePolished()
        screenFraction = 0.8
        screensize = qtw.QDesktopWidget().availableGeometry(self).size()
        sH = screenFraction*screensize
        return sH

    def adjustWindowSize(self, size=None):
        # Updates the size of the central widget, window and menus.
        if size is None:
            size = self.sizeHint()
        # the next line will need to replace the one after
        # once I make the window re-sizable
        #self.resize(size)
        self.setFixedSize(size)
        self.layout().activate() # this rebuilds the window layout, and forces the centralWidget to fill the space

    def resizeEvent(self, event):
        super().resizeEvent(event)
        self.adjustWindowSize(size=event.size())

    def initRealAndLEED(self, active):
        if active:
            # An acceptable LEED input has been loaded.
            # Update all controls with the correct stuff
            self.real = RealSpace(self.leedParams)
            self.leed = LEEDPattern(self.leedParams)
            for ctrl in self.allCtrls:
                self.updateCtrlAfterFileLoad(ctrl)
            self.doms.initPopup()
            # initially set focus to a widget that does not respond to
            # wheelEvent.
            self.doms.toggle.setFocus()

        self.connectControlEvents(active)

    def connectControlEvents(self, active):
        self.buts = [self.enWidg.enUp, self.enWidg.enDown,
                     self.rotWidg.cw, self.rotWidg.ccw,
                     self.rotWidg.h10, self.rotWidg.v10,
                     self.rotWidg.h01, self.rotWidg.v01,
                     self.doms.toggle]
        # First disconnect all controls, as this prevents multiple connections  # LOOK AT Qt::UniqueConnection
        # from being established every time a new file is opened
        self.disconnectAll()

        if active:
            self.enWidg.text.textModified.connect(self.energyChanged)
            self.rotWidg.text.textModified.connect(self.rotationChanged)
            [but.clicked.connect(self.buttonPressed) for but in self.buts]

    def disconnectAll(self):
        for (widg, slt) in zip([self.enWidg.text, self.rotWidg.text],
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

    def updateCtrlAfterFileLoad(self, ctrl):
        if ctrl not in self.allCtrls:
            raise ValueError('Unknown control')

        leed = self.leed
        real = self.real

        if ctrl == self.rotWidg.text:
            ctrl.setText('0.0')
        elif ctrl == self.enWidg.text:
            ctrl.setText(f'{leed.max_energy/2:.1f}')
            ctrl.setLimits(10, leed.max_energy)
        elif ctrl == self.enWidg.limits:
            ctrl.setText('Min = 10 eV\n'
                         f'Max = {leed.max_energy:.1f} eV')
        elif ctrl == self.doms:
            ctrl.updateText(f'{leed.n_domains} inequivalent domain(s)')
            ctrl.setTips(text='Click to see all the superlattice '
                              'matrices that generate the domains.')
            ctrl.hide = False
            if leed.n_domains == 1:
                ctrl.toggle.setEnabled(False)
            else:
                ctrl.toggle.setEnabled(True)
        elif ctrl == self.cellShapes:
            self.insertText(ctrl,
                            f"Slab: {real.surf.cell_shape}, {real.surf.group}. "
                            f"Bulk: {real.bulk.cell_shape}, {real.bulk.group}",
                            'center')
            ctrl.setStatusTip('Cell shapes and plane groups '
                              'for the whole slab and its bulk.')
        elif ctrl == self.realSpace:
            ctrl.plotLattices()
        elif ctrl == self.recSpace:
            ctrl.plotLEED(self.doms.hide)
        else:
            pass
            # Nothing should be done for the others
            # (except enabling and connecting events, which is done elsewhere)

        self.statusBar().showMessage('Ready')

    def get_energy(self): #is it possible to get rid of this?
        return float(self.enWidg.text.text())

    def get_angle(self): #is it possible to get rid of this?
        return float(self.rotWidg.text.text())

    def closeEvent(self, event):
        reply = qtw.QMessageBox.question(self,
                                     'Message',
                                     'Are you sure to quit?',
                                     qtw.QMessageBox.Yes | qtw.QMessageBox.No,
                                     qtw.QMessageBox.No)

        if reply == qtw.QMessageBox.Yes:
            super().closeEvent(event)
        else:
            event.ignore()

    def mouseReleaseEvent(self, event):
        """
        Re-implement mouseReleaseEvent to hide the pop-up with matrices when
        clicking anywhere on the window
        """
        doms = self.doms
        if (doms.isEnabled()
            and doms.matricesPopup.isVisible()
            and not doms.text.underMouse()):
            # the reason for not acting when doms.text is underMouse() is that
            # the event handling in this case is done within DomsBlock
            doms.matricesPopup.hide()
        # TESTING
        # for annot in self.recSpace.annots:
            # annot.show()
        # END TESTING
        super().mouseReleaseEvent(event)

    # The next event handlers take care of user-induced changes
    # on the controls
    def energyChanged(self, eOld, eNew):
        self.recSpace.plotLEED(self.doms.hide)

    def rotationChanged(self, rotOld, rotNew):
        self.realSpace.plotLattices()
        self.recSpace.plotLEED(self.doms.hide)

    def buttonPressed(self, event):
        pressed = self.sender()
        if pressed in self.enWidg.subWidgs:
            self.enWidg.on_buttonPressed(pressed)
        elif pressed in self.rotWidg.subWidgs:
            self.rotWidg.on_buttonPressed(pressed)
        elif pressed == self.doms.toggle:
            self.doms.togglePressed()
            self.recSpace.plotLEED(self.doms.hide)

    # The next two event handlers are for drag-drop of files
    def dragEnterEvent(self, event):
        mimeData = event.mimeData()
        if mimeData.hasUrls():   # Check that the user is dropping a file
            fname = [url.toLocalFile() for url in mimeData.urls()]

            # now handle differently the cases in which the user has
            # dropped multiple files:
            if len(fname) > 1:
                # for now, just take the first one. Then there will be a
                # popup asking which one to load
                fname = fname[0]
                event.mimeData().setUrls(fname)
            # check that the file contains some readable content
            loadedParams = self.parseLEEDFile(fname, checkOnly=True)
            if loadedParams is not None:
                event.acceptProposedAction()

    def dropEvent(self, event):
        loadedParams = self.parseLEEDFile([url.toLocalFile()
                                           for url
                                           in event.mimeData().urls()])
        if loadedParams:
            event.acceptProposedAction()
            self.loadLEEDInput(loadedParams, True)
            self.isSaved = True
        else:
            self.fileReadError()

    def fileReadError(self, fname=None, msg=None):
        if msg is None:
            msg = self.fileErr
        if msg != '':
            self.statusBar().showMessage(msg)

    # def showEvent(self, event):
        # # print("begin showEvent")
        # super().showEvent(event)

        # # print("super().showEvent finished")
        # tmpscreencenter = self.recSpace.geometry().center()
        # tmpLayCtr = self.recSpace.mplLayout.geometry().center()
        # ctmp = self.window().centralWidget().mapToGlobal(
        # # ctmp = self.window().mapToGlobal(
               # self.recSpace.geometry().center())
        # # print('\n\n', tmpscreencenter, tmpLayCtr, ctmp, self.recSpace.parentWidget().mapToGlobal(
               # # self.recSpace.geometry().center()), self.recSpace.parentWidget().mapToGlobal(qtc.QPoint(0,0)))
        # # print('\n\n', self.pos(), tmpscreencenter,
              # # ctmp, self.recSpace.parentWidget()==winWidg, '\n')
        # for tmp in self.tmps:
            # # print(f"before: {tmp.center}", end=' ')
            # tmp.center = ctmp
            # tmp.show()
            # # print(f"after: {tmp.center}")

    def moveEvent(self, event):
        super().moveEvent(event)
        delta = event.pos() - event.oldPos()

        for annot in self.recSpace.annots:
            # NB: one cannot use += for in-place assignment since this redefines
            # the .center attribute that is a @property!
            annot.center = qtc.QPoint(*annot.center.round()) + delta
            if annot.isVisible():
                # move only if visible to handle the situation in which the
                # window is moved while the annotations are displayed: this way
                # annotations move together with the window. If one would have
                # this call even when not visible, the move is issued before the
                # next show, and it will mess the position of the annotation
                annot.move(annot.pos() + delta)

