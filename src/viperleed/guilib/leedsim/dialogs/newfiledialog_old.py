"""
===============================================
      ViPErLEED Graphical User Interface
===============================================
 *** module guilib.leedsim.NewFileDialog ***

Created: 2020-01-11
Author: Michele Riva

"""

import re

import numpy as np
import PyQt5.QtCore as qtc
import PyQt5.QtGui as qtg
import PyQt5.QtWidgets as qtw

from viperleed.guilib.classes.lattice2d import Lattice2D as Lattice
from viperleed.guilib.leedsim.classes.woods import Woods                        # TODO: maybe the old one?
from viperleed.guilib.leedsim.dialogs.dialogbulk3dsym import Bulk3DSymDialog
from viperleed.guilib.widgetslib import AllGUIFonts

angstrom = ' \u212b'
degrees = '\u00b0'

class NewFileDialog(qtw.QDialog):
    dialogEdited = qtc.pyqtSignal(dict) # a dictionary of the LEED parameters
    
    def __init__(self, parent=None, oldParams=dict()):
        # oldParams is used to restore the old pattern if the users presses a 
        # 'cancel' button, and to initialize the values of the current 
        # parameters
        #super(NewFileDialog, self).__init__(parent)
        super().__init__(parent)
        self.setWindowModality(qtc.Qt.WindowModal)
        flags = self.windowFlags()
        flags &= ~qtc.Qt.WindowCloseButtonHint #disable close button
        self.setWindowFlags(flags)
        self.setWindowTitle('Create new LEED input')
        
        self.compose()
        self.woodL = Woods()
        self.bulk_3d_sym_dialog = Bulk3DSymDialog(self)
        self.open(oldParams)
    
    def open(self, params=dict()):
        #---- update some quantities for this session ----#
        self.oldParams = params  # LEED parameters when opened
        
        self.initLattices()
        self.initControlValues()

        super().open()
    
    def initLattices(self):
        if self.oldParams:  # dictionary not empty
            self.supMatrix = self.oldParams['SUPERLATTICE']
            surfBasis = self.oldParams['surfBasis']
            surfGroup = self.oldParams['surfGroup']
            bulkGroup = self.oldParams['bulkGroup']
        else:
            self.supMatrix = np.array([[1, 0], [0, 1]])
            a = [2, 0]
            b = [3*np.cos(np.radians(100)), 3*np.sin(np.radians(100))]
            surfBasis = np.array([a, b])
            surfGroup = 'p1'
            bulkGroup = surfGroup
        
        bulkBasis = np.dot(np.linalg.inv(self.supMatrix), surfBasis)
        
        self.surfLatt = Lattice(surfBasis, group=surfGroup)
        self.bulkLatt = Lattice(bulkBasis, group=bulkGroup)
        
        if self.oldParams:
            self.bulkLatt.group.screws_glides = (self.oldParams['bulk3Dsym'],
                                                 self.bulkLatt.cell_shape)

    def compose(self):
        labFont = AllGUIFonts().labelFont
        bigLabFont = qtg.QFont(AllGUIFonts().largeTextFont)
        bigLabFont.setPointSize(11)
        smallText = AllGUIFonts().smallTextFont
        
        #-------- Initialize the widgets --------#
        # here only the immutable values are set #
        
        # Input type selection
        lattInputLabel = qtw.QLabel('Lattice unit cells')
        lattInputLabel.setFont(bigLabFont)
        allWidgs = [lattInputLabel]
        
        # radio buttons for type selection
        self.bulkInput = qtw.QRadioButton('Bulk')
        self.bulkInput.setFont(smallText)
        self.surfInput = qtw.QRadioButton('Surface')
        self.surfInput.setFont(smallText)
        allWidgs.extend([self.bulkInput, self.surfInput])
        
        # lattice shape drop down menus
        latShapeLabel = qtw.QLabel('Shape')
        latShapeLabel.setFont(labFont)
        self.bulkShape = qtw.QComboBox()
        self.surfShape = qtw.QComboBox()
        for combo in [self.bulkShape, self.surfShape]:
            combo.setFont(smallText)
            combo.addItems(['Oblique', 'Rectangular', 'Square',
                            'Rhombic', 'Hexagonal'])
            combo.setSizePolicy(qtw.QSizePolicy.Fixed,
                                qtw.QSizePolicy.Preferred)
            combo.adjustSize()
        allWidgs.extend([self.bulkShape, self.surfShape])
        
        # lattice parameters
        aLab = qtw.QLabel('a = ')
        bLab = qtw.QLabel('b = ')
        alphaLab = qtw.QLabel('\u03b1 = ')
        aLab.setFont(labFont)
        bLab.setFont(labFont)
        alphaLab.setFont(labFont)
        lattLabels = [aLab, bLab, alphaLab, latShapeLabel]
        [lab.adjustSize() for lab in lattLabels]
        w = max([lab.width() for lab in lattLabels])
        [lab.setMaximumWidth(w) for lab in lattLabels]
        
        self.aB = qtw.QLineEdit('')
        self.bB = qtw.QLineEdit('')
        self.alphaB = qtw.QLineEdit('')
        
        self.aS = qtw.QLineEdit('')
        self.bS = qtw.QLineEdit('')
        self.alphaS = qtw.QLineEdit('')
        self.bulkParams = [self.aB, self.bB, self.alphaB]
        self.surfParams = [self.aS, self.bS, self.alphaS]
        lattParBoxes = [*self.bulkParams, *self.surfParams]
        for tbox in lattParBoxes:
            tbox.setFont(labFont)
            tbox.setSizePolicy(qtw.QSizePolicy.Fixed, qtw.QSizePolicy.Fixed)
            tbox.setMaximumWidth(self.bulkShape.width())
        
        allWidgs.extend([*lattLabels, *lattParBoxes])
        
        # groups
        groupLab = qtw.QLabel('Group')
        groupLab.setFont(labFont)
        self.bulkGroup = qtw.QComboBox()
        self.surfGroup = qtw.QComboBox()
        for combo in [self.bulkGroup, self.surfGroup]:
            combo.setFont(smallText)
            combo.setSizePolicy(qtw.QSizePolicy.Fixed,
                                qtw.QSizePolicy.Preferred)
            combo.setMinimumWidth(self.aB.width())
        
        allWidgs.extend([self.bulkGroup, self.surfGroup])
        
        # Button to input extra symmetry operations for bulk
        self.bulk_3d_sym = qtw.QPushButton('Extra bulk\nsymmetry')
        self.bulk_3d_sym.setFont(AllGUIFonts().buttonFont)
        self.bulk_3d_sym.setSizePolicy(qtw.QSizePolicy.Fixed,
                                           qtw.QSizePolicy.Preferred)
        self.bulk_3d_sym.setMinimumWidth(self.aB.width())
        allWidgs.append(self.bulk_3d_sym)
        
        
        # Text and button for reduction to highest symmetry                     # TODO: this does not really work nicely. Also, its layout is not inserted (commented out below)
        self.highSymTxt = qtw.QLabel('The lattices could have higher symmetry')
        self.highSymBut = qtw.QPushButton('Reduce to higher symmetry')
        for widg in [self.highSymTxt, self.highSymBut]:
            widg.setFont(smallText)
            widg.setSizePolicy(qtw.QSizePolicy.Fixed, qtw.QSizePolicy.Preferred)
        self.highSymBut.setMinimumWidth(self.highSymBut.width()*0.4)
        self.highSymTxt.setMinimumWidth(self.highSymBut.width()*0.4)
        
        # Superstructure input: dropdown and 2x2 matrix
        superstructLab = qtw.QLabel('Reconstruction periodicity')
        superstructLab.setFont(labFont)
        woodsLab = qtw.QLabel("Wood's:")
        woodsLab.setFont(smallText)
        woodsLab.setSizePolicy(qtw.QSizePolicy.Fixed, qtw.QSizePolicy.Preferred)
        self.woods = qtw.QComboBox()
        self.woods.setEditable(True)
        self.woods.setInsertPolicy(qtw.QComboBox.InsertAlphabetically)
        self.woods.setFont(smallText)
        self.woods.setSizePolicy(qtw.QSizePolicy.Fixed,
                                 qtw.QSizePolicy.Preferred)
        self.woods.setSizeAdjustPolicy(qtw.QComboBox.AdjustToContents)
        self.woods.setMinimumContentsLength(12)
        
        superlatticeLab = qtw.QLabel('M = ')
        superlatticeLab.setFont(labFont)
        parFont = qtg.QFont(labFont)
        parFont.setPointSize(30)
        superlatticeLeftP = qtw.QLabel('(')
        superlatticeLeftP.setFont(parFont)
        superlatticeRightP = qtw.QLabel(')')
        superlatticeRightP.setFont(parFont)
        self.superlattice = np.array([[qtw.QLineEdit(), qtw.QLineEdit()],
                                      [qtw.QLineEdit(), qtw.QLineEdit()]])
        intValidator = qtg.QIntValidator()
        for mij in self.superlattice.ravel():
            mij.setFont(labFont)
            mij.setMaximumWidth(35)
            mij.setSizePolicy(qtw.QSizePolicy.Fixed, qtw.QSizePolicy.Preferred)
            mij.setValidator(intValidator)
        
        allWidgs.extend([superstructLab, self.woods, superlatticeLab,
                         *self.superlattice.ravel()])
        
        # maximum LEED energy
        emaxLab = qtw.QLabel('Max. Energy:')
        self.eMax = qtw.QLineEdit()
        for ctrl in [emaxLab, self.eMax]:
            ctrl.setSizePolicy(qtw.QSizePolicy.Fixed, qtw.QSizePolicy.Preferred)
            ctrl.setFont(labFont)
        self.eMax.setMaximumWidth(70)
        allWidgs.extend([emaxLab,self.eMax])
        
        # 'Done' and 'Cancel' buttons
        self.doneBut = qtw.QPushButton('&Done')
        self.cancelBut = qtw.QPushButton('&Cancel')
        for but in [self.doneBut, self.cancelBut]:
            but.setFont(AllGUIFonts().buttonFont)
            but.setSizePolicy(qtw.QSizePolicy.Fixed, qtw.QSizePolicy.Preferred)
        allWidgs.extend([self.doneBut, self.cancelBut])
        
        # status bar
        self.sBar = qtw.QStatusBar()
        self.sBar.setSizeGripEnabled(False)
        
        #------- Make sure all styles are up to date ------#
        [widg.ensurePolished() for widg in allWidgs]
        
        #------- Place controls in layouts -------#
        lattsLay = qtw.QGridLayout()
        lattsLay.setSpacing(5)
        lattsLay.setContentsMargins(5, 5, 5, 5)  # left, top, right, bottom
        lattsLay.setSizeConstraint(qtw.QLayout.SetFixedSize)
        lattsLay.addWidget(lattInputLabel, 0, 1, 1, 2)
        lattsLay.addWidget(self.bulkInput, 1, 1, 1, 1)
        lattsLay.addWidget(self.surfInput, 1, 2, 1, 1)
        lattsLay.addWidget(latShapeLabel, 2, 0, 1, 1)
        lattsLay.addWidget(self.bulkShape, 2, 1, 1, 1)
        lattsLay.addWidget(self.surfShape, 2, 2, 1, 1)
        lattsLay.addWidget(aLab, 3, 0, 1, 1)
        lattsLay.addWidget(self.aB, 3, 1, 1, 1)
        lattsLay.addWidget(self.aS, 3, 2, 1, 1)
        lattsLay.addWidget(bLab, 4, 0, 1, 1)
        lattsLay.addWidget(self.bB, 4, 1, 1, 1)
        lattsLay.addWidget(self.bS, 4, 2, 1, 1)
        lattsLay.addWidget(alphaLab, 5, 0, 1, 1)
        lattsLay.addWidget(self.alphaB, 5, 1, 1, 1)
        lattsLay.addWidget(self.alphaS, 5, 2, 1, 1)
        lattsLay.addWidget(groupLab, 6, 0, 1, 1)
        lattsLay.addWidget(self.bulkGroup, 6, 1, 1, 1)
        lattsLay.addWidget(self.surfGroup, 6, 2, 1, 1)
        lattsLay.addWidget(self.bulk_3d_sym, 7, 1, 1, 1)
        
        hSymLay = qtw.QVBoxLayout()
        hSymLay.setSpacing(4)
        hSymLay.setContentsMargins(0, 10, 0, 10)
        hSymLay.addWidget(self.highSymTxt)
        hSymLay.addWidget(self.highSymBut)
        hSymLay.addStretch(1)
        
        lattsLay.addLayout(hSymLay, 8, 0, 1, 3)
        
        for lab in [lattInputLabel, self.bulkInput, self.surfInput]:
            lattsLay.setAlignment(lab, qtc.Qt.AlignCenter)
        for lab in [*lattLabels, groupLab]:
            lattsLay.setAlignment(lab, qtc.Qt.AlignRight | qtc.Qt.AlignVCenter)
        
        matrixLay = qtw.QGridLayout()
        matrixLay.setSpacing(2)
        matrixLay.setContentsMargins(0, 0, 0, 0)
        [matrixLay.addWidget(self.superlattice[i,j], i, j+1)
         for i in range(2) for j in range(2)]
        outerMatrixLay = qtw.QHBoxLayout()
        outerMatrixLay.addWidget(superlatticeLab)
        outerMatrixLay.addWidget(superlatticeLeftP)
        outerMatrixLay.addLayout(matrixLay)
        outerMatrixLay.addWidget(superlatticeRightP)
        for txt in [superlatticeLab, superlatticeLeftP, superlatticeRightP]:
            outerMatrixLay.setAlignment(txt, qtc.Qt.AlignCenter)
        outerMatrixLay.addStretch(1)
        
        woodsLay = qtw.QHBoxLayout()
        woodsLay.addWidget(woodsLab)
        woodsLay.addWidget(self.woods)
        woodsLay.addStretch(1)
        
        superstructLay = qtw.QVBoxLayout()
        superstructLay.addWidget(superstructLab)
        superstructLay.addLayout(woodsLay)
        superstructLay.addLayout(outerMatrixLay)
        superstructLay.addStretch(1)
        
        eLay = qtw.QHBoxLayout()
        eLay.setSpacing(4)
        eLay.setContentsMargins(5, 5, 5, 5)
        eLay.addWidget(emaxLab)
        eLay.addWidget(self.eMax)
        eLay.addStretch(1)
        eLay.setAlignment(emaxLab, qtc.Qt.AlignRight | qtc.Qt.AlignVCenter)
        
        botButsLay = qtw.QHBoxLayout()
        botButsLay.addWidget(self.sBar)
        botButsLay.addWidget(self.doneBut)
        botButsLay.addWidget(self.cancelBut)
        
        # Place the layouts in parent layouts, and assemble 
        # the layout of the QDialog
        latticeInputsLay = qtw.QVBoxLayout()
        latticeInputsLay.addLayout(lattsLay)
        latticeInputsLay.addStretch(1)
        
        diagLay = qtw.QGridLayout()
        diagLay.addLayout(latticeInputsLay, 0, 0, 2, 1)
        diagLay.addLayout(superstructLay, 0, 2, 1, 1)
        #diagLay.addLayout(hSymLay, 1, 2, 1, 1)
        diagLay.addLayout(eLay, 1, 2, 1, 1)
        diagLay.addLayout(botButsLay, 2, 0, 1, 3)
        diagLay.setColumnMinimumWidth(1, 15)
        diagLay.setAlignment(eLay, qtc.Qt.AlignLeft | qtc.Qt.AlignVCenter)
        diagLay.setSizeConstraint(qtw.QLayout.SetFixedSize)
        
        self.setLayout(diagLay)
        
        #------ Connect signals ------#
        self.connectSignals()
    
    def connectSignals(self):
        self.bulkInput.toggled.connect(self.bulkToggled)
        for shape in [self.bulkShape, self.surfShape]:
            shape.activated.connect(self.lattShapeChanged)
        for group in [self.bulkGroup, self.surfGroup]:
            group.activated.connect(self.groupChanged)
        self.woods.activated.connect(self.woodsSelected)
        self.woods.lineEdit().editingFinished.connect(self.woodsSelected)
        for tbox in [*self.bulkParams, *self.surfParams]:
            tbox.textEdited.connect(self.lattParamsChanged)
            tbox.editingFinished.connect(self.updateLatticeParameters)
        for el in self.superlattice.ravel():
            el.textEdited.connect(self.matrixChanged)
        for but in [self.doneBut, self.cancelBut]:
            but.clicked.connect(self.onButtonPressed)
            # The next two lines prevent the dialog to be closed automatically
            # by firing one of the buttons if the user presses enter/return
            # while focus is on widgets other than the buttons. The
            # reimplementation of keyPressEvent takes care of the behavior on 
            # pressing enter/return
            but.setAutoDefault(False)
            but.setDefault(False)
        self.bulk_3d_sym.clicked.connect(self.on_bulk_3d_pressed)
        self.highSymBut.clicked.connect(self.onReducePressed)
    
    # The next 7 functions are the handler of control changes
    
    def bulkToggled(self, checked):
        bulk = [*self.bulkParams, self.bulkShape]
        surf = [*self.surfParams, self.surfShape]
        if checked:
            enable = bulk
            disable = surf
        else:
            enable = surf
            disable = bulk
        [box.setEnabled(True) for box in enable]
        [box.setEnabled(False) for box in disable]
        
        self.updateLatticeRestrictions()
    
    def lattShapeChanged(self, signal=None):
        #--- Set the appropriate relations depending on the selection ---#
        if self.bulkInput.isChecked():
            shape = self.bulkShape
            params = self.bulkParams  #[a, b, alpha]
        else:
            shape = self.surfShape
            params = self.surfParams
        
        shape = shape.currentText()
        (a, b, alpha) = tuple(params)
        
        if shape in ['Square', 'Rhombic', 'Hexagonal']:
            b.setEnabled(False)
            b.setText(a.text())
        else:
            b.setEnabled(True)
        
        if shape in ['Square', 'Rectangular', 'Hexagonal']:
            alpha.setEnabled(False)
            if shape in ['Square', 'Rectangular']:
                alpha.setText('90' + degrees)
            else:
                alpha.setText('120' + degrees)
        else:
            alpha.setEnabled(True)
        
        #--- react on the changes of lattice parameters ---#
        self.lattParamsChanged()
        
        #--- update the lattice shape of the other ---#
        if self.bulkInput.isChecked():
            other = 'surf'
            self.surfShape.setCurrentText(self.surfLatt.cell_shape)
        else:
            other = 'bulk'
            self.bulkShape.setCurrentText(self.bulkLatt.cell_shape)
        
        #--- and update the dependent controls ---#
        if signal is not None:
            # method has been called as a result of a user action
            other = 'both'
        self.updateLatticeParameters(other)
        self.updateWoodsList()
        self.updateGroups()
    
    def lattParamsChanged(self, newText=None):
        if newText is not None:
            # The function is called as a result of a user interaction
            #
            # check that the input is valid (the user could be typing text)
            # ADD MATH PARSER HERE
            reNumber = re.compile(r'^\s*\d+([.]\d+)?\s*[\u212b\u00b0]?$')
            number = reNumber.match(newText)
            if number is not None:  # input is a number (with unit)
                self.updateLatticeRestrictions()
        
        if self.bulkInput.isChecked():
            type = 'bulk'
        else:
            type = 'surf'
        
        basis = self.getBasis(type)
        if basis is not None:  # edit value is OK
            self.updateLattices(basis)
            self.updateSymmReduction()
    
    def matrixChanged(self, newEl=None):
        if newEl is None:
            newEl = self.superlattice[0, 0].text()
        #--- check that the input is an integer ---#
        v = self.superlattice[0, 0].validator()
        valid = v.validate(newEl, 1)[0] == qtg.QValidator.Acceptable
        
        if valid: # update
            texts = [el.text() for el in self.superlattice.ravel()]
            m = np.array([int(text) for text in texts]).reshape(2, 2)
            #--- check that the input gives a non-singular matrix ---#
            if np.linalg.det(m) != 0:
                self.supMatrix = m.reshape(2, 2)
            
                self.updateLattices()
                self.updateLatticeParameters()
                self.updateLatticeRestrictions()
                self.updateWoods()
                self.sBar.clearMessage()
            else:
                self.sBar.showMessage('Matrix is singular!', 2000)
    
    def woodsSelected(self, woodsIdx=None):
        if woodsIdx is None:
            # the call to the function came from the editingFinished of the
            # QLineEdit
            woods = self.woods.lineEdit().text()
        else:
            # the call came by a user click on a QComboBox item
            woods = self.woods.itemText(woodsIdx)
        woods = self.fixWoods(woods)
        if woods is not None:
            m = self.woodL.woodsToMatrix(woods, self.bulkLatt.basis)
            if m is None:  # woods gives an incommensurate lattice
                self.sBar.showMessage(woods + ' invalid: Lattice is '
                                      + 'incommensurate.', 7000)
                self.woods.removeItem(self.woods.currentIndex())
                self.updateWoods()
            else:
                self.updateSuperlattice(m)
                self.matrixChanged()
        else:
            self.sBar.showMessage("Invalid Wood's syntax.", 1000)
    
    def groupChanged(self, newGroup=None):
        if newGroup is not None:
            self.updateLatticeGroups()
    
    def onButtonPressed(self, checked):
        btn = self.sender()
        if btn == self.doneBut:
            # Automatically make lattices high symmetry before exiting          # TODO: this screws up for the extra bulk operations. The high-symmetry button should rather appear after edits.
            # self.highSymBut.click()
            params = self.packParameters()
            if params is not None:
                self.dialogEdited.emit(params)
                self.accept()
        elif btn == self.cancelBut:
            self.dialogEdited.emit(self.oldParams)
            self.reject()
    
    def onReducePressed(self, checked):
        tBulk = self.bulkLatt.make_high_symmetry()
        tSurf = self.surfLatt.make_high_symmetry()
        m = np.dot(np.dot(tSurf, self.supMatrix), np.linalg.inv(tBulk))
        #  now make entries integer
        m = m.round().astype(int)
        
        self.updateLatticeParameters()
        self.updateSuperlattice(m)
        self.updateLatticeRestrictions()
        self.updateWoodsList()
        self.updateGroups()
    
    def on_bulk_3d_pressed(self, checked):
        self.bulk_3d_sym_dialog.update_operations(self.bulkLatt)
        if self.bulk_3d_sym_dialog.exec() == qtw.QDialog.Accepted:
            operations = self.bulk_3d_sym_dialog.extra_operations()
            self.bulkLatt.group.screws_glides = (operations,
                                                 self.bulkLatt.cell_shape)
    
    # The next 8 function are explicitly called to handle changes of controls
    
    def updateLattices(self, basis = None):
        #------- Update lattices and their lattice parameters -------#
        # If bulk is selected, the surface parameters will be changed
        # at constant bulk parameters. Otherwise the other way around.
        
        if self.bulkInput.isChecked():
            if basis is None:
                basis = self.bulkLatt.basis
            else:
                bulk = self.bulkLatt
                bulk.basis = basis
            latt = self.surfLatt
            m = self.supMatrix
        else:
            if basis is None:
                basis = self.surfLatt.basis
            else:
                surf = self.surfLatt
                surf.basis = basis
            latt = self.bulkLatt
            m = np.linalg.inv(self.supMatrix)
        
        latt.basis = np.dot(m, basis)
        
    def updateLatticeGroups(self):
        bGroup = self.bulkGroup.currentText()
        sGroup = self.surfGroup.currentText()
        allGroups = self.bulkLatt.group.allGroups
        for (latt, group) in zip([self.bulkLatt, self.surfLatt],
                                 [bGroup, sGroup]):
            if group in allGroups:
                bulk_3d = latt.group.screws_glides
                latt.group = group
                latt.group.screws_glides = (bulk_3d, latt.cell_shape)
        
        # Also, update the options for extra bulk operations
        self.bulk_3d_sym_dialog.update_operations(self.bulkLatt)
        # and deactivate button in case there is nothing to add
        self.bulk_3d_sym.setEnabled(self.bulk_3d_sym_dialog.n_extra_ops)
    
    def updateLatticeParameters(self, which='both'):
        if which == 'bulk':
            params = [self.bulkParams]
            latts = [self.bulkLatt]
        elif which == 'surf':
            params = [self.surfParams]
            latts = [self.surfLatt]
        else:
            params = [self.bulkParams, self.surfParams]
            latts = [self.bulkLatt, self.surfLatt]
        for (ctrls,latt) in zip(params, latts):
            a, b, alpha = latt.lattice_parameters
            ctrls[0].setText('{:.4f}{}'.format(a, angstrom))
            ctrls[1].setText('{:.4f}{}'.format(b, angstrom))
            ctrls[2].setText('{:.1f}{}'.format(alpha, degrees))
    
    def updateSuperlattice(self, m=None):
        if m is None:
            m = self.supMatrix
        m = np.array(m)
        for (el, mij) in zip(self.superlattice.ravel(), m.ravel()):
            el.setText(str(int(mij)))
        self.supMatrix = m
    
    def updateLatticeRestrictions(self):
        self.lattShapeChanged()
    
    def updateWoodsList(self):
        # loads the example list of woods notations
        # for the selected bulk shape
        shape = self.bulkLatt.cell_shape
        self.woods.clear()
        self.woods.addItems(sorted(self.woodL.examples[shape]))
        self.updateWoods()
    
    def updateWoods(self):
        woods = Woods().matrixToWoods(self.supMatrix, self.bulkLatt.basis)
        if woods is None:
            woods = 'None'
        idx = self.woods.findText(woods, flags = qtc.Qt.MatchExactly)
        if idx < 0:  # item is not already present
            self.woods.addItem(woods)
            self.woodL.examples[self.bulkLatt.cell_shape] |= {woods}
        self.woods.setCurrentText(woods)
    
    def updateGroups(self):
        for (ctrl, latt) in zip([self.bulkGroup, self.surfGroup],
                                [self.bulkLatt, self.surfLatt]):
            shape = latt.cell_shape
            group = latt.group
            groups = group.groupsForShape[shape]
            ctrl.clear()
            ctrl.addItems(groups)
            ctrl.setCurrentText(group.group)
        # update the groups in the lattices, in case the 
        # current text and the group do not match
        self.updateLatticeGroups()
    
    def updateSymmReduction(self):
        tBulk = self.bulkLatt.high_symm_transform()
        tSurf = self.surfLatt.high_symm_transform()
        id = [[1, 0], [0, 1]]
        isHighSym = np.array_equal(tBulk, id) and np.array_equal(tSurf, id)
        for widg in [self.highSymBut, self.highSymTxt]:
            if isHighSym:
                widg.hide()
            else:
                widg.show()
    
    # and the next 4 functions are utilities
    
    def packParameters(self):
        reUnits = re.compile(r'(\d+(\.\d+)?)\s*[\u212b\u00b0eV]*') 
        # angstrom, degrees, eV
        eMax = reUnits.match(self.eMax.text())
        if eMax is None:
            return None
        eMax = float(eMax.group(1))
        
        params = {'eMax': eMax,
                  'SUPERLATTICE': self.supMatrix,
                  'surfBasis': self.surfLatt.basis,
                  'surfGroup': self.surfLatt.group,
                  'bulkGroup': self.bulkLatt.group,
                  'bulk3Dsym': self.bulkLatt.group.screws_glides}
        
        return params
    
    def getBasis(self, which = 'bulk'):
        reUnits = re.compile(r'(\d+(\.\d+)?)\s*[\u212b\u00b0]?')
        if which == 'bulk':
            params = [self.aB.text(), self.bB.text(), self.alphaB.text()]
        else:
            params = [self.aS.text(), self.bS.text(), self.alphaS.text()]
        params = [reUnits.match(param) for param in params]
        if None in params:
            return None
        params = [param.group(1) for param in params]
        
        for param in params:
            try:
                float(param)
            except:
                return None
        
        # All entries are floats
        (a, b, alpha) = (float(param) for param in params)
        alpha = np.radians(alpha)
        if a > 0 and b > 0:
            return np.array([[a, 0],[b*np.cos(alpha), b*np.sin(alpha)]])
        return None
    
    def fixWoods(self, woods):
        if woods is None:
            return None
        
        reWoods = re.compile(
            r'''
            ^(?P<prefix> ((?![rR])[\sa-zA-Z])*)     # * Multi-letter
                                                    #   prefix -> will turn 
                                                    #   into p or c
            \(?                                     # * Optional open 
                                                    #   parenthesis
            (?P<g1int>\d+)?[*]?                     # * Dir1, integer part,
                                                    #   optional (+ optional 
                                                    #   asterisk for 
                                                    #   multiplication)
            (?:[sqSQ]{0,2}[rR\u221a][tT]?\(?(?P<g1rt>\d+)\)?)?
                                                    # * Dir1, sqrt (and 
                                                    #   similar) part
            [x\u00d7]                               # * One times character
            (?P<g2int>\d+)?[*]?                     # * Dir2, integer part,
                                                    #   optional (+ optional 
                                                    #   asterisk for 
                                                    #   multiplication)
            (?:[sqSQ]{0,2}[rR\u221a][tT]?\(?(?P<g2rt>\d+)\)?)?
                                                    # * Dir2, sqrt (and 
                                                    #   similar) part
            \)?                                     # * Optional closing
                                                    #   parenthesis
            ([rR](?P<alpha>\d+(\.\d+)?)\u00b0?)?$   # * rotation block
            ''', re.VERBOSE)
        
        m = reWoods.match(woods)
        if m is None:
            return None
        
        groups = m.groupdict()
        if 'c' in groups['prefix']:  # force prefix to be p or c
            groups['prefix'] = 'c'
        else:
            groups['prefix'] = 'p'
        
        for key in ['g1int' ,'g2int']:
            # replace the empty integer parts with 1
            if groups[key] is None:
                groups[key] = 1
            else:
                groups[key] = float(groups[key])
        
        toFormat = []
        dirs = [(groups['g1int'], groups['g1rt']),
                (groups['g2int'], groups['g2rt'])]
        for (integral, rt) in dirs:
            # Now combine integer and square root parts
            if rt is not None:
                sq = float(rt)
                (intSq, rad) = self.woodL.squareToProdOfSquares(sq)
                integral *= np.round(np.sqrt(intSq))
            else:
                rad = 1
            
            dirStr = ''
            if rad > 1:
                dirStr = '\u221a' + str(int(rad))
            if not dirStr:
                dirStr = str(int(integral))
            else:
                if integral > 1:
                    dirStr = str(int(integral)) + dirStr
            toFormat.append(dirStr)
        
        woodsFixed = '{}({})'.format(groups['prefix'],'\u00d7'.join(toFormat))
        
        #finally handle the angle
        if groups['alpha'] is not None:
            alpha = float(groups['alpha'])
            cA = np.abs(np.cos(np.radians(alpha)))
            if cA > 1e-3 and 1-cA > 1e-3:
                # angle is neither 0, 90, nor 180
                woodsFixed += 'R{:.1f}'.format(alpha) + degrees
        
        return woodsFixed
    
    def initControlValues(self):
        #---- Explicitly set some values ----#
        self.bulkInput.setChecked(True)
        # always set bulk checked at session startup
        # -> might remember choice later
        self.bulkShape.setCurrentText(self.bulkLatt.cell_shape)
        self.surfShape.setCurrentText(self.surfLatt.cell_shape)
        if 'eMax' in self.oldParams:
            e = self.oldParams['eMax']
        else:
            e = 700
        self.eMax.setText(str(e) + ' eV')
        for widg in [self.highSymTxt, self.highSymBut]:
            widg.hide()
        
        #---- Call some updating functions ----#
        self.updateLatticeParameters()
        self.updateSuperlattice()
        self.updateLatticeRestrictions()
        self.updateWoodsList()
        self.updateGroups()
    
    # Re-implement keyPressEvent so that pressing enter or return when Done or
    # Cancel have keyboard focus fires them
    
    def keyPressEvent(self, event):
        if event.key() in [qtc.Qt.Key_Enter or qtc.Qt.Key_Return]:
            focusWidg = self.focusWidget()
            if isinstance(focusWidg, qtw.QLineEdit):
                return
            elif (isinstance(focusWidg, qtw.QComboBox)
                  and focusWidg.isEditable()):
                return
            elif focusWidg == self.cancelBut:
                self.cancelBut.click()
                return
            self.doneBut.click()
            return
        #super(NewFileDialog, self).keyPressEvent(event)
        super().keyPressEvent(event)
