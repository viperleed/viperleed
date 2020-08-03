"""
===============================================
      ViPErLEED Graphical User Interface
===============================================
 *** module guilib.leedsim.ExportCSVDialog ***

Created: 2020-01-11
Author: Michele Riva

"""

from fractions import Fraction

import numpy as np
import PyQt5.QtCore as qtc
import PyQt5.QtGui as qtg
import PyQt5.QtWidgets as qtw

import guilib as gl


class ExportCSVDialog(qtw.QDialog):
    exportSelected = qtc.pyqtSignal(list) # the list of lines formatted for export
    
    def __init__(self, leed, parent=None):
        if not isinstance(leed, gl.LEEDPattern):
            raise
        #super(ExportCSVDialog, self).__init__(parent)
        super().__init__(parent)
        self.leed = leed
        self.setWindowModality(qtc.Qt.WindowModal)
        flags = self.windowFlags()
        flags &= ~qtc.Qt.WindowCloseButtonHint #disable close button
        self.setWindowFlags(flags)
        self.setWindowTitle('Select export data')
        
        self.compose()
        self.open()
    
    def compose(self):
        font = gl.AllGUIFonts().labelFont
        if self.leed.nDoms > 1:  # more than one domain
            domColors = self.leed.domColors
        else:
            domColors = [(0, 0, 0)]  # black
        
        # textbox for optional name to place in the header
        stuctNameLab = qtw.QLabel('Structure name (optional):')
        self.structName = qtw.QLineEdit('')
        self.structName.setFont(font)
        self.structName.setToolTip('This name will be included in the header'
                                   ' of the exported file')
        
        nameLay = qtw.QVBoxLayout()
        nameLay.addWidget(stuctNameLab)
        nameLay.addWidget(self.structName)
        nameLay.addStretch(1)
        
        txt = qtw.QLabel('Select which domains to export')
        txt.setFont(font)
        
        # Radio buttons to select whether all domains should be exported
        # or only some specific ones
        all = qtw.QRadioButton('All')
        visible = qtw.QRadioButton('Visible domains')
        selection = qtw.QRadioButton('Other')
        self.exportRadio = (all, visible, selection)
        
        for radio in self.exportRadio:
            radio.setFont(font)
        self.exportRadio[0].setChecked(True)
        
        # And 'Export' and 'Cancel' buttons
        self.doneBut = qtw.QPushButton('Export')
        self.cancelBut = qtw.QPushButton('Cancel')
        for but in [self.doneBut, self.cancelBut]:
            but.setFont(font)
            but.setSizePolicy(qtw.QSizePolicy.Fixed, qtw.QSizePolicy.Fixed)
        butsLay = qtw.QHBoxLayout()
        butsLay.addStretch(1)
        butsLay.addWidget(self.doneBut)
        butsLay.addWidget(self.cancelBut)
        
        # Prepare also as many tick boxes as there are domains
        self.domTicks = [qtw.QCheckBox() for dom in range(self.leed.nDoms)]
        for dom, (color, tick) in enumerate(zip(domColors, self.domTicks)):
            tick.setFont(font)
            tick.setText('Dom. %d' % (dom + 1))
            p = tick.palette()  # for coloring text the same as the domains
            p.setColor(qtg.QPalette.WindowText, qtg.QColor.fromRgbF(*color))
            tick.setPalette(p)
            tick.setSizePolicy(qtw.QSizePolicy.Fixed, qtw.QSizePolicy.Preferred)
            tick.hide()
        
        # and arrange them in a QGridLayout, in a similar manner as the 
        # matrices in the MatricesPopup
        nDoms = len(self.domTicks)  # this is 1, 2, 3, 4, 6, 8, or 12
        if nDoms in range(1, 4):  # 1 -- 3
            nRows = 1
        elif nDoms in range(4, 10):  # 4 -- 8
            nRows = 2
        else:
            nRows = 3
        nCols = int(nDoms/nRows)
        
        ticksLay = qtw.QGridLayout()
        [ticksLay.addWidget(tick, index//nCols, index % nCols)
         for (index, tick) in enumerate(self.domTicks)]
        
        # And build the dialog by putting the widgets in a layout
        diagLay = qtw.QVBoxLayout()
        diagLay.addLayout(nameLay)
        diagLay.addWidget(txt)
        [diagLay.addWidget(radio) for radio in self.exportRadio]
        diagLay.addLayout(ticksLay)
        diagLay.addStretch(1)
        diagLay.addLayout(butsLay)
        # diagLay.setSizeConstraint(QLayout.SetFixedSize)
        
        self.setLayout(diagLay)
        
        # In case there is only one domain, there's no reason to show all these
        # controls, except for the optional name and the buttons
        if nDoms == 1:
            children = gl.get_all_children_widgets(diagLay)
            keepVisible = (gl.get_all_children_widgets(nameLay)
                           | gl.get_all_children_widgets(butsLay))
            for child in children - keepVisible:
                child.hide()
        
        # Finally connect some signals
        self.connectControls()
    
    def connectControls(self):
        self.exportRadio[2].toggled.connect(self.exportSomeTriggered)
        for but in [self.doneBut, self.cancelBut]:
            but.clicked.connect(self.onButtonPressed)
    
    def exportSomeTriggered(self, checked=None):
        if checked is None:
            return
        
        if checked:
            [tick.show() for tick in self.domTicks]
        else:
            [tick.hide() for tick in self.domTicks]
    
    def onButtonPressed(self, checked):
        btn = self.sender()
        if btn == self.doneBut:
            if self.exportRadio[0].isChecked():
                # export all domains
                [tick.setChecked(True) for tick in self.domTicks]
            elif self.exportRadio[1].isChecked():
                # export visible domains
                # TODO: read in which ones are visible. Probably an attribute 
                # in self.leed
                # for now behaves the same as the previous
                [tick.setChecked(True) for tick in self.domTicks]
            formatted = self.formatBeams()
            self.exportSelected.emit(formatted)
            self.accept()
        elif btn == self.cancelBut:
            self.reject()
    
    def formatBeams(self):
        domains = [i for (i, tick) in enumerate(self.domTicks)
                   if tick.isChecked()]
        
        beams = self.leed.get_equivalentSpots(domains)
        fractions, groups, overlaps, extinctDoms = zip(*beams)
        
        # Process the data to get them in a better format for packing each
        # line; also find the correct number of characters for each column
        hk = []
        gg = []
        dd = []
        # lengths is a dictionary that will contain the maximum
        # number of characters appearing in each part of the columns and
        # it is used to figure out the correct padding
        # (int*** are the lengths of the integer parts of floats)
        lengths = {'numL': [], 'denL': [], 'fractL': [], 'intHKL': [],
                   'intGL': [], 'domsL': []}
        
        for fract in fractions:
            ind = fract.partition(', ')
            h = Fraction(ind[0])
            k = Fraction(ind[2])
            hk.append((h, k))
            g = np.dot((h, k), self.leed.bulkR.basis)
            gg.append(g)
            
            lengths['numL'].extend([len(str(h.numerator)), 
                                    len(str(k.numerator))]
                                   )
            lengths['denL'].extend([len(str(h.denominator)),
                                    len(str(k.denominator))]
                                    )
            lengths['fractL'].extend([len(str(h)), len(str(k))])
            lengths['intHKL'].extend([len(str(int(float(h)//1))),
                                      len(str(int(float(k)//1)))])
            lengths['intGL'].extend([len(str(int(g[0]//1))),
                                     len(str(int(g[1]//1)))])
        
        for (overlap, extinct) in zip(overlaps, extinctDoms):
            overlapStr = [str(dom) for dom in overlap]
            for e, (s, ext) in enumerate(zip(overlapStr, extinct)):
                if ext:
                    overlapStr[e] = '(%s)'%(s)
            newdd = '+'.join(dom for dom in overlapStr)
            dd.append(newdd)
            lengths['domsL'].append(len(newdd))
        
        for key in lengths.keys():
            lengths[key] = max(lengths[key])
        lengths['groupL'] = max(len(str(np.max(groups))), len('group'))
        lengths['domsL'] = max(lengths['domsL'], len('domain(s)'))
        
        toExport = [*self.formatHeader(lengths)]
        
        for (fract, vecR, group, doms) in zip(hk, gg, groups, dd):
            line = []
            h = fract[0]
            k = fract[1]
            fractHK = []
            for hk in [h, k]:
                lPad = lengths['numL'] - len(str(hk.numerator))
                rPad = lengths['fractL'] - lPad - len(str(hk))
                fractHK.append(' '*lPad + str(hk) + ' '*rPad)
            line.append('({})'.format('|'.join(fractHK)))
            # 2) floating indices
            h = float(h)
            k = float(k)
            for hk in [h, k]:
                lPad = lengths['intHKL'] - len(str(hk).partition('.')[0])
                line.append(' '*lPad + '{:.5f}'.format(hk))
            # 3) gx gy
            for gg in vecR:
                lPad = lengths['intGL'] - len(str(gg).partition('.')[0])
                line.append(' '*lPad + '{:.5f}'.format(gg))
            # 4) group index
            line.append(' '*(lengths['groupL'] - len(str(group))) + str(group))
            # 5) overlapping domains
            line.append(' '*(lengths['domsL'] - len(doms)) + doms)
            
            toExport.append(','.join(line) + ',')
        
        return toExport
    
    def formatHeader(self, lengths):
        # 1st Header (taken care of by LEED_GUI)
        # * String saying where was this exported from
        # * include filename and path of files used for exporting, if saved
        #
        # 2nd Header:
        # * optional structure name
        # * max energy
        # * bulk shape (group); lattice parameters and basis
        # * surf shape (group); lattice parameters
        # * total number of domains, and number of domains exported
        # * for each exported domain: basis and superlattice matrix
        # * finally an uncommented header line for the columns
        
        header = []
        txt = self.structName.text()
        if txt:
            header.append('# Structure: ' + txt + '\n#')
        header.append(
            '# Max. LEED Energy: {:.1f} eV'.format(self.leed.maxEnergy),
            )
        
        for (txt, lattR) in zip(['Bulk', 'Surface'],
                                [self.leed.bulkR, self.leed.surfR]):
            basis = lattR.reciprocal_basis()
            shape = lattR.type
            group = lattR.group.group
            latticeText = self.formatLatticeBasis(basis, shape)
            forHeader = '# {} lattice: {} ({}); {}'.format(txt, shape, group,
                                                           latticeText)
            header.append(forHeader)
            if txt == 'Bulk':
                header.append('#' + ' '*8 + self.formatBasisVectors(basis))
        
        exportDoms = [dom for (dom, tick) in enumerate(self.domTicks)
                          if tick.isChecked()]
        matrices = self.leed.domSuperlattices[exportDoms]
        
        header.append('#\n# Exporting beams from %d domain(s) out of %d '
                      'symmetry-equivalent ones'
                      % (len(matrices), self.leed.nDoms))
        
        basis = self.leed.bulkR.reciprocal_basis()
        for (dom, m) in zip(exportDoms, matrices):
            supBas = np.dot(m, basis)
            domtxt = ' Domain %d - ' % (dom + 1)
            header.append('#' + domtxt + self.formatBasisVectors(supBas))
            header.append('#' + ' '*len(domtxt) + 'Superlattice: '
                          + np.array2string(m, 
                                            separator=' ').replace('\n',''))
        
        # Prepare the headers of the columns
        # 0) Some description
        header.append(
              '#\n'
              '# * h and k are the surface Miller indices of each LEED spot;'
              ' both fractional'
              '\n#   and floating-point versions are provided.'
              '\n# * gx and gy are the horizontal and vertical components of'
              ' reciprocal lattice'
              '\n#   vectors in AA^(-1), and include a factor of 2*pi.'
              '\n# * Beams with the same absolute value of "group" are'
              ' symmetry equivalent (at'
              '\n#   normal incidence); extinct spots have a negative "group"'
              ' index.'
              '\n# * The domains contributing to each spot are listed in'
              ' "domain(s)". Domains'
              '\n#   contributing with glide-extinct spots are reported in'
              ' parentheses.')
        
        # 1) string version of the (h k) indices
        numL = lengths['numL']
        denL = lengths['denL']
        fractL = lengths['fractL']
        intHKL = lengths['intHKL']
        intGL = lengths['intGL']
        groupL = lengths['groupL']
        
        lPad = numL
        if numL == fractL:  # only integer spots
            lPad -= 1
        rPad = fractL - lPad -1
        
        colHeaders = ['\n(' + ' '*lPad + 'h' + ' '*rPad
                      + '|' + ' '*lPad + 'k' + ' '*rPad + ')']
        
        # 2) float version of the (h k) indices
        for hk in ['h', 'k']:
            colHeaders.append(' '*(1 + (intHKL + 5)//2)
                              + hk
                              + ' '*(4 + intHKL - (intHKL + 5)//2))
        
        # 3) float gx gy
        intL = len(str(np.ceil(self.leed.doms[0].max())))
        for gg in ['gx', 'gy']:
            colHeaders.append(' '*((intGL + 5)//2)
                              + gg
                              + ' '*(4 + intGL - (intGL + 5)//2))
        
        # 4) domain group and overlapping domains
        colHeaders.extend([' '*max(0, groupL - 5) + 'group',
                           ' '*max(0, 2*len(matrices) - 9) + 'domain(s)'])
        header.append(','.join(colHeaders) + ',')
        
        return header
    
    def formatLatticeBasis(self, basis, shape, skipMatrix=False):
        (a, b) = tuple(np.linalg.norm(basis, axis=1))
        alpha = np.degrees(np.arccos(np.dot(basis[0], basis[1])/(a*b)))
        
        txt = 'a = {:.4f} AA'.format(a)
        txtLst = [txt]
        
        if shape not in ['Square', 'Hexagonal', 'Rhombic']:
            txtLst.append('b = {:.4f} AA'.format(b))
        if shape not in ['Square', 'Rectangular']:
            txtLst.append('alpha = {:.2f} deg'.format(alpha))
        
        return '; '.join(txtLst)
    
    def formatBasisVectors(self, basis):
        txt = ('Basis: '
               + 'a = {}'.format(np.array2string(basis[0], precision=4,
                                                 floatmode='fixed', 
                                                 separator=' ',
                                                 suppress_small=True))
               + '; '
               + 'b = {}'.format(np.array2string(basis[1], precision=4,
                                                 floatmode='fixed', 
                                                 separator=' ',
                                                 suppress_small=True))
               )
        return txt
