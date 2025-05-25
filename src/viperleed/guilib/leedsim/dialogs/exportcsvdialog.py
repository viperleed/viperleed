"""
===============================================
      ViPErLEED Graphical User Interface
===============================================
 *** module guilib.leedsim.ExportCSVDialog ***

Created: 2020-01-11
Author: Michele Riva

"""

from quicktions import Fraction

import numpy as np
import PyQt5.QtCore as qtc
import PyQt5.QtGui as qtg
import PyQt5.QtWidgets as qtw

from viperleed.guilib.leedsim.classes.leedpattern import LEEDPattern            # TODO: maybe the old one?
from viperleed.guilib.widgetslib import AllGUIFonts
from viperleed.guilib.widgetslib import get_all_children_widgets


class ExportCSVDialog(qtw.QDialog):
    # this signal contains the arguments to be passed on to export_csv
    exportSelected = qtc.pyqtSignal(dict)


    def __init__(self, leed, parent=None):
        if not isinstance(leed, LEEDPattern):
            raise
        #super(ExportCSVDialog, self).__init__(parent)
        super().__init__(parent)
        self.leed = leed
        self.setWindowModality(qtc.Qt.WindowModal)
        flags = self.windowFlags()
        flags &= ~qtc.Qt.WindowCloseButtonHint  # disable close button
        self.setWindowFlags(flags)
        self.setWindowTitle('Select export data')

        self.compose()
        self.connectControls()  # connect some signals
        self.open()

    def compose(self):
        font = AllGUIFonts().labelFont
        if self.leed.n_domains > 1:  # more than one domain
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
        self.domTicks = [qtw.QCheckBox() for dom in range(self.leed.n_domains)]
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
            children = get_all_children_widgets(diagLay)
            keepVisible = (get_all_children_widgets(nameLay)
                           | get_all_children_widgets(butsLay))
            for child in children - keepVisible:
                child.hide()

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
            params = self.pack_export_params()
            self.exportSelected.emit(params)
            self.accept()
        elif btn == self.cancelBut:
            self.reject()

    def pack_export_params(self):
        """
        Prepares the parameters that define what is to be exported
        """
        params = dict()
        params['domains'] = ([i for (i, tick) in enumerate(self.domTicks)
                             if tick.isChecked()],)
        params['name'] = self.structName.text()

        # Other parameters (leed, source, file names) are taken care of
        # by LEEDPatternSimulator

        return params

