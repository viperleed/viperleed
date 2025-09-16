"""Module rotationblock of viperleed.gui.leedsim.widgets.

Defines the RotationBlock widget, the user input for the visual
orientation of the displayed LEED pattern, as well as its
corresponding real-space version.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2020-01-12'
__license__ = 'GPLv3+'

import numpy as np
import PyQt5.QtCore as qtc
import PyQt5.QtWidgets as qtw

from viperleed.gui.widgets.lib import AllGUIFonts
from viperleed.gui.widgets.spinboxes import TextBoxWithButtons


class RotationBlock(TextBoxWithButtons):
    def __init__(self, parent=None):
        params = {
            'labelText': '&Rotation (\u00b0)',
            'textBoxText': '\u2014',
            'botButText': '\u21bb', #CW
            'topButText': '\u21ba', #CCW
            'parent': parent,
            'textBoxTip': 'Change rotation of lattices and pattern. '\
                          'Positive angles are counterclockwise',
            'topButTip': 'Rotate lattices and pattern 10° clockwise'\
                         'Hold Ctrl down for finer control',
            'botButTip': 'Rotate lattices and pattern 10° counterclockwise'\
                         'Hold Ctrl down for finer control',
            'textBoxWidth': 85
        }

        super().__init__(**params)
        self.setTips(**params)
        self.text.setLimits(-180, 180, 'cyclic')
        self.text.setStep(10, 'add')

        #set aliases for easier identification
        self.ccw = self.topBut
        self.cw = self.botBut

        #and build the bottom part
        self.makeBottomWidget()

    def makeBottomWidget(self):
        # bottomWidget looks like this:
        #
        #    'Align bulk'
        #     [1 0] hor ver
        #     [0 1] hor ver
        #
        # will be composed with a QGridLayout, that is then assigned to self.bottomWidget

        bWLay = qtw.QGridLayout()
        bWLay.setSpacing(0)
        bWLay.setContentsMargins(0, 0, 0, 0)

        # prepare the sub-widgets
        rotBulkLab = qtw.QLabel('Align bulk:')
        oneZeroLab = qtw.QLabel("[1 0]")
        zeroOneLab = qtw.QLabel("[0 1]")
        self.h10 = qtw.QPushButton("Hor.")
        self.v10 = qtw.QPushButton("Ver.")
        self.h01 = qtw.QPushButton("Hor.")
        self.v01 = qtw.QPushButton("Ver.")
        horVerButs = ((self.h10, self.v10), (self.h01, self.v01))

        # add the relevant ones to the list of sub-widgets
        self.subWidgs.extend([self.h10, self.v10, self.h01, self.v01,
                             *np.ravel(horVerButs)])

        # set fonts
        [lab.setFont(AllGUIFonts().smallTextFont)
         for lab in [rotBulkLab,oneZeroLab,zeroOneLab]]
        [but.setFont(AllGUIFonts().smallButtonFont)
         for but in np.ravel(horVerButs)]

        # set size policies and adjust the sizes to account for the new font
        for lab in [rotBulkLab, oneZeroLab, zeroOneLab]:
            lab.setSizePolicy(qtw.QSizePolicy.Fixed, qtw.QSizePolicy.Preferred)
            lab.adjustSize()
        for but in np.ravel(horVerButs):
            but.setSizePolicy(qtw.QSizePolicy.Preferred, qtw.QSizePolicy.Fixed)
            but.setMaximumHeight(
                round(qtw.QPushButton().sizeHint().height()*.84)
                )

        # set tooltips
        tips = ['Rotate lattices and pattern to bring the [1 0] '
                    + 'bulk vector horizontal',
                'Rotate lattices and pattern to bring the [1 0] '
                    + 'bulk vector vertical',
                'Rotate lattices and pattern to bring the [0 1] '
                    + 'bulk vector horizontal',
                'Rotate lattices and pattern to bring the [0 1] '
                    + 'bulk vector vertical']
        for (but, tip) in zip(np.ravel(horVerButs), tips):
            but.setStatusTip(tip)
            but.setToolTip(tip)

        # now add the widgets to the layout
        bWLay.addWidget(rotBulkLab, 0, 0, 1, 3, qtc.Qt.AlignLeft)
        bWLay.addWidget(oneZeroLab, 1, 0, 1, 1,
                        qtc.Qt.AlignHCenter | qtc.Qt.AlignLeft)
        bWLay.addWidget(zeroOneLab, 2, 0, 1, 1,
                        qtc.Qt.AlignHCenter | qtc.Qt.AlignLeft)
        [bWLay.addWidget(horVerButs[i][j], i+1, j+1, 1, 1, qtc.Qt.AlignCenter)
         for i in range(2) for j in range(2)]

        # set the minimum heights of rows and columns
        bWLay.setColumnMinimumWidth(0, round(oneZeroLab.width()*1.2))

        # Add the layout to the bottomWidget
        self.bottomWidget.setLayout(bWLay)
        self.bottomWidget.layout().activate()

    def on_botWidgPressed(self, pressed):
        direction = 0
        vert = False
        if pressed in (self.h01, self.v01):
            direction = 1
        if pressed in (self.v10, self.v01):
            vert = True
        angle = self.window().real.angle_for_horizontal_bulk(direction)
        if vert:
            angle += 90
        if qtw.qApp.keyboardModifiers() == qtc.Qt.ControlModifier:
            angle += 180
        self.text.updateText(angle)
