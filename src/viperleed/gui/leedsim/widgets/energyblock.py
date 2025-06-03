"""Module energyblock of viperleed.gui.leedsim.widgets.

Defines the EnergyBlock widget, the control for setting the LEED energy
in a simulated LEED pattern.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2020-01-12'
__license__ = 'GPLv3+'

import PyQt5.QtCore as qtc
import PyQt5.QtWidgets as qtw

from viperleed.gui.basewidgets import TextBoxWithButtons
from viperleed.gui.widgetslib import AllGUIFonts


class EnergyBlock(TextBoxWithButtons):
    def __init__(self, parent=None):
        params = {'labelText': '&Energy (eV)',
                  'textBoxText': '\u2014',
                  'topButText': '\u25b2',  # energy up
                  'botButText': '\u25bc',  # energy down
                  'parent': parent,
                  'textBoxTip': 'Change primary energy of the electron beam',
                  'topButTip': 'Increase energy by 20%. ' \
                               'Hold Ctrl down for a 5% increase',
                  'botButTip': 'Decrease energy by 20%. ' \
                               'Hold Ctrl down for a 5% decrease',
                  'textBoxWidth': 85
                  }

        super().__init__(**params)
        self.setTips(**params)
        self.text.setStep(1.2, 'scale')

        #set aliases for easier identification
        self.enUp = self.topBut
        self.enDown = self.botBut

        # for reasons unknown the down arrow is typeset weirdly
        # (larger than it should) -> use smaller size
        self.enDown.setFont(AllGUIFonts().smallButtonFont)

        self.makeBottomWidget()

    def makeBottomWidget(self):
        # The bottomWidget in this case contains one QlineEdit with
        # two text lines that state the minimum and maximum energy
        #
        self.limits = qtw.QLabel('Min = 10 eV\nMax = \u2014')
        self.limits.setFont(AllGUIFonts().smallTextFont)
        self.limits.adjustSize()

        bwLay = qtw.QHBoxLayout()
        bwLay.setSpacing(0)
        bwLay.setContentsMargins(0,0,0,0)
        bwLay.addWidget(self.limits, qtc.Qt.AlignLeft | qtc.Qt.AlignTop)

        self.bottomWidget.setLayout(bwLay)
        self.subWidgs.extend([self.limits])

        self.layout().setRowMinimumHeight(5, self.limits.height())
