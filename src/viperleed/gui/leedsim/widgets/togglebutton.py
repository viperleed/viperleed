"""Module togglebutton of viperleed.gui.leedsim.widgets.

Defines the ToggleButton class, a QPushButton with a larger-than-normal
size.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2020-01-12'
__license__ = 'GPLv3+'

import PyQt5.QtWidgets as qtw


class ToggleButton(qtw.QPushButton):

    def sizeHint(self):
        self.ensurePolished()
        refbutton = qtw.QPushButton(self.text())
        refbutton.setFont(self.font())
        return refbutton.sizeHint()*1.2

    def minimumSizeHint(self):
        return self.sizeHint()
