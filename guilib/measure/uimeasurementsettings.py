"""Module uimeasurementsettings of viperleed.guilib.measure

===============================================
      ViPErLEED Graphical User Interface
===============================================

Created: 2021-10-29
Author: Michele Riva
Author: Florian Doerr

Defines the SettingsEditor class.
"""
import configparser
from pathlib import Path

import PyQt5.QtCore as qtc
import PyQt5.QtWidgets as qtw

# ViPErLEED modules
from viperleed import guilib as gl
from viperleed.guilib.basewidgets import QDoubleValidatorNoDot
from viperleed.guilib.measure import uimeasurement

TITLE = 'Measurement Settings'


class SettingsEditor(qtw.QDialog):
    """Settings editor window."""

    def __init__(self, parent=None):
        """Initialize window."""
        super().__init__(parent)
        # Keep references to controls, dialogs, and some globals
        self._ctrls = {
            'start_energy': qtw.QLineEdit(''),
            'end_energy': qtw.QLineEdit(''),
            'delta_energy': qtw.QLineEdit(''),
            'measurement_time' : qtw.QLineEdit(''),
            'limit_continuous' : qtw.QLineEdit(''),
            'cycle_time' : qtw.QLineEdit(''),
            'measure_this' : qtw.QLineEdit(''),
            'save': qtw.QPushButton("Apply changes"),
            'undo': qtw.QPushButton("Undo"),
            }

        self._dialogs = {}
        self._glob = {}

        self.__para_validator = ('start_energy', 'end_energy', 'delta_energy',
                             'measurement_time', 'limit_continuous',
                             'cycle_time')
        self.__para_text = ('measure_this',)
        self.__file_name = (uimeasurement.DEFAULT_CONFIG_PATH
                            / 'viperleed_config.ini')
        self.__settings = None
        self.__read_settings()
        self.__compose()
        # Set window properties
        self.setWindowTitle(TITLE)
        self.setAcceptDrops(False)

    def __compose(self):
        """Prepare settings editor."""
        self.setLayout(qtw.QGridLayout())

        self._ctrls['save'].setFont(gl.AllGUIFonts().buttonFont)
        self._ctrls['save'].ensurePolished()
        self._ctrls['save'].clicked.connect(self.__on_save_pressed)
        self._ctrls['save'].setEnabled(True)

        self._ctrls['undo'].setFont(gl.AllGUIFonts().buttonFont)
        self._ctrls['undo'].ensurePolished()
        self._ctrls['undo'].clicked.connect(self.__on_undo_pressed)
        self._ctrls['undo'].setEnabled(True)

        layout = self.layout()

        for key in self.__para_validator:
            self._ctrls[key].setFont(gl.AllGUIFonts().labelFont)
            self._ctrls[key].ensurePolished()
            self._ctrls[key].setValidator(QDoubleValidatorNoDot())
            self._ctrls[key].validator().setLocale(qtc.QLocale.c())
            text = self.__settings.get('measurement_settings', key)
            self._ctrls[key].setText(text)

        for key in self.__para_text:
            self._ctrls[key].setFont(gl.AllGUIFonts().labelFont)
            self._ctrls[key].ensurePolished()
            text = self.__settings.get('measurement_settings', key)
            self._ctrls[key].setText(text)

        layout.addWidget(qtw.QLabel('Start energy ='), 1, 1, 1, 1)
        layout.addWidget(self._ctrls['start_energy'], 1, 2, 1, 1)
        layout.addWidget(qtw.QLabel('End energy ='), 2, 1, 1, 1)
        layout.addWidget(self._ctrls['end_energy'], 2, 2, 1, 1)
        layout.addWidget(qtw.QLabel('Delta energy ='), 3, 1, 1, 1)
        layout.addWidget(self._ctrls['delta_energy'], 3, 2, 1, 1)
        layout.addWidget(qtw.QLabel('measurement_time ='), 4, 1, 1, 1)
        layout.addWidget(self._ctrls['measurement_time'], 4, 2, 1, 1)
        layout.addWidget(qtw.QLabel('limit_continuous ='), 5, 1, 1, 1)
        layout.addWidget(self._ctrls['limit_continuous'], 5, 2, 1, 1)
        layout.addWidget(qtw.QLabel('cycle_time ='), 6, 1, 1, 1)
        layout.addWidget(self._ctrls['cycle_time'], 6, 2, 1, 1)
        layout.addWidget(qtw.QLabel('important quantity ='), 7, 1, 1, 1)
        layout.addWidget(self._ctrls['measure_this'], 7, 2, 1, 1)
        layout.addWidget(self._ctrls['undo'], 8, 1, 1, 1)
        layout.addWidget(self._ctrls['save'], 8, 2, 1, 1)

    def __read_settings(self):
        """Read configuration file."""
        self.__settings = configparser.ConfigParser(comment_prefixes='/',
                                                    allow_no_value=True,
                                                    strict=False)
        self.__settings.read(self.__file_name)

    def __on_undo_pressed(self):
        """Undo changes in QLineEdits."""
        for key in (*self.__para_validator, *self.__para_text):
            text = self.__settings.get('measurement_settings', key)
            self._ctrls[key].setText(text)

    def __on_save_pressed(self):
        """Save changes to settings."""
        for key in (*self.__para_validator, *self.__para_text):
            text = self._ctrls[key].displayText()
            self.__settings.set('measurement_settings', key, text)
        with open(self.__file_name, 'w') as configfile:
            self.__settings.write(configfile)
