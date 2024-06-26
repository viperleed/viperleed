"""Module measurementdialog of viperleed.guilib.measure.dialogs.

========================================
   ViPErLEED Graphical User Interface
========================================

Created: 2024-06-24
Author: Florian Dörr, Michele Riva

Defines the MeasurementDialog that handles selecting the measurement
type and its settings.
"""

from pathlib import Path
from time import localtime, strftime
import shutil

from PyQt5 import QtCore as qtc
from PyQt5 import QtWidgets as qtw

from viperleed.guilib.measure.classes.settings import MissingSettingsFileError
from viperleed.guilib.measure.classes.settings import ViPErLEEDSettings
from viperleed.guilib.measure.hardwarebase import DEFAULTS_PATH
from viperleed.guilib.measure.measurement import ALL_MEASUREMENTS
from viperleed.guilib.measure.widgets.pathselector import PathSelector
from viperleed.guilib.widgets.basewidgets import QNoDefaultDialogButtonBox


class MeasurementDialog(qtw.QDialog):
    """Dialog that handles selecting measurements."""

    measurement_selected = qtc.pyqtSignal(object, ViPErLEEDSettings)
    settings_not_found = qtc.pyqtSignal(Path, Exception)

    def __init__(self, parent=None, **kwargs):
        """Initialize dialog."""
        super().__init__(parent=parent, **kwargs)
        self._ctrls = {
            'type_selection': qtw.QComboBox(),
            'settings_folder': PathSelector(select_file=False),
            'settings_file': qtw.QComboBox(),
            }
        self._cfg_dir = Path()
        self._compose_and_connect()

    @property
    def cfg_dir(self):
        """Return the default config directory for user settings."""
        return self._cfg_dir

    @cfg_dir.setter
    def cfg_dir(self, new_cfg_dir):
        """Set the default config directory for user settings."""
        self._cfg_dir = Path(new_cfg_dir)
        had_no_path = True
        if self._ctrls['settings_folder'].path:
            had_path = False
        self._ctrls['settings_folder'].path = self._cfg_dir
        if had_no_path:
            self._ctrls['settings_folder'].path_changed.emit(self._cfg_dir)

    def _compose_and_connect(self):
        """Place children widgets and connect signals."""
        layout = qtw.QVBoxLayout()
        for name, cls in ALL_MEASUREMENTS.items():
            self._ctrls['type_selection'].addItem(name, userData=cls)
        _bbox = QNoDefaultDialogButtonBox
        buttons = _bbox(_bbox.Ok | _bbox.Cancel)

        layout.addWidget(self._ctrls['type_selection'])
        layout.addWidget(self._ctrls['settings_folder'])
        layout.addWidget(self._ctrls['settings_file'])
        layout.addWidget(buttons)
        self._ctrls['settings_folder'].path_changed.connect(
            self._find_appropriate_settings
            )
        self._ctrls['type_selection'].currentIndexChanged.connect(
            self._find_appropriate_settings
            )
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        self.setLayout(layout)

    @qtc.pyqtSlot()
    def _find_appropriate_settings(self, *_):
        """Find appropriate settings for the selected measurement type."""
        cls = self._ctrls['type_selection'].currentData()
        settings_folder = self._ctrls['settings_folder'].path
        self._ctrls['settings_file'].clear()
        if not settings_folder or settings_folder == Path():
            return
        matching_settings = cls.find_matching_settings_files(
                                None, settings_folder, False, False
                                )
        for settings in matching_settings:
            self._ctrls['settings_file'].addItem(settings.stem, settings)

    def accept(self):
        """Emit selected measurement type and settings path and close."""
        cls = self._ctrls['type_selection'].currentData()
        settings_path = self._ctrls['settings_file'].currentData()

        if not settings_path:
            default_path = cls.find_matching_settings_files(
                None, DEFAULTS_PATH, False, True
                )[0]
            current_time = strftime("_%Y-%m-%d_%H-%M-%S", localtime())
            name = cls.__name__
            settings_path = self.cfg_dir / (name + current_time + '.ini')
            shutil.copy2(default_path, settings_path)

        config = ViPErLEEDSettings()
        try:
            config.read(settings_path)
        except MissingSettingsFileError as err:
            self.settings_not_found.emit(settings_path, err)
            super().reject()
            return
        else:
            self.measurement_selected.emit(cls, config)

        super().accept()
