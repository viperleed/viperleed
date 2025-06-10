"""Module new_measurement_dialog of viperleed.guilib.measure.dialogs.

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


class SelectNewMeasurementDialog(qtw.QDialog):
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
            'clone_settings': qtw.QCheckBox(),
            }
        self._cfg_dir = Path()
        self.setWindowTitle('Select measurement type')
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
            had_no_path = False
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
        self._ctrls['clone_settings'].setText('Create new settings'
                                              ' from old settings.')
        self._ctrls['clone_settings'].setLayoutDirection(qtc.Qt.RightToLeft)
        layout.addWidget(self._ctrls['clone_settings'])
        layout.addWidget(buttons)
        self._ctrls['settings_folder'].path_changed.connect(
            self._find_appropriate_settings
            )
        self._ctrls['type_selection'].currentIndexChanged.connect(
            self._find_appropriate_settings
            )
        self._ctrls['settings_file'].currentTextChanged.connect(
            self._switch_clone_settings
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
        self._ctrls['settings_file'].addItem(
            'Create new settigns from default settings', 'default'
            )
        if not settings_folder or settings_folder == Path():
            return
        matching_settings = cls.find_matching_settings_files(
                                None, settings_folder, False,
                                )
        for settings in matching_settings:
            self._ctrls['settings_file'].addItem(settings.stem, settings)

    def _clone_and_return_settings(self, cls, source_path):
        """Take from source path and clone to a path.

        Parameters
        ----------
        cls : MeasureEnergyCalibration or TimeResolved or IVVideo
            The measurement class for which settings have been selected.
        source_path : Path
            The path to the source settings.

        Returns
        -------
        settings_path : Path
            The path to the new cloned settings.
        """
        current_time = strftime("_%Y-%m-%d_%H-%M-%S", localtime())
        try:
            name = next(n for n, c in ALL_MEASUREMENTS.items()
                        if c is cls)
        except StopIteration:
            name = cls.__name__  # or should we raise some sensible exception?
        settings_path = self.cfg_dir / (name + current_time + '.ini')
        shutil.copy2(source_path, settings_path)
        return settings_path

    @qtc.pyqtSlot()
    def _switch_clone_settings(self, *_):
        """Disable/enable clone settings choice."""
        current_choice = self._ctrls['settings_file'].currentData()
        enable = not current_choice == 'default'
        self._ctrls['clone_settings'].setEnabled(enable)

    @qtc.pyqtSlot()
    def accept(self):
        """Emit selected measurement class and settings and close.

        Emtis
        -----
        measurement_selected
            If the settings was found and successfuly read.
        settings_not_found
            If the settings was not found. Carries the
            settings path and the MissingSettingsFileError.
        """
        cls = self._ctrls['type_selection'].currentData()
        settings_path = self._ctrls['settings_file'].currentData()
        default_path = cls.find_matching_settings_files(None, DEFAULTS_PATH,
                                                        False,)[0]

        if not settings_path or settings_path in ('default', default_path):
            settings_path = self._clone_and_return_settings(cls, default_path)
        elif self._ctrls['clone_settings'].isChecked():
            settings_path = self._clone_and_return_settings(cls, settings_path)

        config = ViPErLEEDSettings()
        try:
            config.read(settings_path)
        except MissingSettingsFileError as err:
            self.settings_not_found.emit(settings_path, err)
            super().reject()
            return

        self.measurement_selected.emit(cls, config)
        super().accept()

    def showEvent(self, event):          # pylint: disable=invalid-name
        """Show self."""
        self._find_appropriate_settings()
        super().showEvent(event)
