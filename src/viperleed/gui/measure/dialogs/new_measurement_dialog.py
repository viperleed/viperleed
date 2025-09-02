"""Module new_measurement_dialog of viperleed.gui.measure.dialogs.

Defines the SelectNewMeasurementDialog that handles selecting the
measurement type and its settings.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    'Florian Dörr (@FlorianDoerr)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2024-06-24'
__license__ = 'GPLv3+'

from pathlib import Path
from time import localtime
from time import strftime
import shutil

from PyQt5 import QtCore as qtc
from PyQt5 import QtWidgets as qtw

from viperleed.gui.measure.classes.settings import MissingSettingsFileError
from viperleed.gui.measure.classes.settings import ViPErLEEDSettings
from viperleed.gui.measure.hardwarebase import get_default_path
from viperleed.gui.measure.measurement import ALL_MEASUREMENTS
from viperleed.gui.measure.widgets.pathselector import PathSelector
from viperleed.gui.widgets.buttons import QNoDefaultDialogButtonBox


default = object()


class SelectNewMeasurementDialog(qtw.QDialog):
    """Dialog that handles selecting measurements.

    Signals
    -------
    measurement_selected : type, ViPErLEEDSettings
        Emitted when a measurement has been selected. Carries
        the measurement class and the requested settings.
    settings_not_found : Path, str
        Emitted when a settings file could not be read. Carries the
        path to the missing settings and the error message as a str.
    """

    measurement_selected = qtc.pyqtSignal(type, ViPErLEEDSettings)
    settings_not_found = qtc.pyqtSignal(Path, str)

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
        had_no_path = not self._ctrls['settings_folder'].path
        self._ctrls['settings_folder'].path = self._cfg_dir
        if had_no_path:
            self._ctrls['settings_folder'].path_changed.emit(self._cfg_dir)

    @property
    def selected_type(self):
        """Return the selected measurement type."""
        return self._ctrls['type_selection'].currentData()

    @property
    def selected_file(self):
        """Return the selected settings file."""
        return self._ctrls['settings_file'].currentData()

    def _compose_and_connect(self):
        """Place children widgets and connect signals."""
        layout = qtw.QFormLayout()
        for name, cls in ALL_MEASUREMENTS.items():
            self._ctrls['type_selection'].addItem(name, userData=cls)
        _bbox = QNoDefaultDialogButtonBox
        buttons = _bbox(_bbox.Ok | _bbox.Cancel)
        self._ctrls['clone_settings'].setText('Create new settings'
                                              ' from old settings')
        self._ctrls['clone_settings'].setLayoutDirection(qtc.Qt.RightToLeft)

        layout.addRow('Measurement type:', self._ctrls['type_selection'])
        layout.addRow('Settings location:', self._ctrls['settings_folder'])
        layout.addRow('Settings file:', self._ctrls['settings_file'])
        layout.addRow(self._ctrls['clone_settings'])
        layout.addRow(buttons)

        self._ctrls['settings_folder'].path_changed.connect(
            self._find_appropriate_settings
            )
        self._ctrls['type_selection'].currentIndexChanged.connect(
            self._find_appropriate_settings
            )
        self._ctrls['settings_file'].currentTextChanged.connect(
            self._update_clone_settings_enabled
            )
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        self.setLayout(layout)

    @qtc.pyqtSlot()
    def _find_appropriate_settings(self, *_):
        """Find appropriate settings for the selected measurement type."""
        settings_folder = self._ctrls['settings_folder'].path
        self._ctrls['settings_file'].clear()
        self._ctrls['settings_file'].addItem(
            'Create new settings from defaults', default
            )
        if not settings_folder:
            return
        matching_settings = self.selected_type.find_matching_settings_files(
                                directory=settings_folder, match_exactly=False,
                                )
        for settings in matching_settings:
            self._ctrls['settings_file'].addItem(settings.stem, settings)

    def _duplicate_settings_file(self, cls, source_path):
        """Create a copy of the settings file at `source_path`.

        A clone of the chosen settings is created in the default config
        directory for user settings. The name consists of the class
        __name__ the file is associated with and the currrent date and
        time at creation.

        Parameters
        ----------
        cls : type
            The measurement class for which settings have been selected.
        source_path : Path
            The path to the settings file to be duplicated.

        Returns
        -------
        settings_path : Path
            The path to the new cloned settings.
        """
        current_time = strftime("%Y-%m-%d_%H-%M-%S", localtime())
        name = cls.display_name.replace(' ', '_')
        settings_path = self.cfg_dir / f'{name}_{current_time}.ini'
        shutil.copy2(source_path, settings_path)
        return settings_path

    @qtc.pyqtSlot()
    def _update_clone_settings_enabled(self, *_):
        """Disable/enable clone settings choice."""
        enable = self.selected_file is not default
        self._ctrls['clone_settings'].setEnabled(enable)

    @qtc.pyqtSlot()
    def accept(self):
        """Emit selected measurement class and settings and close.

        Emits
        -----
        measurement_selected
            If the settings was found and successfuly read.
        settings_not_found
            If the settings was not found. Carries the
            settings path and the MissingSettingsFileError.
        """
        cls = self.selected_type
        settings_path = self.selected_file
        default_path, *_ = cls.find_matching_settings_files(
            directory=get_default_path(), match_exactly=False,
            )

        if not settings_path or settings_path in (default, default_path):
            settings_path = self._duplicate_settings_file(cls, default_path)
        elif self._ctrls['clone_settings'].isChecked():
            settings_path = self._duplicate_settings_file(cls, settings_path)

        config = ViPErLEEDSettings()
        try:
            config.read(settings_path)
        except MissingSettingsFileError as err:
            self.settings_not_found.emit(settings_path, str(err))
            super().reject()
            return

        self.measurement_selected.emit(cls, config)
        super().accept()

    def showEvent(self, event):          # pylint: disable=invalid-name
        """Find known measurement settings, then show this dialog."""
        self._find_appropriate_settings()
        super().showEvent(event)
