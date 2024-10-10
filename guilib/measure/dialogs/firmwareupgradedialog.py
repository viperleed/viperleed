"""Module firmwareupgradedialog of viperleed.guilib.measure.dialogs.

========================================
   ViPErLEED Graphical User Interface
========================================

Created: 2024-01-30
Author: Florian Dörr, Michele Riva

Defines the FirmwareUpgradeDialog class that handles user interaction
while upgrading firmware. The Arduino command-line interface tool that
allows compiling and uploading the viper-ino.ino sketch to the hardware
is integrated into this dialog.
The Arduino CLI is available from
https://github.com/arduino/arduino-cli/releases
"""

from pathlib import Path
import re
from zipfile import ZipFile

from PyQt5 import QtCore as qtc
from PyQt5 import QtWidgets as qtw

from viperleed.guilib.measure import hardwarebase as base
from viperleed.guilib.measure.classes import settings
from viperleed.guilib.measure.classes.ioarduinocli import ArduinoCLIInstaller
from viperleed.guilib.measure.classes.ioarduinocli import (
    FirmwareArchiveUploader
    )
from viperleed.guilib.measure.classes.ioarduinocli import FirmwareUploader
from viperleed.guilib.measure.classes.ioarduinocli import FirmwareVersionInfo
from viperleed.guilib.measure.classes.ioarduinocli import NOT_SET
from viperleed.guilib.measure.widgets.pathselector import PathSelector


_INVOKE = qtc.QMetaObject.invokeMethod
# For development use only. If GET_ARCHIVED_FIRMWARE is true, firmware
# will be taken from .zip archives. If false, regular folders will be
# searched. Useful for uploading firmware in development.
GET_ARCHIVED_FIRMWARE = True


def without_cpp_comments(file):
    """Yield non-empty lines from file that are not comments."""
    in_comment_block = False
    for line in file:
        line, *_ = line.split('//')
        line = line.strip()
        if not line:
            continue
        if in_comment_block and line.endswith('*/'):
            in_comment_block = False
            continue   # The next line may be a good one
        if line.startswith('/*'):
            in_comment_block = True
        if in_comment_block:
            continue
        yield line


def get_firmware_version_from_ino_file(file_name):
    """Return a Version from the contents of file_name."""
    major, minor = None, None
    with file_name.open(encoding='utf-8') as file:
        for line in without_cpp_comments(file):
            if 'FIRMWARE_VERSION_MAJOR' in line:
                try:
                    major = int(re.sub(r'\D', '', line))
                except ValueError:
                    continue
            if 'FIRMWARE_VERSION_MINOR' in line:
                try:
                    minor = int(re.sub(r'\D', '', line))
                except ValueError:
                    continue
            if major is not None and minor is not None:
                return base.Version(major, minor)
    raise ValueError(f'No firmware version in {file_name}')


def get_version_from_folder_name(folder_name):
    """Return a Version from a folder name."""
    *_, version_str = folder_name.split('_')
    try:
        return base.Version(version_str)
    except ValueError:  # Not a valid folder name
        pass
    return None


class FirmwareUpgradeDialog(qtw.QDialog):
    """Dialog to handle user interaction when upgrading firmware."""

    error_occurred = qtc.pyqtSignal(tuple)

    def __init__(self, parent=None):
        """Initialize dialog."""
        super().__init__(parent=parent)
        self._children = {
            'ctrls': {
                'controllers': qtw.QComboBox(),
                'firmware_path': PathSelector(select_file=False),
                'firmware_versions': qtw.QComboBox(),
                },
            'buttons' : {
                'refresh': qtw.QPushButton('&Refresh'),
                'upload': qtw.QPushButton('&Upload firmware'),
                'done': qtw.QPushButton('&Done'),
                'upgrade_cli': qtw.QPushButton('Upgrade &Arduino CLI'),
                },
            'labels': {
                'firmware_version': qtw.QLabel('Installed firmware '
                                               f'version: {NOT_SET}'),
                'highest_version': qtw.QLabel('Most recent firmware '
                                              f'version: {NOT_SET}'),
                },
            }

        self._uploader = (FirmwareArchiveUploader() if GET_ARCHIVED_FIRMWARE
                          else FirmwareUploader())
        self._upload_thread = qtc.QThread()
        self._uploader.moveToThread(self._upload_thread)
        self._upload_thread.start()

        self._downloader = ArduinoCLIInstaller()
        self._download_thread = qtc.QThread()
        self._downloader.moveToThread(self._download_thread)
        self._download_thread.start()

        self._progress_bar = qtw.QProgressBar()
        self._progress_bar.setValue(0)

        self.setWindowTitle('Upgrade ViPErLEED box firmware')
        # This path is set to the firmware path given in system
        # settings if present when opening the dialog.
        self.controls['firmware_path'].path = Path().resolve()
        self._compose()
        self._connect()

    @property
    def buttons(self):
        """Return the buttons of this dialog."""
        return self._children['buttons']

    @property
    def controls(self):
        """Return the controls of this dialog."""
        return self._children['ctrls']

    @property
    def labels(self):
        """Return the labels of this dialog."""
        return self._children['labels']

    def _clean_up(self):
        """Clean up before closing dialog."""
        self._progress_bar.setValue(0)

    def _compose(self):
        """Place children widgets."""
        self.setWindowFlags(self.windowFlags()
                            & ~qtc.Qt.WindowContextHelpButtonHint)

        for btn in self.buttons.values():
            try:
                btn.setAutoDefault(False)
            except AttributeError:
                pass

        layout = qtw.QVBoxLayout()
        layout.setSpacing(layout.spacing() + 10)
        layout.addLayout(self._compose_controller_selection())
        layout.addLayout(self._compose_info_section())
        layout.addLayout(self._compose_path_selection())
        layout.addLayout(self._compose_firmware_selection())
        layout.addLayout(self._compose_progress_bar())
        layout.addLayout(self._compose_upgrade_and_done_button())
        self.setLayout(layout)

    def _compose_controller_selection(self):
        """Return a layout of the controller dropdown and refresh button."""
        layout = qtw.QHBoxLayout()
        layout.addWidget(qtw.QLabel('Select controller:'))
        layout.addWidget(self.controls['controllers'], stretch=1)
        layout.addWidget(self.buttons['refresh'])
        return layout

    def _compose_firmware_selection(self):
        """Return a layout of the firmware dropdown and upload button."""
        layout = qtw.QHBoxLayout()
        layout.addWidget(qtw.QLabel('Select firmware version:'))
        layout.addWidget(self.controls['firmware_versions'], stretch=1)
        layout.addWidget(self.buttons['upload'])
        self.buttons['upload'].setEnabled(False)
        return layout

    def _compose_info_section(self):
        """Return a layout of the controller info."""
        layout = qtw.QHBoxLayout()
        for label in self.labels.values():
            layout.addWidget(label)
        return layout

    def _compose_path_selection(self):
        """Return a layout of the PathSelector."""
        layout = qtw.QHBoxLayout()
        layout.addWidget(qtw.QLabel('Select firmware folder:'))
        layout.addWidget(self.controls['firmware_path'])
        return layout

    def _compose_progress_bar(self):
        """Return a layout of the QProgressBar."""
        layout = qtw.QHBoxLayout()
        layout.addWidget(qtw.QLabel('Progess:'))
        layout.addWidget(self._progress_bar)
        return layout

    def _compose_upgrade_and_done_button(self):
        """Return a layout of the upgrade and the done buttons."""
        layout = qtw.QHBoxLayout()
        layout.addWidget(self.buttons['upgrade_cli'])
        layout.addStretch(1)
        layout.addWidget(self.buttons['done'])
        return layout

    def _connect(self):
        """Connect children signals."""
        self.buttons['done'].clicked.connect(self.accept)
        self.buttons['refresh'].clicked.connect(self._detect_controllers)
        self.buttons['upload'].clicked.connect(self._upload)
        self.buttons['upgrade_cli'].clicked.connect(
            self._upgrade_arduino_cli_and_cores
            )
        self.controls['firmware_path'].path_changed.connect(
            self._find_firmware_versions
            )
        self.controls['controllers'].currentTextChanged.connect(
            self._update_ctrl_labels
            )
        self._downloader.error_occurred.connect(self.error_occurred)
        self._uploader.error_occurred.connect(self.error_occurred)
        self._uploader.controllers_detected.connect(
            self._on_controllers_detected
            )
        self._uploader.upload_finished.connect(self._on_upload_finished)
        self._downloader.cli_installation_finished.connect(self._on_cli_done)
        self._uploader.cli_failed.connect(self._on_cli_done)
        self._downloader.progress_occurred.connect(self._progress_bar.setValue)
        self._uploader.progress_occurred.connect(self._progress_bar.setValue)

    @qtc.pyqtSlot(bool, bool)
    def _continue_open(self, is_installed, is_outdated):
        """Open FirmwareUpgradeDialog if the Arduino CLI is_installed.

        If the Arduino CLI is installed this will open up the
        FirmwareUpgradeDialog. If the CLI is not installed, the dialog
        asking the user to install the CLI will be opened.

        Parameters
        ----------
        is_installed : bool
            True if the Arduino CLI is installed on the PC.
        is_outdated : bool
            True if the CLI needs to be updated to be compatible.

        Returns
        -------
        None.
        """
        base.safe_disconnect(self._downloader.cli_found, self._continue_open)
        # The error_occurred signal is reconnected because it is
        # disconnected in the open() method of the dialog.
        self._downloader.error_occurred.connect(self.error_occurred)
        if is_installed and not is_outdated:
            super().open()
            return
        if is_installed:
            disclaimer_text = (
                '<p>We detected an outdated Arduino CLI version on '
                'your system. The CLI needs to be updated as this '
                'version is incompatible with the tool. </p>'
                )
        else:
            disclaimer_text = (
                '<p>We could not find the Arduino CLI on your system. '
                'If it is installed, you can set its location in the '
                'System Settings menu.</p>'
                )
        self._make_cli_install_disclaimer(disclaimer_text)

    def _ctrl_enable(self, enable):
        """Enable/disable all widgets.

        Parameters
        ----------
        enable : bool
            Whether widgets should be enabled.

        Returns
        -------
        None.
        """
        self._enable_buttons(enable)
        for widget in self.controls.values():
            widget.setEnabled(enable)

    @qtc.pyqtSlot()
    def _detect_controllers(self):
        """Detect connected ViPErLEED controllers."""
        self._enable_buttons(False)
        detect_viperino = True
        _INVOKE(self._uploader, 'get_viperleed_hardware',
                qtc.Q_ARG(bool, detect_viperino))

    def _enable_buttons(self, enable):
        """Enable/disable buttons.

        Parameters
        ----------
        enable : bool
            Whether buttons should be enabled.

        Returns
        -------
        None.
        """
        for button in self.buttons.values():
            button.setEnabled(enable)

    def _enable_upload_button(self):
        """Enable upload if controller and firmware are selected."""
        ctrl = self.controls['controllers'].currentData()
        version = self.controls['firmware_versions'].currentData()
        can_upload = bool(ctrl and version)
        self.buttons['upload'].setEnabled(can_upload)

    @qtc.pyqtSlot()
    def _find_firmware_versions(self, *_):
        """Search for available firmware versions."""
        firmware_dict = {}
        file_path = self.controls['firmware_path'].path

        if not file_path or file_path == Path():  # User did not select path.
            self._update_combo_box('firmware_versions', firmware_dict)
            return

        _get_firmware = (self._get_archived_firmware if GET_ARCHIVED_FIRMWARE
                         else self._get_firmware)
        firmware_dict = _get_firmware(file_path)

        self._update_combo_box('firmware_versions', firmware_dict)
        self._find_most_recent_firmware_version()

    def _find_most_recent_firmware_version(self):
        """Detect most recent firmware suitable for controller."""
        controller = self.controls['controllers'].currentData()
        nr_versions = self.controls['firmware_versions'].count()
        if not controller:
            # Note that we return here before resetting the controller
            # labels as this is already done in _update_ctrl_labels.
            return
        versions = []
        for i in range(nr_versions):
            firmware = self.controls['firmware_versions'].itemData(i)
            if controller['name'] in firmware.folder_name:
                versions.append(firmware.version)
        try:
            max_version = max(versions)
        except ValueError:  # No firmware available
            max_version = NOT_SET

        self.labels['highest_version'].setText(
            f'Most recent firmware version: {max_version}'
            )

    def _get_archived_firmware(self, file_path):
        """Get firmware versions that are in a .zip archive.

        Parameters
        ----------
        file_path : Path
            The location where to look for firmware.

        Returns
        -------
        firmware_dict : dict
            A dict containing a FirmwareVersionInfo for each firmware
            folder. The keys are the text that is later on displayed
            in the QComboBox.
        """
        firmware_dict = {}
        for file in file_path.glob('*.zip'):
            with ZipFile(file, mode='r') as archive:
                folders = (n for n in archive.namelist() if '/' in n)
                top_level = (folder.split('/')[0] for folder in folders)
                folder_and_versions = (
                    (f, get_version_from_folder_name(f))
                    for f in top_level
                    )
                try:
                    folder_name, version = next(folder_and_versions)
                except StopIteration:  # Some non-firmware zip file
                    continue
                firmware_dict[folder_name] = FirmwareVersionInfo(folder_name,
                                                                 version, file)
        return firmware_dict

    def _get_firmware(self, file_path):
        """Get firmware versions that are in a regular folder.

        The directory file_path is searched recursively for supported
        .ino files. This search can potentially take very long if
        performed on a large directory.

        Parameters
        ----------
        file_path : Path
            The location where to look for firmware.

        Returns
        -------
        firmware_dict : dict
            A dict containing a FirmwareVersionInfo for each firmware
            folder. The keys are the text that is later on displayed
            in the QComboBox.
        """
        firmware_dict = {}
        for file_name in file_path.rglob('viper-ino.ino'):
            try:
                version = get_firmware_version_from_ino_file(file_name)
            except ValueError:  # No firmware version in file_name
                continue
            folder_name = file_name.parents[1].name
            folder = str(file_name.parents[1].relative_to(file_path))
            firmware_dict[folder] = FirmwareVersionInfo(
                folder_name, version, file_name.parents[2]
                )
        return firmware_dict

    def _make_cli_install_disclaimer(self, reason):
        """Create the dialog asking the user to install Arduino CLI.

        Parameters
        ----------
        reason : str
            The reason why the CLI needs to be updated/installed as str.

        Returns
        -------
        None.
        """
        accept_btn_text = 'Agree and install Arduino CLI'
        disclaimer = qtw.QMessageBox(parent=self)
        disclaimer.setWindowTitle('Arduino CLI not found')
        disclaimer.setTextFormat(qtc.Qt.RichText)
        disclaimer.setText(
            'The firmware-upgrade tool requires the Arduino command-line '
            'interface (CLI).<p>This is a third-party software that is <b>'
            'not part of ViPErLEED</b>.</p><p>The Arduino CLI is released '
            'under the <a href=https://www.gnu.org/licenses/gpl-3.0.html>'
            'GNU GPL-v3</a> license. You can find the source code and more '
            'information at <a href=https://github.com/arduino/arduino-cli>'
            'github.com/arduino/arduino-cli</a>.</p>'
            f'{reason}'
            '<p>You can also download and install it automatically. Please '
            'make sure to select the desired install location first in the '
            'Settings menu, otherwise the Arduino CLI will be installed '
            'in the default location.</p>'
            f'<p>By clicking on {accept_btn_text!r} you consent to:'
            f'<ul><li>Downloading and installing the Arduino CLI in '
            f'{self._downloader.base_path}</li>'
            '<li>The <a href=https://github.com/arduino/arduino-cli/blob/'
            'master/LICENSE.txt>terms and conditions</a> for usage of the '
            'Arduino CLI</li></ul></p>'
            '<p>Make sure you are connected to the internet before '
            'proceeding.</p>'
            '<p>Note that a window requesting administrator permissions '
            'may pop up if you are not running the program as an '
            'administrator. It may also ask you for permission to install '
            'software from Adafruit Industries LLC Ports and the Arduino USB '
            'driver. You have to accept all of the above for the installation '
            'to succeed.</p>'
            )
        disclaimer.addButton(qtw.QPushButton('Cancel'), disclaimer.RejectRole)
        accept = qtw.QPushButton(accept_btn_text)
        disclaimer.addButton(accept, disclaimer.AcceptRole)
        disclaimer.exec_()
        button = disclaimer.clickedButton()
        if button is accept:
            self._on_cli_done(False)
            _INVOKE(self._downloader, 'get_arduino_cli_from_git')
            super().open()

    @qtc.pyqtSlot()
    @qtc.pyqtSlot(bool)
    def _on_cli_done(self, successful=False):
        """Enable/disable widgets while keeping done button enabled.

        Parameters
        ----------
        successful : bool
            Whether installation of the Arduino CLI was successful.

        Returns
        -------
        None.
        """
        self._ctrl_enable(successful)
        if successful:
            self._enable_upload_button()
        else:
            self.buttons['upgrade_cli'].setEnabled(True)
        self.buttons['done'].setEnabled(True)

    @qtc.pyqtSlot(dict)
    def _on_controllers_detected(self, data_dict):
        """Display detected controllers and enable buttons.

        Parameters
        ----------
        data_dict : dict of dicts
            A dict of available controllers. Each key represents a
            controller and the associated value is a dict with more
            information about the controller. The keys are directly
            used to display the controllers. The expected {key: value}
            pairs in each controller dict are:
            'port': str
                COM port address
            'name': str
                Board name
            'fqbn': str
                Fully qualified board name
            'version': hardwarebase.Version or str
                Firmware version of the device, if the information is
                available, otherwise str

        Returns
        -------
        None.
        """
        self._enable_buttons(True)
        self._update_combo_box('controllers', data_dict)

    @qtc.pyqtSlot()
    def _on_upload_finished(self):
        """Enable buttons after upload or failed upload."""
        self._ctrl_enable(True)

    def _fetch_firmware_path_from_settings(self):
        """Set firmware_path to the path given in system settings."""
        # Note that we do not store the path in the system settings.
        # The user will initially always be directed to the folder
        # specified in the system settings.
        path_getter = settings.SystemSettings()
        firmware_path = path_getter.get('PATHS', 'firmware', fallback=None)
        if firmware_path:
            self.controls['firmware_path'].path = Path(firmware_path)

    @qtc.pyqtSlot()
    def _update_ctrl_labels(self, *_):
        """Display controller firmware version."""
        ctrl = self.controls['controllers'].currentData() or {}
        version = ctrl.get('version', NOT_SET)
        self.labels['firmware_version'].setText(
            f'Installed firmware version: {version}'
            )
        self.labels['highest_version'].setText(
            f'Most recent firmware version: {NOT_SET}'
            )
        self._find_most_recent_firmware_version()

    def _update_combo_box(self, which_combo, data_dict):
        """Replace displayed firmware/controllers with the detected ones.

        Parameters
        ----------
        which_combo : str
            A str determining which QComboBox should be updated.
        data_dict : dict
            A dict of available controllers or firmware. The keys must
            be strings as they are displayed as the item text in the
            QComboBox. Values are stored as the item data.

        Returns
        -------
        None.
        """
        self.controls[which_combo].clear()
        for item_text, item_data in data_dict.items():
            self.controls[which_combo].addItem(item_text, userData=item_data)
        self._enable_upload_button()

    @qtc.pyqtSlot()
    def _upgrade_arduino_cli_and_cores(self):
        """Upgrade the Arduino CLI."""
        self._on_cli_done(False)
        _INVOKE(self._downloader, 'get_arduino_cli_from_git')

    @qtc.pyqtSlot()
    def _upload(self):                                                          # TODO: uploading firmware should add settings file to defaults if missing and update the existing settings file for that specific controller
        """Upload selected firmware to selected controller."""
        selected_ctrl = self.controls['controllers'].currentData()
        firmware = self.controls['firmware_versions'].currentData()
        if not selected_ctrl or not firmware:
            return
        warning = qtw.QMessageBox(parent=self)
        warning.setWindowTitle('About to upload new firmware')
        warning.addButton(qtw.QPushButton('Abort'), warning.RejectRole)
        accept = qtw.QPushButton('Upload')
        warning.addButton(accept, warning.AcceptRole)
        warning.setText(
            'Uploading firmware to '
            f'{self.controls["controllers"].currentText()}. '
            'This will delete all firmware that is currently installed '
            'on the device. Do not disconnect the controller during the '
            f'upload. To proceed press {accept.text()!r}.'
            )
        warning.exec_()
        button = warning.clickedButton()
        if button is accept:
            _INVOKE(self._uploader, 'compile', qtc.Q_ARG(dict, selected_ctrl),
                    qtc.Q_ARG(FirmwareVersionInfo, firmware))
            self._ctrl_enable(False)

    @qtc.pyqtSlot()
    def accept(self):
        """Clean up, then accept."""
        self._clean_up()
        super().accept()

    @qtc.pyqtSlot()
    def close(self):
        """Stop threads, then quit."""
        self._upload_thread.quit()
        self._download_thread.quit()
        super().close()

    @qtc.pyqtSlot()
    def open(self):
        """Check if the Arduino CLI is installed before opening dialog.

        This method interrupts the regular opening of the dialog and
        will stall until the FirmwareUploader has finished the check if
        the Arduino CLI is installled. The firmware_path is also set
        to the path specified in the system settings before opening.

        Returns
        -------
        None.
        """
        self._fetch_firmware_path_from_settings()
        self._downloader.cli_found.connect(self._continue_open)
        base.safe_disconnect(self._downloader.error_occurred,
                             self.error_occurred)
        self._progress_bar.setValue(0)
        _INVOKE(self._downloader, 'update_cli_path_from_settings')
        _INVOKE(self._uploader, 'update_cli_path_from_settings')
        _INVOKE(self._downloader, 'is_cli_installed')
