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

# NOTES: best way to go to also rename Arduino Micro to ViPErLEED is
# to (1) create a copy of the boards.txt file that comes with the CLI
# AFTER installing the avr core (TODO: check where it is, see also the
# notex.txt file) -- call it boards.txt.ViPErLEED; (2) in boards.txt.ViPErLEED
# replace 'Arduino Micro' with 'ViPErLEED' (unless the file was already there)
# (3) rename boards.txt to boards.txt.bak, and boards.txt.ViPErLEED to
# boards.txt; (4) compile and upload; (5) undo no.3.

from collections import namedtuple
import json
from pathlib import Path
import shutil
import subprocess
import sys
from zipfile import ZipFile

from PyQt5 import (QtWidgets as qtw,
                   QtCore as qtc,
                   QtNetwork as qtn)

from viperleed.guilib.measure import hardwarebase as base
from viperleed.guilib.measure.widgets.pathselector import PathSelector


_INVOKE = qtc.QMetaObject.invokeMethod
NOT_SET = '\u2014'


FirmwareVersionInfo = namedtuple('FirmwareVersionInfo',
                                 ('folder_name', 'version', 'path'))


class ViPErLEEDFirmwareError(base.ViPErLEEDErrorEnum):
    """Errors related to the FirmwareUpgradeDialog."""

    ERROR_CONTROLLER_NOT_FOUND = (
        500,
        'Controller at port {} is no longer connected. Try re-plugging it.'
        )
    ERROR_ARDUINO_CLI_NOT_FOUND = (
        501,
        'The Arduino command-line interface was not found in {}. Make sure the '
        'Arduino CLI is installed before trying to detect any devices or '
        'uploading firmware.'
        )
    ERROR_INSTALL_FAILED = (
        502,
        'Failed to download and install the latest Arduino CLI. Make sure the '
        'PC is connected to the internet.'
        )
    ERROR_NO_SUITABLE_CLI = (
        503,
        'Could not find a suitable precompiled version of the Arduino CLI. '
        'You may have to compile from the source at '
        'https://github.com/arduino/arduino-cli'
        )
    ERROR_ARDUINO_CLI_FAILED = (
        504,
        'The Arduino CLI failed. Try to reinstall it. '
        'Return code={}. The error was: {}'
        )


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
                'firmware_version': qtw.QComboBox(),
                },
            'buttons' : {
                'refresh': qtw.QPushButton('&Refresh'),
                'upload': qtw.QPushButton('&Upload firmware'),
                'done': qtw.QPushButton('&Done'),
                'upgrade': qtw.QPushButton('Upgrade &Arduino CLI')
                },
            'labels': {
                'ctrl_type': qtw.QLabel(f'Controller type: {NOT_SET}'),
                'firmware_version': qtw.QLabel('Installed firmware '
                                               f'version: {NOT_SET}'),
                'highest_version': qtw.QLabel('Most recent firmware '
                                              f'version: {NOT_SET}'),
                },
            'timers': {
                #TODO: maybe add timers (trigger detect devices instead of button)
                },
            }

        self._uploader = FirmwareUploader()
        self._upload_thread = qtc.QThread()
        self._uploader.moveToThread(self._upload_thread)
        self._upload_thread.start()

        self._progress_bar = qtw.QProgressBar()
        self._progress_bar.setValue(0)

        self.setWindowTitle('Upgrade ViPErLEED box firmware')
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

    def _make_cli_install_disclaimer(self):
        """Create the dialog asking the user to install Arduino CLI."""
        disclaimer = qtw.QMessageBox(parent=self)
        disclaimer.setText(
            'In order to use the firmware upgrade tool you need to install '
            'the Arduino command-line interface. This software can be found '
            'under https://github.com/arduino/arduino-cli and is published '
            'under the GPL-3.0 license. By clicking accept, you consent to '
            'downloading and installing the Arduino CLI in the folder you '
            'specified in the settings menu. Make sure you are connected '
            'to the internet before clicking accept.'
            )
        disclaimer.addButton(qtw.QPushButton('Cancel'),
                             disclaimer.RejectRole)
        accept = qtw.QPushButton('Accept')
        disclaimer.addButton(accept, disclaimer.AcceptRole)
        accept.clicked.connect(self._download_and_install_arduino_cli)
        disclaimer.open()

    def _compose_controller_selection(self):
        """Compose controller QComboBox and refresh button."""
        layout = qtw.QHBoxLayout()
        layout.addWidget(qtw.QLabel('Select controller:'))
        layout.addWidget(self.controls['controllers'], stretch=1)
        layout.addWidget(self.buttons['refresh'])
        return layout

    def _compose_upgrade_and_done_button(self):
        """Compose the upgrade and the done buttons."""
        layout = qtw.QHBoxLayout()
        layout.addWidget(self.buttons['upgrade'])
        layout.addStretch(1)
        layout.addWidget(self.buttons['done'])
        return layout

    def _compose_firmware_selection(self):
        """Compose firmware QComboBox and upload button."""
        layout = qtw.QHBoxLayout()
        layout.addWidget(qtw.QLabel('Select firmware version:'))
        layout.addWidget(self.controls['firmware_version'], stretch=1)
        layout.addWidget(self.buttons['upload'])
        self.buttons['upload'].setEnabled(False)
        return layout

    def _compose_info_section(self):
        """Compose controller info."""
        layout = qtw.QHBoxLayout()
        for label in self.labels.values():
            layout.addWidget(label)
        return layout

    def _compose_path_selection(self):
        """Compose PathSelector."""
        layout = qtw.QHBoxLayout()
        layout.addWidget(qtw.QLabel('Select firmware folder:'))
        layout.addWidget(self.controls['firmware_path'])
        return layout

    def _compose_progress_bar(self):
        """Compose QProgressBar."""
        layout = qtw.QHBoxLayout()
        layout.addWidget(qtw.QLabel('Progess:'))
        layout.addWidget(self._progress_bar)
        return layout

    def _connect(self):
        """Connect children signals."""
        self.buttons['done'].clicked.connect(self.accept)
        self.buttons['refresh'].clicked.connect(self._detect_controllers)
        self.buttons['upload'].clicked.connect(self._upload)
        self.buttons['upgrade'].clicked.connect(self._upgrade_arduino_cli)
        self.controls['firmware_path'].path_changed.connect(
            self._find_firmware_versions
            )
        self.controls['controllers'].currentTextChanged.connect(
            self._update_ctrl_labels
            )
        self._uploader.error_occurred.connect(self.error_occurred)
        self._uploader.controllers_detected.connect(
            self._on_controllers_detected
            )
        self._uploader.upload_finished.connect(self._on_upload_finished)
        self._uploader.cli_installation_finished.connect(
            self._on_cli_install_done
            )
        self._uploader.progress_occurred.connect(self._progress_bar.setValue)

    @qtc.pyqtSlot(bool)
    def _continue_open(self, is_installed):
        """Open FirmwareUpgradeDialog if the Arduino CLI is_installed.

        If the Arduino CLI is installed this will open up the
        FirmwareUpgradeDialog. If the CLI is not installed, the dialog
        asking the user to install the CLI will be opened.

        Parameters
        ----------
        is_installed : bool
            True if the Arduino CLI is installed on the PC.

        Returns
        -------
        None.
        """
        base.safe_disconnect(self._uploader.cli_found, self._continue_open)
        # The error_occurred signal is reconnected because it is
        # disconnected in the open() method of the dialog.
        self._uploader.error_occurred.connect(self.error_occurred)
        if is_installed:
            super().open()
            return
        self._make_cli_install_disclaimer()

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
        widgets = (*self.buttons.values(), *self.controls.values())
        for widget in widgets:
            widget.setEnabled(enable)

    @qtc.pyqtSlot()
    def _detect_controllers(self):
        """Detect connected ViPErLEED controllers."""
        self._enable_buttons(False)
        detect_viperino = True
        _INVOKE(self._uploader, 'get_viperleed_hardware',
                qtc.Q_ARG(bool, detect_viperino))

    @qtc.pyqtSlot()
    def _download_and_install_arduino_cli(self):
        """Install the Arduino CLI.

        This will trigger the download and installation of the Arduino
        CLI. The FirmwareUpgradeDialog will be opened while the CLI
        is being installed.

        Returns
        -------
        None.
        """
        self._on_cli_install_done(False)
        _INVOKE(self._uploader, 'get_arduino_cli_from_git')
        super().open()

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
        version = self.controls['firmware_version'].currentData()
        can_upload = bool(ctrl and version)
        self.buttons['upload'].setEnabled(can_upload)

    @qtc.pyqtSlot()
    def _find_firmware_versions(self, *args):
        """Search for available firmware versions."""
        firmware_dict = {}
        f_path = self.controls['firmware_path'].path

        if not f_path or f_path == Path():  # User did not select path.
            self._update_combo_box('firmware_version', firmware_dict)
            return

        for file in f_path.glob('*.zip'):
            with ZipFile(file, mode='r') as archive:
                # !!! We are assuming here that the first folder in the
                # .zip file matches the name of the firmware version.
                folder_name, *_ = archive.namelist()[0].split('/')
                *_, version_str = folder_name.split('_')
                try:
                    version = base.Version(version_str)
                except ValueError:
                    # Some non-firmware zip file
                    continue
                firmware_dict[folder_name] = FirmwareVersionInfo(folder_name,
                                                                 version, file)
        self._update_combo_box('firmware_version', firmware_dict)
        self._find_most_recent_firmware_version()

    def _find_most_recent_firmware_version(self):
        """Detect most recent firmware suitable for controller."""
        controller = self.controls['controllers'].currentData()
        nr_versions = self.controls['firmware_version'].count()

        if not controller or not controller['name_raw']:
            return

        versions = []
        for i in range(nr_versions):
            firmware = self.controls['firmware_version'].itemData(i)
            if controller['name_raw'] in firmware.folder_name:
                versions.append(firmware.version)
        try:
            max_version = max(versions)
        except ValueError:  # No firmware available
            max_version = NOT_SET

        self.labels['highest_version'].setText(
            f'Most recent firmware version: {max_version}'
            )

    @qtc.pyqtSlot(bool)
    def _on_cli_install_done(self, enable):
        """Enable/disable widgets while keeping done button enabled.

        Parameters
        ----------
        enable : bool
            Whether widgets should be enabled.

        Returns
        -------
        None.
        """
        self._ctrl_enable(enable)
        if enable:
            self._enable_upload_button()
        self.buttons['done'].setEnabled(True)

    @qtc.pyqtSlot(dict)
    def _on_controllers_detected(self, data_dict):
        """Display detected controllers and enable buttons.

        Parameters
        ----------
        data_dict : dict
            A dict of available controllers.

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

    @qtc.pyqtSlot()
    def _update_ctrl_labels(self, *args):
        """Display controller firmware version."""
        ctrl = self.controls['controllers'].currentData() or {}
        version = ctrl.get('version', NOT_SET)
        box_id = ctrl.get('box_id', NOT_SET)
        self.labels['firmware_version'].setText(
            f'Installed firmware version: {version}'
            )
        self.labels['ctrl_type'].setText(
            f'Controller type: {box_id}'
            )
        self.labels['highest_version'].setText(
            f'Most recent firmware version: {NOT_SET}'
            )
        self._find_most_recent_firmware_version()

    @qtc.pyqtSlot()
    def _upgrade_arduino_cli(self):
        """Upgrade the Arduino CLI."""
        self._ctrl_enable(False)
        _INVOKE(self._uploader, 'get_arduino_cli_from_git')

    def _update_combo_box(self, which_combo, data_dict):
        """Replace displayed firmware/controllers with the detected ones.
        
        Parameters
        ----------
        which_combo : str
            A str determining which QComboBox should be updated.
        data_dict : dict
            A dict of available controllers or firmware.

        Returns
        -------
        None.
        """
        self.controls[which_combo].clear()
        for item_text, item_data in data_dict.items():
            self.controls[which_combo].addItem(item_text, userData=item_data)
        self._enable_upload_button()

    @qtc.pyqtSlot()
    def _upload(self):
        """Upload selected firmware to selected controller."""
        if not self.controls['controllers'].currentData():
            return
        if not self.controls['firmware_version'].currentData():
            return

        selected_ctrl = self.controls['controllers'].currentData()
        firmware = self.controls['firmware_version'].currentData()
        tmp_path = self.controls['firmware_path'].path / 'tmp_'
        upload = True
        _INVOKE(self._uploader, 'compile', qtc.Q_ARG(dict, selected_ctrl),
                qtc.Q_ARG(FirmwareVersionInfo, firmware),
                qtc.Q_ARG(Path, tmp_path), qtc.Q_ARG(bool, upload))
        self._ctrl_enable(False)

    @qtc.pyqtSlot()
    def accept(self):
        """Clean up, then accept."""
        self._clean_up()
        super().accept()

    @qtc.pyqtSlot()
    def open(self):
        """Check if the Arduino CLI is installed before opening dialog.

        This method interrupts the regular opening of the dialog and
        will stall until the FirmwareUploader has finished the check if
        the Arduino CLI is installled.
        """
        self._uploader.cli_found.connect(self._continue_open)
        base.safe_disconnect(self._uploader.error_occurred,
                             self.error_occurred)
        _INVOKE(self._uploader, 'is_cli_installed')


class FirmwareUploader(qtc.QObject):
    """Worker that handles firmware uploads in another thread."""
    error_occurred = qtc.pyqtSignal(tuple)

    # Emitted when the upload of firmware to the controller is finished
    # or when it has failed.
    upload_finished = qtc.pyqtSignal()

    # Emitted when downloading and installing the Arduino CLI, cores and
    # libraries. Sent boolean is true when the full installation
    # succeeded or when only the upgrading of existing cores has failed.
    # The upgrading of existing cores is allowed to fail because the
    # old cores can still be used to operate the FirmwareUploader.
    cli_installation_finished = qtc.pyqtSignal(bool)

    # Emitted after connected controllers have been detected. Carries
    # with it the connected ViPErLEED controllers.
    controllers_detected = qtc.pyqtSignal(dict)

    # Emitted after the check if the Arduino CLI is installed.
    cli_found = qtc.pyqtSignal(bool)

    # Emitted with a value that is indicative of how much of
    # the currently running process has been completed.
    progress_occurred = qtc.pyqtSignal(int)

    def __init__(self, parent=None):
        """Initialise the firmware uploader."""
        super().__init__(parent=parent)
        self.base_path = Path().resolve() / 'hardware/arduino/arduino-cli'
        self.network = qtn.QNetworkAccessManager()
        self.network.setTransferTimeout(timeout=2000)
        # _cli_os_name is the name of the OS specific Arduino CLI
        # version that has to be downloaded from github for
        # installation. It is determined automatically before
        # downloading the CLI.
        self._cli_os_name = ''

    def _get_arduino_cores(self):
        """Return a dictionary of the cores currently installed.

        Returns
        -------
        dict
            The following {key: value} pairs are only present if the
            core detection does not fail. The returned dict will be
            empty if the detection fails.

            'id': str
                qualified name of the core
            'installed': str
                Version currently installed
            'latest': str
                Latest version available
            'name': str
                Descriptive name of core
            'maintainer': str
                Maintainer of the core code
            'website': str
                Url of the maintainer
            'email': str
                Email address of the maintainer
            'boards': dict
                Dictionary of Arduino boards that use this core
        """
        cli = self._get_arduino_cli()
        if not cli:
            return {}
        cli = subprocess.run([cli, 'core', 'list', '--format', 'json'],
                              capture_output=True)
        try:
            cli.check_returncode()
        except subprocess.CalledProcessError:
            base.emit_error(self,
                ViPErLEEDFirmwareError.ERROR_ARDUINO_CLI_FAILED,
                cli.returncode,
                cli.stderr
                )
            self.cli_installation_finished.emit(False)
            return {}
        return json.loads(cli.stdout)

    def _get_installed_cli_version(self):
        """Detect version of installed Arduino CLI.

        Returns
        -------
        Version string : str
            The version of the Arduino CLI.
        """
        ver = None

        try:
            ver = subprocess.run(
                ['arduino-cli', 'version', '--format', 'json'],
                capture_output=True
                )
        except FileNotFoundError:
            pass

        if not ver:
            arduino_cli = 'arduino-cli'
            if 'win' in sys.platform and not 'darwin' in sys.platform:
                arduino_cli += '.exe'
            arduino_cli = self.base_path.joinpath(arduino_cli)
            try:
                ver = subprocess.run(
                    [arduino_cli, 'version', '--format', 'json'],
                    capture_output=True
                    )
            except FileNotFoundError:
                pass

        if ver:
            ver = json.loads(ver.stdout)
            return ver['VersionString']
        return '0.0.0'

    @qtc.pyqtSlot(qtn.QNetworkReply)
    def _download_newest_cli(self, reply):
        """Download newest version of Arduino CLI from github.

        Download newest version of the CLI for the OS if it doesn't
        match the currently installed version of it. Downloading the CLI
        is skipped if it is already installed.

        Parameters
        ----------
        reply : QNetworkReply
            Contains the newest CLI version names and download
            urls for these CLIs for various platforms.

        Returns
        -------
        None.
        """
        base.safe_disconnect(self.network.finished, self._download_newest_cli)

        # Check if connection failed.
        if reply.error():
            base.emit_error(self,
                ViPErLEEDFirmwareError.ERROR_INSTALL_FAILED
                )
            self.cli_installation_finished.emit(False)
            return
        self.progress_occurred.emit(10)

        latest = json.loads(reply.readAll().data())

        platform = sys.platform
        if 'darwin' in platform:
            self._cli_os_name = 'macOS'
        elif 'win' in platform and 'cyg' not in platform:
            self._cli_os_name = 'Windows'
        elif 'linux' in platform:
            self._cli_os_name = 'Linux'
        else:
            base.emit_error(self,
                ViPErLEEDFirmwareError.ERROR_NO_SUITABLE_CLI
                )
            self.cli_installation_finished.emit(False)
            return

        url_latest = ''
        for asset in latest['assets']:
            # This should always pick the 32 bit version if
            # present, as it comes alphabetically earlier
            if self._cli_os_name in asset['name']:
                newest_version = asset['name'].split('_')[1]
                url_latest = asset['browser_download_url']
                self._cli_os_name = asset['name']
                break
        else:
            base.emit_error(self,
                ViPErLEEDFirmwareError.ERROR_NO_SUITABLE_CLI
                )
            self.cli_installation_finished.emit(False)
            return
        self.progress_occurred.emit(15)

        installed_ver = self._get_installed_cli_version()
        if newest_version == installed_ver:
            self.progress_occurred.emit(35)
            self._install_and_upgrade_cores()
            return

        request = qtn.QNetworkRequest(qtc.QUrl(url_latest))
        # Since Qt does not come with SSL support the RedirectAttribute
        # must be set in order to get the file via http.
        request.setAttribute(qtn.QNetworkRequest.FollowRedirectsAttribute, True)
        self.network.finished.connect(self._install_arduino_cli,
                                      type=qtc.Qt.UniqueConnection)
        self.network.get(request)
        self.progress_occurred.emit(25)

    def _install_and_upgrade_cores(self):
        """Download AVR core and upgrade existing cores and libraries."""
        self.progress_occurred.emit(60)
        self._install_arduino_core('arduino:avr')
        self.progress_occurred.emit(80)
        self._upgrade_arduino_cores_and_libraries()
        self.progress_occurred.emit(100)
        self.cli_installation_finished.emit(True)

    @qtc.pyqtSlot(qtn.QNetworkReply)
    def _install_arduino_cli(self, reply):
        """Extract downloaded Arduino CLI from .zip archive.

        Parameters
        ----------
        reply : QNetworkReply
            Contains the Arduino CLI in a .zip archive

        Returns
        -------
        None.
        """
        base.safe_disconnect(self.network.finished,
                             self._install_arduino_cli)
        self.progress_occurred.emit(40)
        with open(self.base_path / self._cli_os_name, 'wb') as archive:
            archive.write(reply.readAll())
        shutil.unpack_archive(self.base_path / self._cli_os_name, self.base_path)
        self.progress_occurred.emit(50)
        self._install_and_upgrade_cores()

    def _install_arduino_core(self, core_name):                                            # TODO: add progress bar
        """Install a given core to the Arduino CLI.

        Parameters
        ----------
        core_name : str
            Should be one of the Arduino core names (e.g., arduino:avr).
            It's easier to get this information from a call to
            get_viperleed_hardware()

        Returns
        -------
        None.
        """
        cli = self._get_arduino_cli()
        if not cli:
            return
        cli = subprocess.run([cli, 'core', 'install', core_name],
                              capture_output=True)
        try:
            cli.check_returncode()
        except subprocess.CalledProcessError:
            base.emit_error(self,
                ViPErLEEDFirmwareError.ERROR_ARDUINO_CLI_FAILED,
                cli.returncode,
                cli.stderr
                )
            self.cli_installation_finished.emit(False)

    def _get_arduino_cli(self, get_from_git=False):                                        # TODO: add progress bar
        """Pick the correct Arduino CLI tool.

        The choice is based on the current operating system.
        If the tool is not available in the arduino-cli folder
        it can be downloaded from github.

        Parameters
        ----------
        get_from_git : bool
            Download the latest version of the Arduino CLI
            from GitHub if no version is available locally

        Returns
        -------
        arduino_cli : Path
            Path to correct executable
        """
        # See if, by chance, it is available system-wide
        try:
            subprocess.run(['arduino-cli', '-h'])
        except FileNotFoundError:
            pass
        else:
            return Path('arduino-cli')

        # Look for the executable in the arduino-cli subfolder
        try:
            self.base_path.mkdir()
        except FileExistsError:
            # folder is there, it may contain the executable
            pass
        else:
            if get_from_git:
                # Get the executables from github
                self.get_arduino_cli_from_git()
            else:
                base.emit_error(self,
                    ViPErLEEDFirmwareError.ERROR_ARDUINO_CLI_NOT_FOUND,
                    self.base_path
                    )
                self.cli_installation_finished.emit(False)
                return

        # See if there is the correct executable in arduino-cli
        arduino_cli = self.base_path / 'arduino-cli'
        if 'win' in sys.platform and not 'darwin' in sys.platform:
            arduino_cli = arduino_cli.with_suffix('.exe')
        if not arduino_cli.is_file():
            base.emit_error(self,
                ViPErLEEDFirmwareError.ERROR_ARDUINO_CLI_NOT_FOUND,
                self.base_path
                )
            self.cli_installation_finished.emit(False)
            return
        return arduino_cli

    def _get_boards(self):
        """Get a list of the available Arduino boards.

        Returns
        -------
        list
        """
        cli = self._get_arduino_cli()
        if not cli:
            return []
        cli = subprocess.run([cli, 'board', 'list', '--format', 'json'],
                              capture_output=True)
        try:
            cli.check_returncode()
        except subprocess.CalledProcessError:
            base.emit_error(self,
                ViPErLEEDFirmwareError.ERROR_ARDUINO_CLI_FAILED,
                cli.returncode,
                cli.stderr
                )
            self.cli_installation_finished.emit(False)
            return []
        boards = json.loads(cli.stdout)
        return [b for b in boards if 'matching_boards' in b]

    @qtc.pyqtSlot()
    def _upgrade_arduino_cores_and_libraries(self):                                                      # TODO: add progress bar
        """Upgrade the outdated cores and libraries."""
        cli = self._get_arduino_cli()
        if not cli:
            return

        cli_update = subprocess.run([cli, 'update'], capture_output=True)
        try:
            cli_update.check_returncode()
        except subprocess.CalledProcessError:
            base.emit_error(self,
                ViPErLEEDFirmwareError.ERROR_INSTALL_FAILED
                )
            self.cli_installation_finished.emit(True)
            return

        self.progress_occurred.emit(90)

        cli_upgrade = subprocess.run([cli, 'upgrade'], capture_output=True)
        try:
            cli_upgrade.check_returncode()
        except subprocess.CalledProcessError:
            base.emit_error(self,
                ViPErLEEDFirmwareError.ERROR_ARDUINO_CLI_FAILED,
                cli_upgrade.returncode,
                cli_upgrade.stderr
                )
            self.cli_installation_finished.emit(True)

    @qtc.pyqtSlot(dict, FirmwareVersionInfo, Path, bool)
    def compile(self, selected_ctrl, firmware, tmp_path, upload):
        """Compile viper-ino for the specified board.

        Parameters
        ----------
        selected_ctrl : dict
            A dict representing the selected controller. It contains the
            controller name, COM port, and fully qualified board name.
        firmware : FirmwareVersionInfo
            A namedtuple representing the selected firmware that is
            to be uploaded to the selected controller. Contains the
            path from where the firmware is taken.
        tmp_path : Path
            The location where the selected firmware archive is
            to be extracted to. This is not to be confused with
            the path from which the firmware is taken.
        upload : bool
            If true, the extracted firmware will be uploaded to the
            selected controller.

        Returns
        -------
        None.
        """
        self.progress_occurred.emit(0)
        # Check if the selected controller is present at all.
        available_ctrls = self.get_viperleed_hardware(False)

        # Return if Arduino CLI was not found to keep buttons disabled.
        if not available_ctrls:
            return

        if not any(available_ctrls[c]['name'] == selected_ctrl['name'] and
                   available_ctrls[c]['port'] == selected_ctrl['port']
                   for c in available_ctrls):
            base.emit_error(self,
                            ViPErLEEDFirmwareError.ERROR_CONTROLLER_NOT_FOUND,
                            selected_ctrl['port'])
            self.controllers_detected.emit(available_ctrls)
            self.upload_finished.emit()
            return

        cli = self._get_arduino_cli()
        if not cli:
            return

        with ZipFile(firmware.path) as firmware_zip:
            firmware_zip.extractall(tmp_path)
        firmware_extracted = tmp_path / firmware.folder_name / 'viper-ino'
        # Add generic inner folder structure to find .ino file
        argv = ['compile', '--clean', '-b', selected_ctrl['fqbn'], firmware_extracted]
        if upload:
            argv.extend(['-u', '-p', selected_ctrl['port']])
        self.progress_occurred.emit(20)
        cli = subprocess.run([cli, *argv], capture_output=True)
        self.progress_occurred.emit(80)
        # cli = subprocess.run([cli, *argv, '--verbose'], capture_output=True) TODO: track upload progress and adjust progress bar accordingly
        try:
            cli.check_returncode()
        except subprocess.CalledProcessError:
            base.emit_error(self,
                ViPErLEEDFirmwareError.ERROR_ARDUINO_CLI_FAILED,
                cli.returncode,
                cli.stderr.decode()
                )
            self.cli_installation_finished.emit(False)
            return
        self.progress_occurred.emit(85)
        # Remove extracted archive
        shutil.rmtree(tmp_path)
        # Update controller list
        self.progress_occurred.emit(90)
        self.get_viperleed_hardware(True)
        self.progress_occurred.emit(100)
        self.upload_finished.emit()

    @qtc.pyqtSlot()
    def get_arduino_cli_from_git(self):
        """Obtain the latest version of Arduino CLI from GitHub.

        Download and extract the latest version of the Arduino
        command-line interface executables for the current platform
        asynchronously and upgrade libraries/cores. For simplicity, the
        32bit version is downloaded for both Linux and Windows. First
        the required CLI version is detected, then it is downloaded from
        github. When the installation is finished, or has failed, the
        cli_installation_finished signal is emitted.

        Returns
        -------
        None.
        """
        self.progress_occurred.emit(0)
        request = qtn.QNetworkRequest(qtc.QUrl(
            'https://api.github.com/repos/arduino/arduino-cli/releases/latest'
            ))
        self.network.finished.connect(self._download_newest_cli,
                                      type=qtc.Qt.UniqueConnection)
        self.network.get(request)

    @qtc.pyqtSlot(bool)
    def get_viperleed_hardware(self, detect_viperino):
        """Return a dict of the ViPErLEED Arduino boards.

        In fact, it returns all the ViPErLEED as well as Arduino Micro
        boards. After calling this function, one should decide whether
        to keep only those boards whose 'name' contains 'ViPErLEED'.
        
        Parameters
        ----------
        detect_viperino : bool
            ViPErLEED controller types will be detected if true.

        Returns
        -------
        viper_boards : dict
            A dict of the detected Arduino Micro boards containing dicts
            with information about the controllers. The following
            {key: value} pairs are present in controller dicts.

            'port': str
                COM port address
            'name': str
                Board name
            'fqbn': str
                Fully qualified board name
            'version': hardwarebase.Version or str
                Firmware version of the device, if the information is
                available, otherwise str
            'box_id': int or str
                int if box_id is set, str otherwise
            'name_raw': str
                Controller type name
        """
        viperleed_names = ('ViPErLEED', 'Arduino Micro')
        boards = self._get_boards()

        # Return if Arduino CLI was not found to keep buttons disabled.
        if not boards:
            return {}

        # Get all available Arduino Micro controllers.
        board_names = [b['matching_boards'][0]['name'] for b in boards]
        viper_boards = []
        for name in viperleed_names:
            for board, board_name in zip(boards, board_names):
                if name in board_name:
                    viper_boards.append(board)

        # Extract data from controller list.
        ctrl_dict = {}
        for board in viper_boards:
            port = board['port']['address']
            board = board['matching_boards'][0]
            ctrl = f'{port} {board["name"]}'
            ctrl_dict[ctrl] = {
                'port': port,
                'name': board['name'],
                'fqbn': board['fqbn'],
                'version': NOT_SET,
                'box_id': NOT_SET,
                # name_raw is used for firmware detection.
                'name_raw': '',
                }

        if not detect_viperino:
            return ctrl_dict

        # Detect ViPErLEED controllers.
        ctrl_with_id = []
        for name, (_, info) in base.get_devices('controller').items():
            port = name.split(' (')[1][:-1]
            ctrl = next((c for c in ctrl_dict
                        if ctrl_dict[c]['port'] == port), None)
            if not ctrl:
                continue
            ctrl_dict[ctrl]['version'] = info.get('firmware', NOT_SET)
            box_id = info.get('box_id')
            # box_id of the measuring ViPErinoController is 0! Check
            # must be != None because of that.
            if box_id != None:
                ctrl_dict[ctrl]['box_id'] = box_id
                ctrl_with_id.append(ctrl)

        # Adjust name for ViPErLEED controllers.
        for ctrl in ctrl_with_id:
            new_ctrl_name = ctrl_dict[ctrl]['port'] + ' '
            if ctrl_dict[ctrl]['box_id'] == 0:
                new_ctrl_name += 'ViPErLEED Data Acquisition'
                ctrl_dict[ctrl]['name_raw'] = 'ViPErLEED_Data_Acquisition'
            ctrl_dict[new_ctrl_name] = ctrl_dict.pop(ctrl)

        self.controllers_detected.emit(ctrl_dict)
        return ctrl_dict

    @qtc.pyqtSlot()
    def is_cli_installed(self):
        """Check if Arduino CLI is installed."""
        cli = self._get_arduino_cli()
        self.cli_found.emit(bool(cli))

    def moveToThread(self, thread):
        """Move self and self.network to thread."""
        super().moveToThread(thread)
        self.network.moveToThread(thread)
