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
# latest = requests.get('https://api.github.com/repos/arduino'
                      # '/arduino-cli/releases/latest').json()
# Use print(latest['tag_name']) to get the latestest version string of the
# arduino-cli. Get arduino-cli version --format json -> convert to dict and
# and compare versions -> update if they don't match. arduino-cli update
# and upgrade only update the cores

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

from viperleed.guilib.measure.widgets.pathselector import PathSelector
from viperleed.guilib.measure import hardwarebase as base


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
        'The Arduino CLI failed. Something might be wrong with the '
        'installation. Return code={}. The error was: {}'
        )


class FirmwareUpgradeDialog(qtw.QDialog):
    """Dialog to handle user interaction when upgrading firmware."""

    error_occurred = qtc.pyqtSignal(tuple)

    def __init__(self, parent=None):
        """Initialize dialog."""
        super().__init__(parent=parent)
        self.__children = {
            'ctrls': {
                'controllers': qtw.QComboBox(),
                'firmware_path': PathSelector(select_file=False),
                'firmware_version': qtw.QComboBox(),
                },
            'buttons' : {
                'refresh': qtw.QPushButton('&Refresh'),
                'upload': qtw.QPushButton('&Upload firmware'),
                'done': qtw.QPushButton('&Done'),
                'upgrade': qtw.QPushButton('&Upgrade CLI')
                },
            'controller_info': {
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

        self.__uploader = FirmwareUploader()
        self.__upload_thread = qtc.QThread()
        self.__uploader.moveToThread(self.__upload_thread)
        self.__upload_thread.start()

        self.__progress_bar = qtw.QProgressBar()
        self.__progress_bar.setValue(0)

        self.__cli_info = qtw.QMessageBox(parent=self)

        self.__OS_name = ''

        self.setWindowTitle('Upgrade ViPErLEED box firmware')
        self.__children['ctrls']['firmware_path'].path = Path().resolve()
        self.__compose()
        self.__connect()

    def __clean_up(self):
        """Clean up before closing dialog."""
        self.__progress_bar.setValue(0)

    def __compose(self):
        """Place children widgets."""
        self.setWindowFlags(self.windowFlags()
                            & ~qtc.Qt.WindowContextHelpButtonHint)

        for btn in self.__children['buttons'].values():
            try:
                btn.setAutoDefault(False)
            except AttributeError:
                pass

        layout = qtw.QVBoxLayout()
        layout.setSpacing(layout.spacing() + 10)
        layout.addLayout(self.__compose_controller_selection())
        layout.addLayout(self.__compose_info_section())
        layout.addLayout(self.__compose_path_selection())
        layout.addLayout(self.__compose_firmware_selection())
        layout.addLayout(self.__compose_progress_bar())
        layout.addLayout(self.__compose_done_button())
        self.setLayout(layout)
        self.__compose_cli_info()

    def __compose_cli_info(self):
        """Creates the dialog asking the user to install Arduino CLI."""
        info = self.__cli_info
        info.setText(
            'In order to use the firmware upgrade tool you need to install '
            'the Arduino command-line interface. This software can be found '
            'under https://github.com/arduino/arduino-cli and is published '
            'under the GPL-3.0 license. By clicking accept, you consent to '
            'downloading and installing the Arduino CLI in the folder you '
            'specified in the settings menu. Make sure you are connected '
            'to the internet before clicking accept.'
            )
        info.addButton(qtw.QPushButton('Cancel'), info.RejectRole)
        accept = qtw.QPushButton('Accept')
        info.addButton(accept, info.AcceptRole)
        accept.clicked.connect(self.__download_and_install_arduino_cli)

    def __compose_controller_selection(self):
        layout = qtw.QHBoxLayout()
        layout.addWidget(qtw.QLabel('Select controller:'))
        layout.addWidget(self.__children['ctrls']['controllers'], stretch=1)
        layout.addWidget(self.__children['buttons']['refresh'])
        return layout

    def __compose_done_button(self):
        layout = qtw.QHBoxLayout()
        layout.addWidget(self.__children['buttons']['upgrade'])
        layout.addStretch(1)
        layout.addWidget(self.__children['buttons']['done'])
        return layout

    def __compose_firmware_selection(self):
        layout = qtw.QHBoxLayout()
        layout.addWidget(qtw.QLabel('Select firmware version:'))
        layout.addWidget(self.__children['ctrls']['firmware_version'],
                         stretch=1)
        layout.addWidget(self.__children['buttons']['upload'])
        self.__children['buttons']['upload'].setEnabled(False)
        return layout

    def __compose_info_section(self):
        layout = qtw.QHBoxLayout()
        for label in self.__children['controller_info'].values():
            layout.addWidget(label)
        return layout

    def __compose_path_selection(self):
        layout = qtw.QHBoxLayout()
        layout.addWidget(qtw.QLabel('Select firmware folder:'))
        layout.addWidget(self.__children['ctrls']['firmware_path'])
        return layout

    def __compose_progress_bar(self):
        layout = qtw.QHBoxLayout()
        layout.addWidget(qtw.QLabel('Progess:'))
        layout.addWidget(self.__progress_bar)
        return layout

    def __connect(self):
        """Connect children signals."""
        self.__children['buttons']['done'].clicked.connect(self.accept)
        self.__children['buttons']['refresh'].clicked.connect(self.__detect)
        self.__children['buttons']['upload'].clicked.connect(self.__upload)
        self.__children['buttons']['upgrade'].clicked.connect(
            self.__upgrade_arduino_cli
            )
        self.__children['ctrls']['firmware_path'].path_changed.connect(
            self.__get_firmware_versions
            )
        self.__children['ctrls']['controllers'].currentTextChanged.connect(
                    self.__update_ctrl_labels
            )
        self.__uploader.error_occurred.connect(self.error_occurred)
        self.__uploader.controllers_detected.connect(self.__update_combo_box)
        self.__uploader.controllers_detected.connect(self.__enable_buttons)
        self.__uploader.upload_finished.connect(self.__ctrl_enable)
        self.__uploader.cli_installation_finished.connect(
            self.__ctrl_enable_with_enabled_done
            )
        self.__uploader.process_progress.connect(self.__progress_bar.setValue)

    @qtc.pyqtSlot(bool)
    def __continue_open(self, is_installed):
        """Continue open() after check is finished.

        If the Arduino CLI is installed this will open up the
        FirmwareUpgradeDialog. If the CLI is not installed, the dialog
        asking the user to install the CLI will be opened.
        """
        base.safe_disconnect(self.__uploader.cli_is_installed,
                             self.__continue_open)
        self.__uploader.error_occurred.connect(self.error_occurred)
        if is_installed:
            super().open()
            return
        self.__cli_info.open()

    @qtc.pyqtSlot(bool)
    def __ctrl_enable_with_enabled_done(self, enable):
        """Enable/disable controls while keeping done enabled."""
        self.__ctrl_enable(enable)
        self.__children['buttons']['done'].setEnabled(True)

    @qtc.pyqtSlot(bool)
    def __ctrl_enable(self, enable):
        """Enable/disable all controls."""
        for btn in self.__children['buttons']:
            self.__children['buttons'][btn].setEnabled(enable)
        for slct in self.__children['ctrls']:
            self.__children['ctrls'][slct].setEnabled(enable)

    @qtc.pyqtSlot()
    def __detect(self):
        """Detect connected ViPErLEED controllers."""
        emit_controllers = True
        for btn in self.__children['buttons']:
            self.__children['buttons'][btn].setEnabled(False)
        _INVOKE(self.__uploader, 'get_viperleed_hardware',
                qtc.Q_ARG(bool, emit_controllers))

    @qtc.pyqtSlot()
    def __download_and_install_arduino_cli(self):
        """Install the Arduino CLI.

        This will trigger the download and installation of the Arduino
        CLI. The FirmwareUpgradeDialog will be opened while the CLI
        is being installed.
        """
        self.__ctrl_enable_with_enabled_done(False)
        _INVOKE(self.__uploader, 'get_arduino_cli_from_git')
        super().open()

    @qtc.pyqtSlot(str, dict)
    def __enable_buttons(self, *args):
        """Enable buttons again."""
        for btn in self.__children['buttons']:
            self.__children['buttons'][btn].setEnabled(True)

    @qtc.pyqtSlot()
    def __get_firmware_versions(self, *args):
        """Search for available firmware versions."""
        firmware_dict = {}
        f_path = self.__children['ctrls']['firmware_path'].path

        if not f_path or f_path == Path():
            self.__update_combo_box('firmware_version', firmware_dict)
            return

        for file in f_path.glob('*.zip'):
            with ZipFile(file, mode='r') as zf:
                # !!! We are assuming here that the first folder in the
                # .zip file matches the name of the firmware version.
                folder_name, *_ = zf.namelist()[0].split('/')
                *_, version_str = folder_name.split('_')
                try:
                    version = base.Version(version_str)
                except ValueError:
                    # Some non-firmware zip file
                    continue
                firmware_dict[folder_name] = FirmwareVersionInfo(folder_name,
                                                                 version, file)
        self.__update_combo_box('firmware_version', firmware_dict)
        self.__get_most_recent_firmware_version()

    def __get_most_recent_firmware_version(self, *args):
        """Detect most recent firmware suitable for controller."""
        ctrl = self.__children['ctrls']['controllers'].currentData()
        nr_versions = self.__children['ctrls']['firmware_version'].count()
        max_version = base.Version('0.0')

        if not ctrl or not ctrl['name_raw']:
            return

        i = 0
        while i < nr_versions:
            firm_ver = self.__children['ctrls']['firmware_version'].itemData(i)
            if ctrl['name_raw'] in firm_ver.folder_name:
                if firm_ver.version > max_version:
                    max_version = firm_ver.version
            i += 1

        if max_version == base.Version('0.0'):
            max_version = NOT_SET

        self.__children['controller_info']['highest_version'].setText(
            f'Most recent firmware version: {max_version}'
            )

    @qtc.pyqtSlot()
    def __update_ctrl_labels(self, *args):
        """Display controller firmware version."""
        ctrl = self.__children['ctrls']['controllers'].currentData() or {}
        version = ctrl.get('version', NOT_SET)
        box_id = ctrl.get('box_id', NOT_SET)
        self.__children['controller_info']['firmware_version'].setText(
            f'Installed firmware version: {version}'
            )
        self.__children['controller_info']['ctrl_type'].setText(
            f'Controller type: {box_id}'
            )
        self.__children['controller_info']['highest_version'].setText(
            f'Most recent firmware version: {NOT_SET}'
            )
        self.__get_most_recent_firmware_version()

    @qtc.pyqtSlot()
    def __upgrade_arduino_cli(self):
        """Upgrade the Arduino CLI."""
        self.__ctrl_enable(False)
        _INVOKE(self.__uploader, 'get_arduino_cli_from_git')

    @qtc.pyqtSlot(str, dict)
    def __update_combo_box(self, which_combo, data_dict):
        """Replace displayed firmware/controllers with the detected ones."""
        self.__children['ctrls'][which_combo].clear()
        for key, value in data_dict.items():
            self.__children['ctrls'][which_combo].addItem(key, userData=value)
        ctrl = self.__children['ctrls']['controllers'].currentData()
        version = self.__children['ctrls']['firmware_version'].currentData()
        if not ctrl or not version:
            self.__children['buttons']['upload'].setEnabled(False)
        else:
            self.__children['buttons']['upload'].setEnabled(True)

    @qtc.pyqtSlot()
    def __upload(self):
        """Upload selected firmware to selected controller."""
        if not self.__children['ctrls']['controllers'].currentData():
            return
        if not self.__children['ctrls']['firmware_version'].currentData():
            return

        selected_ctrl = self.__children['ctrls']['controllers'].currentData()
        firmware = self.__children['ctrls']['firmware_version'].currentData()
        tmp_path = self.__children['ctrls']['firmware_path'].path / 'tmp_'
        upload = True
        _INVOKE(self.__uploader, 'compile', qtc.Q_ARG(dict, selected_ctrl),
                qtc.Q_ARG(FirmwareVersionInfo, firmware),
                qtc.Q_ARG(Path, tmp_path), qtc.Q_ARG(bool, upload))
        self.__ctrl_enable(False)

    @qtc.pyqtSlot()
    def accept(self):
        """Clean up, then accept."""
        self.__clean_up()
        super().accept()

    def open(self):
        """Check if the Arduino CLI is installed before opening dialog.

        This method interrupts the regular opening of the dialog and
        will stall until the FirmwareUploader has finished the check if
        the Arduino CLI is installled.
        """
        self.__uploader.cli_is_installed.connect(self.__continue_open)
        base.safe_disconnect(self.__uploader.error_occurred,
                             self.error_occurred)
        _INVOKE(self.__uploader, 'is_cli_installed')

class FirmwareUploader(qtc.QObject):
    """Worker that handles firmware uploads in another thread."""
    error_occurred = qtc.pyqtSignal(tuple)

    # Emitted when the upload of firmware to the controller is finished
    # or when it has failed. In both cases the signal is emitted with a
    # True boolean.
    upload_finished = qtc.pyqtSignal(bool)

    # Emitted when downloading and installing the Arduino CLI, cores and
    # libraries. Sent boolean is true when the full installation
    # succeeded or when only the upgrading of existing cores has failed.
    # The upgrading of existing cores is allowed to fail because the
    # old cores can still be used to operate the FirmwareUploader.
    cli_installation_finished = qtc.pyqtSignal(bool)

    # Emitted after connected controllers have been detected. Carries
    # with it the connected ViPErLEED controllers.
    controllers_detected = qtc.pyqtSignal(str, dict)

    # Emitted after the check if the Arduino CLI is installed.
    cli_is_installed = qtc.pyqtSignal(bool)

    # Emitted with a value that is indicative of how much of
    # the currently running process has been completed.
    process_progress = qtc.pyqtSignal(int)

    def __init__(self, parent=None):
        """Initialise the firmware uploader."""
        self.base_path = Path().resolve() / 'hardware/arduino/arduino-cli'
        self.manager = qtn.QNetworkAccessManager()
        self.manager.setTransferTimeout(timeout=2000)
        super().__init__(parent=parent)

    def __get_arduino_cores(self):
        """Return a dictionary of the cores currently installed.

        Returns
        -------
        dict
            The following {key: value} pairs are present

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
        cli = self.__get_arduino_cli()
        if not cli:
            return
        cli = subprocess.run([cli, 'core', 'list', '--format', 'json'],
                              capture_output=True)
        try:
            cli.check_returncode()
        except subprocess.CalledProcessError as err:
            base.emit_error(self,
                ViPErLEEDFirmwareError.ERROR_ARDUINO_CLI_FAILED,
                cli.returncode,
                cli.stderr
                )
            self.cli_installation_finished.emit(False)
            return
        return json.loads(cli.stdout)

    def __get_installed_cli_version(self):
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
    def __get_newest_arduino_cli(self, reply):
        """Download newest version of Arduino CLI from github.

        Download newest version of the CLI for the OS if it doesn't
        match the currently installed version of it. Downloading the CLI
        is skipped if it is already installed.

        Parameters
        ----------
        reply : QNetworkReply
            Contains the newest CLI version names and download
            urls for these CLIs for various platforms.
        """

        base.safe_disconnect(self.manager.finished,
                             self.__get_newest_arduino_cli)

        # Check if connection failed.
        if reply.error() != 0:
            base.emit_error(self,
                ViPErLEEDFirmwareError.ERROR_INSTALL_FAILED
                )
            self.cli_installation_finished.emit(False)
            return
        self.process_progress.emit(10)

        latest = json.loads(reply.readAll().data())

        platform = sys.platform
        if 'darwin' in platform:
            self.__OS_name = 'macOS'
        elif 'win' in platform and 'cyg' not in platform:
            self.__OS_name = 'Windows'
        elif 'linux' in platform:
            self.__OS_name = 'Linux'
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
            if self.__OS_name in asset['name']:
                newest_version = asset['name'].split('_')[1]
                url_latest = asset['browser_download_url']
                self.__OS_name = asset['name']
                break
        else:
            base.emit_error(self,
                ViPErLEEDFirmwareError.ERROR_NO_SUITABLE_CLI
                )
            self.cli_installation_finished.emit(False)
            return
        self.process_progress.emit(15)

        installed_ver = self.__get_installed_cli_version()
        if newest_version == installed_ver:
            self.process_progress.emit(35)
            self.__install_and_upgrade_cores()
            return

        request = qtn.QNetworkRequest(qtc.QUrl(url_latest))
        # Since Qt does not come with SSL support the RedirectAttribute
        # must be set in order to get the file via http.
        request.setAttribute(qtn.QNetworkRequest.FollowRedirectsAttribute, 21)
        self.manager.finished.connect(self.__install_arduino_cli,
                                      type=qtc.Qt.UniqueConnection)
        self.manager.get(request)
        self.process_progress.emit(25)

    def __install_and_upgrade_cores(self):
        """Download AVR core and upgrade existing cores and libraries."""
        self.process_progress.emit(60)
        self.__install_arduino_core('arduino:avr')
        self.process_progress.emit(80)
        self.__upgrade_arduino_cores_and_libraries()
        self.process_progress.emit(100)
        self.cli_installation_finished.emit(True)

    @qtc.pyqtSlot(qtn.QNetworkReply)
    def __install_arduino_cli(self, reply):
        """Extract downloaded Arduino CLI from .zip archive.

        Parameters
        ----------
        reply : QNetworkReply
            Contains the Arduino CLI in a .zip archive
        """
        base.safe_disconnect(self.manager.finished,
                             self.__install_arduino_cli)
        self.process_progress.emit(40)
        with open(self.base_path / self.__OS_name, 'wb') as archive:
            archive.write(reply.readAll())
        shutil.unpack_archive(self.base_path / self.__OS_name, self.base_path)
        self.process_progress.emit(50)
        self.__install_and_upgrade_cores()

    def __install_arduino_core(self, core_name):                                            # TODO: add progress bar
        """Install a given core to the Arduino CLI.

        Parameters
        ----------
        core_name : str
            Should be one of the Arduino core names (e.g., arduino:avr).
            It's easier to get this information from a call to
            get_viperleed_hardware()
        """
        cli = self.__get_arduino_cli()
        if not cli:
            return
        cli = subprocess.run([cli, 'core', 'install', core_name],
                              capture_output=True)
        try:
            cli.check_returncode()
        except subprocess.CalledProcessError as err:
            base.emit_error(self,
                ViPErLEEDFirmwareError.ERROR_ARDUINO_CLI_FAILED,
                cli.returncode,
                cli.stderr
                )
            self.cli_installation_finished.emit(False)

    def __get_arduino_cli(self, get_from_git=False):                                        # TODO: add progress bar
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
        arduino_cli = 'arduino-cli'
        if 'win' in sys.platform and not 'darwin' in sys.platform:
            arduino_cli += '.exe'
        arduino_cli = self.base_path.joinpath(arduino_cli)
        if not arduino_cli.is_file():
            base.emit_error(self,
                ViPErLEEDFirmwareError.ERROR_ARDUINO_CLI_NOT_FOUND,
                self.base_path
                )
            self.cli_installation_finished.emit(False)
            return
        return arduino_cli

    def __get_boards(self):
        """Get a dict of the available Arduino boards.

        Returns
        -------
        dict
        """
        cli = self.__get_arduino_cli()
        if not cli:
            return False
        cli = subprocess.run([cli, 'board', 'list', '--format', 'json'],
                              capture_output=True)
        try:
            cli.check_returncode()
        except subprocess.CalledProcessError as err:
            base.emit_error(self,
                ViPErLEEDFirmwareError.ERROR_ARDUINO_CLI_FAILED,
                cli.returncode,
                cli.stderr
                )
            self.cli_installation_finished.emit(False)
            return False
        boards = json.loads(cli.stdout)
        return [b for b in boards if 'matching_boards' in b]

    @qtc.pyqtSlot()
    def __upgrade_arduino_cores_and_libraries(self):                                                      # TODO: add progress bar
        """Upgrade the outdated cores and libraries."""
        cli = self.__get_arduino_cli()
        if not cli:
            return

        cli2 = subprocess.run([cli, 'update'], capture_output=True)
        try:
            cli2.check_returncode()
        except subprocess.CalledProcessError as err:
            base.emit_error(self,
                ViPErLEEDFirmwareError.ERROR_INSTALL_FAILED
                )
            self.cli_installation_finished.emit(True)
            return

        self.process_progress.emit(90)

        cli2 = subprocess.run([cli, 'upgrade'], capture_output=True)
        try:
            cli2.check_returncode()
        except subprocess.CalledProcessError as err:
            base.emit_error(self,
                ViPErLEEDFirmwareError.ERROR_ARDUINO_CLI_FAILED,
                cli2.returncode,
                cli2.stderr
                )
            self.cli_installation_finished.emit(True)

    @qtc.pyqtSlot(dict, FirmwareVersionInfo, Path, bool)
    def compile(self, s_ctrl, firmware, tmp_path, upload):
        """Compile viper-ino for the specified board.

        Parameters
        ----------
        s_ctrl : dict
            A dict representing the selected controller.
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
        """
        self.process_progress.emit(0)
        # Check if the selected controller is present at all.
        available_ctrls = self.get_viperleed_hardware(False)

        # Return if Arduino CLI was not found to keep buttons disabled.
        if available_ctrls == False:
            return False

        if not any(available_ctrls[c]['name'] == s_ctrl['name'] and
                   available_ctrls[c]['port'] == s_ctrl['port']
                   for c in available_ctrls):
            base.emit_error(self,
                            ViPErLEEDFirmwareError.ERROR_CONTROLLER_NOT_FOUND,
                            s_ctrl['port'])
            self.controllers_detected.emit('controllers', available_ctrls)
            self.upload_finished.emit(True)
            return

        cli = self.__get_arduino_cli()
        if not cli:
            return

        with ZipFile(firmware.path) as firmware_zip:
            firmware_zip.extractall(tmp_path)
        firmware_extracted = tmp_path / firmware.folder_name / 'viper-ino'
        # Add generic inner folder structure to find .ino file
        argv = ['compile', '--clean', '-b', s_ctrl['fqbn'], firmware_extracted]
        if upload:
            argv.extend(['-u', '-p', s_ctrl['port']])
        self.process_progress.emit(20)
        cli = subprocess.run([cli, *argv], capture_output=True)
        self.process_progress.emit(80)
        # cli = subprocess.run([cli, *argv, '--verbose'], capture_output=True)
        try:
            cli.check_returncode()
        except subprocess.CalledProcessError as err:
            base.emit_error(self,
                ViPErLEEDFirmwareError.ERROR_ARDUINO_CLI_FAILED,
                cli.returncode,
                cli.stderr.decode()
                )
            self.cli_installation_finished.emit(False)
            return
        self.process_progress.emit(85)
        # Remove extracted archive
        shutil.rmtree(tmp_path)
        # Update controller list
        self.process_progress.emit(90)
        self.get_viperleed_hardware(True)
        self.process_progress.emit(100)
        self.upload_finished.emit(True)

    @qtc.pyqtSlot()
    def get_arduino_cli_from_git(self):
        """Obtain the latest version of Arduino CLI from GitHub.

        Download and extract the latest version of the Arduino
        command-line interface executables for the current platform.
        For simplicity, the 32bit version is downloaded for both Linux
        and Windows.

        Returns
        -------
        None.
        """
        self.process_progress.emit(0)
        request = qtn.QNetworkRequest(qtc.QUrl('https://api.github.com/repos/'
                                      'arduino/arduino-cli/releases/latest'))
        self.manager.finished.connect(self.__get_newest_arduino_cli,
                                      type=qtc.Qt.UniqueConnection)
        self.manager.get(request)

    @qtc.pyqtSlot(bool)
    def get_viperleed_hardware(self, emit_controllers):
        """Return a list of the ViPErLEED Arduino boards.

        In fact, it returns all the ViPErLEED as well as Arduino Micro
        boards. After calling this function, one should decide whether
        to keep only those boards whose 'name' contains 'ViPErLEED'.

        Returns
        -------
        viper_boards : dict
            The detected Arduino Micro boards.
        """
        viperleed_names = ('ViPErLEED', 'Arduino Micro')
        boards = self.__get_boards()

        # Return if Arduino CLI was not found to keep buttons disabled.
        if boards == False:
            return False

        # Get all available Arduino Micro controllers.
        board_names = [b['matching_boards'][0]['name'] for b in boards]
        viper_boards = []
        for name in viperleed_names:
            for i, board_name in enumerate(board_names):
                if name in board_name:
                    viper_boards.append(boards[i])

        # Extract data from controller list.
        ctrl_dict = {}
        for b in viper_boards:
            ctrl = b['port']['address'] + ' ' + b['matching_boards'][0]['name']
            ctrl_dict[ctrl] = {}
            ctrl_dict[ctrl]['port'] = b['port']['address']
            ctrl_dict[ctrl]['name'] = b['matching_boards'][0]['name']
            ctrl_dict[ctrl]['fqbn'] = b['matching_boards'][0]['fqbn']
            ctrl_dict[ctrl]['version'] = NOT_SET
            ctrl_dict[ctrl]['box_id'] = NOT_SET
            # name_raw is used for firmware detection.
            ctrl_dict[ctrl]['name_raw'] = None

        if not emit_controllers:
            return ctrl_dict

        # Detect ViPErLEED controllers.
        ctrl_with_id = []
        for name, (cls, info) in base.get_devices('controller').items():
            for ctrl in ctrl_dict:
                if ctrl_dict[ctrl]['port'] == name.split(' (')[1][:-1]:
                    if info.get('firmware'):
                        ctrl_dict[ctrl]['version'] = info['firmware']
                    if info.get('box_id') != None:
                        ctrl_dict[ctrl]['box_id'] = info['box_id']
                        ctrl_with_id.append(ctrl)

        # Adjust name for ViPErLEED controllers.
        for ctrl in ctrl_with_id:
            new_ctrl_name = ctrl_dict[ctrl]['port'] + ' '
            if ctrl_dict[ctrl]['box_id'] == 0:
                new_ctrl_name += 'ViPErLEED Data Acquisition'
                ctrl_dict[ctrl]['name_raw'] = 'ViPErLEED_Data_Acquisition'
            ctrl_dict[new_ctrl_name] = ctrl_dict.pop(ctrl)

        self.controllers_detected.emit('controllers', ctrl_dict)

    @qtc.pyqtSlot()
    def is_cli_installed(self):
        """Check if Arduino CLI is installed.
        
        Returns
        -------
        None
        """
        cli = self.__get_arduino_cli()
        if cli:
            self.cli_is_installed.emit(True)
        else:
            self.cli_is_installed.emit(False)

    def moveToThread(self, thread):
        """Overloaded moveToThread.
        
        Moves the QNetworkAccessManager to the same thread
        as the FirmwareUploader.
        
        Returns
        -------
        None
        """
        super().moveToThread(thread)
        self.manager.moveToThread(thread)
