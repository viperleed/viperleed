"""Module ioarduinocli of viperleed.guilib.measure.classes.

========================================
   ViPErLEED Graphical User Interface
========================================

Created: 2024-04-03
Author: Florian Dörr, Michele Riva

Defines the ArduinoCLI, ArduinoCLIInstaller, FirmwareArchiveUploader,
and FirmwareUploader classes that handle Arduino CLI and Arduino core
downloads and Arduino firmware upgrades.
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

from PyQt5 import QtCore as qtc
from PyQt5 import QtNetwork as qtn

from viperleed.guilib.measure import hardwarebase as base
from viperleed.guilib.measure.classes import settings


NOT_SET = '\u2014'
_MIN_CLI_VERSION = '0.19.0'
# Arduino cores that are required for operation. They can be
# installed via the FirmwareUpgradeDialog.
REQUIRED_CORES = ('arduino:avr', )


# FirmwareVersionInfo represents the selected firmware that
# is to be uploaded to the selected controller. It is used
# when detecting firmware in the FirmwareUpgradeDialog and
# when uploading firmware with the FirmwareUploader
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
    ERROR_CORE_NOT_FOUND = (
        505,
        'The installed Arduino CLI version is missing the required core(s) {}. '
        'Press "Upgrade Arduino CLI" to install the required core.'
        )

class ArduinoCLI(qtc.QObject):
    """Base class that looks for the Arduino CLI in the file system."""

    error_occurred = qtc.pyqtSignal(tuple)

    # Emitted after the check whether the Arduino CLI is installed.
    # Two bool values, the first one is true if a CLI version is
    # installed, the second one is true if the CLI version is too old.
    cli_found = qtc.pyqtSignal(bool, bool)

    def __init__(self, parent=None):
        """Initialise the Arduino CLI getter."""
        super().__init__(parent=parent)
        self.base_path = Path().resolve() / 'hardware/arduino/arduino-cli'
        self.update_cli_path_from_settings()

    def get_arduino_cli(self):
        """Pick the correct Arduino CLI tool.

        The choice is based on the current operating system.

        Returns
        -------
        arduino_cli : Path
            Path to correct executable

        Raises
        ------
        FileNotFoundError
            If the Arduino CLI was not found.
        """
        # See if, by chance, it is available system-wide
        try:
            subprocess.run(['arduino-cli', '-h'], check=True)
        except (FileNotFoundError, subprocess.CalledProcessError):
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
            raise FileNotFoundError('Arduino CLI not found')

        # See if there is the correct executable in arduino-cli
        arduino_cli = self.base_path / 'arduino-cli'
        if 'win' in sys.platform and 'darwin' not in sys.platform:
            arduino_cli = arduino_cli.with_suffix('.exe')
        if not arduino_cli.is_file():
            raise FileNotFoundError('Arduino CLI not found')
        return arduino_cli

    def get_installed_cli_version(self):
        """Detect version of installed Arduino CLI.

        Returns
        -------
        version : str
            The version of the Arduino CLI.
        """
        try:
            cli = self.get_arduino_cli()
        except FileNotFoundError:
            return '0.0.0'

        try:
            ver_json = subprocess.run([cli, 'version', '--format', 'json'],
                                      capture_output=True, check=True)
        except subprocess.CalledProcessError:
            return '0.0.0'
        ver = json.loads(ver_json.stdout)
        return ver['VersionString']

    def get_installed_cores(self):
        """Detect installed Arduino CLI cores.

        Returns
        -------
        cores : list of dict
            The installed Arduino CLI cores. Each core is represented
            by a dict. The key 'id' of each dict returns the core name.
            The following {key: value} pairs are only present for each
            dict if the core detection does not fail.
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
            'boards': list of dict
                List of Arduino boards that use this core.
                Each board is represented as a dict.

        Raises
        ------
        FileNotFoundError
            If the Arduino CLI was not found.
        subprocess.CalledProcessError
            If the Arduino CLI failed to detect installed cores.

        Emits
        -----
        error_occurred(ViPErLEEDFirmwareError.ERROR_ARDUINO_CLI_NOT_FOUND)
            If the Arduino CLI was not found.
        error_occurred(ViPErLEEDFirmwareError.ERROR_ARDUINO_CLI_FAILED)
            If the Arduino CLI failed to detect installed cores.
        """
        try:
            cli = self.get_arduino_cli()
        except FileNotFoundError:
            base.emit_error(
                self,
                ViPErLEEDFirmwareError.ERROR_ARDUINO_CLI_NOT_FOUND,
                self.base_path,
                )
            raise

        try:
            cores_json = subprocess.run(
                [cli, 'core', 'list', '--format', 'json'],
                capture_output=True, check=True
                )
        except subprocess.CalledProcessError as err:
            self.on_arduino_cli_failed(err)
            raise

        return json.loads(cores_json.stdout)

    @qtc.pyqtSlot()
    def is_cli_installed(self):
        """Check if Arduino CLI is installed.

        Emits
        -----
        cli_found(bool, bool)
            Whether the Arduino CLI was found on the PC and whether
            the CLI needs to be updated to be compatible with the
            package.
        """
        try:
            self.get_arduino_cli()
        except FileNotFoundError:
            self.cli_found.emit(False, True)
            return
        version = self.get_installed_cli_version()
        # Check if the version is at least 0.19.0, as older versions
        # are not compatible with the FirmwareUploader. We need the
        # keyword matching_boards in the dictionary containing the
        # detected boards.
        self.cli_found.emit(True, version < _MIN_CLI_VERSION)

    def on_arduino_cli_failed(self, err):
        """Report an Arduino CLI process failure."""
        base.emit_error(
            self,
            ViPErLEEDFirmwareError.ERROR_ARDUINO_CLI_FAILED,
            err.returncode,
            err.stderr.decode()
            )

    @qtc.pyqtSlot()
    def update_cli_path_from_settings(self):
        """Read the base path for the CLI from system settings."""
        path_getter = settings.SystemSettings()
        cli_path = path_getter.get('PATHS', 'arduino_cli', fallback=None)
        if cli_path:
            self.base_path = Path(cli_path)


class ArduinoCLIInstaller(ArduinoCLI):
    """Worker that handles downloads in another thread."""

    # Emitted when downloading and installing the Arduino CLI, cores and
    # libraries. Sent boolean is true when the full installation
    # succeeded.
    cli_installation_finished = qtc.pyqtSignal(bool)

    # Emitted with a value that is indicative of how much of
    # the currently running process has been completed.
    progress_occurred = qtc.pyqtSignal(int)

    def __init__(self, parent=None):
        """Initialise the Arduino CLI downloader."""
        super().__init__(parent=parent)
        # Notice that we do not give self as a parent for self.network.
        # That's because moving a QNetworkAccessManager automatically
        # to a new thread together with its parent (at is normally the
        # case with parent-child relationships) seems to be broken in
        # Qt5. We move self.network explicitly in self.moveToThread.
        # Assigning a parent here would break the overridden
        # self.moveToThread, as children with a parent cannot be moved
        # 'independently'.
        self.network = qtn.QNetworkAccessManager()
        self.network.setTransferTimeout(timeout=2000)
        # _archive_name is the name of the OS specific Arduino
        # CLI version that has to be downloaded from github for
        # installation. It is determined automatically before
        # downloading the CLI. It is also the file name of the
        # zipped CLI.
        self._archive_name = ''

        # Whenever an error occurs perform _on_install_failed.
        self.error_occurred.connect(self._on_install_failed)

    @qtc.pyqtSlot(qtn.QNetworkReply)
    def _download_newest_cli(self, reply):
        """Download newest version of Arduino CLI from github.

        Download newest version of the CLI for the OS if it doesn't
        match the currently installed version. Downloading the CLI
        is skipped if it is already installed.

        Parameters
        ----------
        reply : QNetworkReply
            Contains the newest CLI version names and download
            urls for these CLIs for various platforms.

        Emits
        -----
        progress_occurred
            Emitted after downloading version data, extracting version
            data, and getting installed version. Last emit is either
            done immediately if installed version matches the newest
            available version, or after initiating download of newest
            version if installed version is outdated.
        error_occurred(ViPErLEEDFirmwareError.ERROR_INSTALL_FAILED)
            If the QNetworkReply from the QNetworkAccessManager
            contains an error.
        error_occurred(ViPErLEEDFirmwareError.ERROR_NO_SUITABLE_CLI)
            If the OS is not supported by any precompiled
            Arduino CLI version.
        cli_installation_finished(False)
            If any error_occurred.
        """
        base.safe_disconnect(self.network.finished, self._download_newest_cli)

        # Check if connection failed.
        if reply.error():
            base.emit_error(
                self,
                ViPErLEEDFirmwareError.ERROR_INSTALL_FAILED
                )
            return
        self.progress_occurred.emit(10)

        latest = json.loads(reply.readAll().data())

        platform = sys.platform
        if 'darwin' in platform:
            os_name = 'macOS'
        elif 'win' in platform and 'cyg' not in platform:
            os_name = 'Windows'
        elif 'linux' in platform:
            os_name = 'Linux'
        else:
            base.emit_error(
                self,
                ViPErLEEDFirmwareError.ERROR_NO_SUITABLE_CLI
                )
            return

        url_latest = ''
        for asset in latest['assets']:
            # This should always pick the 32 bit version if
            # present, as it comes alphabetically earlier
            if os_name in asset['name']:
                newest_version = asset['name'].split('_')[1]
                url_latest = asset['browser_download_url']
                self._archive_name = asset['name']
                break
        else:
            base.emit_error(
                self,
                ViPErLEEDFirmwareError.ERROR_NO_SUITABLE_CLI
                )
            return
        self.progress_occurred.emit(15)

        installed_ver = self.get_installed_cli_version()
        if newest_version == installed_ver:
            self.progress_occurred.emit(60)
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
        """Download AVR core and upgrade existing cores and libraries.

        Emits
        -----
        progress_occurred
            Emitted after the Arduino cores have been installed and
            after the cores and libraries have been upgraded.
        cli_installation_finished(True)
            If the installation succeeded.
        """
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
            Contains the Arduino CLI in a .zip archive.

        Emits
        -----
        progress_occurred
            Emitted after newest CLI version has been downloaded
            and after it has been installed.
        """
        base.safe_disconnect(self.network.finished,
                             self._install_arduino_cli)
        self.progress_occurred.emit(40)
        with open(self.base_path / self._archive_name, 'wb') as archive:
            archive.write(reply.readAll())
        shutil.unpack_archive(self.base_path / self._archive_name,
                              self.base_path)
        # self._archive_name is no longer needed.
        self._archive_name = ''
        self.progress_occurred.emit(60)
        self._install_and_upgrade_cores()

    def _install_arduino_core(self, core_name):                                # TODO: add progress bar
        """Install a given core to the Arduino CLI.

        Parameters
        ----------
        core_name : str
            Should be one of the Arduino core names (e.g., arduino:avr).

        Emits
        -----
        error_occurred(ViPErLEEDFirmwareError.ERROR_ARDUINO_CLI_NOT_FOUND)
            If the Arduino CLI was not found.
        error_occurred(ViPErLEEDFirmwareError.ERROR_ARDUINO_CLI_FAILED)
            If installing the Arduino cores failed. This might
            mean that core_name is not a known core.
        cli_installation_finished(False)
            If any error_occurred.
        """
        try:
            cli = self.get_arduino_cli()
        except FileNotFoundError:
            base.emit_error(
                self,
                ViPErLEEDFirmwareError.ERROR_ARDUINO_CLI_NOT_FOUND,
                self.base_path
                )
            return

        try:
            subprocess.run([cli, 'core', 'install', core_name],
                           capture_output=True, check=True)
        except subprocess.CalledProcessError as err:
            self.on_arduino_cli_failed(err)

    @qtc.pyqtSlot(tuple)
    def _on_install_failed(self, _):
        """Installation of Arduino CLI failed."""
        self.cli_installation_finished.emit(False)

    @qtc.pyqtSlot()
    def _upgrade_arduino_cores_and_libraries(self):
        """Upgrade the outdated cores and libraries.

        Emits
        -----
        progress_occurred
            Emitted after cores and libraries have been updated.
        error_occurred(ViPErLEEDFirmwareError.ERROR_ARDUINO_CLI_NOT_FOUND)
            If the Arduino CLI was not found.
        error_occurred(ViPErLEEDFirmwareError.ERROR_INSTALL_FAILED)
            If the update method of the CLI failed. This may mean the
            internet connection may have an issue.
        error_occurred(ViPErLEEDFirmwareError.ERROR_ARDUINO_CLI_FAILED)
            If the Arduino CLI failed at upgrading the Arduino CLI.
            In this case the issue should not be the internet as the
            update succeeded. The Arduino CLI might be broken.
        cli_installation_finished(False)
            If any error_occurred.
        """
        try:
            cli = self.get_arduino_cli()
        except FileNotFoundError:
            base.emit_error(
                self,
                ViPErLEEDFirmwareError.ERROR_ARDUINO_CLI_NOT_FOUND,
                self.base_path
                )
            return

        try:
            subprocess.run([cli, 'update'], capture_output=True, check=True)
        except subprocess.CalledProcessError:
            base.emit_error(
                self,
                ViPErLEEDFirmwareError.ERROR_INSTALL_FAILED
                )
            return

        self.progress_occurred.emit(90)

        try:
            subprocess.run([cli, 'upgrade'], capture_output=True, check=True)
        except subprocess.CalledProcessError as err:
            self.on_arduino_cli_failed(err)

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

        Emits
        -----
        progress_occurred
            Emitted at the beginning to signal the start of the
            Arduino CLI installation/upgrade process.
        """
        self.progress_occurred.emit(0)
        request = qtn.QNetworkRequest(qtc.QUrl(
            'https://api.github.com/repos/arduino/arduino-cli/releases/latest'
            ))
        self.network.finished.connect(self._download_newest_cli,
                                      type=qtc.Qt.UniqueConnection)
        self.network.get(request)

    def moveToThread(self, thread):     # pylint: disable=invalid-name
        """Move self and self.network to thread."""
        super().moveToThread(thread)
        self.network.moveToThread(thread)


class FirmwareUploader(ArduinoCLI):
    """Worker that handles firmware uploads in another thread."""

    # Emitted when the upload of firmware to the controller is finished
    # or when it has failed.
    upload_finished = qtc.pyqtSignal()

    # Emitted when the Arduino CLI fails.
    cli_failed = qtc.pyqtSignal()

    # Emitted after connected controllers have been detected. Carries
    # with it the connected ViPErLEED controllers.
    controllers_detected = qtc.pyqtSignal(dict)

    # Emitted with a value that is indicative of how much of
    # the currently running process has been completed.
    progress_occurred = qtc.pyqtSignal(int)

    def _check_missing_cores(self):
        """Return the necessary cores that are missing.

        Returns
        -------
        missing : list of str
            List of the missing cores.
        """
        cores = self.get_installed_cores()

        missing_cores = []
        for required_core in REQUIRED_CORES:
            if not any(core['id'] == required_core for core in cores):
                missing_cores.append(required_core)
        return missing_cores

    def _controller_missing(self, port):
        """Return whether there is no controller on the given port.

        Returns
        -------
        missing : bool
            True if the controller is not present.

        Emits
        -----
        error_occurred(ViPErLEEDFirmwareError.ERROR_CONTROLLER_NOT_FOUND)
            If there was no available controller on the given port.
        controllers_detected(available_ctrls)
            If the controller is indeed missing, the remaining
            controllers are emitted.
        """
        available_ctrls = self.get_viperleed_hardware(False)
        selected_ctrl_exists = any(self.ctrls_with_port(available_ctrls, port))
        if selected_ctrl_exists:
            return False
        base.emit_error(self,
                        ViPErLEEDFirmwareError.ERROR_CONTROLLER_NOT_FOUND,
                        port)
        self.controllers_detected.emit(available_ctrls)
        return True

    def _get_boards(self):
        """Get a list of the available Arduino boards.

        Returns
        -------
        boards : list of dicts
            A list that contains each matching board as a dict.
            Each dict representing a controller has the following
            {key: value} pairs:
            'matching_boards' : list of dict
                Holds the board name and the fully qualified board name.
            'port' : dict
                Holds information about the port, most importantly the
                address.

        Emits
        -----
        error_occurred(ViPErLEEDFirmwareError.ERROR_ARDUINO_CLI_NOT_FOUND)
            If the Arduino CLI was not found.
        error_occurred(ViPErLEEDFirmwareError.ERROR_ARDUINO_CLI_FAILED)
            If the Arduino CLI failed to detect the connected boards.
        cli_failed()
            If the Arduino CLI process failed.
        """
        try:
            cli = self.get_arduino_cli()
        except FileNotFoundError:
            base.emit_error(
                self,
                ViPErLEEDFirmwareError.ERROR_ARDUINO_CLI_NOT_FOUND,
                self.base_path
                )
            self.cli_failed.emit()
            return []

        try:
            boards_json = subprocess.run(
                [cli, 'board', 'list', '--format', 'json'],
                capture_output=True, check=True
                )
        except subprocess.CalledProcessError as err:
            self.on_arduino_cli_failed(err)
            self.cli_failed.emit()
            return []
        boards = json.loads(boards_json.stdout)
        return [b for b in boards if 'matching_boards' in b]

    @qtc.pyqtSlot(dict, FirmwareVersionInfo)
    def compile(self, selected_ctrl, firmware):
        """Compile firmware and upload it to the specified board.

        Parameters
        ----------
        selected_ctrl : dict
            A dict representing the selected controller. It contains the
            controller name, COM port, and fully qualified board name.
        firmware : FirmwareVersionInfo
            A namedtuple representing the selected firmware that is
            to be uploaded to the selected controller. Contains the
            path from where the firmware is taken.

        Emits
        -----
        progress_occurred
            Emitted at the start, after checking if all required Arduino
            cores are installed, after checking if the selected
            controller is present, after uploading the firmware, and
            at the end.
        error_occurred(ViPErLEEDFirmwareError.ERROR_ARDUINO_CLI_NOT_FOUND)
            If the Arduino CLI was not found.
        error_occurred(ViPErLEEDFirmwareError.ERROR_CORE_NOT_FOUND)
            If a required Arduino core was missing.
        error_occurred(ViPErLEEDFirmwareError.ERROR_CONTROLLER_NOT_FOUND)
            If there was no available controller on the given port.
        error_occurred(ViPErLEEDFirmwareError.ERROR_ARDUINO_CLI_FAILED)
            If the Arduino CLI failed to upload the selected
            firmware to the controller.
        controllers_detected(available_ctrls)
            If the requested controller was no longer present.
            Contains the detected controllers.
        cli_failed()
            If any Arduino CLI process failed, or if a required Arduino
            core was missing.
        upload_finished()
            Emitted if the selected controller was no longer connected
            when trying upload the firmware, or if the upload was
            successful.
        """
        self.progress_occurred.emit(0)
        try:
            cli = self.get_arduino_cli()
        except FileNotFoundError:
            base.emit_error(
                self,
                ViPErLEEDFirmwareError.ERROR_ARDUINO_CLI_NOT_FOUND,
                self.base_path
                )
            self.cli_failed.emit()
            return

        # Check if the required cores are among the installed cores.
        try:
            missing_cores = self._check_missing_cores()
        except(subprocess.CalledProcessError, FileNotFoundError):
            self.cli_failed.emit()
            return
        if missing_cores:
            base.emit_error(self,
                            ViPErLEEDFirmwareError.ERROR_CORE_NOT_FOUND,
                            missing_cores)
            self.cli_failed.emit()
            return
        self.progress_occurred.emit(10)

        # Check if the selected controller is present at all.
        if self._controller_missing(selected_ctrl['port']):
            self.upload_finished.emit()
            return
        self.progress_occurred.emit(20)

        firmware_file = firmware.path / firmware.folder_name / 'viper-ino'

        argv = ('compile',
                '--clean',  # Remove compiled binary when done
                '--fqbn', selected_ctrl['fqbn'],
                firmware_file,
                # Also upload to board after compilation
                '--upload', '--port', selected_ctrl['port'])

        # subprocess.run([cli, *argv, '--verbose'], capture_output=True, check=True) TODO: track upload progress and adjust progress bar accordingly
        try:
            subprocess.run([cli, *argv], capture_output=True, check=True)
        except subprocess.CalledProcessError as err:
            self.on_arduino_cli_failed(err)
            # No return here to ensure that buttons are enabled and
            # tmp folder (in FirmwareArchiveUploader) is removed.

        # Update controller list
        self.progress_occurred.emit(80)
        self.get_viperleed_hardware(True)
        self.progress_occurred.emit(100)
        self.upload_finished.emit()

    @qtc.pyqtSlot(bool)
    def get_viperleed_hardware(self, detect_viperino):
        """Return a dict of the ViPErLEED Arduino boards.

        In fact, it returns all the ViPErLEED as well as Arduino Micro
        boards.

        Parameters
        ----------
        detect_viperino : bool
            If True, a check, whether the detected boards already have
            ViPErLEED firmware installed, will be performed. viperino
            controllers are renamed to their respective unique name
            in the returned dict and their firmware version is set. If
            detect_viperino is False, then the controllers_detected
            signal will not be emitted after a successful controller
            detection. If no controllers were detected, the signal is
            emitted with an empty dict to force a reset of the
            displayed controllers on the dialog side.

        Returns
        -------
        ctrl_dict : dict
            A dict of the detected Arduino Micro boards containing dicts
            with information about the controllers. keys are unique
            names of the detected controllers, including their address,
            with format '<controller name> (<address>)'. Values are
            dictionaries with the following {key: value} pairs:
            'port': str
                COM port address
            'name': str
                Board name
            'fqbn': str
                Fully qualified board name
            'version': hardwarebase.Version or str
                Firmware version of the device, if the information is
                available, otherwise str

        Emits
        -----
        controllers_detected(ctrl_dict)
            Contains the detected ViPErLEED controllers. Emitted if
            detect_viperino is True or if no boards have been detected.
        """
        viperleed_names = ('ViPErLEED', 'Arduino Micro')
        boards = self._get_boards()
        ctrl_dict = {}

        if not boards:
            self.controllers_detected.emit(ctrl_dict)
            return ctrl_dict

        # Get all available Arduino Micro controllers.
        board_names = [b['matching_boards'][0]['name'] for b in boards]
        viper_boards = []
        for name in viperleed_names:
            for board, board_name in zip(boards, board_names):
                if name in board_name:
                    viper_boards.append(board)

        # Extract data from controller list.
        for board in viper_boards:
            port = board['port']['address']
            board = board['matching_boards'][0]
            ctrl = f'{board["name"]} ({port})'
            ctrl_dict[ctrl] = {
                'port': port,
                'name': board['name'],
                'fqbn': board['fqbn'],
                'version': NOT_SET,
                }

        if not detect_viperino:
            return ctrl_dict

        # Detect ViPErLEED controllers.
        for name, (cls, info) in base.get_devices('controller').items():
            port = info.get('address')
            ctrl = next(self.ctrls_with_port(ctrl_dict, port), None)
            if not ctrl:
                continue
            ctrl_dict[ctrl]['version'] = info.get('firmware', NOT_SET)
            if not getattr(cls, 'box_id', None):
                continue
            # Notice the -2: the last two entries in name are the
            # serial number and the '(COM<port>)' bits, which we
            # don't need.
            ctrl_dict[ctrl]['name'] = '_'.join(name.split()[:-2])
            ctrl_dict[name] = ctrl_dict.pop(ctrl)

        self.controllers_detected.emit(ctrl_dict)
        return ctrl_dict

    @staticmethod
    def ctrls_with_port(ctrls, port):
        """Return a generator selecting the controllers with given port.

        Parameters
        ----------
        ctrls : dict
            A dictionary of controllers.
        port : str
            The desired COM port.

        Returns
        -------
        generator
            The controllers at the desired port.
        """
        return (ctrl for ctrl, ctrl_dict in ctrls.items()
                if ctrl_dict.get('port') == port)


class FirmwareArchiveUploader(FirmwareUploader):
    """Worker that handles archived firmware uploads in another thread."""

    @qtc.pyqtSlot(dict, FirmwareVersionInfo)
    def compile(self, selected_ctrl, firmware):
        """Extract and compile viper-ino for the specified board.

        Parameters
        ----------
        selected_ctrl : dict
            A dict representing the selected controller. It contains the
            controller name, COM port, and fully qualified board name.
        firmware : FirmwareVersionInfo
            A namedtuple representing the selected firmware that is
            to be uploaded to the selected controller. Contains the
            path from where the firmware is taken.

        Emits
        -----
        progress_occurred
            Emitted at the start, after extracting the firmware, after
            checking if all required Arduino cores are installed, after
            checking if the selected controller is present, after
            uploading the firmware, and at the end.
        error_occurred(ViPErLEEDFirmwareError.ERROR_ARDUINO_CLI_NOT_FOUND)
            If the Arduino CLI was not found.
        error_occurred(ViPErLEEDFirmwareError.ERROR_CORE_NOT_FOUND)
            If a required Arduino core was missing.
        error_occurred(ViPErLEEDFirmwareError.ERROR_CONTROLLER_NOT_FOUND)
            If there was no available controller on the given port.
        error_occurred(ViPErLEEDFirmwareError.ERROR_ARDUINO_CLI_FAILED)
            If the Arduino CLI failed to upload the selected
            firmware to the controller.
        controllers_detected(available_ctrls)
            If the requested controller was no longer present.
            Contains the detected controllers.
        cli_failed()
            If the Arduino CLI process failed, or if a required Arduino
            core was missing.
        upload_finished()
            Emitted if the selected controller was no longer connected
            when trying upload the firmware, or if the upload was
            successful.
        """
        # The folder to which the selected archive is extracted.
        tmp_path = firmware.path.parent / 'tmp_'
        with ZipFile(firmware.path) as firmware_zip:
            firmware_zip.extractall(tmp_path)
        self.progress_occurred.emit(5)
        firmware = firmware._replace(path=tmp_path)
        super().compile(selected_ctrl, firmware)
        shutil.rmtree(tmp_path)
