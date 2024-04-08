"""Module ioarduinocli of viperleed.guilib.measure.classes.

========================================
   ViPErLEED Graphical User Interface
========================================

Created: 2024-04-03
Author: Florian Dörr, Michele Riva

Defines the ArduinoCLI, ArduinoCLIInstaller, and FirmwareUploader
classes that handle Arduino CLI and Arduino core downloads and Arduino
firmware upgrades.
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


class ArduinoCLI(qtc.QObject):
    """Base class that can get the Arduino CLI."""

    def __init__(self, parent=None):
        """Initialise the Arduino CLI getter."""
        super().__init__(parent=parent)
        self.base_path = Path().resolve() / 'hardware/arduino/arduino-cli'
        self.update_cli_path()

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

    @qtc.pyqtSlot()
    def update_cli_path(self):
        """Update the Arduino CLI path."""
        path_getter = settings.SystemSettings()
        cli_path = path_getter.get('PATHS', 'arduino_cli', fallback=None)
        if cli_path:
            self.base_path = Path(cli_path)


class ArduinoCLIInstaller(ArduinoCLI):
    """Worker that handles downloads in another thread."""
    error_occurred = qtc.pyqtSignal(tuple)

    # Emitted when downloading and installing the Arduino CLI, cores and
    # libraries. Sent boolean is true when the full installation
    # succeeded.
    cli_installation_finished = qtc.pyqtSignal(bool)

    # Emitted after the check if the Arduino CLI is installed.
    cli_found = qtc.pyqtSignal(bool)

    # Emitted with a value that is indicative of how much of
    # the currently running process has been completed.
    progress_occurred = qtc.pyqtSignal(int)

    def __init__(self, parent=None):
        """Initialise the Arduino CLI downloader."""
        super().__init__(parent=parent)
        self.network = qtn.QNetworkAccessManager()
        self.network.setTransferTimeout(timeout=2000)
        # _archive_name is the name of the OS specific Arduino
        # CLI version that has to be downloaded from github for
        # installation. It is determined automatically before
        # downloading the CLI. It is file name of the zipped CLI.
        self._archive_name = ''

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

        Emits
        -----
        progress_occurred
            If progress occurred.
        error_occurred(ViPErLEEDFirmwareError.ERROR_INSTALL_FAILED)
            If the QNetworkReply from the QNetworkAccessManager
            contains an error.
        error_occurred(ViPErLEEDFirmwareError.ERROR_NO_SUITABLE_CLI)
            If the OS is not supported by any precompiled
            Arduino CLI version.
        cli_installation_finished(False)
            If the installation failed.
        """
        base.safe_disconnect(self.network.finished, self._download_newest_cli)

        # Check if connection failed.
        if reply.error():
            base.emit_error(self,
                ViPErLEEDFirmwareError.ERROR_INSTALL_FAILED
                )
            self._install_failed()
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
            base.emit_error(self,
                ViPErLEEDFirmwareError.ERROR_NO_SUITABLE_CLI
                )
            self._install_failed()
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
            base.emit_error(self,
                ViPErLEEDFirmwareError.ERROR_NO_SUITABLE_CLI
                )
            self._install_failed()
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

        Emits
        -----
        error_occurred(ViPErLEEDFirmwareError.ERROR_ARDUINO_CLI_NOT_FOUND)
            If the Arduino CLI was not found.
        error_occurred(ViPErLEEDFirmwareError.ERROR_ARDUINO_CLI_FAILED)
            If the Arduino CLI failed.
        cli_installation_finished(False)
            If the installation failed.
        """
        try:
            cli = self.get_arduino_cli()
        except FileNotFoundError:
            base.emit_error(self,
                ViPErLEEDFirmwareError.ERROR_ARDUINO_CLI_NOT_FOUND,
                self.base_path
                )
            self._install_failed()
            return {}

        try:
            cores = subprocess.run([cli, 'core', 'list', '--format', 'json'],
                                   capture_output=True, check=True)
        except subprocess.CalledProcessError as err:
            base.emit_error(self,
                ViPErLEEDFirmwareError.ERROR_ARDUINO_CLI_FAILED,
                err.returncode,
                err.stderr
                )
            self._install_failed()
            return {}
        return json.loads(cores.stdout)

    def _get_installed_cli_version(self):
        """Detect version of installed Arduino CLI.

        Returns
        -------
        Version string : str
            The version of the Arduino CLI.
        """
        try:
            cli = self.get_arduino_cli()
        except FileNotFoundError:
            return '0.0.0'

        ver = subprocess.run([cli, 'version', '--format', 'json'],
                             capture_output=True, check=True)
        ver = json.loads(ver.stdout)
        return ver['VersionString']

    def _install_and_upgrade_cores(self):
        """Download AVR core and upgrade existing cores and libraries.

        Returns
        -------
        None.

        Emits
        -----
        progress_occurred
            If progress occurred.
        cli_installation_finished(True)
            If the installation succeeded.
        """
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

        Emits
        -----
        progress_occurred
            If progress occurred.
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
        self.progress_occurred.emit(50)
        self._install_and_upgrade_cores()

    def _install_arduino_core(self, core_name):                                # TODO: add progress bar
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

        Emits
        -----
        error_occurred(ViPErLEEDFirmwareError.ERROR_ARDUINO_CLI_NOT_FOUND)
            If the Arduino CLI was not found.
        error_occurred(ViPErLEEDFirmwareError.ERROR_ARDUINO_CLI_FAILED)
            If the Arduino CLI failed.
        cli_installation_finished(False)
            If the installation failed.
        """
        try:
            cli = self.get_arduino_cli()
        except FileNotFoundError:
            base.emit_error(self,
                ViPErLEEDFirmwareError.ERROR_ARDUINO_CLI_NOT_FOUND,
                self.base_path
                )
            self._install_failed()
            return

        try:
            subprocess.run([cli, 'core', 'install', core_name],
                           capture_output=True, check=True)
        except subprocess.CalledProcessError as err:
            base.emit_error(self,
                ViPErLEEDFirmwareError.ERROR_ARDUINO_CLI_FAILED,
                err.returncode,
                err.stderr
                )
            self._install_failed()

    def _install_failed(self):
        """Installation of Arduino CLI failed."""
        self.cli_installation_finished.emit(False)

    @qtc.pyqtSlot()
    def _upgrade_arduino_cores_and_libraries(self):
        """Upgrade the outdated cores and libraries.

        Returns
        -------
        None.

        Emits
        -----
        progress_occurred
            If progress occurred.
        error_occurred(ViPErLEEDFirmwareError.ERROR_ARDUINO_CLI_NOT_FOUND)
            If the Arduino CLI was not found.
        error_occurred(ViPErLEEDFirmwareError.ERROR_INSTALL_FAILED)
            If the update method of the CLI failed. This may mean the
            internet connection may have an issue.
        error_occurred(ViPErLEEDFirmwareError.ERROR_ARDUINO_CLI_FAILED)
            If the Arduino CLI failed at upgrading the Arduino CLI.
            In this case the issue should not be the internet as the
            update must have went through.
        cli_installation_finished(False)
            If the installation failed.
        """
        try:
            cli = self.get_arduino_cli()
        except FileNotFoundError:
            base.emit_error(self,
                ViPErLEEDFirmwareError.ERROR_ARDUINO_CLI_NOT_FOUND,
                self.base_path
                )
            self._install_failed()
            return

        try:
            subprocess.run([cli, 'update'], capture_output=True, check=True)
        except subprocess.CalledProcessError:
            base.emit_error(self,
                ViPErLEEDFirmwareError.ERROR_INSTALL_FAILED
                )
            self._install_failed()
            return

        self.progress_occurred.emit(90)

        try:
            subprocess.run([cli, 'upgrade'], capture_output=True, check=True)
        except subprocess.CalledProcessError as err:
            base.emit_error(self,
                ViPErLEEDFirmwareError.ERROR_ARDUINO_CLI_FAILED,
                err.returncode,
                err.stderr
                )
            self._install_failed()

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
            If progress occurred.
        """
        self.progress_occurred.emit(0)
        request = qtn.QNetworkRequest(qtc.QUrl(
            'https://api.github.com/repos/arduino/arduino-cli/releases/latest'
            ))
        self.network.finished.connect(self._download_newest_cli,
                                      type=qtc.Qt.UniqueConnection)
        self.network.get(request)

    @qtc.pyqtSlot()
    def is_cli_installed(self):
        """Check if Arduino CLI is installed.

        Emits
        -----
        cli_found(bool)
            Whether the Arduino CLI was found on the PC.
        """
        try:
            self.get_arduino_cli()
        except FileNotFoundError:
            self.cli_found.emit(False)
            return
        self.cli_found.emit(True)

    def moveToThread(self, thread):     # pylint: disable=invalid-name
        """Move self and self.network to thread."""
        super().moveToThread(thread)
        self.network.moveToThread(thread)


class FirmwareUploader(ArduinoCLI):
    """Worker that handles firmware uploads in another thread."""
    error_occurred = qtc.pyqtSignal(tuple)

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

    def _get_boards(self):
        """Get a list of the available Arduino boards.

        Returns
        -------
        list

        Emits
        -----
        error_occurred(ViPErLEEDFirmwareError.ERROR_ARDUINO_CLI_NOT_FOUND)
            If the Arduino CLI was not found.
        error_occurred(ViPErLEEDFirmwareError.ERROR_ARDUINO_CLI_FAILED)
            If the Arduino CLI failed.
        cli_failed()
            If the Arduino CLI process failed.
        """
        try:
            cli = self.get_arduino_cli()
        except FileNotFoundError:
            base.emit_error(self,
                ViPErLEEDFirmwareError.ERROR_ARDUINO_CLI_NOT_FOUND,
                self.base_path
                )
            self.cli_failed.emit()
            return []

        try:
            boards = subprocess.run([cli, 'board', 'list', '--format', 'json'],
                                    capture_output=True, check=True)
        except subprocess.CalledProcessError as err:
            base.emit_error(self,
                ViPErLEEDFirmwareError.ERROR_ARDUINO_CLI_FAILED,
                err.returncode,
                err.stderr
                )
            self.cli_failed.emit()
            return []
        boards = json.loads(boards.stdout)
        return [b for b in boards if 'matching_boards' in b]

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

        Emits
        -----
        progress_occurred
            If progress occurred.
        error_occurred(ViPErLEEDFirmwareError.ERROR_CONTROLLER_NOT_FOUND)
            If there was no available controller the given port.
        error_occurred(ViPErLEEDFirmwareError.ERROR_ARDUINO_CLI_NOT_FOUND)
            If the Arduino CLI was not found.
        error_occurred(ViPErLEEDFirmwareError.ERROR_ARDUINO_CLI_FAILED)
            If the Arduino CLI failed.
        controllers_detected(available_ctrls)
            If the requested controller was no longer present.
            Contains the detected controllers.
        cli_failed()
            If the Arduino CLI process failed.
        upload_finished()
            If the upload either finished or because the selected
            controller was no longer present.
        """
        self.progress_occurred.emit(0)
        # Check if the selected controller is present at all.
        available_ctrls = self.get_viperleed_hardware(False)
        selected_ctrl_exists = any(
            self.ctrls_with_port(available_ctrls, selected_ctrl['port']))
        if not selected_ctrl_exists:
            base.emit_error(self,
                            ViPErLEEDFirmwareError.ERROR_CONTROLLER_NOT_FOUND,
                            selected_ctrl['port'])
            self.controllers_detected.emit(available_ctrls)
            self.upload_finished.emit()
            return

        try:
            cli = self.get_arduino_cli()
        except FileNotFoundError:
            base.emit_error(self,
                ViPErLEEDFirmwareError.ERROR_ARDUINO_CLI_NOT_FOUND,
                self.base_path
                )
            self.cli_failed.emit()
            return

        with ZipFile(firmware.path) as firmware_zip:
            firmware_zip.extractall(tmp_path)
        firmware_extracted = tmp_path / firmware.folder_name / 'viper-ino'     # TODO: Add generic inner folder structure to find .ino file
        argv = ['compile', '--clean', '-b', selected_ctrl['fqbn'], firmware_extracted]
        if upload:
            argv.extend(['-u', '-p', selected_ctrl['port']])
        self.progress_occurred.emit(20)
        # subprocess.run([cli, *argv, '--verbose'], capture_output=True, check=True) TODO: track upload progress and adjust progress bar accordingly
        try:
            subprocess.run([cli, *argv], capture_output=True, check=True)
        except subprocess.CalledProcessError as err:
            base.emit_error(self,
                ViPErLEEDFirmwareError.ERROR_ARDUINO_CLI_FAILED,
                err.returncode,
                err.stderr.decode()
                )
            # No return here to ensure that tmp folder is
            # removed and buttons are enabled.
        self.progress_occurred.emit(85)
        # Remove extracted archive
        shutil.rmtree(tmp_path)
        # Update controller list
        self.progress_occurred.emit(90)
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

        Emits
        -----
        controllers_detected(ctrl_dict)
            Contains the detected ViPErLEED controllers.
        """
        viperleed_names = ('ViPErLEED', 'Arduino Micro')
        boards = self._get_boards()

        if not boards:
            self.controllers_detected.emit({})
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
            ctrl = f'{board["name"]} ({port})'
            ctrl_dict[ctrl] = {
                'port': port,
                'name': board['name'],
                'fqbn': board['fqbn'],
                'version': NOT_SET,
                'box_id': NOT_SET,
                }

        if not detect_viperino:
            return ctrl_dict

        # Detect ViPErLEED controllers.
        for name, (_, info) in base.get_devices('controller').items():
            port = info.get('address')
            ctrl = next(self.ctrls_with_port(ctrl_dict, port), None)
            if not ctrl:
                continue
            ctrl_dict[ctrl]['version'] = info.get('firmware', NOT_SET)
            box_id = info.get('box_id')
            # box_id of the measuring ViPErinoController is 0! Check
            # must be is not None because of that.
            if box_id is not None:
                ctrl_dict[ctrl]['box_id'] = box_id
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
        return (ctrl for ctrl in ctrls if ctrls[ctrl]['port'] == port)
