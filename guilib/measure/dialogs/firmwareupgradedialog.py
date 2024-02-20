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
                   QtCore as qtc)
import requests             # Look at https://wiki.qt.io/Download_Data_from_URL

from viperleed.guilib.measure.widgets.pathselector import PathSelector
from viperleed.guilib.measure import hardwarebase as base


_INVOKE = qtc.QMetaObject.invokeMethod
NOT_SET = '\u2014'


FirmwareVersionInfo = namedtuple('FirmwareVersionInfo',
                                 ('folder_name', 'version', 'path'))


class FirmwareUpgradeDialog(qtw.QDialog):
    """Dialog to handle user interaction when upgrading firmware."""

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
                'upload': qtw.QPushButton('&Upload'),
                'done': qtw.QPushButton('&Done'),
                },
            'controller_info': {
                'ctrl_type': qtw.QLabel(f'Controller type: {NOT_SET}'),
                'firmware_version': qtw.QLabel(
                    f'Installed firmware version: {NOT_SET}'),
                'highest_version': qtw.QLabel(
                    f'Most recent firmware version: {NOT_SET}'),
                },
            'timers': {
                #TODO: maybe add timers (trigger detect devices instead of button)
                },
            }

        self.__uploader = FirmwareUploader()
        self.__upload_thread = qtc.QThread()
        self.__uploader.moveToThread(self.__upload_thread)
        self.__upload_thread.start()

        self.setWindowTitle('Upgrade ViPErLEED box firmware')
        self.__children['ctrls']['firmware_path'].path = Path().resolve()
        self.__compose()
        self.__connect()

    def __compose(self):
        """Place children widgets."""
        # Remove the '?' button from the title bar
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
        layout.addLayout(self.__compose_done_button())
        self.setLayout(layout)

    def __compose_info_section(self):
        layout = qtw.QHBoxLayout()
        for lbl in self.__children['controller_info'].values():
            layout.addWidget(lbl)
        return layout

    def __compose_controller_selection(self):
        layout = qtw.QHBoxLayout()
        layout.addWidget(qtw.QLabel('Select controller:'))
        layout.addWidget(self.__children['ctrls']['controllers'], stretch=1)
        layout.addWidget(self.__children['buttons']['refresh'])
        return layout

    def __compose_path_selection(self):
        layout = qtw.QHBoxLayout()
        layout.addWidget(qtw.QLabel('Select firmware folder:'))
        layout.addWidget(self.__children['ctrls']['firmware_path'])
        return layout

    def __compose_firmware_selection(self):
        layout = qtw.QHBoxLayout()
        layout.addWidget(qtw.QLabel('Select firmware version:'))
        layout.addWidget(self.__children['ctrls']['firmware_version'],
                         stretch=1)
        layout.addWidget(self.__children['buttons']['upload'])
        self.__children['buttons']['upload'].setEnabled(False)
        return layout

    def __compose_done_button(self):
        layout = qtw.QHBoxLayout()
        layout.addStretch(1)
        layout.addWidget(self.__children['buttons']['done'])
        return layout

    def __connect(self):
        """Connect children signals."""
        self.__children['buttons']['done'].clicked.connect(self.accept)
        self.__children['buttons']['refresh'].clicked.connect(self.__detect)
        self.__children['buttons']['upload'].clicked.connect(self.__upload)
        self.__children['ctrls']['firmware_path'].path_changed.connect(
            self.__get_firmware_versions
            )
        self.__children['ctrls']['controllers'].currentTextChanged.connect(
                    self.__update_ctrl_labels
            )
        self.__uploader.controllers_detected.connect(self.__update_combo_box)
        self.__uploader.controllers_detected.connect(self.__allow_refresh)
        self.__uploader.upload_finished.connect(self.__ctrl_enable)

    def accept(self):
        """Clean up, then accept."""
        self.__clean_up()
        super().accept()

    def __clean_up(self):
        """Clean up before closing dialog."""
        # No-op for now

    @qtc.pyqtSlot(str, dict)
    def __allow_refresh(self, *args):
        """Allow refreshing and uploading again."""
        for btn in self.__children['buttons']:
            self.__children['buttons'][btn].setEnabled(True)

    @qtc.pyqtSlot(bool)
    def __ctrl_enable(self, enable):
        """Enable/disable controls."""
        for btn in self.__children['buttons']:
            self.__children['buttons'][btn].setEnabled(enable)
        for slct in self.__children['ctrls']:
            self.__children['ctrls'][slct].setEnabled(enable)

    def __detect(self):
        """Detect connected ViPErLEED controllers."""
        emit_controllers = True
        for btn in self.__children['buttons']:
            self.__children['buttons'][btn].setEnabled(False)
        _INVOKE(self.__uploader, 'get_viperleed_hardware',
                qtc.Q_ARG(bool, emit_controllers))

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

    @qtc.pyqtSlot(str, dict)
    def __update_combo_box(self, which_combo, data_dict):
        """Replace displayed items with the detected ones."""
        self.__children['ctrls'][which_combo].clear()
        for key, value in data_dict.items():
            self.__children['ctrls'][which_combo].addItem(key, userData=value)
        ctrl = self.__children['ctrls']['controllers'].currentData()
        version = self.__children['ctrls']['firmware_version'].currentData()
        if not ctrl or not version:
            self.__children['buttons']['upload'].setEnabled(False)
        else:
            self.__children['buttons']['upload'].setEnabled(True)

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

        _INVOKE(self.__uploader, '__compile', qtc.Q_ARG(dict, selected_ctrl),
                qtc.Q_ARG(FirmwareVersionInfo, firmware),
                qtc.Q_ARG(Path, tmp_path), qtc.Q_ARG(bool, upload))
        self.__ctrl_enable(False)


class FirmwareUploader(qtc.QObject):
    """Worker that handles firmware uploads in another thread."""
    upload_finished = qtc.pyqtSignal(bool)

    controllers_detected = qtc.pyqtSignal(str, dict)

    def __init__(self, parent=None):
        """Initialise the firmware uploader."""
        self.base_path = Path().resolve() / 'hardware/arduino/arduino-cli'
        super().__init__(parent=parent)

    @qtc.pyqtSlot(dict, FirmwareVersionInfo, Path, bool)
    def __compile(self, s_ctrl, firmware, tmp_path, upload):
        """Compile viper-ino for the specified board.

        Parameters
        ----------
        """
        # Check if the selected controller is present at all.
        available_ctrls = self.get_viperleed_hardware(False)

        if not any(available_ctrls[c]['name'] == s_ctrl['name'] and
                   available_ctrls[c]['port'] == s_ctrl['port']
                   for c in available_ctrls):
            print('Selected controller is no longer present.')
            self.controllers_detected.emit('controllers', available_ctrls)
            self.upload_finished.emit(True)
            return

        cli = self.get_arduino_cli()
        with ZipFile(firmware.path) as firmware_zip:
            firmware_zip.extractall(tmp_path)
        firmware_extracted = tmp_path / firmware.folder_name / 'viper-ino'
        # Add generic inner folder structure to find .ino file
        argv = ['compile', '--clean', '-b', s_ctrl['fqbn'], firmware_extracted]
        if upload:
            argv.extend(['-u', '-p', s_ctrl['port']])

        cli = subprocess.run([cli, *argv], capture_output=True)
        try:
            cli.check_returncode()
        except subprocess.CalledProcessError as err:
            raise RuntimeError('Arduino CLI failed. Return code='
                               f'{cli.returncode}. The error was:\n'
                               + cli.stderr.decode())
        # Remove extracted archive
        shutil.rmtree(tmp_path)
        # Update controller list
        self.get_viperleed_hardware(True)
        self.upload_finished.emit(True)

    def upgrade_arduino_cli(self):                                                      # TODO: add progress bar
        """Upgrade the outdated cores and libraries."""
        cli = self.get_arduino_cli()
        cli = subprocess.run([cli, 'upgrade'], capture_output=True)
        try:
            cli.check_returncode()
        except subprocess.CalledProcessError as err:
            raise RuntimeError('Arduino CLI failed. Return code='
                               f'{cli.returncode}. The error was:\n{cli.stderr}')

    def get_arduino_cores(self):
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
        cli = self.get_arduino_cli()
        cli = subprocess.run([cli, 'core', 'list', '--format', 'json'],
                              capture_output=True)
        try:
            cli.check_returncode()
        except subprocess.CalledProcessError as err:
            raise RuntimeError('Arduino CLI failed. Return code='
                               f'{cli.returncode}. The error was:\n{cli.stderr}')
        return json.loads(cli.stdout)

    def install_arduino_core(self, core_name):                                            # TODO: add progress bar
        """Install a given core to the Arduino CLI.

        Parameters
        ----------
        core_name : str
            Should be one of the Arduino core names (e.g., arduino:avr).
            It's easier to get this information from a call to
            get_viperleed_hardware()
        """
        cli = self.get_arduino_cli()
        cli = subprocess.run([cli, 'core', 'install', core_name],
                              capture_output=True)
        try:
            cli.check_returncode()
        except subprocess.CalledProcessError as err:
            raise RuntimeError(f'Failed to install {core_name}. Return code='
                               f'{cli.returncode}. The error was:\n{cli.stderr}')

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
        latest = requests.get('https://api.github.com/repos/arduino'
                              '/arduino-cli/releases/latest').json()
        platform = sys.platform
        if 'darwin' in platform:
            correct_name = 'macOS'
        elif 'win' in platform and 'cyg' not in platform:
            correct_name = 'Windows'
        elif 'linux' in platform:
            correct_name = 'Linux'
        else:
            raise RuntimeError('Could not find a suitable precompiled '
                               'version of the Arduino CLI. You may have '
                               'to compile from the source at '
                               'https://github.com/arduino/arduino-cli')
        url_latest = ''
        for asset in latest['assets']:
            # This should always pick the 32 bit version if
            # present, as it comes alphabetically earlier
            if correct_name in asset['name']:
                url_latest = asset['browser_download_url']
                correct_name = asset['name']
                break
        else:
            raise RuntimeError('Could not find a suitable precompiled '
                               'version of the Arduino CLI. You may have '
                               'to compile from the source at '
                               'https://github.com/arduino/arduino-cli')

        with open(self.base_path / correct_name, 'wb') as archive:
            archive.write(requests.get(url_latest).content)
        shutil.unpack_archive(self.base_path / correct_name, self.base_path)

    def get_arduino_cli(self, get_from_git=False):                                        # TODO: add progress bar
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
                raise RuntimeError('Arduino command-line interface not found. '
                                   f'It should be available in {self.base_path}, '
                                   'but the directory does not exists.')

        # See if there is the correct executable in arduino-cli
        arduino_cli = 'arduino-cli'
        if 'win' in sys.platform and not 'darwin' in sys.platform:
            arduino_cli += '.exe'
        arduino_cli = self.base_path.joinpath(arduino_cli)
        if not arduino_cli.is_file():
            raise RuntimeError('Arduino CLI not found. You can download it '
                               'from https://github.com/arduino/arduino-cli '
                               f'and unpack it inside {self.base_path}')
        return arduino_cli

    def get_boards(self):
        """Get a dict of the available Arduino boards.

        Returns
        -------
        dict
        """
        cli = self.get_arduino_cli()
        cli = subprocess.run([cli, 'board', 'list', '--format', 'json'],
                              capture_output=True)
        try:
            cli.check_returncode()
        except subprocess.CalledProcessError as err:
            raise RuntimeError(
                'Arduino CLI failed with return code '
                f'{cli.returncode}. The error was:\n{cli.stderr}'
                ) from err
        boards = json.loads(cli.stdout)
        return [b for b in boards if 'matching_boards' in b]

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
        boards = self.get_boards()

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





