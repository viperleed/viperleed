"""Upload ViPErLEED sketch via Arduino CLI.

Provides some functionality to interface with the Arduino command-line
interface tool that allows compiling and uploading the viper-ino.ino
sketch to the hardware.

The Arduino CLI is available from https://github.com/arduino/arduino-cli/releases
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__created__ = '2021-05-09'

import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import warnings

import requests


def get_arduino_cli_from_git():
    """Obtain the latest version of Arduino CLI from GitHub.

    Download and extract the latest version of the Arduino
    command-line interface executables for the current platform.
    For simplicity, the 32bit version is downloaded for both Linux
    and Windows.

    Returns
    -------
    None.
    """
    latest = requests.get("https://api.github.com/repos/arduino"
                          "/arduino-cli/releases/latest").json()
    if 'nt' in os.name:
        correct_name = 'Windows'
    elif 'mac' in os.name:
        correct_name = 'macOS'
    elif 'posix' in os.name:
        correct_name = 'Linux'
    else:
        raise RuntimeError("Could not find a suitable precompiled "
                           "version of the Arduino CLI. You may have"
                           "to compile from the source at "
                           "https://github.com/arduino/arduino-cli")
    url_latest = ''
    for asset in latest['assets']:
        # This should always pick the 32 bit version if
        # present, as it comes alphabetically earlier
        if correct_name in asset['name']:
            url_latest = asset['browser_download_url']
            correct_name = asset['name']
            break
    if not url_latest:
        raise RuntimeError("Could not find a suitable precompiled "
                           "version of the Arduino CLI. You may have "
                           "to compile from the source at "
                           "https://github.com/arduino/arduino-cli")

    base_path = Path(__file__).parent.resolve() / 'arduino-cli'
    with open(base_path / correct_name, 'wb') as archive:
        archive.write(requests.get(url_latest).content)
    shutil.unpack_archive(base_path / correct_name, base_path)


def get_arduino_cli(get_from_git=False):
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
    base_path = Path(__file__).parent.resolve() / 'arduino-cli'
    try:
        base_path.mkdir()
    except FileExistsError:
        # folder is there, it may contain the executable
        pass
    else:
        if get_from_git:
            # Get the executables from github
            get_arduino_cli_from_git()
        else:
            raise RuntimeError("Arduino command-line interface not found. "
                               f"It should be available in {base_path}, "
                               "but the directory does not exists.")

    # See if there is the correct executable in arduino-cli
    arduino_cli = 'arduino-cli'
    if 'win' in sys.platform:
        arduino_cli += '.exe'
    arduino_cli = base_path.joinpath(arduino_cli)
    if not arduino_cli.is_file():
        raise RuntimeError("Arduino CLI not found. You can download it "
                           "from https://github.com/arduino/arduino-cli "
                           f"and unpack it inside {base_path}")
    return arduino_cli


def get_boards():
    """Get a dict of the available Arduino boards.

    Returns
    -------
    dict
        TODO: keys and values. I have to try out with a connected board
    """
    cli = get_arduino_cli()
    cli = subprocess.run([cli, 'board', 'list', '--format', 'json'],
                          capture_output=True)
    try:
        cli.check_returncode()
    except subprocess.CalledProcessError as err:
        raise RuntimeError("Arduino CLI failed with return code "
                           f"{cli.returncode}. The error was:\n{cli.stderr}")
    return json.loads(cli.stdout)


def install_arduino_core(core_name):
    """Install a given core to the Arduino CLI.

    Parameters
    ----------
    core_name : str
        Should be one of the Arduino core names (e.g., arduino:avr).
        It's easier to get this information from a call to
        get_viperleed_hardware()
    """
    cli = get_arduino_cli()
    cli = subprocess.run([cli, 'core', 'install', core_name],
                          capture_output=True)
    try:
        cli.check_returncode()
    except subprocess.CalledProcessError as err:
        raise RuntimeError(f"Failed to install {core_name}. Return code="
                           f"{cli.returncode}. The error was:\n{cli.stderr}")


def get_arduino_cores():
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
    cli = get_arduino_cli()
    cli = subprocess.run([cli, 'core', 'list', '--format', 'json'],
                          capture_output=True)
    try:
        cli.check_returncode()
    except subprocess.CalledProcessError as err:
        raise RuntimeError("Arduino CLI failed. Return code="
                           f"{cli.returncode}. The error was:\n{cli.stderr}")
    return json.loads(cli.stdout)


def upgrade_arduino_cli():
    """Upgrade the outdated cores and libraries."""
    cli = get_arduino_cli()
    cli = subprocess.run([cli, 'upgrade'], capture_output=True)
    try:
        cli.check_returncode()
    except subprocess.CalledProcessError as err:
        raise RuntimeError("Arduino CLI failed. Return code="
                           f"{cli.returncode}. The error was:\n{cli.stderr}")


def get_viperleed_hardware():
    """Return a list of the ViPErLEED Arduino boards.

    In fact, it returns all the ViPErLEED as well as Arduino Micro
    boards. After calling this function, one should decide whether
    to keep only those boards whose 'name' contains 'ViPErLEED'.
    """
    viperleed_names = ('ViPErLEED', 'Arduino Micro')
    boards = get_boards()
    if not boards:
        raise RuntimeError("No Arduino boards connected.")
    board_names = [board['name'] for board in boards]
    viper_boards = []
    for name in viperleed_names:
        for i, board_name in enumerate(board_names):
            if name in board_name:
                viper_boards.append(boards[i])

    if not viper_boards:
        raise RuntimeError("Could not find ViPErLEED hardware. "
                           "Check that it is connected.")
    return viper_boards


def compile(for_board, upload=False):
    """Compile viper-ino for the specified board.

    Parameters
    ----------
    for_board : dict??
        The board for which we want to compile the firmware
    upload : bool
        Decide whether the sketch should also be uploaded
        after successful compilation
    """
    cli = get_arduino_cli()
    viperino = Path(__file__).parent.resolve() / 'viper-ino.ino'

    argv = ['compile', '--clean', '-b', for_board['fqbn'], viperino]
    if upload:
        argv.append(['-u', '-p', for_board['port']]

    cli = subprocess.run([cli, *argv], capture_output=True)
    try:
        cli.check_returncode()
    except subprocess.CalledProcessError as err:
        raise RuntimeError("Arduino CLI failed. Return code="
                           f"{cli.returncode}. The error was:\n{cli.stderr}")


if __name__ == '__main__':
    print(get_boards())
    install_arduino_core('arduino:avr')
    print(list_arduino_cores())
