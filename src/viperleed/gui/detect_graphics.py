"""Module detect_graphics of viperleed.gui.

Defines functionality to detect whether the current system is capable
of displaying widgets.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-05-22'
__license__ = 'GPLv3+'

import os
from pathlib import Path
from multiprocessing import Process
import subprocess
import sys

try:
    from PyQt5.QtWidgets import QApplication
except ImportError:
    _HAS_PYQT = False
else:
    import PyQt5
    from PyQt5 import QtCore as qtc
    _HAS_PYQT = True

# How long (in seconds) to wait before deciding there's no graphics?
_TIMEOUT = 1
_UNIX_RUNTIME_DIR = os.environ.get('XDG_RUNTIME_DIR', '')


def _init_dummy_qapp():
    """Initialize a dummy QApplication instance.

    This function is only meant to be used for checking whether the
    current platform has graphical capabilities. Calling this will
    block forever on a system that has no graphical capability.

    Returns
    -------
    None.
    """
    suppress_file_permission_warnings()
    _ = QApplication([])
    sys.exit(0)


def has_graphics():
    """Return whether the system has graphical capabilities.

    This function will block for one second on systems that do not
    have graphical capabilities.

    Returns
    -------
    bool
    """
    if not _HAS_PYQT:
        return False
    proc = Process(target=_init_dummy_qapp)
    proc.start()
    proc.join(_TIMEOUT)
    if proc.is_alive():
        # The QApplication is still trying to start up. This means we
        # have no display.
        proc.terminate()
        return False
    return True


def has_pyqt():
    """Return whether PyQt5 is installed."""
    return _HAS_PYQT


def find_missing_qt_dependencies():
    """Return a dict of missing PyQt dependencies."""
    if not has_pyqt():
        raise NotImplementedError('Cannot check Qt dependencies without PyQt')
    if sys.platform.startswith('win'):
        # No evidence of missing dependencies on Windows
        return {}
    try:
        return Qt5DependencyFinder().find_missing_dependencies()
    except NotImplementedError:
        return {}


def suppress_file_permission_warnings():
    """Silence warnings concerning file permissions by QStandardPaths."""
    # The source of the warning is a too-permissive file mode on
    # the default runtime-created (once per reboot cycle) directory,
    # stored in the XDG_RUNTIME_DIR environment variable. The warning
    # is emitted since Qt 5.15.3 in the private checkXdgRuntimeDir of
    # qtbase/src/corelib/io/qstandarpaths_unix.cpp, called in
    # QStandardPaths.writableLocation(QStandardPaths.RuntimeLocation).
    # We can't really use that static method to get the path that
    # will be used, as the return value depends on whether the mode
    # is correct: a temporary directory is created in case of invalid
    # permissions. While we could reproduce the C++ code, we would
    # need to do this in a Qt-version-dependent manner. Instead,
    # it is simpler to explicitly create a ViPErLEED-GUI-dedicated
    # directory with the expected permissions (0o0700) and set
    # the XDG_RUNTIME_DIR environment variable accordingly.
    if not _UNIX_RUNTIME_DIR:
        # Not a UNIX platform. Don't bother, all the other
        # implementations don't care about file permissions
        return
    dirname = 'viperleed-gui'
    new_runtime_dir = Path(_UNIX_RUNTIME_DIR).resolve() / dirname
    try:
        new_runtime_dir.mkdir(mode=0o0700, parents=True, exist_ok=True)
    except OSError:
        # Failed to create. Probably a PermissionError.
        # Try another location that should be writable.
        standard_path = qtc.QStandardPaths
        tmp = Path(standard_path.writableLocation(standard_path.TempLocation))
        new_runtime_dir = (tmp / dirname).resolve()
        new_runtime_dir.mkdir(mode=0o0700, parents=True, exist_ok=True)
    os.environ['XDG_RUNTIME_DIR'] = str(new_runtime_dir)


class Qt5DependencyFinder:
    """Helper for finding missing dependencies of Qt5 platform plugins."""

    def __init__(self):
        """Initialize an instance."""
        self._root = (
            Path(PyQt5.__file__).resolve().parent  # The PyQt5 folder
            / 'Qt5/plugins/platforms'
            )
        assert self._root.exists()

    @classmethod
    def find_install_for_libs(cls, missing):
        """Return a suggestion on how to install `missing` libraries."""
        finder_name = f'_find_install_for_libs_{sys.platform}'
        try:
            finder = getattr(cls, finder_name)
        except AttributeError:
            return ''
        return finder(missing)

    def find_missing_dependencies(self):
        """Return a dictionary of missing dependencies for this platform."""
        finder_name = f'_find_missing_deps_{sys.platform}'
        try:
            finder = getattr(self, finder_name)
        except AttributeError:
            raise NotImplementedError('Cannot find missing Qt dependencies on '
                                      f'platform {sys.platform!r}') from None
        try:
            return finder()
        except (FileNotFoundError, subprocess.SubprocessError):
            return {}

    @classmethod
    def _find_install_for_libs_linux(cls, missing):
        """Return a suggestion on how to install missing libraries on Linux."""
        install_libs = cls._list_install_for_libs_linux(missing)
        if not install_libs:
            return ''
        return 'sudo apt install ' + ' '.join(install_libs)

    def _find_missing_deps_linux(self):
        """Return missing dependencies of `lib_names` on a Linux machine."""
        lib_names = _QtUnixPlatformsGetter.get()
        missing = {}
        for lib_name in lib_names:
            for so_file in self._root.rglob(f'*{lib_name}.so*'):
                cmd = ['ldd', str(so_file)]
                res = subprocess.run(cmd, capture_output=True, check=True)
                missing_this = [
                    line.split('=>')[0].strip()
                    for line in res.stdout.decode().splitlines()
                    # ldd always writes "not found" for missing stuff
                    # pylint: disable-next=magic-value-comparison
                    if 'not found' in line
                    ]
                if not missing_this:
                    continue
                this_so = str(so_file.relative_to(self._root).as_posix())
                missing_this = {Path(m.split()[0]).stem.replace('.so', '')
                                for m in missing_this}
                missing[this_so] = sorted(missing_this)
        return missing

    @staticmethod
    def _list_install_for_libs_linux(missing):
        """Return library names that make `missing` ones available on Linux."""
        # Build a regular-expression pattern for apt-cache
        rgx = '|'.join(missing)
        try:
            res = subprocess.run(['apt-cache', 'search', rgx],
                                 capture_output=True,
                                 timeout=2,  # seconds
                                 check=True)
        except (FileNotFoundError, subprocess.SubprocessError):
            return []
        libs = [line.split(' - ')[0]
                for line in res.stdout.decode().splitlines()]
        install_libs = []
        for missing_lib in missing:
            available_libs = [lib for lib in libs
                              if lib.startswith(missing_lib)]
            if not available_libs:
                return []  # Only return install_libs if all are found
            available_libs.sort(reverse=True)  # More recent first
            try:           # Prefer a non '-dev' version if it exists
                install_libs.append(next(
                    lib for lib in available_libs
                    # pylint: disable-next=magic-value-comparison
                    if '-dev' not in lib
                    ))
            except StopIteration:
                install_libs.append(next(iter(available_libs)))
        return install_libs


# This class is really just a container
# pylint: disable-next=too-few-public-methods
class _QtUnixPlatformsGetter:
    """A helper for finding the names of Qt platform plugins on UNIX."""

    @classmethod
    def get(cls):
        """Return a set of platform names for the current Qt version."""
        version = qtc.QT_VERSION_STR.replace('.', '_')
        try:
            getter = getattr(cls, f'_get_{version}')
        except AttributeError:
            raise NotImplementedError('Unsupported Qt version '
                                      + qtc.QT_VERSION_STR) from None
        return getter()

    @staticmethod
    def _get_5_15():
        """Return the known platforms used by Qt 5.15 by default."""
        # In order of importance used by Qt, according to the
        # implementation in qguiapplication.cpp, specifically
        # QGuiApplicationPrivate::createPlatformIntegration
        default = 'xcb'
        env_default = os.environ.get('QT_QPA_DEFAULT_PLATFORM', default)
        env_session = os.environ.get('XDG_SESSION_TYPE')
        env_desktop = (os.environ.get('XDG_CURRENT_DESKTOP')
                       or os.environ.get('XDG_SESSION_DESKTOP'))
        env_platform = os.environ.get('QT_QPA_PLATFORM')

        if env_platform:     # QT_QPA_PLATFORM always takes precedence
            return {env_platform}
        if not env_session:  # No XDG_SESSION_TYPE defined. Use default
            return {env_default}
        # Since Qt does direct string comparisons, it doesn't make
        # much sense to do anything more fancy here. Hence, disable
        # the magic-value-comparisons (R2004) below
        if env_session == 'x11':      # pylint: disable=R2004
            return {default, env_default}   # Always XCB on X11
        if env_session == 'wayland':  # pylint: disable=R2004
             # pylint: disable-next=R2004
            is_gnome = env_desktop and 'gnome' in env_desktop
            return {env_default} if is_gnome else {env_session, env_default}
        return {env_session, env_default}

    _get_5_15_17 = _get_5_15
    _get_5_15_16 = _get_5_15
    _get_5_15_15 = _get_5_15
    _get_5_15_14 = _get_5_15
    _get_5_15_13 = _get_5_15
    _get_5_15_12 = _get_5_15
    _get_5_15_11 = _get_5_15
    _get_5_15_10 = _get_5_15
    _get_5_15_9  = _get_5_15
    _get_5_15_8  = _get_5_15
    _get_5_15_7  = _get_5_15
    _get_5_15_6  = _get_5_15
    _get_5_15_5  = _get_5_15
    _get_5_15_4  = _get_5_15
    _get_5_15_3  = _get_5_15
    _get_5_15_2  = _get_5_15
    _get_5_15_1  = _get_5_15
    _get_5_15_0  = _get_5_15
