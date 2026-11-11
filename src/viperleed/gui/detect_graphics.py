"""Module detect_graphics of viperleed.gui.

Defines functionality to detect whether the current system is capable
of displaying widgets, and which Qt binding (PyQt5 or PyQt6) can be
used for this. The bindings are tried in a subprocess, as a broken
Qt installation may crash the interpreter. The candidates, and their
order, are the same as for viperleed.gui.qt_compat (see the binding
module there). The first one that works is then selected for use by
viperleed.gui.qt_compat.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2025-05-22'
__license__ = 'GPLv3+'

from enum import IntEnum
import importlib
import importlib.util
from multiprocessing import Pipe
from multiprocessing import Process
from multiprocessing.connection import wait
import os
from pathlib import Path
import subprocess
import sys

from viperleed.gui.qt_compat.binding import get_binding
from viperleed.gui.qt_compat.binding import get_candidates
from viperleed.gui.qt_compat.binding import set_binding

# How long (in seconds) to wait for the checker subprocess to start
# and import a Qt binding. This can take several seconds, e.g., on the first
# execution on macOS or with a cold disk cache. It is not related to
# the presence of a display, hence the generous value.
_STARTUP_TIMEOUT = 60
# How long (in seconds) to wait for a QApplication to be created,
# after the binding was imported, before deciding there's no graphics.
_DISPLAY_TIMEOUT = 5
# The checker subprocess exits with _EXIT_CODE_OFFSET + PyQtSanity.
# The offset keeps the exit codes away from those used by Python
# itself (e.g., 1 for an unhandled exception, also while the
# subprocess is being started), which then count as crashes.
_EXIT_CODE_OFFSET = 100
_UNIX_RUNTIME_DIR = os.environ.get('XDG_RUNTIME_DIR', '')


class PyQtSanity(IntEnum):
    """Status of the current installation of PyQt5 or PyQt6."""

    OK = 0
    NOT_FOUND = 1      # ModuleNotFoundError
    IMPORT_ERROR = 2   # ImportError: perhaps a Qt version mismatch
    RUNTIME_CRASH = 3  # "Aborted (core dumped)", typical ABI mismatch
    NO_DISPLAY = 4     # OK, but graphics support is missing
    UNRESPONSIVE = 5   # Importing the binding did not finish in time


class PyQtSanityChecker:
    """A callable that checks whether a Qt binding is usable.

    This class will check the environment only once. If you suspect
    the environment to have changed, call .forget_status() first to
    force a new explicit check. If a binding can be imported, it is
    selected for use by viperleed.gui.qt_compat.
    """

    _status = None
    _binding = None

    def __call__(self):
        """Check whether a Qt binding is usable.

        Returns
        -------
        PyQtSanity
            The status of the Qt installation.

        Raises
        ------
        QtBindingNotFoundError
            If the QT_API environment variable has an invalid value.
        """
        cls = type(self)
        if cls._status is None:
            cls._status, cls._binding = self._check_pyqt_sanity()
            if cls._binding:
                set_binding(cls._binding)
        return cls._status

    @property
    def binding(self):
        """Return the name of the binding that could be imported, or None.

        This is None also before checking, or if the check was skipped
        via assume_ok.
        """
        return type(self)._binding

    @classmethod
    def assume_ok(cls):
        """Take Qt as usable without checking the environment."""
        cls._status = PyQtSanity.OK
        cls._binding = None

    @classmethod
    def forget_status(cls):
        """Forget any previous checks of the environment."""
        cls._status = None
        cls._binding = None

    @staticmethod
    def _check_pyqt_sanity():
        """Check whether a Qt binding is usable.

        The check runs in a subprocess, not to crash the current
        interpreter. Two time limits apply. The subprocess has up
        to _STARTUP_TIMEOUT seconds to start and import one of the
        candidate Qt bindings. It then has _DISPLAY_TIMEOUT seconds
        to create a QApplication. Creating a QApplication may hang
        on a system that has no working display (e.g., if the DISPLAY
        environment variable points to an X server that does not
        answer).

        Returns
        -------
        sanity : PyQtSanity
            The status of the Qt installation.
        binding : str or None
            The name of the binding that could be imported. None if
            none could, or if the subprocess did not report it.
        """
        imports_done, child_end = Pipe(duplex=False)
        proc = Process(target=_import_dummy_pyqt_and_run_app,
                       args=(child_end, get_candidates()))
        proc.start()
        # Close our copy of the child's end, so that imports_done
        # becomes ready (end of file) also if the subprocess dies
        child_end.close()
        try:
            ready, binding = _wait_for_imports(imports_done, proc.sentinel)
        finally:
            imports_done.close()
        if not ready:
            proc.terminate()
            return PyQtSanity.UNRESPONSIVE, None
        proc.join(_DISPLAY_TIMEOUT)
        if proc.is_alive():
            # The QApplication did not crash the interpreter, but it is
            # still trying to start up. This means we have no display.
            proc.terminate()
            return PyQtSanity.NO_DISPLAY, binding
        sanity = _sanity_from_exit_code(proc.exitcode)
        usable = sanity in {PyQtSanity.OK, PyQtSanity.NO_DISPLAY}
        return sanity, binding if usable else None


check_pyqt_sanity = PyQtSanityChecker()


def has_graphics():
    """Return whether the system has graphical capabilities.

    On systems that do not have graphical capabilities, this function
    may block for a few seconds (see PyQtSanityChecker).

    Returns
    -------
    bool
    """
    if not has_pyqt():
        return False
    return check_pyqt_sanity() is not PyQtSanity.NO_DISPLAY


def has_pyqt():
    """Return whether PyQt5 or PyQt6 is installed and importable."""
    return check_pyqt_sanity() in {PyQtSanity.OK, PyQtSanity.NO_DISPLAY}


def find_missing_qt_dependencies():
    """Return a dict of missing PyQt dependencies."""
    if check_pyqt_sanity() is PyQtSanity.NOT_FOUND:
        raise NotImplementedError('Cannot check Qt dependencies without PyQt')
    if sys.platform.startswith('win'):
        # No evidence of missing dependencies on Windows
        return {}
    binding = check_pyqt_sanity.binding or _find_installed_binding()
    if not binding:
        return {}
    try:
        return QtDependencyFinder(binding).find_missing_dependencies()
    except NotImplementedError:
        return {}


def _find_installed_binding():
    """Return the first candidate binding that is installed, or None.

    The binding is not imported. If the check failed, this is the
    binding whose (missing) dependencies are most likely relevant.
    """
    return next((binding for binding in get_candidates()
                 if importlib.util.find_spec(binding) is not None),
                None)


def suppress_file_permission_warnings(binding=None):
    """Silence warnings concerning file permissions by QStandardPaths.

    Parameters
    ----------
    binding : str, optional
        The name of the Qt binding to be used. If not given, the one
        selected by viperleed.gui.qt_compat. It is imported only on
        UNIX, and only if the runtime directory cannot be created.
    """
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
        qtc = importlib.import_module(f'{binding or get_binding()}.QtCore')
        standard_path = qtc.QStandardPaths
        tmp = Path(standard_path.writableLocation(
            standard_path.StandardLocation.TempLocation
            ))
        new_runtime_dir = (tmp / dirname).resolve()
        new_runtime_dir.mkdir(mode=0o0700, parents=True, exist_ok=True)
    os.environ['XDG_RUNTIME_DIR'] = str(new_runtime_dir)


class QtDependencyFinder:
    """Helper for finding missing dependencies of Qt platform plugins."""

    def __init__(self, binding):
        """Initialize an instance.

        Parameters
        ----------
        binding : str
            The name of the Qt binding (e.g., 'PyQt5') whose
            dependencies should be checked. Only the package itself
            is imported, not its Qt modules.
        """
        pyqt = importlib.import_module(binding)
        self._binding = binding
        qt_version = binding[-1]  # '5' for PyQt5, '6' for PyQt6
        self._root = (
            Path(pyqt.__file__).resolve().parent  # The PyQt folder
            / f'Qt{qt_version}/plugins/platforms'
            )
        assert self._root.exists()

    @classmethod
    def find_install_for_libs(cls, missing):
        """Return a suggestion on how to install `missing` libraries."""
        if not missing:
            return ''
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
        lib_names = _QtUnixPlatformsGetter.get(self._binding)
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
        assert missing
        # Build a regular-expression pattern for apt-cache
        rgx = '|'.join(missing)
        try:
            res = subprocess.run(['apt-cache', 'search', rgx],
                                 capture_output=True,
                                 timeout=3,  # seconds
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


def _sanity_from_exit_code(exit_code):
    """Return the PyQtSanity of the checker subprocess from `exit_code`."""
    try:
        return PyQtSanity(exit_code - _EXIT_CODE_OFFSET)
    except ValueError:
        # Assume there was a crash. One could in principle
        # check whether the exit code is negative, or one
        # of the signals, but this keeps this checker as
        # OS-independent as possible.
        return PyQtSanity.RUNTIME_CRASH


def _wait_for_imports(imports_done, sentinel):
    """Wait until a binding was imported, or the subprocess exited.

    Parameters
    ----------
    imports_done : multiprocessing.connection.Connection
        The receiving end of the pipe, to which the checker subprocess
        sends the name of the binding once it was imported.
    sentinel : int
        The sentinel of the checker subprocess.

    Returns
    -------
    ready : list
        Those among `imports_done` and `sentinel` that are ready
        within _STARTUP_TIMEOUT seconds.
    binding : str or None
        The name of the binding that was imported, if it was sent.
    """
    ready = wait((imports_done, sentinel), _STARTUP_TIMEOUT)
    if imports_done not in ready:
        return ready, None
    try:
        return ready, imports_done.recv()
    except (EOFError, OSError):  # The subprocess ended without sending
        return ready, None


def _import_dummy_pyqt_and_run_app(imports_done, candidates):
    """Attempt importing a Qt binding and its modules.

    This function is meant to be used in a child process to prevent
    crashing the current interpreter if the Qt installation is broken.
    Calling this may block forever on a system that has no graphical
    capability. It will also exit the interpreter if importing all
    of the `candidates` fails. In this case, the exit code gives
    information on the cause (see _exit_with): IMPORT_ERROR if any
    binding was found but could not be imported, NOT_FOUND otherwise.

    Parameters
    ----------
    imports_done : multiprocessing.connection.Connection
        The sending end of a pipe. The name of the binding is sent,
        and the connection closed, once all its modules were imported.
    candidates : Sequence of str
        The names of the bindings to be tried, in order.

    Returns
    -------
    None.
    """
    # Right now, we only exit the interpreter with little extra
    # information on what may have caused the problem. One could
    # populate a dict of info in the failure cases, storing the
    # exception message and the traceback, and write it to a
    # temporary (e.g., json) file for better context. We can consider
    # doing so in the future if it can help debugging user problems.
    failures = set()
    binding = None
    for candidate in candidates:
        sanity = _try_importing(candidate)
        if sanity is PyQtSanity.OK:
            binding = candidate
            break
        failures.add(sanity)
    if binding is None:
        _exit_with(PyQtSanity.IMPORT_ERROR
                   if PyQtSanity.IMPORT_ERROR in failures
                   else PyQtSanity.NOT_FOUND)
    imports_done.send(binding)
    imports_done.close()
    _init_dummy_qapp(binding)


def _try_importing(binding):
    """Return the PyQtSanity of importing the modules of `binding`."""
    try:
        importlib.import_module(binding)
    except ModuleNotFoundError:
        return PyQtSanity.NOT_FOUND
    except ImportError:
        return PyQtSanity.IMPORT_ERROR
    sub_modules = (
        f'{binding}.QtCore',
        f'{binding}.QtWidgets',  # Also loads QtGui
        )
    for module in sub_modules:
        try:
            importlib.import_module(module)
        except ModuleNotFoundError as exc:
            # The binding may be a namespace package without its
            # modules, e.g., if only PyQt5-sip and PyQt5-Qt5 are left
            # after uninstalling PyQt5. Any other missing module (e.g.,
            # PyQt5.sip) means that the installation is broken.
            not_found = exc.name == module
            return (PyQtSanity.NOT_FOUND if not_found
                    else PyQtSanity.IMPORT_ERROR)
        except ImportError:
            return PyQtSanity.IMPORT_ERROR
    return PyQtSanity.OK


def _exit_with(sanity):
    """Exit the checker subprocess with an exit code for `sanity`."""
    sys.exit(_EXIT_CODE_OFFSET + sanity.value)


def _init_dummy_qapp(binding):
    """Initialize a dummy QApplication instance.

    This function is only meant to be used for checking whether the
    current platform has graphical capabilities. Calling this may
    block forever on a system that has no graphical capability. It
    always exits the interpreter otherwise, with exit code NO_DISPLAY
    if the QApplication has no screen, OK if it has one.

    Parameters
    ----------
    binding : str
        The name of the Qt binding to be used.

    Returns
    -------
    None.
    """
    suppress_file_permission_warnings(binding)
    qtw = importlib.import_module(f'{binding}.QtWidgets')
    app = qtw.QApplication([])
    has_screen = app.primaryScreen() is not None
    sanity = PyQtSanity.OK if has_screen else PyQtSanity.NO_DISPLAY
    _exit_with(sanity)


# This class is really just a container
# pylint: disable-next=too-few-public-methods
class _QtUnixPlatformsGetter:
    """A helper for finding the names of Qt platform plugins on UNIX.

    This class imports a Qt binding dynamically and may crash the
    interpreter if check_pyqt_sanity() is PyQtSanity.RUNTIME_CRASH.
    """

    # Mapping {version-str: platform-getter-name} for known versions
    _version_to_getter_name = {
        **{f'5.15.{minor}': '_get_5_15' for minor in range(18)},
        }

    @classmethod
    def get(cls, binding):
        """Return a set of platform names for the Qt version of `binding`."""
        try:
            qtc = importlib.import_module(f'{binding}.QtCore')
        except ImportError:
            return set()
        try:
            getter_name = cls._version_to_getter_name[qtc.QT_VERSION_STR]
        except KeyError:
            version = qtc.QT_VERSION_STR.replace('.', '_')
            getter_name = f'_get_{version}'
        try:
            getter = getattr(cls, getter_name)
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
