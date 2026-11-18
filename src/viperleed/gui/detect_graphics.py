"""Module detect_graphics of viperleed.gui.

Defines functionality to detect whether the current system is capable
of displaying widgets, and which Qt binding (PyQt5 or PyQt6) can be
used for this. The bindings are tried one at a time, each in its own
subprocess, as a broken Qt installation may crash the interpreter. The
candidates, and their order, are the same as for viperleed.gui.qt_compat
(see the binding module there). The first one that works is then
selected for use by viperleed.gui.qt_compat.
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
from typing import NamedTuple
from typing import Optional

from viperleed.gui.qt_compat.binding import PYQT5
from viperleed.gui.qt_compat.binding import PYQT6
from viperleed.gui.qt_compat.binding import QtApiError
from viperleed.gui.qt_compat.binding import QtBindingNotFoundError
from viperleed.gui.qt_compat.binding import QtBindingTooOldError
from viperleed.gui.qt_compat.binding import get_binding
from viperleed.gui.qt_compat.binding import get_candidates
from viperleed.gui.qt_compat.binding import import_binding
from viperleed.gui.qt_compat.binding import set_binding

# How long (in seconds) to wait for the checker subprocess to start
# and import a Qt binding. This can take several seconds, e.g., on
# the first execution on macOS or with a cold disk cache. It is not
# related to the presence of a display, hence the generous value.
_STARTUP_TIMEOUT = 60
# How long (in seconds) to wait for a QApplication to be created,
# after the binding was imported, before deciding there's no graphics.
_DISPLAY_TIMEOUT = 5
# The checker subprocess exits with _EXIT_CODE_OFFSET + PyQtSanity.
# The offset keeps the exit codes away from those used by Python
# itself (e.g., 1 for an unhandled exception, also while the
# subprocess is being started), which then count as crashes.
_EXIT_CODE_OFFSET = 100
# The modules of a Qt binding that the GUI needs. The checker imports
# them, and their compiled extensions (e.g., QtGui.abi3.so) are looked
# at for missing system libraries.
_QT_MODULES = ('QtCore', 'QtGui', 'QtWidgets')
_UNIX_RUNTIME_DIR = os.environ.get('XDG_RUNTIME_DIR', '')


class PyQtSanity(IntEnum):
    """Status of the current installation of PyQt5 or PyQt6."""

    OK = 0
    NOT_FOUND = 1      # ModuleNotFoundError
    IMPORT_ERROR = 2   # ImportError: perhaps a Qt version mismatch
    RUNTIME_CRASH = 3  # "Aborted (core dumped)", typical ABI mismatch
    NO_DISPLAY = 4     # OK, but graphics support is missing
    UNRESPONSIVE = 5   # Importing the binding did not finish in time
    INVALID_QT_API = 6  # The QT_API environment variable is invalid
    TOO_OLD = 7         # Older than in qt_compat.binding.MIN_VERSIONS

    @property
    def allows_fallback(self):
        """Return whether another binding may work after this status.

        Another binding may work if this one is missing, broken, too
        old, or crashes. A missing display, or a system so busy that
        importing does not finish in time, would affect any binding.
        """
        return self in {PyQtSanity.NOT_FOUND,
                        PyQtSanity.IMPORT_ERROR,
                        PyQtSanity.RUNTIME_CRASH,
                        PyQtSanity.TOO_OLD}

    @property
    def importable(self):
        """Return whether this status means that Qt can be imported."""
        return self in {PyQtSanity.OK, PyQtSanity.NO_DISPLAY}


class BindingCheck(NamedTuple):
    """The outcome of checking a Qt binding."""

    binding: Optional[str]   # None if no binding could be checked
    sanity: PyQtSanity
    # As QT_VERSION_STR. None if unknown (e.g., QtCore failed to load)
    qt_version: Optional[str] = None


class PyQtSanityChecker:
    """A callable that checks whether a Qt binding is usable.

    This class will check the environment only once. If you suspect
    the environment to have changed, call .forget_status() first to
    force a new explicit check. If a binding can be imported, it is
    selected for use by viperleed.gui.qt_compat.

    The candidate bindings (see qt_compat.binding.get_candidates) are
    checked one at a time, each in its own subprocess, in order. The
    check stops at the first binding that can be used, or after a
    status that would affect any binding (see allows_fallback of
    PyQtSanity).
    """

    _results = ()
    _summary = None  # The BindingCheck that the status is about

    @property
    def binding(self):
        """Return the name of the binding the status is about, or None.

        For a status that allows importing Qt, this is the binding
        that was imported (and selected for use by qt_compat). If the
        check of the last binding stopped the checks (e.g., because
        there is no display), it is that binding. Otherwise, it is the
        first binding that is installed but unusable (e.g., broken).
        It is None if no binding was found, before checking, or if the
        check was skipped via assume_ok.
        """
        summary = type(self)._summary
        return summary.binding if summary else None

    @property
    def qt_version(self):
        """Return the Qt version of binding, as QT_VERSION_STR, or None.

        It is None also if the binding was not checked, or if it could
        not load QtCore.
        """
        summary = type(self)._summary
        return summary.qt_version if summary else None

    @property
    def results(self):
        """Return a BindingCheck for each binding checked.

        The bindings are in the order they were checked. Those that
        are not installed are included, with sanity NOT_FOUND. It is
        a single BindingCheck(None, INVALID_QT_API) if the QT_API
        environment variable prevents selecting any binding, and empty
        before checking, or if the check was skipped via assume_ok.
        """
        return type(self)._results

    def __call__(self):
        """Check whether a Qt binding is usable.

        Returns
        -------
        PyQtSanity
            The status of the Qt installation.
        """
        cls = type(self)
        if cls._summary is None:
            cls._results = self._check_bindings()
            cls._summary = self._summarize(cls._results)
            if cls._summary.binding and cls._summary.sanity.importable:
                set_binding(cls._summary.binding)
        return cls._summary.sanity

    @classmethod
    def assume_ok(cls):
        """Take Qt as usable without checking the environment."""
        cls._summary = BindingCheck(None, PyQtSanity.OK)
        cls._results = ()

    @classmethod
    def forget_status(cls):
        """Forget any previous checks of the environment."""
        cls._summary = None
        cls._results = ()

    @classmethod
    def _check_binding(cls, binding):
        """Return the BindingCheck of `binding`, done in a subprocess.

        The check runs in a subprocess, not to crash the current
        interpreter. Two time limits apply. The subprocess has up
        to _STARTUP_TIMEOUT seconds to start and import `binding`.
        It then has _DISPLAY_TIMEOUT seconds to create a QApplication.
        Creating a QApplication may hang on a system that has no
        working display (e.g., if the DISPLAY environment variable
        points to an X server that does not answer).

        Parameters
        ----------
        binding : str
            The name of the binding to be checked.

        Returns
        -------
        BindingCheck
            The status of `binding`, and its Qt version if known.
        """
        imports_done, child_end = Pipe(duplex=False)
        proc = Process(target=_QtProbe(binding), args=(child_end,))
        proc.start()
        # Close our copy of the child's end, so that imports_done
        # becomes ready (end of file) also if the subprocess dies
        child_end.close()
        try:
            ready, qt_version = cls._wait_for_imports(imports_done,
                                                      proc.sentinel)
        finally:
            imports_done.close()
        if not ready:
            cls._stop(proc)
            return BindingCheck(binding, PyQtSanity.UNRESPONSIVE)
        proc.join(_DISPLAY_TIMEOUT)
        if proc.is_alive():
            # The QApplication did not crash the interpreter, but it is
            # still trying to start up. This means we have no display.
            cls._stop(proc)
            return BindingCheck(binding, PyQtSanity.NO_DISPLAY, qt_version)
        sanity = cls._sanity_from_exit_code(proc.exitcode)
        return BindingCheck(binding, sanity, qt_version)

    @classmethod
    def _check_bindings(cls):
        """Return the BindingCheck of the bindings checked.

        Returns
        -------
        results : tuple
            See the results property.
        """
        try:
            candidates = get_candidates()
        except QtApiError:
            return (BindingCheck(None, PyQtSanity.INVALID_QT_API),)
        results = []
        for binding in candidates:
            # A binding that is not installed needs no subprocess
            installed = importlib.util.find_spec(binding) is not None
            result = (cls._check_binding(binding) if installed
                      else BindingCheck(binding, PyQtSanity.NOT_FOUND))
            results.append(result)
            if not result.sanity.allows_fallback:
                break
        return tuple(results)

    @staticmethod
    def _sanity_from_exit_code(exit_code):
        """Return the PyQtSanity of the checker from its `exit_code`."""
        try:
            return PyQtSanity(exit_code - _EXIT_CODE_OFFSET)
        except ValueError:
            # Assume there was a crash. One could in principle
            # check whether the exit code is negative, or one
            # of the signals, but this keeps this checker as
            # OS-independent as possible.
            return PyQtSanity.RUNTIME_CRASH

    @staticmethod
    def _stop(proc):
        """Terminate the checker subprocess `proc`, and wait for it."""
        proc.terminate()
        # Do not leave a zombie process behind on UNIX
        proc.join(_DISPLAY_TIMEOUT)

    @staticmethod
    def _summarize(results):
        """Return the BindingCheck that summarizes all `results`.

        Parameters
        ----------
        results : Sequence of BindingCheck
            In the order they were checked. There is at least one.

        Returns
        -------
        summary : BindingCheck
            That of the last binding checked, if it stopped the
            checks. Otherwise, that of the first binding that is
            installed but unusable, or BindingCheck(None, NOT_FOUND)
            if none is installed.
        """
        *_, last = results
        if not last.sanity.allows_fallback:
            return last
        unusable = (result for result in results
                    if result.sanity is not PyQtSanity.NOT_FOUND)
        return next(unusable, BindingCheck(None, PyQtSanity.NOT_FOUND))

    @staticmethod
    def _wait_for_imports(imports_done, sentinel):
        """Wait until the subprocess is done importing, or exited.

        Parameters
        ----------
        imports_done : multiprocessing.connection.Connection
            The receiving end of the pipe, to which the checker
            subprocess sends the Qt version once done importing.
        sentinel : int
            The sentinel of the checker subprocess.

        Returns
        -------
        ready : list
            Those among `imports_done` and `sentinel` that are ready
            within _STARTUP_TIMEOUT seconds.
        qt_version : str or None
            The Qt version sent by the subprocess, if any.
        """
        ready = wait((imports_done, sentinel), _STARTUP_TIMEOUT)
        if imports_done not in ready:
            return ready, None
        try:
            return ready, imports_done.recv()
        except (EOFError, OSError):  # Exited without sending
            return ready, None


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
    return check_pyqt_sanity().importable


def find_missing_qt_dependencies():
    """Return a dict of missing PyQt dependencies."""
    sanity = check_pyqt_sanity()
    if sanity is PyQtSanity.NOT_FOUND:
        raise NotImplementedError('Cannot check Qt dependencies without PyQt')
    if sanity is PyQtSanity.INVALID_QT_API:
        raise NotImplementedError('Cannot check Qt dependencies with an '
                                  'invalid QT_API environment variable')
    if sanity is PyQtSanity.TOO_OLD:
        raise NotImplementedError('Cannot check Qt dependencies of a Qt '
                                  'binding that is too old')
    if sys.platform.startswith('win'):
        # No evidence of missing dependencies on Windows
        return {}
    try:
        # The binding is None if none could be imported
        finder = QtDependencyFinder(check_pyqt_sanity.binding)
    except QtBindingNotFoundError:  # None is installed at all
        return {}
    try:
        return finder.find_missing_dependencies()
    except NotImplementedError:  # Not supported on this platform
        return {}


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
    """Helper for finding missing system libraries needed by Qt."""

    def __init__(self, binding=None):
        """Initialize an instance.

        Parameters
        ----------
        binding : str, optional
            The name of the Qt binding (e.g., 'PyQt5') whose
            dependencies should be checked. If not given, the first
            installed candidate (see qt_compat.binding.get_candidates).
            It is not imported here.

        Raises
        ------
        QtBindingNotFoundError
            If `binding` is not given, and no candidate is installed.
        """
        self._binding = binding or self._find_installed_binding()

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

    @staticmethod
    def _find_installed_binding():
        """Return the first candidate binding that is installed.

        The binding is not imported. If checking Qt failed, this is
        the binding whose (missing) dependencies are most likely
        relevant.

        Raises
        ------
        QtBindingNotFoundError
            If no candidate binding is installed.
        """
        candidates = get_candidates()
        try:
            return next(binding for binding in candidates
                        if importlib.util.find_spec(binding) is not None)
        except StopIteration:
            raise QtBindingNotFoundError(
                f'None of the Qt bindings {candidates} is installed'
                ) from None

    def _find_missing_deps_linux(self):
        """Return missing dependencies of Qt on a Linux machine.

        Returns
        -------
        missing : dict
            Keys are the paths of the files that miss libraries:
            the compiled Qt modules of the binding, relative to its
            package folder, and the platform plugins, relative to
            their folder. Values are the sorted names of the missing
            libraries.
        """
        missing = self._find_missing_deps_of_modules_linux()
        try:
            missing.update(self._find_missing_deps_of_platforms_linux())
        except NotImplementedError:  # No platform rules for the binding
            pass
        return missing

    def _find_missing_deps_of_modules_linux(self):
        """Return missing dependencies of the binding's Qt modules.

        The modules are not imported: libraries missing for them
        make importing the binding fail (IMPORT_ERROR).
        """
        missing = {}
        for root in self._get_package_folders():
            for module in _QT_MODULES:
                for so_file in root.glob(f'{module}.*so'):
                    missing_libs = self._find_missing_libs_linux(so_file)
                    if missing_libs:
                        missing[so_file.name] = missing_libs
        return missing

    def _find_missing_deps_of_platforms_linux(self):
        """Return missing dependencies of the Qt platform plugins.

        Libraries missing for the platform plugins make Qt abort
        when starting (RUNTIME_CRASH). Nothing is imported: nothing
        is looked for if the plugins are not where the wheels of the
        binding have them (see _get_platforms_path).

        Raises
        ------
        NotImplementedError
            If it is not known which platforms the binding uses.
        """
        root = self._get_platforms_path()
        if root is None:
            return {}
        lib_names = _QtUnixPlatformsGetter.get(self._binding)
        missing = {}
        for lib_name in lib_names:
            for so_file in self._find_plugin_files_linux(root, lib_name):
                missing_libs = self._find_missing_libs_linux(so_file)
                if missing_libs:
                    this_so = so_file.relative_to(root).as_posix()
                    missing[this_so] = missing_libs
        return missing

    @staticmethod
    def _find_missing_libs_linux(so_file):
        """Return the sorted names of the libraries `so_file` misses.

        The libraries are found via ldd, which also lists those
        needed by the libraries that `so_file` needs.
        """
        cmd = ['ldd', str(so_file)]
        res = subprocess.run(cmd, capture_output=True, check=True)
        missing = {
            # E.g., 'libGL.so.1 => not found' gives 'libGL'
            Path(line.split('=>')[0].split()[0]).stem.replace('.so', '')
            for line in res.stdout.decode().splitlines()
            # ldd always writes "not found" for missing stuff
            # pylint: disable-next=magic-value-comparison
            if 'not found' in line
            }
        return sorted(missing)

    @staticmethod
    def _find_plugin_files_linux(root, platform):
        """Yield the plugin files for `platform` found in `root`."""
        # The plugin of a platform is libq<platform>.so (e.g.,
        # libqxcb.so), except for wayland before Qt 6.10, whose
        # plugin is libqwayland-generic.so. Any suffix may follow
        # '.so' (e.g., a version number).
        stems = {f'libq{platform}', f'libq{platform}-generic'}
        for so_file in root.rglob(f'libq{platform}*.so*'):
            if so_file.name.split('.so')[0] in stems:
                yield so_file

    def _get_package_folders(self):
        """Yield the folders of the binding's package, as Path objects.

        The binding is not imported. Nothing is yielded if the binding
        is not installed, or not a package.
        """
        spec = importlib.util.find_spec(self._binding)
        roots = (spec.submodule_search_locations or ()) if spec else ()
        yield from (Path(root) for root in roots)

    def _get_platforms_path(self):
        """Return the folder of the Qt platform plugins, or None.

        The folder is found without importing the binding (which may
        crash the interpreter), where the wheels of the binding have
        it: Qt5/plugins/platforms (Qt6/... for PyQt6) in the package
        folder, or Qt/plugins/platforms for PyQt5<=5.15.2 and PyQt6
        6.0. It is None elsewhere, e.g., for a binding that uses the
        Qt libraries of the system.
        """
        major = self._binding[-1]  # E.g., '5' for PyQt5
        for root in self._get_package_folders():
            for qt_folder in (f'Qt{major}', 'Qt'):
                platforms = root / qt_folder / 'plugins' / 'platforms'
                if platforms.is_dir():
                    return platforms
        return None

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


# The only public method is __call__, as this is a subprocess target
# pylint: disable-next=too-few-public-methods
class _QtProbe:
    """The checker of a Qt binding that runs in a subprocess.

    Importing a broken Qt installation may crash the interpreter, and
    creating a QApplication may hang on a system without a display.
    Hence, PyQtSanityChecker runs this in a subprocess, once for each
    binding. An instance is the target of the subprocess: it is pickled
    when the subprocess is spawned (the default on Windows and macOS),
    hence it only stores the name of the binding.
    """

    def __init__(self, binding):
        """Initialize an instance.

        Parameters
        ----------
        binding : str
            The name of the binding to be checked.
        """
        self._binding = binding

    def __call__(self, imports_done):
        """Import the binding and its modules, then exit.

        Calling this may block forever on a system that has no
        graphical capability. Otherwise, it always exits the
        interpreter, with an exit code that gives the outcome (see
        _exit_with): that of _try_importing if importing fails, that
        of _init_dummy_qapp otherwise.

        Parameters
        ----------
        imports_done : multiprocessing.connection.Connection
            The sending end of a pipe. Once done importing, also if
            importing failed, the Qt version (as QT_VERSION_STR, or
            None if QtCore could not be imported) is sent, and the
            connection closed.

        Returns
        -------
        None.
        """
        # Right now, we only exit the interpreter with little extra
        # information on what may have caused the problem. One could
        # populate a dict of info in the failure cases, storing the
        # exception message and the traceback, and write it to a
        # temporary (e.g., json) file for better context. We can
        # consider doing so in the future if it can help debugging
        # user problems.
        sanity, qt_version = self._try_importing(self._binding)
        imports_done.send(qt_version)
        imports_done.close()
        if sanity is not PyQtSanity.OK:
            self._exit_with(sanity)
        self._init_dummy_qapp(self._binding)

    @staticmethod
    def _exit_with(sanity):
        """Exit the subprocess with an exit code for `sanity`."""
        sys.exit(_EXIT_CODE_OFFSET + sanity.value)

    @classmethod
    def _init_dummy_qapp(cls, binding):
        """Initialize a dummy QApplication instance, then exit.

        Calling this may block forever on a system that has no
        graphical capability. It always exits the interpreter
        otherwise, with exit code NO_DISPLAY if the QApplication
        has no screen, OK if it has one.

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
        cls._exit_with(PyQtSanity.OK if has_screen
                       else PyQtSanity.NO_DISPLAY)

    @staticmethod
    def _try_importing(binding):
        """Return the PyQtSanity of importing modules of `binding`.

        Whether the binding is not installed, broken, or too old, is
        decided as when selecting a binding (see import_binding in
        qt_compat.binding).

        Returns
        -------
        sanity : PyQtSanity
            The outcome of importing the modules of `binding`.
        qt_version : str or None
            The Qt version of `binding`, as QT_VERSION_STR, if QtCore
            could be imported (also if other modules could not, or if
            the binding is too old).
        """
        qt_version = None
        # Both may fail alike, but other modules may fail after the
        # version is known. pylint: disable-next=too-many-try-statements
        try:
            qt_version = import_binding(binding).QT_VERSION_STR
            import_binding(binding, _QT_MODULES)
        except QtBindingNotFoundError:
            return PyQtSanity.NOT_FOUND, qt_version
        except QtBindingTooOldError as exc:
            return PyQtSanity.TOO_OLD, exc.qt_version
        except ImportError:
            return PyQtSanity.IMPORT_ERROR, qt_version
        return PyQtSanity.OK, qt_version


# This class is really just a container
# pylint: disable-next=too-few-public-methods
class _QtUnixPlatformsGetter:
    """A helper for finding the names of Qt platform plugins on UNIX.

    The Qt version is deduced from the name of the binding, which is
    not imported (importing it may crash the interpreter).
    """

    # Mapping {binding: platform-getter-name}. The wheels of PyQt5>=5.15
    # (see pyproject.toml) ship Qt 5.15. Qt 6 picks platforms alike from
    # 6.0 to (at least) 6.11.
    _binding_to_getter_name = {
        PYQT5: '_get_5_15',
        PYQT6: '_get_6',
        }

    @classmethod
    def get(cls, binding):
        """Return the platform names that the Qt of `binding` may use.

        Raises
        ------
        NotImplementedError
            If `binding` is not known.
        """
        try:
            getter_name = cls._binding_to_getter_name[binding]
        except KeyError:
            raise NotImplementedError(
                f'Unknown platforms for Qt binding {binding!r}'
                ) from None
        return getattr(cls, getter_name)()

    @classmethod
    def _get_5_15(cls):
        """Return the known platforms used by Qt 5.15 by default."""
        # In order of importance used by Qt, according to the
        # implementation in qguiapplication.cpp, specifically
        # QGuiApplicationPrivate::createPlatformIntegration
        default = 'xcb'
        env_default = os.environ.get('QT_QPA_DEFAULT_PLATFORM', default)
        env_session = os.environ.get('XDG_SESSION_TYPE')
        env_desktop = (os.environ.get('XDG_CURRENT_DESKTOP')
                       or os.environ.get('XDG_SESSION_DESKTOP'))
        env_platforms = cls._get_requested_platforms()

        if env_platforms:    # QT_QPA_PLATFORM always takes precedence
            return env_platforms
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

    @classmethod
    def _get_6(cls):
        """Return the platforms that Qt 6 may use by default."""
        # Qt 6 makes a list of platforms, and uses the first one that
        # loads (see QGuiApplicationPrivate::createPlatformIntegration
        # in qtbase/src/gui/kernel/qguiapplication.cpp). For Linux
        # builds, the list always contains the default platform, xcb.
        # It also contains wayland in a Wayland session: first, except
        # on GNOME before Qt 6.3. The order is irrelevant here.
        env_platforms = cls._get_requested_platforms()
        if env_platforms:  # QT_QPA_PLATFORM always takes precedence
            return env_platforms
        platforms = {'xcb'}
        session = os.environ.get('XDG_SESSION_TYPE')
        # pylint: disable-next=magic-value-comparison  # As in Qt
        if 'WAYLAND_DISPLAY' in os.environ or session == 'wayland':
            platforms.add('wayland')
        return platforms

    @staticmethod
    def _get_requested_platforms():
        """Return the platforms in QT_QPA_PLATFORM, without arguments.

        QT_QPA_PLATFORM may contain several platforms, separated by
        semicolons, each followed by its arguments after a colon
        (e.g., 'wayland;xcb:some-argument').
        """
        env_platform = os.environ.get('QT_QPA_PLATFORM', '')
        platforms = (p.split(':')[0] for p in env_platform.split(';'))
        return {p for p in platforms if p}
