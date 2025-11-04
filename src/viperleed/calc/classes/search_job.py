"""Module search_job of viperleed.calc.classes."""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-07-17'
__license__ = 'GPLv3+'

from abc import ABC
from abc import abstractmethod
from contextlib import contextmanager
import logging
import multiprocessing as mp
import os
from pathlib import Path
import subprocess
import time

import psutil

IS_WINDOWS = os.name == 'nt'  # pylint: disable=magic-value-comparison
# Return code to mimic "command not found" in shells
COMMAND_NOT_FOUND_RETURN_CODE = 9009 if IS_WINDOWS else 127
KILLED_EXIT_CODE = 1

logger = logging.getLogger(__name__)


class SearchJob:
    """Launch and manage the MPI-based search job in a separate process.

    Takes care of directing the input to stdin (if needed), and logging
    stdout and stderr. It also gracefully handles termination of all
    MPI subprocesses in case of a request to terminate the job,
    including KeyBoardInterrupt.

    This is compatible with the "spawn" multiprocessing method.
    """

    def __init__(self, command, input_data=None, log_path=None):
        """Initialize a SearchJob instance.

        Parameters
        ----------
        command : Sequence of str
            The command to execute, for example,
            ['mpirun', '-n', '4', './binary']
        input_data : str, optional
            Input data to send to the subprocess's stdin (rf.info
            contents). No data is sent if not given, None, or an
            empty string. Default is None.
        log_path : str or Path, optional
            File to which stdout and stderr should be redirected. No
            redirection takes place if not given or None. Default is
            None.
        """
        self.command = command
        self.input_data = input_data or ''
        self.log_path = Path(log_path) if log_path else None
        self._mp_proc = None
        self._kill_me_flag = mp.Value('b', False)  # shared bool
        self._return_code = mp.Value('i', -999)    # shared int

    @property
    def returncode(self):
        """Return the exit code of the search subprocess."""
        return self._return_code.value

    def start(self):
        """Start the search subprocess inside a multiprocessing process."""
        command_str = ' '.join(self.command)
        logger.debug(f'Starting search process with command {command_str}.')
        args = (
            self.command,
            self.input_data,
            str(self.log_path.resolve()) if self.log_path else None,
            self._kill_me_flag,
            self._return_code,
            )
        self._mp_proc = mp.Process(target=_run_search_worker, args=args)
        self._mp_proc.start()
        logger.debug(f'Search process started with PID {self._mp_proc.pid}')

    def is_running(self):
        """Return True if the job is still running."""
        return self._mp_proc and self._mp_proc.is_alive()

    def wait(self):
        """Block until the job is complete."""
        if self._mp_proc:
            self._mp_proc.join()

    def terminate(self):
        """Request the job to terminate (including the called subprocess)."""
        self._kill_me_flag.value = True
        if not self._mp_proc:
            return
        self._mp_proc.join(timeout=5)
        if self._mp_proc.is_alive():
            # Hard kill if graceful didn't work
            self._mp_proc.terminate()


class SearchWorkerABC(ABC):
    """Abstract base class for search workers."""

    # pylint: disable-next=too-many-arguments  # Not worth to pack them
    def __init__(self, command, input_data, log_path, kill_flag, return_code):
        """Initialize search woker class."""
        self.command = command
        self.input_data = input_data
        self.log_path = log_path
        self.kill_flag = kill_flag
        self.return_code = return_code
        self._proc = None

    @contextmanager
    def open_log(self):
        """Return a stream for logging."""
        if self.log_path is None:
            yield subprocess.DEVNULL
        else:
            log_path = Path(self.log_path)
            log_path.parent.mkdir(parents=True, exist_ok=True)
            with log_path.open('a', encoding='utf-8') as log_file:
                yield log_file

    def run(self):
        """Run the worker subprocess and handle cleanup."""
        with self.open_log() as log_file:
            try:
                self._run(log_file)
            except FileNotFoundError:
                self.return_code.value = COMMAND_NOT_FOUND_RETURN_CODE
            except KeyboardInterrupt:
                logger.info('KeyboardInterrupt received — killing search job.')
                self.kill()
            except Exception as exc:
                logger.error(f'Error running search job: {exc}', exc_info=True)
                self.kill()
                raise

    @abstractmethod
    def start_subprocess(self, log_file):
        """Start the subprocess and return the Popen object."""

    def kill(self):
        """Forcefully terminate the subprocess and all its children."""
        if self._proc:
            self._terminate_subprocess()
        self.return_code.value = KILLED_EXIT_CODE

    def _cleanup_handles(self):
        """Close platform-specific handles."""

    def _feed_inputs(self):
        """Pipe input_data to STDIN, if necessary."""
        if not self.input_data or not self._proc:
            return
        try:
            # Communicate will close stdin on success.
            # If timeout expire, process continues.
            self._proc.communicate(input=self.input_data, timeout=0.2)
        except subprocess.TimeoutExpired:
            pass
        except (OSError, subprocess.SubprocessError):
            logger.error(
                'Error starting search. Check files SD.TL and rf.info.'
                )
            self.return_code.value = 1
            raise

    def _monitor_and_wait(self):
        """Monitor the process for termination requests."""
        while self._proc.poll() is None:  # not finished yet
            if self.kill_flag.value:
                logger.debug('Termination requested. Killing process tree.')
                self.kill()
                return
            time.sleep(1.0)
        self.return_code.value = self._proc.returncode or 0

    def _run(self, log_file):
        with self.start_subprocess(log_file) as self._proc:
            try:  # pylint: disable=too-many-try-statements
                self._feed_inputs()
                self._monitor_and_wait()
            finally:
                self._cleanup_handles()

    @abstractmethod
    def _terminate_subprocess(self):
        """Forcefully terminate the subprocess and all its children."""


def _psutil_kill_process_tree(pid):
    """Use psutil to terminate a process and all its children."""
    try:
        parent = psutil.Process(pid)
    except psutil.NoSuchProcess:
        return
    to_kill = parent.children(recursive=True) + [parent]
    # Try first gracefully with .terminate
    for proc in to_kill:
        try:
            proc.terminate()
        except psutil.NoSuchProcess:
            pass
    _, to_kill = psutil.wait_procs(to_kill, timeout=2)

    # Then kill hard those that failed gracefully
    for proc in to_kill:
        try:
            proc.kill()
        except psutil.NoSuchProcess:
            pass
    psutil.wait_procs(to_kill, timeout=2)


def _run_search_worker(*args):
    """Top-level entry point for the multiprocessing worker process."""
    worker = SearchWorker(*args)
    worker.run()


# About the disable: we could prevent it by having OS-specific
# concrete implementations in other modules and importing from
# there. Right now this seems unnecessary.
# pylint: disable-next=too-complex
if not IS_WINDOWS:  # Unix implementation

    class SearchWorker(SearchWorkerABC):
        """Unix-specific search worker implementation."""

        def start_subprocess(self, log_file):
            """Start the subprocess and return the Popen object."""
            kwargs = {
                'stdin': subprocess.PIPE,
                'stdout': log_file,
                'stderr': log_file,
                'encoding': 'utf-8',
                'start_new_session': True,
                }
            return subprocess.Popen(self.command, **kwargs)

        def _terminate_subprocess(self):
            """Forcefully terminate the subprocess and all its children."""
            _psutil_kill_process_tree(self._proc.pid)

else:  # Windows implementation

    import ctypes
    import ctypes.wintypes as wt

    # WinAPI constant, not available via ctypes. Value stable since
    # Win 7 (till at least Win 11). Necessary to give enough
    # permissions for job termination.
    PROCESS_ALL_ACCESS = 0x1F0FFF

    def valid_handle(handle, *_):
        """Raise OSError unless handle is a valid handle."""
        if not handle:
            # ctypes.WinError is a factory for OSError
            raise ctypes.WinError(ctypes.get_last_error())
        return handle

    def successful_return(success, *_):
        """Raise OSError unless success indicates a successful call."""
        if not success:
            # ctypes.WinError is a factory for OSError
            raise ctypes.WinError(ctypes.get_last_error())
        return success

    kernel32 = ctypes.WinDLL('kernel32', use_last_error=True)

    kernel32.CreateJobObjectW.argtypes = (wt.LPVOID, wt.LPCWSTR)
    kernel32.CreateJobObjectW.restype = wt.HANDLE
    kernel32.CreateJobObjectW.errcheck = valid_handle

    kernel32.AssignProcessToJobObject.argtypes = (wt.HANDLE, wt.HANDLE)
    kernel32.AssignProcessToJobObject.restype = wt.BOOL
    kernel32.AssignProcessToJobObject.errcheck = successful_return

    kernel32.TerminateJobObject.argtypes = (wt.HANDLE, wt.UINT)
    kernel32.TerminateJobObject.restype = wt.BOOL
    kernel32.TerminateJobObject.errcheck = successful_return

    kernel32.OpenProcess.argtypes = (wt.DWORD, wt.BOOL, wt.DWORD)
    kernel32.OpenProcess.restype = wt.HANDLE
    kernel32.OpenProcess.errcheck = valid_handle

    kernel32.CloseHandle.argtypes = (wt.HANDLE,)
    kernel32.CloseHandle.restype = wt.BOOL
    kernel32.CloseHandle.errcheck = successful_return


    class SearchWorker(SearchWorkerABC):
        """Windows-specific search worker implementation."""

        def __init__(self, *args, **kwargs):
            """Initialize Windows-specific worker object."""
            super().__init__(*args, **kwargs)
            self._job_handle = None   # Job Object HANDLE
            self._proc_handle = None  # process HANDLE from OpenProcess

        def start_subprocess(self, log_file):
            """Start the subprocess and return the Popen object."""
            kwargs = {
                'stdin': subprocess.PIPE,
                'stdout': log_file,
                'stderr': log_file,
                'encoding': 'utf-8',
                'creationflags': subprocess.CREATE_NEW_PROCESS_GROUP,
                }
            # About the disable: We should not use a with here, as
            # otherwise the resources may get closed too early for
            # us to be able to correctly handle the WinAPI handles.
            # The base class takes care of using a with context and
            # cleaning up all handles should any exception occur in
            # here.
            # pylint: disable-next=consider-using-with
            proc = subprocess.Popen(self.command, **kwargs)
            self._create_winapi_handles(proc)
            return proc

        def _cleanup_handles(self):
            """Close WinAPI handles."""
            if self._proc_handle:
                try:
                    kernel32.CloseHandle(self._proc_handle)
                except OSError:
                    pass
            if self._job_handle:
                try:
                    kernel32.CloseHandle(self._job_handle)
                except OSError:
                    pass
            self._proc_handle = None
            self._job_handle = None

        def _create_winapi_handles(self, proc):
            """Assign WinAPI handles to `proc`."""
            # We create a Job Object to wrap any detached processes
            # that may be created by the subprocess call. This is
            # important when using Intel oneAPI, as mpiexec (called
            # by subprocess) spawns a detached hydra_pmi_proxy.exe
            # that manages the various MPI processes. Killing mpiexec
            # with psutil will not kill hydra_pmi_proxy, leaving the
            # search de facto running in the background.
            self._job_handle = kernel32.CreateJobObjectW(None, None)
            self._proc_handle = kernel32.OpenProcess(
                PROCESS_ALL_ACCESS,  # dwDesiredAccess
                False,               # bInheritHandle
                int(proc.pid),       # dwProcessId
                )
            kernel32.AssignProcessToJobObject(self._job_handle,
                                              self._proc_handle)

        def _terminate_subprocess(self):
            """Forcefully terminate the subprocess and all its children."""
            if self._job_handle:
                try:
                    kernel32.TerminateJobObject(self._job_handle,
                                                KILLED_EXIT_CODE)
                finally:
                    self._cleanup_handles()
            # Ensure no stray processes remain via psutil
            if self._proc:
                _psutil_kill_process_tree(self._proc.pid)
