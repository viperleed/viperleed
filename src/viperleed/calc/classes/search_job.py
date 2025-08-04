"""Class SearchJob."""

__authors__ = ("Alexander M. Imre (@amimre)",)
__copyright__ = "Copyright (c) 2019-2025 ViPErLEED developers"
__created__ = "2025-07-17"
__license__ = "GPLv3+"

from contextlib import contextmanager
from pathlib import Path
import logging
import multiprocessing as mp
import psutil
import subprocess
import time

logger = logging.getLogger(__name__)


class SearchJob:
    """Launch and manage a the MPI-based search job in a separate process.

    This takes care of directing the input to stdin, and logging the stdout and
    stderr. It also gracefully handles termination of all MPI-subprocesses in
    case of a request to terminate the job, including KeyBoardInterrupt.

    This is compatible with the "spawn" multiprocessing method.
    """

    def __init__(self, command, input_data, log_path=None):
        """Initialize a SearchJob instance.

        Parameters
        ----------
        command : list of str
            The command to execute, e.g. ['mpirun', '-n', '4', './binary']
        input_data : str
            Input data to send to the subprocess's stdin (rf.info contents).
        log_path : str or Path, optional
            File to which stdout and stderr should be redirected.
        """
        self.command = command
        self.input_data = input_data
        self._mp_proc = None
        self.log_path = Path(log_path) if log_path else None
        self._kill_me_flag = mp.Value("b", False)  # shared bool
        self._return_code = mp.Value("i", -999)  # shared int for return code

    def start(self):
        """Start the search subprocess inside a multiprocessing.Process."""
        logger.debug(
            'Starting search process with command ' f'"{" ".join(self.command)}".'
        )
        self._mp_proc = mp.Process(
            target=_run_search_worker,
            args=(
                self.command,
                self.input_data,
                str(self.log_path),
                self._kill_me_flag,
                self._return_code,
            ),
        )
        self._mp_proc.start()

    def is_running(self):
        """Return True if the job is still running."""
        return self._mp_proc and self._mp_proc.is_alive()

    def wait(self):
        """Block until the job is complete."""
        self._mp_proc.join()

    def terminate(self):
        """Request the job to terminate (including the called subprocess)."""
        self._kill_me_flag.value = True
        if self._mp_proc:
            self._mp_proc.join(timeout=5)
            if self._mp_proc.is_alive():
                self._mp_proc.terminate()  # hard kill if graceful didn't work

    @property
    def returncode(self):
        return self._return_code.value


@contextmanager
def _optional_log_file_path(log_path):
    """Context manager to open a log file or yield subprocess.DEVNULL."""
    if log_path:
        f = open(log_path, "a")
        try:
            yield f
        finally:
            f.close()
    else:
        yield subprocess.DEVNULL


def _run_search_worker(command, input_data, log_path, kill_flag, return_code):
    """Run the search command in a separate process."""

    with _optional_log_file_path(log_path) as log_f:
        try:
            proc = subprocess.Popen(
                command,
                stdin=subprocess.PIPE,
                stdout=log_f,
                stderr=log_f,
                encoding="ascii",
                start_new_session=True,
            )
        except Exception:
            # failed to start the process
            return_code.value = 127  # mimic "command not found"
            return

        # send input data to the process
        try:
            proc.communicate(input=input_data, timeout=0.2)
        except subprocess.TimeoutExpired:
            pass
        except (OSError, subprocess.SubprocessError):
            logger.error("Error starting search. Check files SD.TL and rf.info.")
            return_code.value = 1
            return

        # get the process info using psutil
        try:
            ps_proc = psutil.Process(proc.pid)
        except psutil.Error:  # catch all psutil errors
            # failed to get process info
            return_code.value = 1
            raise

        # Monitoring loop
        def monitor_process():
            """Monitor the process and check for kill requests."""
            while proc.poll() is None:  # not finished yet
                if kill_flag.value:
                    _kill_proc_tree(ps_proc)
                    return
                time.sleep(1.0)

        try:
            monitor_process()
        except KeyboardInterrupt:
            logger.info("Killing process due to Keyboard interrupt.")
            _kill_proc_tree(ps_proc)
            return_code.value = 1
            raise
        except Exception:
            return_code.value = 1
            _kill_proc_tree(ps_proc)

        # pass the return code back to the main process
        return_code.value = proc.returncode


def _kill_proc_tree(ps_proc):
    """Kill the given process and all of its children."""
    for child in ps_proc.children(recursive=True):
        try:
            child.kill()
        except psutil.NoSuchProcess:
            pass
    try:
        ps_proc.kill()
    except psutil.NoSuchProcess:
        pass
