"""Class SearchJob."""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-07-17'
__license__ = 'GPLv3+'

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
        self.log_path = Path(log_path) if log_path else None
        self._mp_proc = None
        self._kill_me_flag = mp.Value("b", False)  # shared bool

    @staticmethod
    def _run_worker(command, input_data, log_path, kill_flag):
        """Run the search command in a separate process."""

        try:
            log_f = open(log_path, "a") if log_path else subprocess.DEVNULL
        except Exception as e:
            return

        try:
            proc = subprocess.Popen(
                command,
                stdin=subprocess.PIPE,
                stdout=log_f,
                stderr=log_f,
                encoding="ascii",
                start_new_session=True,
            )
        except Exception as e:
            return

        try:
            ps_proc = psutil.Process(proc.pid)

            try:
                proc.communicate(input=input_data, timeout=0.2)
            except subprocess.TimeoutExpired:
                pass
            except (OSError, subprocess.SubprocessError):
                logger.error(
                    "Error starting search. Check files SD.TL and rf.info.")

            while proc.poll() is None:
                if kill_flag.value:
                    SearchJob._kill_proc_tree(ps_proc)
                    return
                time.sleep(0.5)

        except KeyboardInterrupt:
            logger.info("Killing process due to Keyboard interrupt.")
            SearchJob._kill_proc_tree(ps_proc)
            raise

        except Exception:
            SearchJob._kill_proc_tree(ps_proc)

        finally:
            print("[SearchJob] Cleaning up.")
            if log_path and log_f != subprocess.DEVNULL:
                log_f.close()

    @staticmethod
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

    def start(self):
        """Start the search subprocess inside a multiprocessing.Process."""
        logger.debug('Starting search process with command '
                     f'"{" ".join(self.command)}".')
        self._mp_proc = mp.Process(
            target=SearchJob._run_worker,
            args=(self.command,
                  self.input_data,
                  str(self.log_path) if self.log_path else None,
                  self._kill_me_flag)
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
        return self._mp_proc.exitcode if self._mp_proc else None
