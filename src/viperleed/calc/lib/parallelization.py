"""Module parallelization of viperleed.calc.lib.

Defines functionality useful for running tasks in multiple processes.

Part of this functionality used to be in lib.leedbase (originally
2019-06-13). Moved here especially to prevent cyclic imports.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Alexandra Mia Imre (@alexmiame)',
)
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2024-03-22'
__license__ = 'GPLv3+'

import logging
import multiprocessing
import os
from pathlib import Path
import time

import psutil

from viperleed.calc.files import parameters


logger = logging.getLogger(__name__)

# Environment variables that limit the number of threads used by
# the BLAS/LAPACK (and OpenMP) libraries of TensErLEED executables
BLAS_THREAD_VARIABLES = (
    'MKL_NUM_THREADS',
    'OMP_NUM_THREADS',
    'OPENBLAS_NUM_THREADS',
    )
# Variables the user may set to choose the number of BLAS threads.
# OMP_NUM_THREADS is shared by all libraries: if the user set it,
# no limit is set automatically. Otherwise, each library is limited
# unless the user set one of its own variables: OpenBLAS and MKL
# prefer their own variables over OMP_NUM_THREADS.
USER_BLAS_THREAD_VARIABLES = (
    'GOTO_NUM_THREADS',
    'MKL_NUM_THREADS',
    'OMP_NUM_THREADS',
    'OPENBLAS_NUM_THREADS',
    )
_SHARED_BLAS_THREAD_VARIABLE = 'OMP_NUM_THREADS'
# The variable set for each library, and those that, if set by the
# user, leave the number of threads of that library to the user
_LIBRARY_BLAS_THREAD_VARIABLES = {
    'MKL_NUM_THREADS': ('MKL_NUM_THREADS',),
    'OPENBLAS_NUM_THREADS': ('GOTO_NUM_THREADS', 'OPENBLAS_NUM_THREADS'),
    }
_DEFAULT_BLAS_THREADS = '1'


def blas_thread_env():
    """Return a copy of os.environ with a BLAS thread limit.

    The returned environment is meant to be passed as the `env` of
    TensErLEED subprocesses that run in parallel, so that a threaded
    BLAS library does not oversubscribe the CPU. If the user set
    OMP_NUM_THREADS to a value that is not blank, the environment
    is returned unchanged (manual mode). Otherwise, OMP_NUM_THREADS
    is set to 1, and so are OPENBLAS_NUM_THREADS (unless the user
    set it or GOTO_NUM_THREADS) and MKL_NUM_THREADS (unless the user
    set it). Values set by the user are never changed. os.environ
    itself is never modified.

    Returns
    -------
    env : dict
        A copy of os.environ, with the BLAS thread limit, if any.
    """
    env = os.environ.copy()
    env.update(_get_blas_thread_limits(env))
    return env


def log_blas_thread_env(calculation):
    """Log the BLAS thread limit that blas_thread_env applies.

    If the user did not set any of USER_BLAS_THREAD_VARIABLES, the
    limit set automatically is logged at debug level. Otherwise, the
    variables set by the user, and the limits that are still set
    automatically, are logged at info level. The message reflects
    the current os.environ, i.e., the same environment that
    blas_thread_env uses. os.environ is not modified.

    Parameters
    ----------
    calculation : str
        The name of the calculation that uses the limit. It is used
        in the message.

    Returns
    -------
    None.
    """
    user_set = _get_user_blas_threads(os.environ)
    limits = _get_blas_thread_limits(os.environ)
    if not user_set:
        logger.debug('BLAS thread limit for %s: %s',
                     calculation, _format_variables(limits))
        return
    automatic = (f'Setting the others automatically: '
                 f'{_format_variables(limits)}.' if limits
                 else 'Not setting a limit automatically.')
    logger.info('Number of BLAS threads for %s set by the user: %s. %s',
                calculation, _format_variables(user_set), automatic)


def _format_variables(variables):
    """Return a string with 'name=value' pairs from a dict."""
    return ', '.join(f'{name}={value}' for name, value in variables.items())


def _get_blas_thread_limits(env):
    """Return the BLAS thread variables to set automatically in `env`.

    Parameters
    ----------
    env : Mapping
        The environment to look into.

    Returns
    -------
    limits : dict
        Names and values of the variables that limit the number of
        BLAS threads and that should be set in `env`. Empty if the
        user set OMP_NUM_THREADS. Never contains a variable that the
        user set.
    """
    user_set = _get_user_blas_threads(env)
    if _SHARED_BLAS_THREAD_VARIABLE in user_set:
        return {}
    limits = dict.fromkeys(BLAS_THREAD_VARIABLES, _DEFAULT_BLAS_THREADS)
    for name, overrides in _LIBRARY_BLAS_THREAD_VARIABLES.items():
        if any(var in user_set for var in overrides):
            del limits[name]
    return limits


def _get_user_blas_threads(env):
    """Return a dict of the USER_BLAS_THREAD_VARIABLES set in `env`.

    Parameters
    ----------
    env : Mapping
        The environment to look into.

    Returns
    -------
    user_set : dict
        Names and values of those USER_BLAS_THREAD_VARIABLES that are
        present in `env` with a value that is not blank.
    """
    return {name: env[name] for name in USER_BLAS_THREAD_VARIABLES
            if env.get(name, '').strip()}


def monitoredPool(rp, poolsize, function, tasks, update_from=Path()):
    """
    The 'function' and 'tasks' arguments are passed on to a multiprocessing
    pool of size 'poolsize' with apply_async. While waiting for the pool to
    finish, the PARAMETERS file is read every second to check whether there is
    a STOP command. If so, the pool is terminated.

    Parameters
    ----------
    rp : Rparams object
        Needed for the parameter update
    poolsize : int
        passed on to multiprocessing.Pool
    function : function
        passed on to multiprocessing.Pool.apply_async
    tasks : list of arguments
        treated like the arguments of pool.map, i.e. each element is passed on
        in a seperate call of 'function' via multiprocessing.Pool.apply_async
    update_from : pathlike
        directory from which PARAMETERS should be read for updates

    Returns
    -------
    None

    """

    def kill_pool(p):
        """Kill the subprocesses, then terminate the pool."""
        for proc in p._pool:
            parent = psutil.Process(proc.pid)
            for child in parent.children(recursive=True):
                child.kill()
        p.terminate()

    def checkPoolResult(r):
        nonlocal pool
        nonlocal killed
        if r != "":
            kill_pool(pool)
            killed = True
        return r

    pool = multiprocessing.Pool(poolsize)
    results = []
    killed = False
    for task in tasks:
        r = pool.apply_async(function, (task,), callback=checkPoolResult)
        results.append(r)
    pool.close()
    try:
        while not all(r.ready() for r in results):
            if killed:
                break
            # See if user wants to STOP
            parameters.update(rp, update_from=update_from)
            if rp.STOP:
                kill_pool(pool)
                logger.info("Stopped by STOP parameter.")
                return
            time.sleep(1)
    except KeyboardInterrupt:
        logger.warning("Stopped by keyboard interrupt.")
        kill_pool(pool)
        raise
    pool.join()
    error = False
    for r in results:
        try:
            v = r.get(timeout=1)
            if v:
                logger.error(v)
                error = True
        except (TimeoutError, multiprocessing.context.TimeoutError):
            logger.error("Failed to get result from execution of {}"
                         .format(function.__name__))
            error = True
    if error:
        raise RuntimeError("Error in parallel execution of {}"
                           .format(function.__name__))
