"""Module parallelization of viperleed.calc.lib.

Defines functionality useful for running tasks in multiple processes.

Part of this functionality used to be in lib.leedbase (originally
2019-06-13). Moved here especially to prevent cyclic imports.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Alexander M. Imre (@amimre)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2019-06-13'
__license__ = 'GPLv3+'

import logging
import multiprocessing
from pathlib import Path
import time

import psutil

from viperleed.calc.files import parameters


logger = logging.getLogger(__name__)


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
