"""Module context of viperleed.calc.lib.

Defines useful context managers used in various bits of calc.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-02-04'
__license__ = 'GPLv3+'

from contextlib import contextmanager
import os
from pathlib import Path


@contextmanager
def execute_in_dir(path, mkdir=False):
    """Safely execute code in a specific directory.

    Parameters
    ----------
    path : str or Path
        The path in which the block of code should be executed
        before returning to the current directory.
    mkdir : bool, optional
        Whether `path` (and its parents) should be created if
        it does not exist yet. Default is False.

    Raises
    ------
    ValueError
        If `path` is not the path to a directory, or if it does
        not exists and `mkdir` was False.
    """
    home = Path.cwd()
    path = Path(path)
    try:
        os.chdir(path)
    except NotADirectoryError as exc:
        raise ValueError(f'{path} is not a directory') from exc
    except FileNotFoundError as exc:
        if not mkdir:
            raise ValueError(f'Directory {path} does not exist. Create it '
                             'beforehand, or call with mkdir=True.') from exc
        path.mkdir(parents=True)
        os.chdir(path)
    try:
        yield
    finally:
        os.chdir(home)
