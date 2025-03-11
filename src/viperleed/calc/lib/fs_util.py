"""Module fs_util of viperleed.calc.lib.

Defines backward-compatible file-system-related functionality. Mostly
functions available in the shutil standard-library module.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-03-11'
__license__ = 'GPLv3+'

import shutil
import sys

PY38 = 3, 8


def copytree_exists_ok(source, destination):
    """Copy the whole tree at the `source` directory to `destination`.

    This is a wrapper around the shutil.copytree function that
    maintains backwards compatibility down to python v3.5.

    Parameters
    ----------
    source : Path
        Base of the directory tree to be copied. Notice that symlinks
        in `source` will NOT be handled correctly for python < 3.8.
    destination : Path
        Path to the directory that will mirror source and its contents.
        It is created if it does not exist yet.

    Returns
    -------
    None.
    """
    if sys.version_info >= PY38:
        # dirs_exist_ok was introduced in python 3.8:
        # https://docs.python.org/3/library/shutil.html#shutil.copytree
        # pylint: disable-next=unexpected-keyword-arg
        shutil.copytree(source, destination, dirs_exist_ok=True)
        return
    # For earlier python versions, we need to do things manually. We
    # use a simplified version of the implementation of copytree from
    # shutil for py3.8. We assume that source and destination are Path
    # objects, and that we don't have anything special like symlinks.
    # The next line will not work in py<3.5 because of exist_ok.
    destination.mkdir(parents=True, exist_ok=True)
    for srcentry in source.glob('*'):
        dstentry = destination / srcentry.name
        if srcentry.is_dir():
            copytree_exists_ok(srcentry, dstentry)
        else:  # file
            shutil.copy2(srcentry, dstentry)
    shutil.copystat(source, destination)
