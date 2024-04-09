"""Handling of TensErLEED Fortran source code files.

Checks for the presence of the TensErLEED source code in a given path,
deals with versioning, and provides the path to the TensErLEED source
code.
"""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-04-09'
__license__ = 'GPLv3+'

import logging
import os
from pathlib import Path

from viperleed import VIPERLEED_TENSORLEED_ENV

# TensErLEED Versions
KNOWN_TL_VERSIONS = (
    '1.6',
    '1.61',
    '1.71',
    '1.72',
    '1.73',  # TODO: use Version when available
    '1.74',
    '1.75',
    '1.76',
    '2.0',
    )

_KNOWN_TENSORLEED_TOP_FOLDERS = {
    'viperleed-tensorleed',
    'tensorleed',  # For backwards compatibility
    }

def get_tensorleed_path(tensorleed_path=None):
    """Return the path to the TensErLEED source code.

    Parameters
    ----------
    tensorleed_path : pathlike, optional
        Path to the viperleed-tensorleed source code, by default None.
        If not given, the $VIPERLEED_TENSORLEED environment variable
        is used, if defined. As a last resort, the current directory
        is used. In all cases, the path given may be the top-level tree
        containing the tensor-LEED code, its parent, or a first-level
        subfolder.

    Returns
    -------
    Path
        Path to the TensErLEED source code, that is, the directory
        that contains, e.g., potentially multiple TensErLEED sources
        as well as the EEASiSSS source code and compiled binaries.

    Raises
    ------
    ValueError
        If neither the tensorleed_path argument nor the
        VIPERLEED_TENSORLEED environment variable are set.
    FileNotFoundError
        If neither tensorleed_path (or VIPERLEED_TENSORLEED), its
        parent, or any of its children, point to one of the known
        tensor-LEED directory structures.
    """
    if tensorleed_path is not None:
        path_ = Path(tensorleed_path)
        return _verify_tensorleed_path(path_)

    # Check environment variable $VIPERLEED_TENSORLEED
    try:
        path_ = Path(os.environ[VIPERLEED_TENSORLEED_ENV])
    except KeyError:  # Environment variable not set. Try CWD below.
        pass
    else:
        return _verify_tensorleed_path(path_)

    # Last resort: try with the current directory
    try:
        return _verify_tensorleed_path(Path.cwd())
    except ValueError:
        raise ValueError(
            'TensErLEED path not specified.\n'
            'Either pass a path to the TensErLEED source code with the '
            '--tensorleed command-line argument, or set the environment '
            f'variable {VIPERLEED_TENSORLEED_ENV}.'
            ) from None


def _verify_tensorleed_path(path_):
    """Return the path to the tensor-LEED folder, starting from `path_`.

    Parameters
    ----------
    path_ : Path
        Path where the tensor-LEED code is looked up.

    Returns
    -------
    verified_path : Path
        Path to the folder that contains tensor-LEED code. The
        following places are checked: `path_`, its containing
        folder, and all of its direct subfolders. `verified_path`
        is identified as a tensor-LEED containing folder if its
        name is one of the known tensor-LEED repository names
        or if it contains any subfolder or zip file whose name starts with
        'TensErLEED'.

    Raises
    ------
    FileNotFoundError
        If `path_` is not the path to an existing directory,
        or if finding a verified tensor-LEED path failed.
    """
    source = path_.resolve()
    if not source.is_dir():
        raise FileNotFoundError(f'Tensor-LEED directory {path_} not found.')

    # Collect various potential source-code places: the given path_,
    # its containing folder, and all of the path_'s direct subfolders
    potential_sources = (
        source,
        source.parent,
        *(d for d in source.iterdir() if d.is_dir())
        )

    # First look for an exact name match with
    # the known top-level tensor-LEED folders
    try:
        source_by_name = next(d for d in potential_sources
                              if d.name in _KNOWN_TENSORLEED_TOP_FOLDERS)
    except StopIteration:  # Not a known folder name
        pass
    else:
        if source_by_name != source:
            logger.warning(
                f'{source_by_name.name!r} directory found in '
                f'{source_by_name.parent}. Using {source_by_name} '
                f'instead of {source}.'
                )
        return source_by_name

    # Not found by name. Try to find if one contains a
    # TensErLEED directory instead. In principle here we
    # could use Rparams.get_tenserleed_directory() but we
    # do not care about versions and globbing is enough.
    try:
        with_tenserleed = (
            d for d in potential_sources
            if any(d.glob(f'{rparams.TENSERLEED_FOLDER_NAME}*'))
            )
    except StopIteration:
        raise FileNotFoundError('Could not find a known tensor-LEED '
                                f'source directory at {source}.') from None
    if with_tenserleed != source:
        logger.warning(
            f'TensErLEED code found in {with_tenserleed}. '
            f'Using the latter instead of {source}.'
            )
    return with_tenserleed
