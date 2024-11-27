"""Module iotensors of viperleed.calc.files.

Defines functionality useful for reading information from Tensor files.

Part of the functionality in this module used to be part of leedbase
(originally 2019-06-13). Moved here to reduce cyclic import issues.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-03-22'
__license__ = 'GPLv3+'

import copy
import logging
import os
from pathlib import Path
import re
import shutil
from zipfile import ZipFile

from viperleed.calc.constants import DEFAULT_TENSORS
from viperleed.calc.files import parameters, poscar, vibrocc
from viperleed.calc.lib.base import copytree_exists_ok

logger = logging.getLogger(__name__)


def getTensors(index, base_dir='', target_dir=''):
    """Fetch Tensors with a given `index` from a folder or an archive.

    Parameters
    ----------
    index : int
        The tensor index to be retrieved.
    base_dir : str or Path, optional
        The folder from which Tensors should be retrieved. It may
        be 'Tensors' or its parent folder. Default is the current
        directory.
    target_dir : str or Path, optional
        The root of the target tree in which the Tensor files should
        be placed. The files are placed into `target_dir`/'Tensors'.
        Default is the current directory.

    Raises
    ------
    RuntimeError
        When no Tensor file is found for `index`.
    OSError
        If any copying/extraction fails.
    """
    base_dir = Path(base_dir).resolve()
    # See if by chance we were given the 'Tensors' folder as base_dir
    # instead of its parent, in which case we navigate one level up.
    if (base_dir.name == DEFAULT_TENSORS
            and not (base_dir/DEFAULT_TENSORS).is_dir()):
        base_dir = base_dir.parent

    tensor_name = f'{DEFAULT_TENSORS}_{index:03d}'
    tensor_folder = base_dir / DEFAULT_TENSORS / tensor_name
    tensor_file = tensor_folder.with_suffix('.zip')
    target_dir = Path(target_dir).resolve()
    unpack_path = target_dir / DEFAULT_TENSORS / tensor_name

    if not tensor_folder.is_dir() and not tensor_file.is_file():
        logger.error(f'{DEFAULT_TENSORS} not found')
        raise RuntimeError(f'No {DEFAULT_TENSORS} folder/zip file '             # TODO: FileNotFoundError?
                           f'for index {index} in {base_dir}')

    if not tensor_folder.is_dir() and tensor_file.is_file():
        logger.info(f'Unpacking {tensor_file.name}...')
        unpack_path.mkdir(parents=True, exist_ok=True)
        try:
            with ZipFile(tensor_file, 'r') as archive:
                archive.extractall(unpack_path)                                 # TODO: maybe it would be nicer to read directly from the zip file
        except OSError:
            logger.error(f'Failed to unpack {tensor_file.name}')
            raise
        return

    assert tensor_folder.is_dir()
    if base_dir != target_dir:
        try:
            copytree_exists_ok(tensor_folder, unpack_path)
        except OSError:
            logger.error(f'Failed to move {DEFAULT_TENSORS} '
                         f'from {tensor_folder.name}')
            raise


def getTensorOriStates(sl, path):
    """Reads POSCAR, PARAMETERS and VIBROCC from the target path, gets the
    original state of the atoms and sites, and stores them in the given
    slab's atom/site oriState variables."""
    path = Path(path).resolve()
    for fn in ["POSCAR", "PARAMETERS", "VIBROCC"]:
        if not (path / fn).is_file():
            logger.error(f"No {fn} file in {path}")
            raise RuntimeError(                                                 # TODO: FileNotFoundError?
                f"Could not check {DEFAULT_TENSORS}: File missing"
                )
    dn = path.parent
    try:
        tsl = poscar.read(path / "POSCAR")
        trp = parameters.read(filename=path/"PARAMETERS")
        parameters.interpret(trp, slab=tsl, silent=True)
        tsl.full_update(trp)
        vibrocc.readVIBROCC(trp, tsl, filename=path/"VIBROCC", silent=True)
        tsl.full_update(trp)
    except Exception as exc:
        logger.error(f"Error checking {DEFAULT_TENSORS}: Error "
                     f"while reading input files in {dn}")
        logger.debug("Exception:", exc_info=True)
        raise RuntimeError(f"Could not check {DEFAULT_TENSORS}: Error "
                           "loading old input files") from exc
    if tsl.n_atoms != sl.n_atoms:
        logger.error(f"POSCAR from {dn} is incompatible with "
                     "current POSCAR.")
        raise RuntimeError(f"{DEFAULT_TENSORS} file incompatible")
    for at in sl:
        tat = tsl.atlist.get(at.num, None)
        if tat is None:
            logger.error(f"POSCAR from {dn} is incompatible with "
                         "current POSCAR.")
            raise RuntimeError(f"{DEFAULT_TENSORS} file incompatible")
        at.copyOriState(tat)
    if len(tsl.sitelist) != len(sl.sitelist):
        logger.error(f"Sites from {dn} input differ from current input.")
        raise RuntimeError(f"{DEFAULT_TENSORS} file incompatible")
    for site in sl.sitelist:
        tsitel = [s for s in tsl.sitelist if site.label == s.label]
        if len(tsitel) != 1:
            logger.error(f"Sites from {dn} input differ from current input.")
            raise RuntimeError(f"{DEFAULT_TENSORS} file incompatible")
        site.oriState = copy.deepcopy(tsitel[0])
