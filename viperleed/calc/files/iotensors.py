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

from viperleed.calc.files import parameters, poscar, vibrocc

logger = logging.getLogger(__name__)


def getTensors(index, base_dir=".", target_dir=".", required=True):
    """Fetches Tensor files from Tensors or archive with specified tensor
    index. If required is set True, an error will be printed if no Tensor
    files are found.
    base_dir is the directory in which the Tensor directory is based.
    target_dir is the directory to which the Tensor files should be moved."""
    dn = "Tensors_"+str(index).zfill(3)
    tensor_dir = (Path(base_dir) / "Tensors").resolve()
    unpack_path = (Path(target_dir) / "Tensors" / dn).resolve()
    zip_path = (tensor_dir / dn).with_suffix(".zip")

    if (os.path.basename(base_dir) == "Tensors"
            and not tensor_dir.is_dir()):
        base_dir = os.path.dirname(base_dir)
    if not (tensor_dir / dn).is_dir():
        if (tensor_dir / dn).with_suffix(".zip").is_file():
            logger.info(f"Unpacking {dn}.zip...")
            os.makedirs(unpack_path, exist_ok=True)
            try:
                with ZipFile(zip_path, 'r') as zip_ref:
                    zip_ref.extractall(unpack_path)                             # TODO: maybe it would be nicer to read directly from the zip file
            except Exception:
                logger.error(f"Failed to unpack {dn}.zip")
                raise
        else:
            logger.error("Tensors not found")
            raise RuntimeError("Tensors not found")
    elif base_dir != target_dir:
        try:
            os.makedirs(unpack_path, exist_ok=True)
            for file in os.listdir(os.path.join(base_dir, "Tensors", dn)):
                shutil.copy2(file, unpack_path)
        except OSError:
            logger.error(f"Failed to move Tensors from {dn}")
            raise


def getTensorOriStates(sl, path):
    """Reads POSCAR, PARAMETERS and VIBROCC from the target path, gets the
    original state of the atoms and sites, and stores them in the given
    slab's atom/site oriState variables."""
    path = Path(path).resolve()
    for fn in ["POSCAR", "PARAMETERS", "VIBROCC"]:
        if not (path / fn).is_file():
            logger.error(f"No {fn} file in {path}")
            raise RuntimeError("Could not check Tensors: File missing")         # TODO: FileNotFoundError?
    dn = path.parent
    try:
        tsl = poscar.read(path / "POSCAR")
        trp = parameters.read(filename=path/"PARAMETERS")
        parameters.interpret(trp, slab=tsl, silent=True)
        tsl.full_update(trp)
        vibrocc.readVIBROCC(trp, tsl, filename=path/"VIBROCC", silent=True)
        tsl.full_update(trp)
    except Exception as exc:
        logger.error("Error checking Tensors: Error while reading "
                     f"input files in {dn}")
        logger.debug("Exception:", exc_info=True)
        raise RuntimeError("Could not check Tensors: Error loading old input "
                           "files") from exc
    if tsl.n_atoms != sl.n_atoms:
        logger.error(f"POSCAR from {dn} is incompatible with "
                     "current POSCAR.")
        raise RuntimeError("Tensors file incompatible")
    for at in sl:
        tat = tsl.atlist.get(at.num, None)
        if tat is None:
            logger.error(f"POSCAR from {dn} is incompatible with "
                         "current POSCAR.")
            raise RuntimeError("Tensors file incompatible")
        at.copyOriState(tat)
    if len(tsl.sitelist) != len(sl.sitelist):
        logger.error(f"Sites from {dn} input differ from current input.")
        raise RuntimeError("Tensors file incompatible")
    for site in sl.sitelist:
        tsitel = [s for s in tsl.sitelist if site.label == s.label]
        if len(tsitel) != 1:
            logger.error(f"Sites from {dn} input differ from current input.")
            raise RuntimeError("Tensors file incompatible")
        site.oriState = copy.deepcopy(tsitel[0])
