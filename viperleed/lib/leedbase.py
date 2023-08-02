# -*- coding: utf-8 -*-
"""
Created on Jun 13 2019

@author: Florian Kraushofer
@author: Alexander M. Imre

Contains LEED- and TLEEDM-specific functions used throughout the tleedm module
"""

import copy
import logging
import multiprocessing
import os
from pathlib import Path
import re
import subprocess
import shutil
import time
from zipfile import ZipFile

import numpy as np
import psutil
from quicktions import Fraction

from viperleed.guilib import get_equivalent_beams
from viperleed.lib.base import cosvec

# The following imports are potentially the cause of ciclic
# imports. They are used exclusively as part of getTensorOriStates
# which could potentially be split off somewhere else 
from viperleed.lib.files import parameters, poscar, vibrocc


# constants for conversion Angstrom and eV <-> atomic units
HARTREE_TO_EV = 27.211396
EV_TO_HARTREE = 1/HARTREE_TO_EV
BOHR_TO_ANGSTROM = 0.529177210903
ANGSTROM_TO_BOHR = 1/BOHR_TO_ANGSTROM

logger = logging.getLogger("tleedm.leedbase")


###############################################
#                FUNCTIONS                    #
###############################################

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
            parameters.updatePARAMETERS(rp, update_from=update_from)
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
    return


def getYfunc(ivfunc, v0i):
    """
    Returns the Y function for a given function. The derivative is
    approximated by the slope between one point to the left and one to the
    right, so the Y function does not have values for the highest and lowest
    energies that are passed.
    Y = L/sqrt(1 + v0i^2 * L^2)
    L = ivfunc' / ivfunc

    Parameters
    ----------
    ivfunc : iterable
        Should be a 2D array, list of lists, or list of tuples containing
        energies in ivfunc[:,0] and values in ivfunc[:,1], both as float.

    v0i : float
        Imaginary part of the inner potential

    Returns
    -------
    yfunc : numpy array
        2D array of the Y-function energies and values

    """
    if len(ivfunc) < 3:
        # cannot calculate derivative, return empty
        return np.array([[]])
    yfunc = None
    ivfunc = np.array(ivfunc)
    for i in range(1, len(ivfunc)-1):
        vals = ivfunc[i-1:i+2, 1]
        lf = ((vals[2] - vals[0]) /
              (2*vals[1] + 1e-100))  # +1e-100 to avoid div by 0
        y = lf / (1 + (lf*v0i)**2)
        try:
            yfunc = np.append(yfunc, [[ivfunc[i, 0], y]], axis=0)
        except ValueError:
            if not yfunc:
                yfunc = np.array([[ivfunc[i, 0], y]])
            else:
                raise
    return yfunc


def _version_from_dirname(dirname):
    try:
        return float(dirname.split('v')[-1])
    except Exception:
        logger.debug("Could not parse version number "
                     f"for directory {dirname}")
        return np.nan


def getTLEEDdir(tensorleed_path, version=None):
    """Finds directories in the 'tensorleed' folder that have names starting
    with 'TensErLEED', then picks the one with the highest version number.
    Returns an absolute path to that directory, eg
    './tensorleed/TensErLEED-v1.6'."""
    if tensorleed_path is None:
        raise RuntimeError("tensorleed_path is None")
    source_dir = tensorleed_path.resolve()
    tl_version_dirs = [dir.resolve() for dir in source_dir.iterdir()
                       if ((source_dir / dir).is_dir()
                       and dir.name.startswith('TensErLEED'))]
    logger.log(1, f"getTLEEDdir: available TensErLEED directories: "
                 f"{[d.name for d in tl_version_dirs]}")
    if not tl_version_dirs:
        raise FileNotFoundError("Could not find any TensErLEED directory.")
    if version:
        logger.log(5, f"getTLEEDdir: Looking for TensErLEED version {version}")
        for tl_dir in tl_version_dirs:
            if np.isclose(version, _version_from_dirname(tl_dir.name)):
                return tl_dir
        # if we get here, we didn't find the requested version
        raise RuntimeError("Could not find requested TensErLEED version "
                           f"{version}.")
    # if no version is specified, return the highest version
    version_numbers = [_version_from_dirname(d.name) for d in tl_version_dirs]
    if all(np.isnan(version_numbers)):
        raise RuntimeError("Could not find any TensErLEED version.")
    highest_tl_version_dir = tl_version_dirs[np.nanargmax(version_numbers)]
    logger.log(1, "getTLEEDdir: highest TensErLEED version is "
               f"{highest_tl_version_dir.name}")
    return highest_tl_version_dir


def getMaxTensorIndex(home=".", zip_only=False):
    """
    Checks the Tensors folder for the highest Tensor index there,
    returns that value, or zero if there is no Tensors folder or no valid
    Tensors zip file. zip_only looks only for zip files, ignoring directories.
    """
    tensor_dir = (Path(home) / "Tensors").resolve()
    if not tensor_dir.is_dir():
        return 0
    indlist = []
    rgx = re.compile(r'Tensors_[0-9]{3}\.zip')
    for f in [f for f in os.listdir(os.path.join(home, "Tensors"))
              if (os.path.isfile(os.path.join(home, "Tensors", f))
                  and rgx.match(f))]:
        m = rgx.match(f)
        if m.span()[1] == 15:  # exact match
            indlist.append(int(m.group(0)[-7:-4]))
    if not zip_only:
        rgx = re.compile(r'Tensors_[0-9]{3}')
        for f in [f for f in os.listdir(os.path.join(home, "Tensors"))
                  if ((tensor_dir / f).is_dir() and rgx.match(f))]:
            m = rgx.match(f)
            if m.span()[1] == 11:  # exact match
                indlist.append(int(m.group(0)[-3:]))
    if indlist:
        return max(indlist)
    return 0


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
            for file in os.path.listdir(os.path.join(base_dir, "Tensors", dn)):
                shutil.copy2(file, unpack_path)
        except Exception:
            logger.error("Failed to move Tensors from {dn}")
            raise
    return None


def getDeltas(index, basedir=".", targetdir=".", required=True):
    """Fetches Delta files from Deltas or archive with specified tensor index.
    If required is set True, an error will be printed if no Delta files are
    found.
    basedir is the directory in which the Delta directory is based.
    targetdir is the directory to which the Tensor files should be moved."""
    dn = "Deltas_"+str(index).zfill(3)
    _basedir, _targetdir = Path(basedir).resolve(), Path(targetdir).resolve()
    zip_path=(_basedir / "Deltas" / dn).with_suffix(".zip")
    if os.path.isdir(_basedir / "Deltas" / dn):
        for f in [f for f in os.listdir(_basedir / "Deltas" / dn)
                  if (os.path.isfile(_basedir / "Deltas" / dn / f)
                      and f.startswith("DEL_"))]:
            try:
                shutil.copy2(_basedir / "Deltas" / dn / f, targetdir)
            except Exception:
                logger.error("Could not copy existing delta files to "
                             "work directory")
                raise
    elif os.path.isfile(zip_path):
        logger.info(f"Unpacking {dn}.zip...")
        try:
            with ZipFile(zip_path, 'r') as zip_ref:
                zip_ref.extractall(_targetdir)                                  # TODO: maybe it would be nicer to read directly from the zip file
        except Exception:
            logger.error(f"Failed to unpack {dn}.zip")
            raise
    elif required:
        logger.error("Deltas not found")
        raise RuntimeError("Deltas not found")
    return None


def getTensorOriStates(sl, path):
    """Reads POSCAR, PARAMETERS and VIBROCC from the target path, gets the
    original state of the atoms and sites, and stores them in the given
    slab's atom/site oriState variables."""
    _path = Path(path).resolve()
    for fn in ["POSCAR", "PARAMETERS", "VIBROCC"]:
        if not (_path / fn).is_file():
            logger.error("File "+fn+" is missing in "+path)
            raise RuntimeError("Could not check Tensors: File missing")
    dn = os.path.basename(path)
    try:
        tsl = poscar.readPOSCAR(_path / "POSCAR")
        trp = parameters.readPARAMETERS(filename=_path / "PARAMETERS")
        parameters.interpretPARAMETERS(trp, slab =tsl, silent=True)
        tsl.fullUpdate(trp)
        vibrocc.readVIBROCC(trp, tsl, filename=_path / "VIBROCC", silent=True)
        tsl.fullUpdate(trp)
    except Exception:
        logger.error("Error checking Tensors: Error while reading "
                     "input files in "+dn)
        logger.debug("Exception:", exc_info=True)
        raise RuntimeError("Could not check Tensors: Error loading old input "
                           "files")
    if len(tsl.atlist) != len(sl.atlist):
        logger.error("POSCAR from "+dn+" is incompatible with "
                     "current POSCAR.")
        raise RuntimeError("Tensors file incompatible")
    for at in sl.atlist:
        tal = [tat for tat in tsl.atlist if at.oriN == tat.oriN]
        if len(tal) != 1:
            logger.error("POSCAR from "+dn+" is incompatible with "
                         "current POSCAR.")
            raise RuntimeError("Tensors file incompatible")
        at.copyOriState(tal[0])
    if len(tsl.sitelist) != len(sl.sitelist):
        logger.error("Sites from "+dn+" input differ from current input.")
        raise RuntimeError("Tensors file incompatible")
    for site in sl.sitelist:
        tsitel = [s for s in tsl.sitelist if site.label == s.label]
        if len(tsitel) != 1:
            logger.error("Sites from "+dn+" input differ from current input.")
            raise RuntimeError("Tensors file incompatible")
        site.oriState = copy.deepcopy(tsitel[0])
    return None


def fortran_compile_batch(tasks, retry=True, logname="fortran-compile.log"):
    """
    Performs a list of fortran compilations.

    Parameters
    ----------
    tasks : iterable of tuples (pre, filename, post)
        Each entry will be passed to fortran_compile.
    retry : bool, optional
        If compilation fails, check whether failure is likely due to missing
        '-mcmodel=medium', and if so, retry. The default is True.
    logname: str, optional
        Name of the log file.

    Returns
    -------
    None on success, else raises RuntimeException.

    """

    r = None
    for (pre, filename, post) in tasks:
        r = fortran_compile(pre=pre, filename=filename, post=post,
                            logname=logname)
        if r:
            break
    if not r:
        return None
    if r == "mcmodel":
        if not retry:
            raise RuntimeError(
                "Compiling fortran code failed, likely due to missing "
                "'-mcmodel=medium'")
        logger.warning("Compiling fortran code failed; retrying with "
                       "'-mcmodel=medium' option.")
        newtasks = [(pre + " -mcmodel=medium", filename, post)
                    for (pre, filename, post) in tasks]
        fortran_compile_batch(newtasks, retry=False, logname=logname)
    else:
        raise RuntimeError("Compiling fortran code failed, unknown return "
                           "value: " + r)
    return None


def fortran_compile(pre="", filename="", post="",
                    logname="fortran-compile.log"):
    """Assembles pre+filename+post to a filename, tries to execute via
    subprocess.run, raises an exception if it fails."""
    fc = pre+" "+filename+" "+post
    fcl = fc.split()
    sep = ""
    if os.path.isfile(logname):
        sep = "\n\n"
    try:
        with open(logname, "a") as log:
            log.write(sep + "############\n# COMPILING: " + fc
                      + "\n############\n\n")
        with open(logname, "a") as log:
            r = subprocess.run(fcl, stdout=log, stderr=log)
    except Exception:
        logger.error("Error compiling "+filename)
        raise
    if r.returncode not in (0, 1):
            raise RuntimeError("Fortran compiler subprocess returned "
                               f"{r.returncode}.")
    if r.returncode == 1 and "-mcmodel=medium" not in fc:
        mcmodel = False
        with open(logname, "r") as log:
            if "relocation truncated to fit:" in log.read():
                mcmodel = True
                logger.warning(f"Compiling file {filename} failed, likely due to "
                               "missing '-mcmodel=medium' flag.")
        if not mcmodel:
            raise RuntimeError("Fortran compiler subprocess returned "
                               f"{r.returncode}.")
        return "mcmodel"
    if r.returncode != 0:
            raise RuntimeError("Fortran compiler subprocess returned "
                               f"{r.returncode}.")
    return None


def checkLattice(ab, eps=1e-3):
    """Takes unit vectors a,b as a 2x2 matrix, returns (lat,t), where lat is
    a string "square", "rectangular", "hexagonal", "rhombic" or "oblique", and
    t is a 2x2 transformation matrix which will transform the cell to obtuse
    if it is rhombic or hexagonal."""
    # Author: Michele Riva; slightly modified for consistency by FK
    t = np.array([[1, 0], [0, 1]])
    c = cosvec(ab[0], ab[1])
    d = (np.linalg.norm(ab[0]) / np.linalg.norm(ab[1])) - 1
    if abs(c) < eps:  # angle is 90°
        if np.abs(d) < eps:
            lat = "square"
        else:
            lat = "rectangular"
    elif np.abs(d) < eps:  # rhombic or hex
        if c > eps:  # angle is acute -> redefine to make it obtuse
            t = np.dot(np.array([[0, -1], [1, 0]]), t)  # keeps the handedness
            ab = np.dot(t, ab)
        c = cosvec(ab[0], ab[1])
        if abs(c + 1/2) < eps:
            lat = "hexagonal"
        else:
            lat = "rhombic"
    else:
        lat = "oblique"
    return (lat, t)


def reduceUnitCell(ab, eps=1e-3):
    """Takes an obtuse unit cell as a (2x2) matrix and reduces it to minimum
    circumference, keeping the area constant. This might reduce oblique unit
    cells to rectangular or hexagonal ones. Returns (ab, t, celltype), where
    ab is the modified unit cell, t is the transformation matrix, and celltype
    is a string describing the unit cell ("square", "rectangular",
    "hexagonal", "rhombic" or "oblique")."""
    # Author: Michele Riva; slightly modified for consistency by FK
    (lat, t) = checkLattice(ab, eps=eps)
    ab = np.dot(t, ab)
    if lat == "oblique":
        # Transform lattice to have the shortest two vectors, with angle
        #  closest to 90°. This might bring it to rect, hex or rhombic.
        # If neither, will anyway transform to have the closest to rect.
        # ALGORITHM for reduction to closest to rect:
        # This is a discrete version of Gram-Schmidt's algorithm to find
        #  orthogonal bases
        # At each iteration:
        # - order vectors by norm, the shortest first
        # - determine the projection of the second on the first, and calculate
        #    the nearest integer kk
        # - subtract from the second the projection calculated above
        # - check whether now the second is the smallest. If yes, repeat,
        #    otherwise finished.
        swap = np.array([[1, 0], [0, 1]])
        # matrix that keeps track of whether a and b are swapped
        while True:  # Swap vectors if needed to get the shortest first
            if np.linalg.norm(ab[0]) > np.linalg.norm(ab[1]):
                t0 = np.array([[0, 1], [1, 0]])
            else:
                t0 = np.array([[1, 0], [0, 1]])
            swap = np.dot(t0, swap)
            t = np.dot(t0, t)
            ab = np.dot(t0, ab)
            kk = int(np.round(np.dot(ab[0], ab[1]) / np.dot(ab[0], ab[0])))
            t0 = np.array([[1, 0], [-kk, 1]])
            t = np.dot(t0, t)
            ab = np.dot(t0, ab)
            if np.linalg.norm(ab[0]) <= np.linalg.norm(ab[1]):
                break
        # Swap vectors back if they were overall swapped
        t = np.dot(swap, t)
        ab = np.dot(swap, ab)
#         END OF ALGORITHM. Now the lattice ab is closest to rectangular. It
#           might be still any shape (square, rect, hex, rhombic, oblique)
        lat, t0 = checkLattice(ab, eps=eps)
        t = np.dot(t0, t)
        ab = np.dot(t0, ab)
#       If ab is still oblique, try to see if it can be transformed to hex or
#         rhombic by choosing "a" not to be the shortest vector of all.
#       If possible, keep the new transformation. Otherwise, stick to the one
#         that makes it closest to rectangular.
        if lat == "oblique":  # lattice is still oblique
            # Re-swapping guarantees that that the matrix has on the first line
            #   the shortest possible vector, and on the second line the second
            #   shortest possible vector.
            # The only possible combinations that can lead to a rhombic/hex are
            #   a'=b+a or a'=b-a, depending on whether the angle is acute or
            #   obtuse, respectively
            t2 = swap
            ab = np.dot(swap, ab)

            c = cosvec(ab[0], ab[1])
            t0 = [[-int(np.sign(c)), 1], [0, 1]]
            t2 = np.dot(t0, t2)
            ab = np.dot(t0, ab)

            lat, t0 = checkLattice(ab, eps=eps)
            t2 = np.dot(t0, t2)

            if lat == "oblique":
                # lattice is still oblique, no transformation is needed (will
                #   keep the one closest to rect)
                t2 = np.array([[1, 0], [0, 1]])
        else:
            t2 = np.array([[1, 0], [0, 1]])
        t = np.dot(t2, t)
    return ab, t, lat


def getLEEDdict(sl, rp):
    """Returns a LEED dict containing information needed by guilib functions"""
    if sl.planegroup == "unknown":
        logger.warning("Generating LEED dictionary for slab with unknown "
                       "plane group!")
    if sl.planegroup in ["pm", "pg", "cm", "rcm", "pmg"]:
        pgstring = sl.planegroup+"[{} {}]".format(sl.orisymplane.par[0],
                                                  sl.orisymplane.par[1])
    else:
        pgstring = sl.planegroup
    if not (abs(np.round(rp.SUPERLATTICE).astype(int) - rp.SUPERLATTICE)
            < 1e-3).all():
        logger.error("getLEEDdict: SUPERLATTICE contains non-integer-valued "
                     "entries.")
        return None
    d = {"eMax": rp.THEO_ENERGIES[1],
         "SUPERLATTICE": rp.SUPERLATTICE.astype(int),
         "surfBasis": sl.ucell[:2, :2].T,
         "surfGroup": pgstring, "bulkGroup": sl.bulkslab.foundplanegroup,
         "bulk3Dsym": sl.bulkslab.getBulk3Dstr(),
         "screenAperture": rp.SCREEN_APERTURE,
         "beamIncidence": (rp.THETA, rp.PHI)}
    # some values can be overwritten via parameters:
    if isinstance(rp.AVERAGE_BEAMS, tuple):
        d["beamIncidence"] = rp.AVERAGE_BEAMS
    # some definitions for bulk symmetry. # TODO: use guilib functions
    allowed_groups = {
        "oblique": ("p1", "p2"),
        "rhombic": ("p1", "p2", "cm", "cmm"),
        "rectangular": (
            "p1", "p2", "pm", "pg", "rcm", "pmm", "pmg", "pgg", "rcmm"),
        "square": (
            "p1", "p2", "pm", "pg", "cm", "cmm", "rcm", "pmm", "pmg", "pgg",
            "rcmm", "p4", "p4m", "p4g"),
        "hexagonal": (
            "p1", "p2", "cm", "cmm", "p3", "p3m1", "p31m", "p6", "p6m")
        }
    allowed_rotations = {
        "oblique": (2,),
        "rhombic": (2,),
        "rectangular": (2,),
        "square": (2, 4),
        "hexagonal": (2, 3, 4)
        }
    allowed_mirrors = {
        "oblique": (),
        "rhombic": ((1, 1), (1, -1)),
        "rectangular": ((1, 0), (0, 1)),
        "square": ((1, 0), (0, 1), (1, 1), (1, -1)),
        "hexagonal": ((1, 0), (0, 1), (1, 1), (1, -1), (1, 2), (2, 1))
        }
    if "group" in rp.SYMMETRY_BULK:
        if (rp.SYMMETRY_BULK["group"].split("[")[0]
                not in allowed_groups[sl.bulkslab.celltype]):
            logger.warning("Group {} given in SYMMETRY_BULK is not allowed "
                           "for bulk cell type '{}'.".format(
                               rp.SYMMETRY_BULK["group"].split("[")[0],
                               sl.bulkslab.celltype))
            rp.setHaltingLevel(2)
        else:
            d["bulkGroup"] = rp.SYMMETRY_BULK["group"]
        bulk_rotations = []
        bulk_mirrors = []
        if "rotation" in rp.SYMMETRY_BULK:
            for order in rp.SYMMETRY_BULK["rotation"]:
                if order not in allowed_rotations[sl.bulkslab.celltype]:
                    logger.warning(
                        "Rotation order {} given in SYMMETRY_BULK is not "
                        "allowed for bulk cell type '{}'.".format(
                            order, sl.bulkslab.celltype))
                    rp.setHaltingLevel(2)
                else:
                    bulk_rotations.append(order)
        if "mirror" in rp.SYMMETRY_BULK:
            for par in rp.SYMMETRY_BULK["mirror"]:
                if par not in allowed_mirrors[sl.bulkslab.celltype]:
                    logger.warning(
                        "Mirror direction {} given in SYMMETRY_BULK is not "
                        "allowed for bulk cell type '{}'.".format(
                            par, sl.bulkslab.celltype))
                    rp.setHaltingLevel(2)
                else:
                    bulk_mirrors.append(par)
        b3ds = ""
        if bulk_rotations:
            b3ds += "r({})".format(", ".join([str(v) for v in bulk_rotations]))
        if bulk_mirrors:
            if b3ds:
                b3ds += ", "
            b3ds += "m({})".format(", ".join([np.array2string(np.array(par),
                                                              separator=",")
                                              for par in bulk_mirrors]))
        if not b3ds:
            d["bulk3Dsym"] = "None"
        else:
            d["bulk3Dsym"] = b3ds
    return d


def getSymEqBeams(sl, rp):
    """Returns a list of tuples ((hf,kf), index), where (hf,kf) are beams and
    index is the group of other beams they are equivalent to"""
    if rp.AVERAGE_BEAMS is False:
        return []
    if not rp.domainParams:
        d = [getLEEDdict(sl, rp)]
    else:
        d = [getLEEDdict(dp.sl, dp.rp) for dp in rp.domainParams]
    if any([v is None for v in d]):
        logger.error("Failed to get beam equivalence list")
        return []
    symeqnames = get_equivalent_beams(*d)
    symeq = []
    rgx = re.compile(r'(?P<h>[-0-9/]+)\s*,\s*(?P<k>[-0-9/]+)')
    for (name, index) in symeqnames:
        m = rgx.match(name)
        if m:
            h, k = m.group("h"), m.group("k")
            try:
                hf, kf = float(Fraction(h)), float(Fraction(k))
            except ValueError:
                logger.error("getBeamCorrespondence: Could not convert beam "
                             "names from getEquivalentBeams to h,k floats")
                raise
            symeq.append(((hf, kf), index))
        else:
            logger.warning("getBeamCorrespondence: Beam name from "
                           "getEquivalentBeams not recognized: "+name)
    return symeq


def getBeamCorrespondence(sl, rp):
    """Compares theoretical and experimental beams, returns a list containing
    a number for each theoretical beam, indicating which experimental beam it
    corresponds to."""
    eps = 1e-5  # for comparing beams
    #   make dictionary {theoretical_beam: experimental_beam}
    beamcorr = {}
    # get symmetry-equivalence list:
    symeq = getSymEqBeams(sl, rp)
    # first, go through experimental beams and assign theoreticals if clear:
    remlist = []
    for eb in rp.expbeams:
        found = False
        for tb in rp.ivbeams:
            if eb.isEqual(tb, eps=eps):
                beamcorr[tb] = eb
                found = True
                break
        if not found:
            # check symmetry equivalent ones to exp beam
            eqbl = []   # hk of equivalent beams
            for (hk, i) in symeq:
                if eb.isEqual_hk(hk, eps=eps):
                    eqbl.extend([hk2 for (hk2, j) in symeq if i == j])
            for hk in eqbl:
                for tb in rp.ivbeams:
                    if tb.isEqual_hk(hk, eps=eps):
                        beamcorr[tb] = eb
                        found = True
                        break
                if found:
                    break
        if not found:
            logger.warning(
                "Experimental beam "+eb.label+" does not have a "
                "theoretical counterpart. Consider adding it to IVBEAMS.")
            rp.setHaltingLevel(1)
            remlist.append(eb)
    for b in remlist:
        rp.expbeams.remove(b)
        logger.warning("Beam "+b.label+" was removed from set of "
                       "experimental beams!")
    # NUMBER OF EXPERIMENTAL BEAMS MUST BE STATIC FROM HERE ON OUT
    # now, for theoretical beams without assignment, see if there is a
    #   symmetry-equivalent beam to assign them to
    notfound = []
    for tb in rp.ivbeams:
        if tb not in beamcorr.keys():
            found = False
            eqbl = []   # hk of equivalent beams
            for (hk, i) in symeq:
                if tb.isEqual_hk(hk, eps=eps):
                    eqbl.extend([hk2 for (hk2, j) in symeq if i == j])
            for hk in eqbl:
                for eb in rp.expbeams:
                    if eb.isEqual_hk(hk, eps=eps):
                        beamcorr[tb] = eb
                        found = True
                        break
                if found:
                    break
            if not found:
                notfound.append(tb.label)
    if notfound:
        logger.debug("No experimental beams found for calculated beams: "
                     + ", ".join(notfound))
    beamcorr_list = [-1]*len(rp.ivbeams)
    for i, tb in enumerate(rp.ivbeams):
        if tb in beamcorr:
            beamcorr_list[i] = rp.expbeams.index(beamcorr[tb])
    return beamcorr_list

# TODO: can eventually become part of compileTask class
def copy_compile_log(rp, logfile, log_name="fortran-compile"):
    """Copy compilation log file to compile_logs (will be moved to SUPP later).

    Parameters
    ----------
    rp : RunParameters
        rp object of the calculation.
    logfile : pathlike or str
        Path to the logfile that should be copied.
    log_name : str, optional
        Name to be used to identify logfile, eg. "refcalc". Default: "fortran-compile"
    """
    _logfile = Path(logfile)
    if rp.compile_logs_dir is None:
        # Compile directory not set. Cannot copy. Do not bother complaining.
        return
    target_path = (rp.compile_logs_dir / log_name).with_suffix(".log")
    try:
        shutil.copy2(_logfile, target_path)
    except OSError as err:
        logger.warning(
            f"Unable to copy compilation log file {str(_logfile)}. Info: {err}"
        )
