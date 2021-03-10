# -*- coding: utf-8 -*-
"""
Created on Jun 13 2019

@author: Florian Kraushofer

Contains LEED- and TLEEDM-specific functions used throughout the tleedm module
"""

import logging
import numpy as np
import re
import subprocess
import os
import shutil
import copy
from fractions import Fraction

from viperleed.guilib import get_equivalent_beams
from viperleed.tleedmlib.base import parseMathSqrt, angle, cosvec
from viperleed.tleedmlib.files.parameters import (
    readPARAMETERS, interpretPARAMETERS)
from viperleed.tleedmlib.files.poscar import readPOSCAR
from viperleed.tleedmlib.files.vibrocc import readVIBROCC

logger = logging.getLogger("tleedm.leedbase")

###############################################
#                 GLOBALS                     #
###############################################
periodic_table = [
    'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na',
    'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr',
    'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',
    'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
    'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',
    'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf',
    'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po',
    'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm',
    'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs',
    'Mt', 'Ds', 'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og']

elementCovalentRadii = {
    "H": 0.31, "He": 0.28, "Li": 1.28, "Be": 0.96,
    "B": 0.84, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57, "Ne": 0.58,
    "Na": 1.66, "Mg": 1.41, "Al": 1.21, "Si": 1.11, "P": 1.07, "S": 1.05,
    "Cl": 1.02, "Ar": 1.06, "K": 2.03, "Ca": 1.76, "Sc": 1.70, "Ti": 1.60,
    "V": 1.53, "Cr": 1.39, "Mn": 1.39, "Fe": 1.32, "Co": 1.26, "Ni": 1.24,
    "Cu": 1.32, "Zn": 1.22, "Ga": 1.22, "Ge": 1.20, "As": 1.19, "Se": 1.20,
    "Br": 1.20, "Kr": 1.16, "Rb": 2.20, "Sr": 1.95, "Y": 1.90, "Zr": 1.75,
    "Nb": 1.64, "Mo": 1.54, "Tc": 1.47, "Ru": 1.46, "Rh": 1.42, "Pd": 1.39,
    "Ag": 1.45, "Cd": 1.44, "In": 1.42, "Sn": 1.39, "Sb": 1.39, "Te": 1.38,
    "I": 1.39, "Xe": 1.40, "Cs": 2.44, "Ba": 2.15, "La": 2.07, "Ce": 2.04,
    "Pr": 2.03, "Nd": 2.01, "Pm": 1.99, "Sm": 1.98, "Eu": 1.98, "Gd": 1.96,
    "Tb": 1.94, "Dy": 1.92, "Ho": 1.92, "Er": 1.89, "Tm": 1.90, "Yb": 1.87,
    "Lu": 1.87, "Hf": 1.75, "Ta": 1.70, "W": 1.62, "Re": 1.51, "Os": 1.44,
    "Ir": 1.41, "Pt": 1.36, "Au": 1.36, "Hg": 1.32, "Tl": 1.45, "Pb": 1.46,
    "Bi": 1.48, "Po": 1.40, "At": 1.50, "Rn": 1.50, "Fr": 2.60, "Ra": 2.21,
    "Ac": 2.15, "Th": 2.06, "Pa": 2.00, "U": 1.96, "Np": 1.90, "Pu": 1.87,
    "Am": 1.80, "Cm": 1.69}
# from Cordero et al., 2008 (DOI: 10.1039/B801115J)

elementAtomicMass = {
    "H": 1.00797, "He": 4.00260, "Li": 6.941, "Be": 9.01218,
    "B": 10.81, "C": 12.011, "N": 14.0067, "O": 15.9994, "F": 18.998403,
    "Ne": 20.179, "Na": 22.98977, "Mg": 24.305, "Al": 26.98154, "Si": 28.0855,
    "P": 30.97376, "S": 32.06, "Cl": 35.453, "K": 39.0983, "Ar": 39.948,
    "Ca": 40.08, "Sc": 44.9559, "Ti": 47.90, "V": 50.9415, "Cr": 51.996,
    "Mn": 54.9380, "Fe": 55.847, "Ni": 58.70, "Co": 58.9332, "Cu": 63.546,
    "Zn": 65.38, "Ga": 69.72, "Ge": 72.59, "As": 74.9216, "Se": 78.96,
    "Br": 79.904, "Kr": 83.80, "Rb": 85.4678, "Sr": 87.62, "Y": 88.9059,
    "Zr": 91.22, "Nb": 92.9064, "Mo": 95.94, "Tc": 98, "Ru": 101.07,
    "Rh": 102.9055, "Pd": 106.4, "Ag": 107.868, "Cd": 112.41, "In": 114.82,
    "Sn": 118.69, "Sb": 121.75, "I": 126.9045, "Te": 127.60, "Xe": 131.30,
    "Cs": 132.9054, "Ba": 137.33, "La": 138.9055, "Ce": 140.12,
    "Pr": 140.9077, "Nd": 144.24, "Pm": 145, "Sm": 150.4, "Eu": 151.96,
    "Gd": 157.25, "Tb": 158.9254, "Dy": 162.50, "Ho": 164.9304, "Er": 167.26,
    "Tm": 168.9342, "Yb": 173.04, "Lu": 174.967, "Hf": 178.49, "Ta": 180.9479,
    "W": 183.85, "Re": 186.207, "Os": 190.2, "Ir": 192.22, "Pt": 195.09,
    "Au": 196.9665, "Hg": 200.59, "Tl": 204.37, "Pb": 207.2, "Bi": 208.9804,
    "Po": 209, "At": 210, "Rn": 222, "Fr": 223, "Ra": 226.0254, "Ac": 227.0278,
    "Pa": 231.0359, "Th": 232.0381, "Np": 237.0482, "U": 238.029}


###############################################
#                FUNCTIONS                    #
###############################################

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


def getTLEEDdir(home="", version=0.):
    """Finds directories in the 'tensorleed' folder that have names starting
    with 'TensErLEED', then picks the one with the highest version number.
    Returnsa relative path to that directory, eg
    './tensorleed/TensErLEED-v1.6'."""
    sd = os.path.join(home, 'tensorleed')
    ls = [dn for dn in os.listdir(sd) if (os.path.isdir(os.path.join(sd, dn))
                                          and dn.startswith('TensErLEED'))]
    highest = 0.0
    founddir = ''
    for dn in ls:
        try:
            f = float(dn.split('v')[-1])
            if f == version:
                founddir = dn
                break
            if f > highest:
                highest = f
                founddir = dn
        except Exception:
            pass
    if founddir != '':
        if version != 0 and f != version:
            logger.error("getTLEEDdir: Could not find requested "
                         "TensErLEED version {}".format(version))
            return ''
        return os.path.join(sd, dn)
    else:
        return ''


def getMaxTensorIndex(home="."):
    """Checks the Tensors folder for the highest Tensor index there,
    returns that value, or zero if there is no Tensors folder or no valid
    Tensors zip file."""
    if not os.path.isdir(os.path.join(home, "Tensors")):
        return 0
    indlist = []
    rgx = re.compile(r'Tensors_[0-9]{3}\.zip')
    for f in [f for f in os.listdir(os.path.join(home, "Tensors"))
              if (os.path.isfile(os.path.join(home, "Tensors", f))
                  and rgx.match(f))]:
        m = rgx.match(f)
        if m.span()[1] == 15:  # exact match
            indlist.append(int(m.group(0)[-7:-4]))
    rgx = re.compile(r'Tensors_[0-9]{3}')
    for f in [f for f in os.listdir(os.path.join(home, "Tensors"))
              if (os.path.isdir(os.path.join(home, "Tensors", f))
                  and rgx.match(f))]:
        m = rgx.match(f)
        if m.span()[1] == 11:  # exact match
            indlist.append(int(m.group(0)[-3:]))
    if indlist:
        return max(indlist)
    return 0


def getTensors(index, basedir=".", targetdir=".", required=True):
    """Fetches Tensor files from Tensors or archive with specified tensor
    index. If required is set True, an error will be printed if no Tensor
    files are found.
    basedir is the directory in which the Tensor directory is based.
    targetdir is the directory to which the Tensor files should be moved."""
    dn = "Tensors_"+str(index).zfill(3)
    if (os.path.basename(basedir) == "Tensors"
            and not os.path.isdir(os.path.join(basedir, "Tensors"))):
        basedir = os.path.dirname(basedir)
    if not os.path.isdir(os.path.join(basedir, "Tensors", dn)):
        if os.path.isfile(os.path.join(basedir, "Tensors", dn+".zip")):
            try:
                logger.info("Unpacking {}.zip...".format(dn))
                os.makedirs(os.path.join(targetdir, "Tensors", dn),
                            exist_ok=True)
                shutil.unpack_archive(os.path.join(basedir, "Tensors",
                                                   dn+".zip"),
                                      os.path.join(targetdir, "Tensors", dn))
            except Exception:
                logger.error("Failed to unpack {}.zip".format(dn))
                raise
        else:
            logger.error("Tensors not found")
            raise RuntimeError("Tensors not found")
    elif basedir != targetdir:
        try:
            os.makedirs(os.path.join(targetdir, "Tensors", dn), exist_ok=True)
            for file in os.path.listdir(os.path.join(basedir, "Tensors", dn)):
                shutil.copy2(file, os.path.join(targetdir, "Tensors", dn))
        except Exception:
            logger.error("Failed to move Tensors from {}".format(dn))
            raise
    return None


def getDeltas(index, basedir=".", targetdir=".", required=True):
    """Fetches Delta files from Deltas or archive with specified tensor index.
    If required is set True, an error will be printed if no Delta files are
    found.
    basedir is the directory in which the Delta directory is based.
    targetdir is the directory to which the Tensor files should be moved."""
    dn = "Deltas_"+str(index).zfill(3)
    if os.path.isdir(os.path.join(basedir, "Deltas", dn)):
        for f in [f for f in os.listdir(os.path.join(basedir, "Deltas", dn))
                  if (os.path.isfile(os.path.join(basedir, "Deltas", dn, f))
                      and f.startswith("DEL_"))]:
            try:
                shutil.copy2(os.path.join(basedir, "Deltas", dn, f), targetdir)
            except Exception:
                logger.error("Could not copy existing delta files to "
                             "work directory")
                raise
    elif os.path.isfile(os.path.join(basedir, "Deltas", dn+".zip")):
        try:
            logger.info("Unpacking {}.zip...".format(dn))
            shutil.unpack_archive(os.path.join(basedir, "Deltas", dn+".zip"),
                                  targetdir)
        except Exception:
            logger.error("Failed to unpack {}.zip".format(dn))
            raise
    elif required:
        logger.error("Deltas not found")
        raise RuntimeError("Deltas not found")
    return None


def getTensorOriStates(sl, path):
    """Reads POSCAR, PARAMETERS and VIBROCC from the target path, gets the
    original state of the atoms and sites, and stores them in the given
    slab's atom/site oriState variables."""
    for fn in ["POSCAR", "PARAMETERS", "VIBROCC"]:
        if not os.path.isfile(os.path.join(path, fn)):
            logger.error("File "+fn+" is missing in "+path)
            raise RuntimeError("Could not check Tensors: File missing")
    dn = os.path.basename(path)
    try:
        tsl = readPOSCAR(os.path.join(path, "POSCAR"))
        trp = readPARAMETERS(filename=os.path.join(path, "PARAMETERS"))
        interpretPARAMETERS(trp, slab=tsl, silent=True)
        tsl.fullUpdate(trp)
        readVIBROCC(trp, tsl, filename=os.path.join(path, "VIBROCC"),
                    silent=True)
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


def fortranCompile(pre="", filename="", post="",
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
    if r.returncode != 0:
        raise RuntimeError("Fortran compiler subprocess returned {}"
                           .format(r.returncode))
    return None


def writeWoodsNotation(ucell):
    """Takes a unit cell (as a (2x2) matrix) and attempts to write it in Woods
    Notation. Returns empty string if no Woods notation is found."""
    # !!! VERY INCOMPLETE, should at least detect simple c(a x b) cases
    # !!! Same functionality exists in guilib; replace at some point
    if ucell[1, 0] == 0 and ucell[0, 1] == 0:
        return("(" + str(int(ucell[0, 0])) + "x" + str(int(ucell[1, 1])) + ")")
    else:
        return ""


def readWoodsNotation(s, ucell):
    """Takes a string that should contain the transformation from the bulk to
    the surface unit cell in Wood notation, as well as a bulk unit cell (from
    which only the surface vectors are read). Returns a 2x2 transformation
    matrix."""
    p = re.compile(r'\s*(?P<type>[PCpc]*)\s*\(\s*(?P<g1>.+)\s*[xX]\s*'
                   + r'(?P<g2>.+)\s*\)\s*[rR]*\s*(?P<alpha>[\d.]*)')
    # this regular expression matches if (any amount of whitespace at any
    #                                                     point is ignored):
    #   - optional: first character is p or c (or P or C) (-> type)
    #   - then there is a '('
    #   - then whatever (at least one character), interpret later (-> g1)
    #   - then 'x' or 'X'
    #   - then whatever (at least one character), interpret later (-> g2)
    #   - then ')'
    #   - then (optional) 'r' or 'R'
    #   - then (optional) an integer or float number (-> alpha)
    m = p.match(s)
    if not m:
        logging.error('Could not read woods notation input '+s)
        return None
    if not m.group('type'):
        t = 'p'
    else:
        t = m.group('type').lower()
    if not m.group('alpha'):
        alpha = 0.0
    else:
        try:
            alpha = float(m.group('alpha'))
        except ValueError:
            logger.error('Could not read Woods notation angle: '
                         + m.group('alpha')+', setting angle to zero')
            alpha = 0.0
    alpha *= np.pi/180
    g1 = parseMathSqrt(m.group('g1'))
    g2 = parseMathSqrt(m.group('g2'))
    # get surface unit cell vectors from bulk unit cell (has surface
    #  periodicity!!):
    if alpha == 0.0 and t == 'p':
        mat = np.array([[g1, 0.], [0., g2]], dtype=float)
    else:
        r = [ucell[:2, 0], ucell[:2, 1]]
        # q = np.linalg.norm(r[1])/np.linalg.norm(r[0])
        # this would be to get from bulk vectors to surface, we have to reverse
        q = 1/(np.linalg.norm(r[1])/np.linalg.norm(r[0]))
        omega = abs(angle(r[0], r[1]))
        # this is always constant in Wood notation, no need to reverse.
        if t == 'p':
            # matrices from: Klaus Hermann; Crystallography and Surface
            #                               Structure (Second Edition, Wiley)
            mat = ((1/np.sin(omega))
                   * np.array([[g1*np.sin(omega-alpha),
                                g1*(1/q)*np.sin(alpha)],
                               [-g2*q*np.sin(alpha),
                                g2*np.sin(omega+alpha)]],
                              dtype=float))
        else:
            mat = ((1/(2*np.sin(omega)))
                   * np.array([[g1*np.sin(omega-alpha)-g2*q*np.sin(alpha),
                                g1*(1/q)*np.sin(alpha)+g2*np.sin(omega+alpha)],
                               [-g1*np.sin(omega-alpha)-g2*q*np.sin(alpha),
                                -g1*(1/q)*np.sin(alpha)+g2*np.sin(omega
                                                                  + alpha)]],
                              dtype=float))
    warn = False
    for i in range(0, 2):
        for j in range(0, 2):
            if abs(mat[i, j] - round(mat[i, j])) < 1e-4:
                mat[i, j] = round(mat[i, j])
            else:
                warn = True
    if warn:
        logger.warning("SUPERLATTICE matrix from Woods notation was "
                       "identified as:\n"+str(mat))
        logger.warning("SUPERLATTICE values do not round to "
                       "integer values. Check SUPERLATTICE parameter.")
    # check whether unit cell follows convention
    abbt = np.dot(np.linalg.inv(mat), np.transpose(ucell[:2, :2]))
    if np.linalg.norm(abbt[0]) > np.linalg.norm(abbt[1]) + 1e-4:
        # swap bulk vectors; -1 to keep handedness
        mat = np.dot(np.array([[0, 1], [-1, 0]]), mat)
    return mat


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
         "screenAperture": rp.SCREEN_APERTURE}
    return d


def getSymEqBeams(sl, rp):
    """Returns a list of tuples ((hf,kf), index), where (hf,kf) are beams and
    index is the group of other beams they are equivalent to"""
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
    #   initialize to "no corresponding beam" for all theoretical beams:
    beamcorr = [-1]*len(rp.ivbeams)
    # get symmetry-equivalence list:
    symeq = getSymEqBeams(sl, rp)
    # first, go through experimental beams and assign theoreticals if clear:
    remlist = []
    for (ne, eb) in enumerate(rp.expbeams):
        found = False
        for (nt, tb) in enumerate(rp.ivbeams):
            if eb.isEqual(tb, eps=eps):
                beamcorr[nt] = ne
                found = True
                break
        if not found:
            # check symmetry equivalent ones to exp beam
            eqbl = []   # hk of equivalent beams
            for (hk, i) in symeq:
                if eb.isEqual_hk(hk, eps=eps):
                    eqbl.extend([hk2 for (hk2, j) in symeq if i == j])
            for hk in eqbl:
                for (nt, tb) in enumerate(rp.ivbeams):
                    if tb.isEqual_hk(hk, eps=eps):
                        beamcorr[nt] = ne
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
    for (nt, tb) in enumerate(rp.ivbeams):
        if beamcorr[nt] == -1:
            found = False
            eqbl = []   # hk of equivalent beams
            for (hk, i) in symeq:
                if tb.isEqual_hk(hk, eps=eps):
                    eqbl.extend([hk2 for (hk2, j) in symeq if i == j])
            for hk in eqbl:
                for (ne, eb) in enumerate(rp.expbeams):
                    if eb.isEqual_hk(hk, eps=eps):
                        beamcorr[nt] = ne
                        found = True
                        break
                if found:
                    break
            if not found:
                logger.debug("No experimental beam found for calculated beam "
                             + tb.label)
    return beamcorr
