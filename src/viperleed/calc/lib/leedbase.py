"""LEED-specific functions used throughout the viperleed calc package."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2019-06-13'
__license__ = 'GPLv3+'

import logging
import os
from pathlib import Path
import re
import shutil
import subprocess
from zipfile import ZipFile

import numpy as np
from quicktions import Fraction

from viperleed.calc.constants import DEFAULT_DELTAS
from viperleed.calc.constants import DEFAULT_TENSORS
from viperleed.calc.lib.math_utils import cosvec
from viperleed.calc.lib.math_utils import lcm
from viperleed.calc.lib.matrix import SingularMatrixError
from viperleed.calc.lib.matrix import ensure_integer_matrix
from viperleed.guilib import get_equivalent_beams


# constants for conversion Angstrom and eV <-> atomic units
HARTREE_TO_EV = 27.211396
EV_TO_HARTREE = 1/HARTREE_TO_EV
BOHR_TO_ANGSTROM = 0.529177210903
ANGSTROM_TO_BOHR = 1/BOHR_TO_ANGSTROM

logger = logging.getLogger(__name__)


###############################################
#                FUNCTIONS                    #
###############################################

def get_superlattice_repetitions(matrix):
    """Return the number of repeats in 2D to cover a superlattice.

    Parameters
    ----------
    matrix : Sequence
        Shape (2, 2). The matrix representing the transformation
        between the "base" lattice and the super-lattice. The
        relation is suprlattice_basis = matrix @ base_basis,
        with the _basis matrices such that in-plane unit vectors
        are rows, i.e., a, b = _basis. Must be non-singular, and
        with all entries (close to) integers.

    Returns
    -------
    repeats : tuple
        Two elements, each corresponds to the number of repetitions
        of the "base" cell along its two unit vectors so as to
        completely cover the superlattice area when back-folded.

    Raises
    ------
    ValueError
        If matrix is singular or an inappropriate shape.
    NonIntegerMatrixError
        If matrix has non-integer elements.
    """
    # Work on a numpy-array copy of the input. Ensure it's float
    # to avoid UFuncTypeError when doing in-place divisions below
    matrix = ensure_integer_matrix(np.copy(matrix).astype(float))
    if matrix.shape != (2, 2):
        raise ValueError(f'Unexpected shape {matrix.shape} for superlattice '
                         'transform matrix. Should be (2, 2)')
    n_repeats = abs(matrix[0, 0] * matrix[1, 1] - matrix[1, 0] * matrix[0, 1])
    n_repeats = round(n_repeats)
    if not n_repeats:
        raise SingularMatrixError('Superlattice transform matrix is singular')

    if n_repeats == 1:  # No area change
        return 1, 1

    # The trick is making the input matrix into a lower-triangular
    # matrix by Gauss elimination. Then we read off the 'number of
    # repeats along the first base vector' from the only non-zero
    # element of the first row. The other direction is just the
    # ratio of n_repeats and the one we already know.
    if not matrix[1, 1]:  # The matrix is upper-triangular. Swap.
        matrix[[0, 1]] = matrix[[1, 0]]
    if matrix[0, 1]:      # The matrix is not triangular. Gauss.
        multiple = lcm(abs(round(matrix[0, 1])),  # lcm needs integers
                       abs(round(matrix[1, 1])))
        first_row, second_row = matrix
        first_row /= first_row[1]
        first_row -= second_row/second_row[1]
        first_row *= multiple  # Ensure integer
    return abs(round(matrix[0, 0])), abs(round(n_repeats / matrix[0, 0]))


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


# TODO: move to iotensors?
def get_tensor_indices(home='', zip_only=False):
    """Yield the indices of all the Tensor files/folders in `home`/Tensors.

    Parameters
    ----------
    home : str or Path
        The base directory, in which a Tensors folder should
        be present. Only files/folders in <home/Tensors> will
        be looked up.
    zip_only : bool, optional
        Search only for (.zip) archives, skipping directories.
        Default is False.

    Yields
    ------
    indices : int
        The unique indices of Tensor files found. Notice that
        sorting of the indices is not guaranteed.
    """
    def get_index(fpath):
        """Return the tensor index from a file path."""
        exists = fpath.is_file() or (None if zip_only else fpath.is_dir())
        if not exists:
            return -1
        *_, ind = fpath.stem.split('_', maxsplit=1)
        try:
            return int(ind)
        except ValueError:
            return -1

    tensors = Path(home, DEFAULT_TENSORS).resolve()
    if not tensors.is_dir():
        return

    _base_pattern = f'{DEFAULT_TENSORS}_[0-9][0-9][0-9]*'
    patterns = (f'{_base_pattern}.zip',)
    if not zip_only:
        patterns += (_base_pattern,)

    likely_files = {f for p in patterns for f in tensors.glob(p)}
    indices = (get_index(f) for f in likely_files)
    yield from (ind for ind in indices if ind > 0)


# TODO: move to iotensors?
def getMaxTensorIndex(home='', zip_only=False):
    """Return the highest index of tensor files/folders in `home`/Tensors.

    Parameters
    ----------
    home : str or Path
        The base directory, in which a Tensors folder should
        be present. Only files/folders in <home/Tensors> will
        be looked up.
    zip_only : bool, optional
        Search only for (.zip) archives, skipping directories.
        Default is False.

    Returns
    -------
    max_index : int
        The largest among the indices found. Zero if no tensor
        file/directories are present.
    """
    try:
        return max(get_tensor_indices(home, zip_only))
    except ValueError:  # No files
        return 0


def getDeltas(index, basedir='', targetdir='', required=True):                  # TODO: some similarities with code in iotensors
    """Fetch delta files with a given `index` from a folder or an archive.

    Parameters
    ----------
    index : int
        The progressive index of the delta-amplitudes file to be
        retrieved. This is identical to the index of the Tensors
        with which the delta-amplitude calculation was performed.
    basedir : str or Path, optional
        The folder from which Deltas should be retrieved. It should
        be the path containing the 'Deltas' folder. Default is the
        current directory.
    targetdir : str or Path, optional
        The path to the directory in which the delta files
        should be placed. Default is the current directory.
    required : bool, optional
        Whether the Deltas_`index` file/folder must be present at
        `basedir`/'Deltas'. Raise RuntimeError if True and the
        file/folder is not found. Default is True.

    Raises
    ------
    RuntimeError
        When no delta file is found for `index`.
    OSError
        If any copying/extraction fails.
    """
    basedir, targetdir = Path(basedir).resolve(), Path(targetdir).resolve()
    delta_folder = basedir / DEFAULT_DELTAS / f'{DEFAULT_DELTAS}_{index:03d}'
    delta_zip = delta_folder.with_suffix('.zip')
    if delta_folder.is_dir():                                                   # TODO: why not copytree? Does it matter that we copy non-files or anything that is non DEL_*?
        for delta_file in delta_folder.glob('DEL_*'):
            if not delta_file.is_file():
                continue
            try:
                shutil.copy2(delta_file, targetdir)
            except OSError:
                logger.error('Could not copy existing delta files to '
                             f'{targetdir.name} directory')
                raise
    elif delta_zip.is_file():
        logger.info(f'Unpacking {delta_zip.name}...')
        try:
            with ZipFile(delta_zip, 'r') as archive:
                archive.extractall(targetdir)                                   # TODO: maybe it would be nicer to read directly from the zip file
        except OSError:
            logger.error(f'Failed to unpack {delta_zip.name}')
            raise
    elif required:
        logger.error(f'{DEFAULT_DELTAS} not found')
        raise RuntimeError(f'No {DEFAULT_DELTAS} folder/zip file '              # TODO: FileNotFoundError
                           f'for index {index} in {basedir}')


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
    c = cosvec(*ab)
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
        c = cosvec(*ab)
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

            c = cosvec(*ab)
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


def reduce_c_vector(c_vec, ab_cell):
    """Compute a Minkowski-reduced version of the third unit-cell vector.

    Parameters
    ----------
    c_vec : numpy.ndarray
        The unit-cell vector to be reduced. Only the first
        two components are taken into consideration. They
        are modified in place.
    ab_cell : numpy.ndarray
        The in-plane unit cell to be used for reducing c.
        Unit vectors are rows, i.e., a, b == ab_cell.

    Returns
    -------
    None.
    """
    # The problem is essentially the closest-vector problem (CVP) that
    # is common in the theory of lattices: we have to find the lattice
    # vector (of the a,b lattice) that is closest to the projection
    # of c on the plane formed by ab_cell. That's the quantity to be
    # removed from c to make it shortest. The trick is that the ab_cell
    # should be Minkowski-reduced beforehand to ensure that also the c
    # vector will be Minkowski-reduced. The theory behind this can be
    # found in doi.org/10.1007/978-3-540-24847-7_26.
    ab_cell, *_ = reduceUnitCell(ab_cell)
    c_frac_ab = c_vec[:2].dot(np.linalg.inv(ab_cell))

    # The closest projection is an integer version of c_frac_ab.
    # It's easiest to take the floor, then consider also -a, -b,
    # and -(a+b), whichever gives the shortest vector. Otherwise
    # we'd have to also consider +a, +b, a-b, b-a, etc...
    c_frac_ab -= np.floor(c_frac_ab)

    increments = (0, 0), (1, 0), (0, 1), (1, 1)
    c_vec[:2] = min(
        ((c_frac_ab - f).dot(ab_cell) for f in increments),
        key=np.linalg.norm
        )


def bulk_3d_string(screws, glides):
    """Return info about bulk screw axes and glide planes as a string.

    Parameters
    ----------
    screws : Sequence
        Items are orders of rotation for screw axes.
    glides : Sequence
        Items are 2-tuples representing fractional coordinates
        of direction vectors of the glide planes.

    Returns
    -------
    bulk_3d_str : str
        Format is 'r(2, 4), m([1,1], [ 1,-1])' if there is
        any screw axes or glide planes, otherwise 'None'.
    """
    b3ds = ''
    if screws:
        screws_str = ', '.join(str(v) for v in screws)
        b3ds += f'r({screws_str})'
    if glides:
        if b3ds:
            b3ds += ', '
        glides_str = ', '.join(np.array2string(direction, separator=',')
                               for direction in glides)
        if glides_str:
            b3ds += f'm({glides_str})'
    return b3ds or 'None'


def getLEEDdict(sl, rp):
    """Return a LEED dict containing information needed by guilib functions."""
    if sl.planegroup == 'unknown':
        logger.warning('Generating LEED dictionary for slab with unknown '
                       'plane group!')
    pgstring = sl.planegroup
    if pgstring in {'pm', 'pg', 'cm', 'rcm', 'pmg'}:
        pgstring += str(sl.orisymplane.par)
    if not (abs(np.round(rp.SUPERLATTICE).astype(int) - rp.SUPERLATTICE)
            < 1e-3).all():
        logger.error('getLEEDdict: SUPERLATTICE contains non-integer-valued '
                     'entries.')
        return None                                                             # TODO: would be better to raise
    # Some values can be overwritten via parameters:
    d = {'eMax': rp.THEO_ENERGIES.max,
         'SUPERLATTICE': rp.SUPERLATTICE.round().astype(int),
         'surfBasis': sl.ab_cell.T,
         'surfGroup': pgstring,
         'bulkGroup': sl.bulkslab.foundplanegroup,
         'bulk3Dsym': sl.bulkslab.get_bulk_3d_str(),
         'screenAperture': rp.SCREEN_APERTURE,
         'beamIncidence': (rp.THETA, rp.PHI)}
    # some values can be overwritten via parameters:
    if isinstance(rp.AVERAGE_BEAMS, tuple):
        d['beamIncidence'] = rp.AVERAGE_BEAMS
    if 'group' not in rp.SYMMETRY_BULK:
        return d

    # Some definitions for bulk symmetry.                                       # TODO: use guilib functions
    allowed_groups = {
        'oblique': ('p1', 'p2'),
        'rhombic': ('p1', 'p2', 'cm', 'cmm'),
        'rectangular': (
            'p1', 'p2', 'pm', 'pg', 'rcm', 'pmm', 'pmg', 'pgg', 'rcmm'),
        'square': (
            'p1', 'p2', 'pm', 'pg', 'cm', 'cmm', 'rcm', 'pmm', 'pmg', 'pgg',
            'rcmm', 'p4', 'p4m', 'p4g'),
        'hexagonal': (
            'p1', 'p2', 'cm', 'cmm', 'p3', 'p3m1', 'p31m', 'p6', 'p6m')
        }
    allowed_rotations = {
        'oblique': (2,),
        'rhombic': (2,),
        'rectangular': (2,),
        'square': (2, 4),
        'hexagonal': (2, 3, 4)
        }
    allowed_mirrors = {
        'oblique': (),
        'rhombic': ((1, 1), (1, -1)),
        'rectangular': ((1, 0), (0, 1)),
        'square': ((1, 0), (0, 1), (1, 1), (1, -1)),
        'hexagonal': ((1, 0), (0, 1), (1, 1), (1, -1), (1, 2), (2, 1))
        }
    if (rp.SYMMETRY_BULK['group'].split('[')[0]
            not in allowed_groups[sl.bulkslab.celltype]):
        hermann, *_ = rp.SYMMETRY_BULK['group'].split('[')
        logger.warning(
            f'Group {hermann} given in SYMMETRY_BULK is not allowed '
            f'for bulk cell shape {sl.bulkslab.celltype!r}.'
            )
        rp.setHaltingLevel(2)
    else:
        d['bulkGroup'] = rp.SYMMETRY_BULK['group']
    bulk_rotations = []
    bulk_mirrors = []
    if 'rotation' in rp.SYMMETRY_BULK:
        for order in rp.SYMMETRY_BULK['rotation']:
            if order not in allowed_rotations[sl.bulkslab.celltype]:
                logger.warning(
                    f'Rotation order {order} given in SYMMETRY_BULK is not '
                    f'allowed for bulk cell shape {sl.bulkslab.celltype!r}.'
                    )
                rp.setHaltingLevel(2)
            else:
                bulk_rotations.append(order)
    if 'mirror' in rp.SYMMETRY_BULK:
        for par in rp.SYMMETRY_BULK['mirror']:
            if par not in allowed_mirrors[sl.bulkslab.celltype]:
                logger.warning(
                    f'Mirror direction {par} given in SYMMETRY_BULK is not '
                    f'allowed for bulk cell shape {sl.bulkslab.celltype!r}.'
                    )
                rp.setHaltingLevel(2)
            else:
                bulk_mirrors.append(par)
    d['bulk3Dsym'] = bulk_3d_string(bulk_rotations, bulk_mirrors)
    return d


def getSymEqBeams(sl, rp):
    """Returns a list of tuples ((hf,kf), index), where (hf,kf) are beams and
    index is the group of other beams they are equivalent to"""
    if rp.AVERAGE_BEAMS is False:
        return []
    if not rp.domainParams:
        d = [getLEEDdict(sl, rp)]
    else:
        d = [getLEEDdict(dp.slab, dp.rpars) for dp in rp.domainParams]
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
def copy_compile_log(rp, logfile, save_as='fortran-compile'):
    """Copy compilation log file to compile_logs (will be moved to SUPP later).

    Parameters
    ----------
    rp : Rparams
        rp object of the calculation.
    logfile : pathlike or str
        Path to the log file that should be copied.
    save_as : str, optional
        Name under which `logfile` should be copied,
        e.g. 'refcalc'. Default is 'fortran-compile'.

    Returns
    -------
    None.
    """
    target_path = (rp.paths.compile_logs / save_as).with_suffix('.log')
    try:
        shutil.copy2(logfile, target_path)
    except OSError as exc:
        logger.warning('Unable to copy compilation log '
                       f'file {logfile}. Info: {exc}')
