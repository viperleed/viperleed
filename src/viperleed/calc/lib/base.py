"""Contains generic functions used in the TensErLEED scripts."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2019-06-13'
__license__ = 'GPLv3+'

import itertools
import logging
import multiprocessing
import os
import re
import shutil
import subprocess
import sys

import numpy as np
import scipy.spatial as sps

logger = logging.getLogger(__name__)
COLLAPSE_EPS = 1e-8  # Default for collapsing fractional coordinates

###############################################
#                EXCEPTIONS                   #
###############################################

class NonIntegerMatrixError(ValueError):                                         # TODO: move somewhere else
    """A matrix that should have integer values does not."""


class SingularMatrixError(ValueError, ZeroDivisionError):
    """A matrix that needs inversion is singular."""


###############################################
#                 CLASSES                     #
###############################################

class BackwardsReader:                                                          # TODO: maybe better to give it __enter__ and __exit__?
    """Simple class for reading a large file in reverse without having to read
    the entire file to memory.
    From http://code.activestate.com/recipes/120686/, adapted for python 3 and
    changed somewhat to fit requirements.
    """

    def readline(self):
        while len(self.data) == 1 and ((self.blkcount * self.blksize)
                                       < self.size):
            self.blkcount = self.blkcount + 1
            line = self.data[0]
            try:
                # read from end of file
                self.f.seek(-self.blksize * self.blkcount, 2)
                self.data = ((self.f.read(self.blksize).decode(self.encoding))
                             + line).splitlines()
            except IOError:    # can't seek before the beginning of the file
                self.f.seek(0)
                self.data = (((self.f.read(self.size - (self.blksize
                                                        * (self.blkcount-1))))
                              .decode(self.encoding) + line).splitlines())

        if len(self.data) == 0:
            return ""

        line = self.data.pop()
        return line + '\n'

    def close(self):
        self.f.close()

    def __init__(self, file, blksize=4096, encoding="utf-8"):
        """initialize the internal structures"""
        # get the file size
        self.size = os.stat(file)[6]
        # how big of a block to read from the file...
        self.blksize = blksize
        # how many blocks we've read
        self.blkcount = 1
        # what encoding to use
        self.encoding = encoding
        self.f = open(file, 'rb')
        # if the file is smaller than the blocksize, read a block,
        # otherwise, read the whole thing...
        if self.size > self.blksize:
            # read from end of file
            self.f.seek(-self.blksize * self.blkcount, 2)
        self.data = ((self.f.read(self.blksize))
                     .decode(self.encoding).splitlines())
        if not self.data:  # File is empty
            return
        # strip the last item if it's empty... a byproduct of the last line
        # having a newline at the end of it
        if not self.data[-1]:
            self.data.pop()




###############################################
#                FUNCTIONS                    #
###############################################

def add_edges_and_corners(cartesian, fractional, releps, ucell, props=None):
    """Add atoms at edges and corners to the Cartesian coordinates given.

    Notice that edges and corners will be added only for dimensions
    up to the smallest among: shape(fractional)[1], len(releps),
    len(ucell). The only requirements are that `fractional` contains
    at most the same number of elements as `cartesian`, and that the
    second axis of `cartesian` and `ucell` are broadcastable. This
    allows, e.g., to pass in 3D coordinates but only add repeats
    for 2D edges and corners.

    Parameters
    ----------
    cartesian : Sequence
        The Cartesian coordinates to which replicas of atoms at
        edges and corners will be appended. Shape (N, 2) or (N, 3).
    fractional : Sequence
        Fractional coordinates corresponding to `cartesian`. Should
        have shape (M, m), with M <= N.
    releps : Sequence
        Tolerance on fractional coordinates to consider an atom from
        `cartesian` to be close to an edge. Shape (2,) or (3,).
    ucell : Sequence
        Unit cell to be used for the displacement of atoms close to
        edges/corners. Shape (2, 2) or (3, 3), and such that unit
        cell vectors are rows. The second dimension should match
        the second one of `cartesian`.
    props : list or None, optional
        If given, it should have the same length as `cartesian`. It
        is interpreted as a list of properties that should also be
        duplicated whenever an atom is duplicated. Default is None.

    Returns
    -------
    cartesian_extended : numpy.ndarray
        Cartesian atomic coordinates to which extra atoms have been
        appended at the end for all those atoms in `cartesian` that
        were at edges/corners.
    props_extended : list or None
        Corresponding properties of the added atoms. None if `props`
        was not given to begin with.
    """
    for i, pos in enumerate(fractional):
        addlist = []
        for posj, epsj, uvec in zip(pos, releps, ucell):
            if abs(posj) < epsj:        # 'left' edge
                addlist.append(cartesian[i] + uvec)
            elif abs(posj - 1) < epsj:  # 'right' edge
                addlist.append(cartesian[i] - uvec)
        if len(addlist) == 2:
            # 2D corner - add the diagonally opposed one
            addlist.append(sum(addlist) - cartesian[i])
        elif len(addlist) == 3:
            # 3D corner - add all diagonally opposed points
            addlist.extend(
                p1 + p2 - cartesian[i]
                for p1, p2 in itertools.combinations(addlist, 2)
                )
            addlist.append(sum(addlist[:3]) - 2 * cartesian[i])

        if not addlist:
            continue
        cartesian = np.concatenate((cartesian, addlist))
        if props is not None:
            props.extend([props[i]]*len(addlist))
    return cartesian, props


def collapse(cartesians, ucell, ucell_inv=None, method='floor'):
    """Collapse Cartesian coordinates to the base cell.

    Parameters
    ----------
    cartesians : numpy.ndarray
        Cartesian coordinates to be collapsed. Shape (N, 2) or (N, 3).
    ucell : numpy.ndarray
        Basis vectors (rows) of the unit cell to be used for collapsing.
    ucell_inv : numpy.ndarray or None, optional
        The inverse of ucell. If not given or None, it is computed
        on the fly. Default is None.
    method : {'floor', 'round'}, optional
        Which method should be used for collapsing. 'floor'
        returns Cartesian coordinates collapsed to the (0, 0)
        cell, i.e., with all fractional coordinates in range
        [0, 1]. 'round' collapses the coordinates so that they
        are as close to zero as possible. Fractional coordinates
        are then in range [-0.5, 0.5].

    Returns
    -------
    collapsed_cartesians : numpy.ndarray
        The Cartesian coordinates, collapsed according to method
    fractional_coordinates : numpy.ndarray
        The corresponding fractional coordinates
    """
    if ucell_inv is None:
        ucell_inv = np.linalg.inv(ucell)
    fractional = collapse_fractional(cartesians.dot(ucell_inv),
                                     method=method)
    return fractional.dot(ucell), fractional


def _floor_eps(eps):
    """Return the floored-int version of values after adding eps."""
    def _floor(values):
        return np.floor(values + eps)
    return _floor


def collapse_fractional(coordinates, method='floor',
                        eps=COLLAPSE_EPS, in_place=False):
    """Collapse fractional coordinates to the base cell.

    Parameters
    ----------
    coordinates : numpy.ndarray
        Fractional coordinates to be collapsed. Shape (N, 2) or (N, 3).
    method : {'floor', 'round'}, optional
        Which method should be used for collapsing.  'round' collapses
        the coordinates so that they are as close to zero as possible.
        Fractional coordinates are then in range [-0.5, 0.5]. 'floor'
        returns fractional coordinates collapsed to the (0, 0) cell,
        that is, with all fractional coordinates in range [-`eps`,
        1-`eps`]. This ensures that collapsing a fractional coordinate
        of (essentially) 1 gives (essentially) zero. Default is 'floor'.
    eps : float or Sequence, optional
        Used only if method == 'floor'. The (fractional) tolerance
        for collapsing. If a sequence, it should have as many items
        as the second axis of coordinates. Default is 1e-8.
    in_place : bool, optional
        If True, the function modifies directly `coordinates`. Can be
        used to save some memory. Default is False.

    Returns
    -------
    collapsed : numpy.ndarray
        The collapsed fractional coordinates. Will be the same object
        passed in if in_place is True-thy.

    Raises
    ------
    ValueError
        If method is not one of the valid methods.
    """
    _methods = {'f': _floor_eps(eps), 'r': np.round}
    try:
        round_ = _methods[method[0]]
    except (TypeError, IndexError, KeyError):
        raise ValueError(f'collapse_fractional: Unknown method={method}')
    if in_place:
        collapsed = coordinates
        collapsed -= round_(coordinates)
    else:
        collapsed = coordinates - round_(coordinates)
    return collapsed


def ensure_integer_matrix(matrix, eps=1e-6):
    """Return a rounded version of matrix. Raise if matrix isn't integer."""
    rounded = np.round(matrix)
    if np.any(abs(matrix - rounded) > eps):
        raise NonIntegerMatrixError(matrix)
    return rounded


def pairwise(iterable):
    """Return a generator of pairs of subsequent elements."""
    try:
        # pairwise is available in py >= 3.10
        return itertools.pairwise(iterable)
    except AttributeError:
        pass
    # Rough equivalent from the python docs
    orig, offset = itertools.tee(iterable)
    next(offset, None)
    return zip(orig, offset)


def rotation_matrix(angle, dim=2):
    """Returns a (2x2) matrix for in-plane rotation of the given rotation
    angle. Set dim=3 to get a 3x3 matrix with rotation in [:2, :2]."""
    if dim < 2:
        raise ValueError('Rotation matrix needs at least dimension 2')
    m = np.eye(dim, dtype=float)
    m[:2, :2] = np.array([[np.cos(angle), -np.sin(angle)],
                          [np.sin(angle), np.cos(angle)]])
    return m


def rotation_matrix_order(order, dim=2):
    """Returns a (2x2) matrix for in-plane rotation of the given rotation
    order. Set dim=3 to get a 3x3 matrix with rotation in [:2, :2]."""
    angle = 2*np.pi/order
    return rotation_matrix(angle, dim=dim)


def fortranContLine(s):
    """Takes a sting that might be too long to fit on a single line of fortran
    code and splits it into continuation lines as necessary. Returns a string
    with ampersands and line breaks."""
    limit = 72  # fortran line length limit
    if len(s) <= limit:
        return s
    o = s[:6]
    s = s[6:]
    while len(s) > (limit-6):   # 6 characters for beginning of line
        o += s[:(limit-6)]
        o += "&\n     &"
        s = s[(limit-6):]
    o += s
    return o


def readVector(s, ucell=None, defRelaltive=False):
    """Takes a string 'xyz[f1 f2 f3]', 'abc[f1 f2 f3]' or just '[f1 f2 f3]'
    and returns the corresponding vector in cartesian coordinates, or None if
    the string cannot be parsed."""
    m = re.match(r'\s*(xyz|abc)?\[\s*(?P<v1>[-0-9.]+)\s+'
                 r'(?P<v2>[-0-9.]+)\s+(?P<v3>[-0-9.]+)\s*\]', s)
    if not m:
        return None
    try:
        v1 = float(m.group("v1"))
        v2 = float(m.group("v2"))
        v3 = float(m.group("v3"))
    except (ValueError, IndexError):
        return None
    if "abc" in s:
        if ucell is None:
            return None
        uct = ucell.T
        return np.dot((v1, v2, v3), ucell.T)
    # xyz
    return np.array([v1, v2, v3])


def cosvec(x, y):
    """
    Returns the cosine of the angle between two vectors

    Parameters
    ----------
    x : numpy.array
        First vector
    y : numpy.array
        Second vector

    Returns
    -------
    float
        Cosine of the angle between the two vectors

    """
    return np.dot(x, y) / (np.linalg.norm(x) * np.linalg.norm(y))


def dict_equal(d1, d2):                                                         # TODO: d1 == d2 works the same
    """
    Checks whether two dictionaries are equal, i.e. contain the same set of
    keys with the same values

    Parameters
    ----------
    d1 : dict
        First dictionary
    d2 : dict
        Second dictionary

    Returns
    -------
    bool
        True if all keys and values match, False otherwise
    """
    if len({k: d1[k] for k in d1 if k in d2 and d1[k] == d2[k]})-len(d1) == 0:
        return True
    return False


def lcm(a, b):
    "Calculate the lowest common multiple of two integers a and b"
    return a * b // np.gcd(a, b)


def max_diff(list_):
    """Find the index in a list such that list[i] - list[i-1] is the largest
    value difference in the list (after sorting). Return i and the difference.
    """
    if len(list_) < 2:
        return None, None
    list_ = sorted(list_)
    maxdiff = max([list_[i] - list_[i-1] for i in range(1, len(list_))])
    m = [i for i in range(1, len(list_)) if list_[i] - list_[i-1] == maxdiff]
    if m:
        return m[0], maxdiff
    return None, maxdiff


def parseMathSqrt(s):
    try:
        f = float(s)
    except ValueError:
        f = 1.0
        if '*' in s:
            sl = s.split('*')
        else:
            sl = [s]
        sl2 = []
        for el in sl:
            if '/' in el:
                sl2.append(el.split('/')[0])
                for s in el.split('/')[1:]:
                    sl2.append(1/parseMathSqrt(s))
            else:
                sl2.append(el)
        for el in sl2:
            try:
                f *= float(el)
            except ValueError:
                if 'sqrt' not in el:
                    logger.error('Could not interpret '+el+' as float value')
                    raise
                else:
                    p = re.compile(r'\s*sqrt\s*\(\s*(?P<value>[\d.]*)\s*\)')
                    m = p.match(el)
                    if m:
                        try:
                            f *= np.sqrt(float(m.group('value')))
                        except ValueError:
                            logger.error('Could not interpret ' +
                                         m.group('value')+' as float value')
                            raise
                    else:
                        logger.error('Could not interpret ' + el
                                     + ' as float value')
                        raise
    return f


def angle(v1, v2):
    """Returns the angle between two 2D vectors"""
    # Use cross product for sine, dot product for cosine
    return np.arctan2(v1[0]*v2[1] - v1[1]*v2[0], v1[0]*v2[0] + v1[1]*v2[1])


def dist_from_line(p1, p2, r):
    """
    Calculates the distance of a point from a line, with the line defined by
    two other points. Works in both 2 and 3 dimensions.

    Parameters
    ----------
    p1, p2 : np.arrays of size 2 or 3
        The points defining the line.
    r : np.array of same size as p1, p2
        The point for which the distance should be determined.

    Returns
    -------
    float
        The distance.
    """
    if len(p1) == 2:
        return (abs((p2[1] - p1[1]) * r[0] - (p2[0] - p1[0]) * r[1]
                    + p2[0] * p1[1] - p2[1] * p1[0])
                / np.sqrt((p2[1] - p1[1])**2 + (p2[0] - p1[0])**2))
    if len(p1) == 3:
        return (np.linalg.norm(np.cross((r - p1), (p2 - p1)))
                / np.linalg.norm(p2 - p1))
    raise ValueError("Vector dimensions have to be either 2 or 3.")


def readToExc(llist):                                                           # TODO: unused; could be an iterator
    """For reading PARAMETERS files; takes a list, returns elements until the
    first one that starts with an exclamation mark."""
    read = True
    newlist = []
    for s in llist:
        if read:
            if s[0] == '!':
                read = False
            else:
                newlist.append(s)
    return newlist


def splitSublists(llist, sep):                                                  # TODO: could be an iterator
    """Takes a list and a separator, splits strings in the list by the
    separator, returns results as list of lists"""
    newlist = []
    sublist = []
    for el in llist:
        if sep not in el:
            sublist.append(el)
        else:
            pl = el.split(sep)
            if pl[0]:
                sublist.append(pl[0])
            newlist.append(sublist)
            sublist = []
            if len(pl) > 1:
                for i in range(1, len(pl)):
                    s = pl[i]
                    if s:
                        sublist.append(s)
                    if i != len(pl) - 1:
                        newlist.append(sublist)
                        sublist = []
    newlist.append(sublist)
    return(newlist)


def splitMaxRight(s, sep):
    """Same as s.split(sep, maxsplit=1), but splitting at the first instance
    from the right."""
    sr = s[::-1]
    L = sr.split(sep, maxsplit=1)
    L.reverse()
    nl = []
    for ns in L:
        nl.append(ns[::-1])
    return nl


def recombineListElements(llist, com):
    """Takes a list, checks in each element whether the first/last characters
    are the given combination character, and if so, combines list elements with
    the list element before/after."""
    i = 0
    newlist = llist[:]
    while i < len(newlist)-1:
        if newlist[i][-1] == com or newlist[i+1][0] == com:
            newlist[i] += newlist.pop(i+1)
        else:
            i += 1
    return newlist


def addUnequalPoints(l1, l2, eps, uniqueLists=False):
    """Adds all points from l1 to l2, if they are not already in l2
    (+- epsilon)."""
    nl2 = l2[:]
    nl1 = l1[:]
    if len(l2) == 0:
        nl2 = nl1
    else:
        if not uniqueLists:
            # first get rid of duplicates in nl1
            tree = sps.KDTree(nl1)
            usepoint = [True] * len(nl1)
            for (i, p) in enumerate(nl1):
                if usepoint[i]:
                    for j in tree.query_ball_point(p, eps)[1:]:
                        usepoint[j] = False
            nl1 = list(itertools.compress(nl1, usepoint))
        # then add remaining elements to l2:
        dl = sps.distance.cdist(np.vstack(tuple(l1)), np.vstack(tuple(l2)),
                                'euclidean')
        for (i, sublist) in enumerate(dl):
            if min(sublist) >= eps:
                nl2.append(l1[i])
    return nl2


def available_cpu_count():
    """ Number of available virtual or physical CPUs on this system, i.e.
    user/real as output by time(1) when called with an optimally scaling
    userspace-only program"""

    # cpuset
    # cpuset may restrict the number of *available* processors
    try:
        with open('/proc/self/status') as f:
            m = re.search(r'(?m)^Cpus_allowed:\s*(.*)$', f.read())
        if m:
            res = bin(int(m.group(1).replace(',', ''), 16)).count('1')
            if res > 0:
                return res
    except IOError:
        pass

    # Python 2.6+
    try:
        return multiprocessing.cpu_count()
    except NotImplementedError:
        pass

    # POSIX
    try:
        res = int(os.sysconf('SC_NPROCESSORS_ONLN'))

        if res > 0:
            return res
    except (AttributeError, ValueError):
        pass

    # Windows
    try:
        res = int(os.environ['NUMBER_OF_PROCESSORS'])
        if res > 0:
            return res
    except (KeyError, ValueError):
        pass

    # BSD
    try:
        sysctl = subprocess.Popen(['sysctl', '-n', 'hw.ncpu'],
                                  stdout=subprocess.PIPE)
        scStdout = sysctl.communicate()[0]
        res = int(scStdout)
        if res > 0:
            return res
    except (OSError, ValueError):
        pass

    # Linux
    try:
        with open('/proc/cpuinfo') as f:
            res = f.read().count('processor\t:')
        if res > 0:
            return res
    except IOError:
        pass

    # Solaris
    try:
        pseudoDevices = os.listdir('/devices/pseudo/')
        res = 0
        for pd in pseudoDevices:
            if re.match(r'^cpuid@[0-9]+$', pd):
                res += 1
        if res > 0:
            return res
    except OSError:
        pass

    # Other UNIXes (heuristic)
    try:
        try:
            with open('/var/run/dmesg.boot') as f:
                dmesg = f.read()
        except IOError:
            dmesgProcess = subprocess.Popen(['dmesg'], stdout=subprocess.PIPE)
            dmesg = dmesgProcess.communicate()[0]
        res = 0
        while '\ncpu' + str(res) + ':' in dmesg:
            res += 1
        if res > 0:
            return res
    except OSError:
        pass

    return -1


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
    if sys.version_info >= (3, 8):
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


def make_unique_list(w_duplicates):                                             # TODO: better function in guilib.helpers
    """Helper function to remove duplicates from list. Does same as creating a set but preservers order.

    Args:
        w_duplicates (iterable): list with duplicates

    Returns:
        list: list with duplictes removed
    """
    unique_list = []
    for item in w_duplicates:
        if item not in unique_list:
            unique_list.append(item)
    return unique_list
