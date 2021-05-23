# -*- coding: utf-8 -*-
"""
Created on Jun 13 2019

@author: Florian Kraushofer

Contains generic functions used in the TensErLEED scripts.
"""

import logging
import numpy as np
import re
import subprocess
import multiprocessing
import os
import scipy.spatial as sps
import itertools

logger = logging.getLogger("tleedm.base")


###############################################
#                 CLASSES                     #
###############################################

class BackwardsReader:
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
                             + line).split("\n")
            except IOError:    # can't seek before the beginning of the file
                self.f.seek(0)
                self.data = (((self.f.read(self.size - (self.blksize
                                                        * (self.blkcount-1))))
                              .decode(self.encoding) + line).split("\n"))

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
                     .decode(self.encoding).split("\n"))
        # strip the last item if it's empty... a byproduct of the last line
        # having a newline at the end of it
        if not self.data[-1]:
            self.data.pop()


###############################################
#                FUNCTIONS                    #
###############################################

def rotMatrix(order):
    """Returns a (2x2) matrix for in-plane rotation of the given rotation
    order."""
    # these explicit definitions are likely useless, but sqrts might be
    #  marginally more accurate than sin/cos
    if order == 2:
        return np.array([[-1, 0], [0, -1]])
    elif order == -3:
        return np.array([[-0.5, -np.sqrt(3)/2], [np.sqrt(3)/2, -0.5]])
    elif order == 3:
        return np.array([[-0.5, np.sqrt(3)/2], [-np.sqrt(3)/2, -0.5]])
    elif order == -4:
        return np.array([[0, 1], [-1, 0]])
    elif order == 4:
        return np.array([[0, -1], [1, 0]])
    elif order == -6:
        return np.array([[0.5, np.sqrt(3)/2], [-np.sqrt(3)/2, 0.5]])
    elif order == 6:
        return np.array([[0.5, -np.sqrt(3)/2], [np.sqrt(3)/2, 0.5]])
    else:
        angle = 2*np.pi/order
        return np.array([[np.cos(angle), np.sin(angle)],
                         [-np.sin(angle), np.cos(angle)]])


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


def readIntRange(s):
    """Takes a string, returns a list of integers. If the string is a single
    int or space-separated ints, the return value is a list containing only
    those int. If the string contains i1-i2 or i1:i2, the list is
    range(i1, i2+1), i.e. contains both i1 and i2."""
    out = []
    sl = s.split()
    for ss in sl:
        try:
            out.append(int(ss))
        except ValueError:
            if "-" in ss:
                spl = ss.split("-")
                try:
                    out.extend(list(range(int(spl[0]), int(spl[1])+1)))
                except (ValueError, IndexError):
                    return []
            elif ":" in ss:
                spl = ss.split(":")
                try:
                    out.extend(list(range(int(spl[0]), int(spl[1])+1)))
                except (ValueError, IndexError):
                    return []
    return list(set(out))


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
        uct = ucell.transpose()
        vec = (v1*uct[0] + v2*uct[1] + v3*uct[2])
    else:  # xyz
        vec = np.array([v1, v2, v3])
    return vec


def readIntLine(line, width=3):
    """
    Reads an (arbitrary length) line of integers with fixed width. Will try
    to interpret everything as integers until the line ends.

    Parameters
    ----------
    line : str
        The line to interpret
    width : integer, optional
        The width of each integer. The default is 3.

    Returns
    -------
    Tuple of integers

    """
    line = line.rstrip()
    out = []
    try:
        while len(line) > 0:
            out.append(int(line[:width]))
            line = line[width:]
    except (ValueError, IndexError):
        raise
    return tuple(out)


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


def dict_equal(d1, d2):
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

# def angle(v1, v2, acute=True):
#     """Returns the angle between two vectors"""
#  !! UNSIGNED!
#     angle = np.arccos(np.dot(v1, v2)
#                       / (np.linalg.norm(v1)*np.linalg.norm(v2)))
#     if acute == True:
#         return angle
#     else:
#         return 2 * np.pi - angle


def angle(v1, v2):
    """Returns the angle between two 2D vectors"""
    # angle = np.arctan2(v2[1],v2[0]) - np.arctan2(v1[1],v1[0])
    # if abs(angle) > np.pi:
    #     angle += -np.sign(angle) * 2*np.pi
    # return angle
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
    elif len(p1) == 3:
        return (np.linalg.norm(np.cross((r - p1), (p2 - p1)))
                / np.linalg.norm(p2 - p1))
    else:
        raise ValueError("Vector dimensions have to be either 2 or 3.")


def readToExc(llist):
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


def splitSublists(llist, sep):
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
