"""Module base of viperleed.calc.lib.

Contains generic functions used in the TensErLEED scripts.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2019-06-13'
__license__ = 'GPLv3+'

from contextlib import contextmanager
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


def parseMathSqrt(s):                                                           # TODO: replace with guilib math parser after refactor
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


def dist_from_line(p1, p2, r):                                                  # TODO: will be removed in #51
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


def readToExc(llist):                                                           # TODO: remove after vibrocc refactor; could be an iterator
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


# TODO: This function is confusing. It is used to 'recombine' elements
# that were space-split before into 'sep'-then-space-split lists. E.g.,
# in OCC_DELTA DISPLACEMENTS:
#    chem start stop step, other_chem start stop step
# --> by space first: ['chem', 'start', 'stop', 'step,',
#                      'other_chem', 'start', 'stop', 'step']
# --> splitSublists:  [['chem', 'start', 'stop', 'step'],
#                      ['other_chem', 'start', 'stop', 'step']]
# It would be way cleaner to do the splitting in the correct order:
# --> split first on comma: ['chem start stop step',
#                            'other_chem start stop step']
# then each item on spaces.
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


# TODO: remove next in #51
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


# TODO: consider using psutil instead. See
# https://psutil.readthedocs.io/en/latest/index.html#psutil.cpu_count
# Adapted from https://stackoverflow.com/questions/1006289/
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
