"""Module for reading and writing Woods notation.

The functionality in this module used to be part of calc.lib.leedbase.
"""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-06-07'
__license__ = 'GPLv3+'

import re
import logging

import numpy as np

from viperleed.calc.lib.base import parseMathSqrt, angle, cosvec

logger = logging.getLogger(__name__)

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
    return mat

def writeWoodsNotation(ucell):
    """Takes a unit cell (as a (2x2) matrix) and attempts to write it in Woods
    Notation. Returns empty string if no Woods notation is found."""
    # !!! VERY INCOMPLETE, should at least detect simple c(a x b) cases
    # !!! Same functionality exists in guilib; replace at some point
    if ucell[1, 0] == 0 and ucell[0, 1] == 0:
        return("(" + str(int(ucell[0, 0])) + "x" + str(int(ucell[1, 1])) + ")")
    else:
        return ""
