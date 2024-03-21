# -*- coding: utf-8 -*-
"""Module utils of viperleed.tleedmlib.classes.slab.

Created on 2023-02-22

@author: Michele Riva (@michele-riva)

Contains functions and classes useful for dealing with BaseSlab objects.
"""

import itertools
import numpy as np


def _cycle(sequence, start=0):                                                  # TODO: could it be useful in other places?
    """Return a generator that cycles though `sequence` beginning at start."""
    if sequence:
        start %= len(sequence)
    _cycled = itertools.cycle(sequence)
    return itertools.islice(_cycled, start, None)


def _left_handed(ab_cell):
    """Return whether a 2D unit cell is left-handed."""
    return np.linalg.det(ab_cell) < 0


def _z_distance(atom1, atom2):                                                  # TODO: could it be useful in other places? If yes, then maybe the right place is .atom?
    """Return the distance along z between two atoms."""
    return abs(atom2.cartpos[2] - atom1.cartpos[2])
