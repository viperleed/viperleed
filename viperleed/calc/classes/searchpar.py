"""Module searchpar of viperleed.calc.classes.

Defines SearchPar, a class containing information about ONE parameter
of the tensor-LEED structural search. This module used to be part of
rparams.py (created on 2019-06-13). Refactored by Michele Riva in Oct
 2023.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-10-23'
__license__ = 'GPLv3+'

import logging

import numpy as np

_LOGGER = logging.getLogger(__name__)


class SearchPar:
    """Stores properties of ONE parameter of the search, i.e. what variation
    of what atom is linked to this parameter."""

    def __init__(self, atom, mode, el, deltaname):
        self.atom = atom
        self.mode = mode
        self.el = el
        self.deltaname = deltaname
        self.steps = 1
        self.edges = (None, None)  # the first and last value in the range
        self.center = 1  # the index closest to "no change" (Fortran index starting at 1)
        self.non_zero = False   # whether the center is truly "unchanged"
        self.restrictTo = None  # None, Index, or other search par
        self.linkedTo = None    # other search par linked via 'atom number'
        self.parabolaFit = {"min": None,
                            "err_co": np.nan, "err_unco": np.nan}
        d = {}
        if mode == "occ":
            el = next(iter(atom.disp_occ.keys()))  # look at any element
            self.steps = len(atom.disp_occ[el])
            self.center = atom.disp_center_index[mode][el] + 1 # (Fortran index starting at 1)
            self.non_zero = (abs(atom.disp_occ[el][self.center-1]
                                 - atom.site.occ[el]) >= 1e-4)
            edges = []
            for ind in (0, -1):
                edges.append(" + ".join("{:.2f} {}".format(
                    atom.disp_occ[e][ind], e) for e in atom.disp_occ
                    if atom.disp_occ[e][ind] > 0.005))
                if edges[-1] == "":
                    edges[-1] = "vac"
            self.edges = tuple(edges)
        else:
            if mode == "geo":
                d = atom.disp_geo
            elif mode == "vib":
                d = atom.disp_vib
        if len(d) > 0 and el != "vac":  # if vac: use defaults
            if el in d:
                k = el
            else:
                k = "all"
            self.steps = len(d[k])
            self.edges = (d[k][0], d[k][-1])
            if k not in atom.disp_center_index[mode]:
                self.center = atom.disp_center_index[mode]["all"] + 1
            else:
                self.center = atom.disp_center_index[mode][k] + 1
            self.non_zero = (np.linalg.norm(d[k][self.center-1]) >= 1e-4)
