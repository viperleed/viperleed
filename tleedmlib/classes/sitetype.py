# -*- coding: utf-8 -*-
"""
Created on Jun 13 2019

@author: Florian Kraushofer

Class storing properties of a site
"""

import logging
import re
import numpy as np

import tleedmlib as tl

logger = logging.getLogger("tleedm.sitetype")


class Sitetype:
    """Site types are identified by (main) element and name, and store
    vibrational amplitude and occupation"""

    def __init__(self, el, name):
        self.el = el.capitalize()
        self.name = name
        self.label = self.el + '_' + self.name
        self.vibamp = {}    # vibrational amplitude per element
        self.occ = {}       # occupation per element

        self.oriState = None    # deep copy of self before a search is applied
        self.mixedEls = []      # stores the relevant rparams.ELEMENT_MIX

    def __str__(self):
        return self.label

    def isEquivalent(self, site2):
        """Checks whether two sites are equivalent, i.e. have the same label
        and the same values for vibrational amplitudes and occupations."""
        if (self.label == site2.label
                and tl.base.dict_equal(self.vibamp, site2.vibamp)
                and tl.base.dict_equal(self.occ, site2.occ)):
            return True
        return False

    def getVibAmp(self, rp, chemel):
        """Calculates a default vibrational amplitude for this site for the
        given chemical element from the element atomic mass, the experimental
        temperature and the material's Debye temperature."""
        if rp.T_DEBYE is None or rp.T_EXPERIMENT is None:
            logger.error("Cannot generate default vibrational amplitudes: "
                         "Temperature or Debye temperature undefined.")
            raise ValueError("Temperature and Debye temperature must be float")
        if chemel in rp.ELEMENT_RENAME:
            el = rp.ELEMENT_RENAME[chemel].capitalize()
        else:
            el = chemel.capitalize()
        if el not in tl.leedbase.periodic_table:
            logger.error(
                "Cannot generate default vibrational amplitude for site "
                + self.label + ": Element " + el + " not recognized.")
            raise ValueError("Element " + el + " not recognized.")
        if el not in tl.leedbase.elementAtomicMass:
            logger.error(
                "Cannot generate default vibrational amplitude for site "
                + self.label + ": Element" + el + " atomic mass unknown.")
            raise NotImplementedError("Element" + el + " atomic mass unknown.")
        scaling = 1.0
        for s in rp.VIBR_AMP_SCALE:
            try:
                label = s.split()[0]
                value = float(s.split()[1])
            except (ValueError, IndexError):
                logger.error("Failed to interpret VIBR_AMP_SCALE parameter "
                             "part, ignoring input: " + s)
                continue
            if value <= 0:
                logger.warning(
                    "VIBR_AMP_SCALE parameter: scaling values have "
                    "to be positive floats, ignoring input: " + s)
                continue
            mstr = re.escape(label)
            mstr = mstr.replace('\\*', '.*')
            m = re.match(mstr, self.label)
            if m:
                if m.end(0) == len(self.label):  # real match
                    scaling = value
        self.vibamp[chemel] = round(
            scaling * (np.sqrt(np.sqrt(1+16*((rp.T_EXPERIMENT
                                              / rp.T_DEBYE)**2))
                               * 109.15 / (tl.leedbase.elementAtomicMass[el]
                                           * rp.T_DEBYE))), 3)
        return
