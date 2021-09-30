# -*- coding: utf-8 -*-
"""
Created on Jun 13 2019

@author: Florian Kraushofer

Class storing properties of a site
"""

import logging
import re
import numpy as np

import viperleed.tleedmlib as tl

logger = logging.getLogger("tleedm.sitetype")

element_symbols = {'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10, 'Na': 11,
                       'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18, 'K': 19, 'Ca': 20, 'Sc': 21,
                       'Ti': 22, 'V': 23, 'Cr': 24, 'Mn': 25, 'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29, 'Zn': 30,
                       'Ga': 31, 'Ge': 32, 'As': 33, 'Se': 34, 'Br': 35, 'Kr': 36, 'Rb': 37, 'Sr': 38, 'Y': 39,
                       'Zr': 40, 'Nb': 41, 'Mo': 42, 'Tc': 43, 'Ru': 44, 'Rh': 45, 'Pd': 46, 'Ag': 47, 'Cd': 48,
                       'In': 49, 'Sn': 50, 'Sb': 51, 'Te': 52, 'I': 53, 'Xe': 54, 'Cs': 55, 'Ba': 56, 'La': 57,
                       'Ce': 58, 'Pr': 59, 'Nd': 60, 'Pm': 61, 'Sm': 62, 'Eu': 63, 'Gd': 64, 'Tb': 65, 'Dy': 66,
                       'Ho': 67, 'Er': 68, 'Tm': 69, 'Yb': 70, 'Lu': 71, 'Hf': 72, 'Ta': 73, 'W': 74, 'Re': 75,
                       'Os': 76, 'Ir': 77, 'Pt': 78, 'Au': 79, 'Hg': 80, 'Tl': 81, 'Pb': 82, 'Bi': 83, 'Po': 84,
                       'At': 85, 'Rn': 86, 'Fr': 87, 'Ra': 88, 'Ac': 89, 'Th': 90, 'Pa': 91, 'U': 92, 'Np': 93,
                       'Pu': 94, 'Am': 95, 'Cm': 96, 'Bk': 97, 'Cf': 98, 'Es': 99, 'Fm': 100, 'Md': 101, 'No': 102,
                       'Lr': 103, 'Rf': 104, 'Db': 105, 'Sg': 106, 'Bh': 107, 'Hs': 108, 'Mt': 109, 'Ds': 110,
                       'Rg': 111, 'Cn': 112, 'Nh': 113, 'Fl': 114, 'Mc': 115, 'Lv': 116, 'Ts': 117, 'Og': 118}

class Sitetype:
    """Site types are identified by (main) element and name, and store
    vibrational amplitude and occupation"""

    def __init__(self, el, name):
        self.el = el.capitalize() # Element
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

class Atom_type(Sitetype):
    """
    Collection of atoms with the same element and same site; inherited fromm Sitetype. Used for EEASiSSS input.
    """
    def __init__(self, el, name, new_bulk):
        super().__init__(el, name)
        self.label = self.el + '_in_' + self.name
        self.rmtmin = 1.0
        self.rmtmax = 2.5
        self.S = 0.32
        self.fxc = 1.0
        self.atoms = []
        self.atomic_number = element_symbols[el]
        self.new_bulk = new_bulk # attribute used in psgen; distinguishes if this atom is part of the auto-added bulk

    @property
    def getVibAmp(self, rp, chemel):
        raise AttributeError("Atom_type class objects are not intended to give a vibration amplitude. Use Site instead.")

    @property
    def isEquivalent(self, site2): #TODO: remove, not needed (I think)
        """Checks whether two sites are equivalent, i.e. have the same label
        and the same values for vibrational amplitudes and occupations."""
        if (self.label == site2.label
                and tl.base.dict_equal(self.vibamp, site2.vibamp)
                and tl.base.dict_equal(self.occ, site2.occ)
                and self.el == site2.el
                and self.new_bulk == site2.new_bulk):
            return True
        return False


    def set_MT_params(self, rmtmin, rmtmax, S):
        self.rmtmin = rmtmin
        self.rmtmax = rmtmax
        self.S = S

    def add_atom(self, atom):
        self.atoms.append(atom)

    def get_atomic_number(self):
        return self.atomic_number