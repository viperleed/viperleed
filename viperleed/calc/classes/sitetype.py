"""Module sitetype of viperleed.calc.classes."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2019-06-13'
__license__ = 'GPLv3+'

import logging
import re

import numpy as np

from viperleed.calc.lib.periodic_table import ATOMIC_MASS
from viperleed.calc.lib.periodic_table import COVALENT_RADIUS
from viperleed.calc.lib.periodic_table import PERIODIC_TABLE

logger = logging.getLogger(__name__)


class Sitetype:
    """Class storing properties of a site.

    Site types are identified by (main) element and name, and store
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
                and self.vibamp == site2.vibamp
                and self.occ == site2.occ):
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
        if el not in PERIODIC_TABLE:
            logger.error(
                "Cannot generate default vibrational amplitude for site "
                + self.label + ": Element " + el + " not recognized.")
            raise ValueError("Element " + el + " not recognized.")
        if el not in ATOMIC_MASS:
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
                               * 109.15 / (ATOMIC_MASS[el]
                                           * rp.T_DEBYE))), 3)
        return


class Atom_type(Sitetype):
    """
    Collection of atoms with the same element and same site; inherited fromm Sitetype. Used for EEASiSSS input.
    """
    def __init__(self, el, name, new_bulk, layer = None):
        super().__init__(el, name)
        self.label = self.el + '_in_' + self.name
        self.covalent_radius = COVALENT_RADIUS[el.capitalize()]
        self.S = 0.32
        self.fxc = 1.0
        self.smallest_NN_dist = 1e10 # smallest NN distance of any atom of this type
        self.atoms = []
        self.bulk_layer = layer

        try:
            self.atomic_number = PERIODIC_TABLE.index(el.capitalize())+1
        except ValueError:
            logger.error('Invalid element symbol during Atom type assignment')
            raise ValueError('Invalid element symbol: ' + el.capitalize())

        self.new_bulk = new_bulk # attribute used in psgen; distinguishes if this atom is part of the auto-added bulk

    @property
    def getVibAmp(self, rp, chemel):
        raise AttributeError("Atom_type class objects are not intended to give a vibration amplitude. Use Site instead.")

    @property
    def isEquivalent(self, site2): #TODO: remove, not needed (I think)
        """Checks whether two sites are equivalent, i.e. have the same label
        and the same values for vibrational amplitudes and occupations."""
        if (self.label == site2.label
                and self.vibamp == site2.vibamp
                and self.occ == site2.occ
                and self.el == site2.el
                and self.new_bulk == site2.new_bulk):
            return True
        return False


    def set_MT_params(self, rmtmin, rmtmax, S):
        self.rmtmin = rmtmin
        self.rmtmax = rmtmax
        self.S = S

    def add_atom(self, atom, NN_dist):
        self.atoms.append(atom)
        if NN_dist < self.smallest_NN_dist:
            self.smallest_NN_dist = NN_dist


    def get_atomic_number(self):
        return self.atomic_number

    def get_type_NN_dist(self):
        return self.smallest_NN_dist

    def get_layer(self):
        return self.bulk_layer
