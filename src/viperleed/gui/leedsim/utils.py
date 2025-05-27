"""Module utils of viperleed.gui.leedsim.

Defines functions used in multiple spots in the leedsim package.
Part of the contents of this module come from the 2025 refactor
of the gui package.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-05-24'
__license__ = 'GPLv3+'


import numpy as np


def screen_radius(energy, aperture):                                            # WILL BE renamed, and docstring needs work
    """
    Returns the radius of the projection of a LEED screen with a given
    aperture that contains LEED spots with 'energy'

    Parameters
    ----------
    energy : float
        primary energy of the electrons
    aperture: float, default=110.0
              degrees of aperture of the solid angle captured by the
              LEED screen. The current (2020-08-17) default value is taken
              from the dimensions of the screen of the ErLEED optics. The
              MCP SpectaLEED from Omicron appears to have an equivalent
              aperture of ~80°
    """
    electron_mass = 9.109e-31  # kg
    electron_charge = 1.60218e-19   # C
    hbar = 1.05457e-34  # J*s
    # The next one is the prefactor in AA^-1 eV^(-1/2)
    rt2me_hbar = np.sqrt(2*electron_mass*electron_charge) / hbar*1e-10

    return rt2me_hbar*np.sqrt(energy)*np.sin(np.radians(aperture)/2)


def sort_hk(index):
    """Key-sorting criterion for (h, k) indices.

    This method can be used as the 'key' optional argument of
    list.sort() or sorted() to sort beams given as (h, k)
    indices in the iterable to be sorted with decreasing h+k
    and, on second pass, decreasing h. E.g., applying this to
    [(-1, 1), (-1, -1), (1, -1), (1, 1)] gives the list
    [(1, 1), (1, -1), (-1, 1), (-1, -1)].

    Parameters
    ----------
    index : iterable
        Two-element iterable (h, k)

    Returns
    -------
    tuple
        Two elements, same types as h and k.
    """
    h_index, k_index = index
    return -(h_index + k_index), -h_index
