"""Module realspace of viperleed.gui.leedsim.classes.

Defines the RealSpace class, containing information useful for
displaying a real-space 2D lattice.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2021-03-13'
__license__ = 'GPLv3+'

import numpy as np

from viperleed.gui.classes.lattice2d import Lattice2D
from viperleed.gui.leedsim.classes.leedparameters import LEEDParameters


class RealSpace():
    """Use to plot a real-space lattice."""

    def __init__(self, params):
        """Initialize RealSpace instance.

        Parameters
        ----------
        params : dict, ConfigParser or LEEDParameters
            The LEED parameters defining this representation
            of real-space bulk and surface lattices.

        Returns
        -------
        None.
        """
        params = LEEDParameters(params)

        self.superlattice = params['SUPERLATTICE']

        # Set up the lattices
        surf_basis = params['surfBasis']

        # self.fov is the real space field of view.
        # contains at least 4 unit cells
        self.fov = 4.2*max(surf_basis.ravel())
        self.surf = Lattice2D(surf_basis,
                              space='real',
                              group=params['surfGroup'],
                              limit=self.fov)
        self.bulk = Lattice2D(self.bulk_basis,
                              space='real',
                              group=params['bulkGroup'],
                              limit=self.fov)
        self.bulk.group.set_screws_glides(params['bulk3Dsym'],
                                          self.bulk.cell_shape)

    @property
    def bulk_basis(self):
        """Return the real-space basis of the bulk lattice."""
        if hasattr(self, 'bulk'):
            return self.bulk.basis
        return np.dot(np.linalg.inv(self.superlattice),
                      self.surf.basis).round(10)

    def angle_for_horizontal_bulk(self, direction):
        """Return the angle that brings a basis vector horizontal.

        Parameters
        ----------
        direction : {0, 1}
            Which of the basis vectors should be considered

        Returns
        -------
        angle : float
            Angle in degrees. Rotating by angle in-plane,  TRUE??
            aligns basis[direction] with the x axis of the
            Cartesian coordinate system.

        Raises
        ------
        ValueError
            If direction is not one between 0 and 1
        """
        if direction not in (0, 1):
            raise ValueError("First positional argument of "
                             "angle_for_horizontal_bulk() "
                             f"should be 0 or 1. Got {direction!r} instead.")

        basis = self.bulk.basis[direction]
        # NB: np.arctan2 takes first the 'y' then the 'x'!
        theta = np.arctan2(basis[1], basis[0])

        return np.degrees(theta)
