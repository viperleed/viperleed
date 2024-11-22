"""
======================================
  ViPErLEED Graphical User Interface
======================================
 *** module guilib.leedsim.classes.realspace ***

Defines the RealSpace class, used for displaying the real-space lattice

Author: Michele Riva
Created: 2021-03-13
"""

import numpy as np

from viperleed import guilib as gl

class RealSpace():
    def __init__(self, params):
        params = gl.LEEDParameters(params)

        self.superlattice = params['SUPERLATTICE']

        # Set up the lattices
        surf_basis = params['surfBasis']

        # self.fov is the real space field of view.
        # contains at least 4 unit cells
        self.fov = 4.2*max(surf_basis.ravel())
        self.surf = gl.Lattice(surf_basis,
                               space='real',
                               group=params['surfGroup'],
                               limit=self.fov)
        self.bulk = gl.Lattice(self.bulk_basis,
                               space='real',
                               group=params['bulkGroup'],
                               limit=self.fov)
        self.bulk.group.screws_glides = (params['bulk3Dsym'],
                                         self.bulk.cell_shape)

    @property
    def bulk_basis(self):
        """
        Real-space basis of the bulk lattice
        """
        if hasattr(self, 'bulk'):
            return self.bulk.basis
        return np.dot(np.linalg.inv(self.superlattice),
                      self.surf.basis).round(10)

    def angle_for_horizontal_bulk(self, direction):
        if direction not in (0, 1):
            raise ValueError("First positional argument of "
                             "angle_for_horizontal_bulk() "
                             f"should be 0 or 1. Got {direction!r} instead.")

        basis = self.bulk.basis[direction]
        # NB: np.arctan2 takes first the 'y' then the 'x'!
        theta = np.arctan2(basis[1], basis[0])

        return np.degrees(theta)

