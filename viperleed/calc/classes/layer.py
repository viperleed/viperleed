# -*- coding: utf-8 -*-
"""Class storing position and atom list of a layer.

Created on Jun 13 2019
"""

import numpy as np

__authors__ = ["Florian Kraushofer (@fkraushofer)",]

class Layer:
    """To be used with Slab; has origin, atoms (a subset of the ones in slab),
    and a number. Sublayers also use the "Layer" class."""

    def __init__(self, slab, num, isBulk=False, sublayer=False):
        self.slab = slab
        self.num = num      # consecutive layer numbering, 0 being highest
        self.isBulk = isBulk    # defined by BLAY in PARAMETERS file
        self.atlist = []    # atoms in this layer
        # cartesian origin: xy from POSCAR with possible displacements,
        #  z is highest atom
        self.cartori = None
        self.cartbotz = None    # z position of lowest atom
        if sublayer:
            # list of candidate positions for rotation / mirror / glide planes
            self.symposlist = []

    def getLayerPos(self):
        """Gets a cartesian origin coordinate for the layer, using z of the
        highest atom. x,y are calculated from the origin of an a,b unit cell at
        that height. Also assigns posInLayer for all atoms in this layer."""
        al = self.atlist[:]     # temporary copy
        al.sort(key=lambda atom: atom.pos[2])
        topat = al[-1]
        botat = al[0]
        oripos = np.array([0., 0., topat.pos[2]])
        self.cartori = np.dot(self.slab.ucell, oripos)
        # this gets x and y correct, but z still in the wrong direction and
        #  with origin as POSCAR
        # -> just take the z from the highest atom.
        self.cartori[2] = topat.cartpos[2]
        self.cartbotz = botat.cartpos[2]
        return
