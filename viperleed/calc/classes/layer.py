# -*- coding: utf-8 -*-
import numpy as np

__authors__ = ["Florian Kraushofer (@fkraushofer)",]
__created__ = "2019-06-13"

class Layer:
    """Class storing position and atom list of a layer.

    To be used with Slab; has origin, atoms (a subset of the ones in slab),
    and a number. Sublayers also use the "Layer" class.

    During the ViPErLEED and TensErLEED calculation, every Slab is be separated
    into layers by cutting the crystal parallel to the surface plane. The way
    this separation is preformed is primarily controlled by the
    :ref:`LAYER_CUTS parameter<ctrunc>` parameter (which is evaluated to to a
    :class:`viperleed.calc.classes.rparams.Rparams` attribute).
    While this separation is in principle arbitrary, it can strongly affect the
    execution time of the calculation.
"""

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
