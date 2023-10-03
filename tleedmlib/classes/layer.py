# -*- coding: utf-8 -*-
"""Module layer of viperleed.tleedmlib.classes.

Created on Jun 13 2019

@author: Florian Kraushofer (@fkraushofer)
@author: Michele Riva (@michele-riva)

Classes storing position and atom list of a Layer, and its SubLayer
subclass. The latter is for atoms of the same chemical species at
the same z position.
"""

import numpy as np


class Layer:
    """A container of atoms residing close to one another along z.

    This is intended to be used with Slab objects. Has origin, atoms
    (a subset of the ones in slab), and a number.

    Attributes
    ----------
    atlist : list of Atom
        The atoms that belong to this layer.
    cartori : numpy.ndarray                                                     @fkraushofer: could you please elaborate a bit more on how (x,y) and z are defined for this one? It's not 100% clear to me.
        The Cartesian origin of this layer. This is derived from
        the position of its highest atom. A call to update_position()
        updates this attribute.
    cartbotz : float
        Z (i.e., out-of-plane) position of the bottom-most atom in
        the layer. A call to update_position() updates this attribute.
    is_bulk : bool
        Whether this layer has bulk character.
    num : int
        A progressive index (zero-based) identifying this layer
        within its slab. Normally, layer.num == 0 for the layer
        closest to the solid/vacuum interface.
    slab : Slab
        The slab to which this layer belongs.
    """

    def __init__(self, slab, num, isBulk=False, sublayer=False):
        """Initialize instance."""
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
