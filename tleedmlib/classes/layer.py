# -*- coding: utf-8 -*-
"""Module layer of viperleed.tleedmlib.classes.

Created on Jun 13 2019

@author: Florian Kraushofer (@fkraushofer)
@author: Michele Riva (@michele-riva)

Classes storing position and atom list of a Layer, and its SubLayer
subclass. The latter is for atoms of the same chemical species at
the same z position.
"""


class LayerError(Exception):
    """Base exception for Layer objects."""


class LayerHasNoAtomsError(LayerError):
    """Operation cannot be performed as the layer is empty."""


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

    def __init__(self, slab, num, is_bulk=False, sublayer=False):
        """Initialize instance."""
        self.slab = slab
        self.num = num
        self.is_bulk = is_bulk
        self.atlist = []
        # cartesian origin: xy from POSCAR with possible displacements,
        #  z is highest atom
        self.cartori = None
        self.cartbotz = None    # z position of lowest atom
        if sublayer:
            # list of candidate positions for rotation / mirror / glide planes
            self.symposlist = []

    def __iter__(self):
        """Return an iterator of Atoms in this Layer."""
        return iter(self.atlist)

    @property
    def n_atoms(self):
        """Return the number of atoms in this layer."""
        return len(self.atlist)

    def update_position(self):
        """Update the Cartesian position of this layer.

        The z of the highest atom is used. x, y are calculated
        from the origin of an a, b unit cell at that height.

        Returns
        -------
        None.
        """
        if not self.atlist:
            raise LayerHasNoAtomsError(
                f'{type(self).__name__} needs atoms to update_position()'
                )
        sorted_atoms = sorted(self, key=lambda atom: atom.pos[2])
        topat = sorted_atoms[-1]
        botat = sorted_atoms[0]

        c_vec = self.slab.ucell.T[2]
        self.cartori = topat.pos[2] * c_vec
        # this gets x and y correct, but z still in the wrong direction and
        #  with origin as POSCAR
        # -> just take the z from the highest atom.
        self.cartori[2] = topat.cartpos[2]
        self.cartbotz = botat.cartpos[2]
