"""Module layer of viperleed.calc.classes.

Classes storing position and atom list of a Layer, and its SubLayer
subclass. The latter is for atoms of the same chemical species at
the same z position.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2019-06-13'
__license__ = 'GPLv3+'

from viperleed.calc.classes.atom_containers import AtomContainer
from viperleed.calc.classes.atom_containers import AtomContainerError


class LayerError(AtomContainerError):
    """Base exception for Layer objects."""


class LayerHasNoAtomsError(LayerError):
    """Operation cannot be performed as the layer is empty."""


class Layer(AtomContainer):                                                     # TODO: modify description of .cartori when flipping .cartpos[2] -- Issue #174
    """A container of atoms residing close to one another along z.

    This is intended to be used with Slab objects. Has origin, atoms
    (a subset of the ones in slab), and a number.

    Attributes
    ----------
    atlist : list of Atom
        The atoms that belong to this layer.
    cartori : numpy.ndarray
        The Cartesian origin of this layer. It is the position where a
        plane passing through the topmost atom in this layer intersects
        the c vector of the slab's unit cell. Notice that `cartori` is
        in the same Cartesian frame as the atoms (i.e., z increases
        going from vacuum into the solid).
        A call to update_position() updates this attribute.
    cartbotz : float
        Z (i.e., out-of-plane) position of the bottom-most atom in
        the layer. A call to update_position() updates this attribute.
    is_bulk : bool
        Whether this layer has bulk character.
    num : int
        A progressive index (zero-based) identifying this layer
        within its slab. Normally, `layer.num == 0` for the layer
        closest to the solid/vacuum interface.
    slab : Slab
        The slab to which this layer belongs.
    """

    def __init__(self, slab, num, is_bulk=False):
        """Initialize instance."""
        self.slab = slab
        self.num = num
        self.is_bulk = is_bulk
        self.atlist = []

        # Note about the next two: they should be kept in memory,
        # rather than re-determined on-the-fly, since the position
        # of the atoms in question may change during a search, but
        # TensErLEED expects atom positions to be referred to the
        # layer position as it was used in the reference calculation
        self.cartori = None
        self.cartbotz = None

    def __contains__(self, value):
        """Return whether this layer contains an atom."""
        return value in self.atlist

    def __iter__(self):
        """Return an iterator of Atoms in this Layer."""
        return iter(self.atlist)

    @property
    def n_atoms(self):
        """Return the number of atoms in this layer."""
        return len(self.atlist)

    @property
    def thickness(self):
        """Return the z thickness of this layer.

        Note
        ----
        This is calculated from the stored positions of this layer,
        and may not be up to date if atoms in this layer have moved
        along z since the last call to update_position().
        """
        return abs(self.cartbotz - self.cartori[2])

    def update_position(self):
        """Update the Cartesian position of this layer from its atoms."""
        if not self.atlist:
            raise LayerHasNoAtomsError(
                f'{type(self).__name__} needs atoms to update_position()'
                )
        sorted_atoms = sorted(self, key=lambda atom: atom.pos[2])
        topat = sorted_atoms[-1]
        botat = sorted_atoms[0]

        self.cartbotz = botat.cartpos[2]
        self.cartori = topat.pos[2] * self.slab.c_vector
        # So far x and y are correct, but z is still in the wrong               # TODO: this will not be necessary when we flip .cartpos[2] -- Issue #174
        # direction and with origin as POSCAR. Take the z directly
        # from the highest atom
        self.cartori[2] = topat.cartpos[2]


class SubLayer(Layer):                                                          # TODO: modify description of .cartori when flipping .cartpos[2] -- Issue #174
    """A Layer with the same chemical element and the same z position.

    Attributes
    ----------
    atlist : list of Atom
        The atoms that belong to this sublayer.
    cartori : numpy.ndarray
        The Cartesian origin of this layer. It is the position where a
        plane passing through the topmost atom in this layer intersects
        the c vector of the slab's unit cell. Notice that `cartori` is
        in the same Cartesian frame as the atoms (i.e., z increases
        going from vacuum into the solid).
        A call to update_position() updates this attribute.
    cartbotz : float
        Z (i.e., out-of-plane) position of the bottom-most atom in
        the layer. A call to update_position() updates this attribute.
    cartpos : numpy.ndarray
        The Cartesian position of this sublayer. Notice that this is
        not an alias for `cartori`. The (x, y) components of the two
        are not the same: `cartori[:2]` is a point along the c vector
        of the unit cell, while `cartpos[:2]` may be a point inside
        (or outside if atoms are not collapsed) the unit cell.
    element : str
        The chemical element of the atoms of this sublayer.
    is_bulk : bool
        Whether this layer has bulk character.
    num : int
        A progressive index (zero-based) identifying this sublayer
        within its slab. Normally, `sublayer.num == 0` for the
        sublayer closest to the solid/vacuum interface.
    pos : numpy.ndarray
        The fractional position of this sublayer. This is the
        fractional position corresponding to `cartpos`, not the
        one corresponding to `cartori`.
    slab : Slab
        The slab to which this sublayer belongs.
    symposlist : list
        Candidate positions for rotation / mirror / glide planes.
    """

    def __init__(self, slab, num, is_bulk=False):
        """Initialize instance."""
        super().__init__(slab, num, is_bulk=is_bulk)
        self.symposlist = []

    def __str__(self):
        """Return a string version of this sublayer."""
        n_atoms = self.n_atoms
        if not n_atoms:
            atoms = '(empty)'
        elif n_atoms == 1:
            atoms = self.atlist[0]
        else:
            atoms = [f'{at}' for at in self]
        return f'Sublayer(num={self.num}, {atoms}, cartbotz={self.cartbotz})'

    @property
    def element(self):
        """Return the chemical element for this sublayer."""
        if not self.atlist:
            raise LayerHasNoAtomsError(
                f'A {type(self).__name__} without atoms has no element'
                )
        return self.atlist[0].el

    @property
    def cartpos(self):
        """Return the Cartesian position of this sublayer."""
        if not self.atlist:
            raise LayerHasNoAtomsError(
                f'A {type(self).__name__} without atoms has no cartpos'
                )
        return self.atlist[0].cartpos

    @property
    def pos(self):
        """Return the fractional position of this sublayer."""
        if not self.atlist:
            raise LayerHasNoAtomsError(
                f'A {type(self).__name__} without atoms has no pos'
                )
        return self.atlist[0].pos
