"""Module errors of viperleed.calc.classes.slab.

Contains exceptions specific to BaseSlab objects.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-02-21'
__license__ = 'GPLv3+'

from viperleed.calc.classes.atom_containers import AtomContainerError


class SlabError(AtomContainerError):
    """Base exception for Slab objects."""


class AlreadyMinimalError(SlabError, RuntimeError):
    """A minimization was requested, but the quantity is already minimal."""


class AtomsTooCloseError(SlabError):
    """At least two atoms are too close, irrespective of their element."""


class EmptySlabError(SlabError):
    """The slab has no atoms."""


class InvalidUnitCellError(SlabError):
    """Exception raised when the unit cell of a slab is inappropriate."""


class MissingBulkSlabError(SlabError, RuntimeError):
    """A bulkslab attribute would be needed but it is not available."""


class MissingElementsError(SlabError):
    """The slab has no elements."""


class MissingLayersError(SlabError):
    """An operation requires layers defined."""


class MissingSublayersError(SlabError):
    """An operation requires sublayers defined."""


class NoBulkRepeatError(SlabError, RuntimeError):
    """Exception raised when failing to find a bulk repeat vector."""


class TooFewLayersError(SlabError):
    """Some layers are present, but it's too few for the operation."""


class VacuumError(SlabError):
    """Something is wrong with the vacuum gap of a slab."""

    def __init__(self, message, fixed_slab):
        """Initialize exception instance.

        Parameters
        ----------
        message : str
            Textual message to be displayed for this exception
        fixed_slab : SurfaceSlab
            A slab that has been modified to fix the problem with
            the vacuum gap.

        Returns
        -------
        None.
        """
        self.message = message
        self.fixed_slab = fixed_slab
        super().__init__(message)


class NotEnoughVacuumError(VacuumError):
    """A slab has too little vacuum."""


class NoVacuumError(NotEnoughVacuumError):
    """A slab has no vacuum gap at all."""


class WrongVacuumPositionError(VacuumError):
    """A slab has enough vacuum, but it's somewhere in the middle."""
