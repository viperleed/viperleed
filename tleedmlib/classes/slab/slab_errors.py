# -*- coding: utf-8 -*-
"""Module slab_errors of viperleed.tleedmlib.classes.slab.

Created on 2023-02-21

@author: Michele Riva (@michele-riva)

Contains exceptions specific to BaseSlab objects.
"""

from viperleed.tleedmlib.classes.atom_containers import AtomContainerError


class SlabError(AtomContainerError):
    """Base exception for Slab objects."""


class AlreadyMinimalError(SlabError, RuntimeError):
    """A minimization was requested, but the quantity is already minimal."""


class InvalidUnitCellError(SlabError):
    """Exception raised when the unit cell of a slab is inappropriate."""


class MissingBulkSlabError(SlabError, RuntimeError):
    """A bulkslab attribute would be needed but it is not available."""


class MissingLayersError(SlabError):
    """An operation requires layers defined."""


class MissingSublayersError(SlabError):
    """An operation requires sublayers defined."""


class NoBulkRepeatError(SlabError, RuntimeError):
    """Exception raised when failing to find a bulk repeat vector."""


class TooFewLayersError(SlabError):
    """Some layers are present, but it's too few for the operation."""
