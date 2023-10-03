# -*- coding: utf-8 -*-
"""Module slab_errors of viperleed.tleedmlib.classes.slab.

Created on 2023-02-21

@author: Michele Riva (@michele-riva)

Contains exceptions specific to BaseSlab objects.
"""

class SlabError(Exception):
    """Base exception for Slab objects."""


class InvalidUnitCellError(SlabError):
    """Exception raised when the unit cell of a slab is inappropriate."""


class NeedsLayersError(SlabError):
    """An operation requires layers defined."""


class NeedsSublayersError(SlabError):
    """An operation requires sublayers defined."""
