"""Module base of viperleed.files.displacements.tokens."""

__authors__ = ("Alexander M. Imre (@amimre)",)
__copyright__ = "Copyright (c) 2019-2025 ViPErLEED developers"
__created__ = "2025-05-12"
__license__ = "GPLv3+"

from abc import ABC, abstractmethod

class TokenParserError(ValueError):
    """Base exception for all token‐parser errors in DISPLACEMENTS."""


class DisplacementsFileToken(ABC):
    """Base class for tokens for the DISPLACEMENTS file."""

    @abstractmethod
    def __eq__(self, other):
        """Compare self to other token object."""

    @abstractmethod
    def __repr__(self):
        """Return representation of token object."""
