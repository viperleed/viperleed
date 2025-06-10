"""Base module for DISPLACEMENT file parsing."""

__authors__ = ('Alexander M. Imre (@amimre)',)
__created__ = '2025-05-12'

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
