"""Module base of viperleed.calc.files.new_displacements.tokens."""

__authors__ = ('Alexandra Mia Imre (@alexmiame)',)
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-05-12'
__license__ = 'GPLv3+'

from abc import ABC, abstractmethod


class TokenParserError(Exception):
    """Base exception for all token‐parser errors in DISPLACEMENTS."""


class DisplacementsFileToken(ABC):
    """Base class for tokens for the DISPLACEMENTS file."""

    @abstractmethod
    def __eq__(self, other):
        """Compare self to other token object."""

    @abstractmethod
    def __str__(self):
        """Return representation of token object."""
