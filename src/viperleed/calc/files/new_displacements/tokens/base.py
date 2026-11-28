"""Module base of viperleed.calc.files.new_displacements.tokens."""

__authors__ = ('Alexandra Mia Imre (@alexmiame)',)
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2025-05-12'
__license__ = 'GPLv3+'

from abc import ABC, abstractmethod
import math
import re

_PLAIN_DECIMAL = re.compile(r'[+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?')


class TokenParserError(Exception):
    """Base exception for all token‐parser errors in DISPLACEMENTS."""


def parse_plain_float(text):
    """Return `text` as a finite float in plain decimal notation.

    Unlike float(), this rejects "nan", "inf", "infinity", and digit
    separators such as "1_0".

    Parameters
    ----------
    text : str
        The text of a single number.

    Returns
    -------
    float
        The value of `text`.

    Raises
    ------
    ValueError
        If `text` is not a plain decimal number, or if it is not
        finite (e.g., "1e999").
    """
    if not _PLAIN_DECIMAL.fullmatch(text):
        msg = f'Not a plain decimal number: "{text}".'
        raise ValueError(msg)
    value = float(text)
    if not math.isfinite(value):
        msg = f'Not a finite number: "{text}".'
        raise ValueError(msg)
    return value


class DisplacementsFileToken(ABC):
    """Base class for tokens for the DISPLACEMENTS file."""

    @abstractmethod
    def __eq__(self, other):
        """Compare self to other token object."""

    @abstractmethod
    def __str__(self):
        """Return representation of token object."""
