"""Module r_factor_type of viperleed.calc.classes.rparams.special.

Defines the RFactorType enumeration, used for parameter R_FACTOR_TYPE.
"""

__authors__ = ('Alexandra M. Imre (@alexmiame)',)
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2025-11-11'
__license__ = 'GPLv3+'

from enum import Enum

from .base import SpecialParameter


class RFactorType(Enum):
    """Enumeration of the R-factor types.

    The value of each member is the integer code that TensErLEED
    uses for the same R-factor type. Members can be created from
    their value, or from a string. Strings may be the code itself
    (e.g., '1'), or a case-insensitive name or synonym (e.g.,
    'pendry', 'Zanazzi-Jona', 'R²').

    Examples
    --------
    >>> RFactorType(1) is RFactorType('rp') is RFactorType.PENDRY
    True
    """

    PENDRY = 1
    R2 = 2
    ZJ = 3
    SMOOTH = 4

    @classmethod
    def _missing_(cls, value):
        """Return the member for a string `value`, or None."""
        if not isinstance(value, str):
            return None
        value = value.strip()
        if value.isdigit():
            try:
                return cls(int(value))
            except ValueError:
                return None
        return _MEMBER_FROM_NAME.get(_norm_name(value))


def _norm_name(name):
    """Return a normalized version of `name` for robust matching."""
    name = name.strip().casefold()
    return (
        name.replace('²', '2')
        .replace('^', '')
        .replace('_', '')
        .replace('-', '')
        .replace(' ', '')
    )


_ALIASES = {
    RFactorType.PENDRY: ('pendry', 'pendry_r', 'rp', 'pend', 'p'),
    RFactorType.R2: ('R2', 'r2', 'r squared'),
    RFactorType.ZJ: ('zj', 'zanazzi jona', 'zanazzi-jona', 'zanazzi_jona',
                     'zannazi', 'z'),
    RFactorType.SMOOTH: ('smooth', 'smoothed', 'smooth pendry',
                         'smooth-pendry', 'schmid', 's'),
    }
_MEMBER_FROM_NAME = {
    _norm_name(alias): member
    for member, aliases in _ALIASES.items()
    for alias in aliases
    }


class _RFactorTypeParameter(SpecialParameter, param='R_FACTOR_TYPE'):
    """Register RFactorType as the class of parameter R_FACTOR_TYPE."""

    @classmethod
    def from_value(cls, value):
        """Return the RFactorType member corresponding to `value`."""
        return RFactorType(value)
