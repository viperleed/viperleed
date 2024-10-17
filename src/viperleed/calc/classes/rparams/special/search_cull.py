"""Module search_cull of viperleed.calc.classes.rparams.special.

Defines the SearchCullType and SearchCull classes. They are convenience
classes for handling user input of the SEARCH_CULL parameter.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-12-14'
__license__ = 'GPLv3+'

from dataclasses import InitVar, dataclass, field
from enum import Enum  # Unfortunately StrEnum was introduced in py3.11
from numbers import Real

from viperleed.calc.lib.dataclass_utils import non_init_field

from ..defaults import NO_VALUE
from .base import SpecialParameter


class SearchCullType(Enum):
    """Enumeration of the implemented culling strategies."""

    GENETIC = 'genetic'
    CLONE = 'clone'
    RANDOM = 'random'
    UNSET = NO_VALUE

    @property
    def is_clone(self):
        """Return whether this is a GENETIC cull type."""
        return self is SearchCullType.CLONE

    @property
    def is_genetic(self):
        """Return whether this is a GENETIC cull type."""
        return self is SearchCullType.GENETIC

    @property
    def is_random(self):
        """Return whether this is a GENETIC cull type."""
        return self is SearchCullType.RANDOM


@dataclass(repr=False)
class SearchCull(SpecialParameter, param='SEARCH_CULL'):
    """A container for information used when culling a search population.

    Attributes
    ----------
    value : int or float
        When a integer, corresponds to the exact number of
        individuals to be culled from a search population.
        When a float, it is a fraction of the population.
        Use the nr_individuals(population) method to always
        retrieve the absolute number of individuals to be
        culled.
    type_ : SearchCullType
        Which culling scheme should be applied. Currently
        supported schemes are:
            CLONE (duplicate some individuals)
            GENETIC (combine traits from pairs of individuals)
            RANDOM (pick random traits)
    """

    _value: InitVar[Real]
    _type: SearchCullType = field(default=NO_VALUE)
    _int: int = non_init_field(default=NO_VALUE)
    _fraction: float = non_init_field(default= NO_VALUE)

    def __post_init__(self, _value):
        """Assign integer or fractional values for this SearchCull."""
        # Most of the functionality in here used to
        # be part of the ParameterInterpreter class
        try:
            cull_float = float(_value)
        except (TypeError, ValueError):
            raise TypeError('SearchCull value should be numeric') from None

        self._assign_int_or_fraction(cull_float)

        # Assign the type_ @property to make sure we correctly attempt
        # to convert the initialization argument to a SearchCullType
        self.type_ = self._type

    def _assign_int_or_fraction(self, cull_float):
        """Assign self._int or self._fraction from cull_float."""
        # See if we have an integer-valued cull
        cull_int = int(cull_float)
        if cull_float >= 1 and abs(cull_float - cull_int) < 1e-6:
            self._int = cull_int
        elif cull_float >= 1:
            raise ValueError('Values greater than one must be integers')
        elif cull_float >= 0:
            self._fraction = cull_float
        else:
            raise ValueError('SearchCull value must be non-negative')

    def __bool__(self):
        """Return whether this SearchCull would cause any culling."""
        return self.value > 0

    def __repr__(self):
        """Return a representation string for this SearchCull."""
        return f'{self.__class__.__name__}({self.value}, type_={self.type_})'

    @property
    def has_type(self):
        """Return whether this SearchCull has a type_ defined."""
        return self.type_ is not SearchCullType.UNSET

    @property
    def is_fractional(self):
        """Return whether this SearchCull is a fraction of the population."""
        return self._fraction is not NO_VALUE

    @property
    def type_(self):
        """Return the type of this SearchCull."""
        return self._type

    @type_.setter
    def type_(self, new_type):
        """Set the type of this SearchCull."""
        try:
            self._type = SearchCullType(new_type)
        except ValueError:
            raise ValueError(f'Unknown SearchCull type {new_type!r}') from None

    @property
    def value(self):
        """Return the (integer or fractional) value for this SearchCull."""
        return self._fraction if self.is_fractional else self._int

    @classmethod
    def from_value(cls, value):
        """Return a SearchCull from a Sequence value."""
        return cls(*value)

    def nr_individuals(self, population):
        """Return the number of individuals to be culled from population."""
        if not isinstance(population, int) or population < 0:
            raise ValueError(f'Invalid population={population!r}. '
                             'Should be a non-negative integer')
        ncull = (round(self.value * population) if self.is_fractional
                 else self.value)
        if ncull > population:
            raise ValueError(f'{self} too large for {population} individuals')
        return ncull
