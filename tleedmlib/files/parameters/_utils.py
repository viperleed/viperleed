# -*- coding: utf-8 -*-
"""Module parameters of viperleed.tleedmlib.

Created on Tue Aug 18 16:56:39 2020

@author: Florian Kraushofer
@author: Alexander M. Imre
@author: Michele Riva

Initial version by Florian Kraushofer in 2020, major rewrite by
Alexander Imre & Michele Riva in June 2023.

Functions for reading from and writing to the PARAMETERS file.
"""

import ast
from collections.abc import Sequence
from dataclasses import dataclass, field
from functools import partialmethod
import logging
from pathlib import Path
import re
import shutil

import numpy as np

from viperleed.tleedmlib.base import (strip_comments, splitSublists,
                                      readVector, readIntRange,
                                      recombineListElements)
from viperleed.tleedmlib.classes import rparams
from viperleed.tleedmlib.files.woods_notation import readWoodsNotation
from viperleed.tleedmlib import periodic_table
from viperleed.tleedmlib.sections._sections import TLEEDMSection as Section

from .errors import (
    ParameterError, ParameterValueError, ParameterParseError,
    ParameterIntConversionError, ParameterFloatConversionError,
    ParameterBooleanConversionError, ParameterNotRecognizedError,
    ParameterNumberOfInputsError, ParameterRangeError,
    ParameterUnknownFlagError, ParameterNeedsFlagError
    )


logger = logging.getLogger('tleedm.files.parameters')


@dataclass(frozen=True)
class NumericBounds:
    """A container for bounds of numeric parameters.

    Attributes
    ----------
    type_ : {int, float}, optional
        Which numeric type these bounds are for.
    range : tuple, optional
        Minimum and maximum value for numeric parameter.
        Either limit can be None, signifying no limitation.
        Default is (None, None).
    accept_limits : tuple
        Whether the min/max endpoints of the range are
        acceptable values. Default is (True, True).
    out_of_range_event : {'fail', 'coerce', 'modulo'}, optional
        How to behave in case a number goes out of bounds.
        Default is 'fail'.
    """
    type_: type = float
    range_: tuple = (None, None)
    accept_limits: tuple = (True, True)
    out_of_range_event: str = 'fail'
    _msg_: str = field(init=False)

    def __post_init__(self):
        """Process assigned attributes.

        Raises
        ------
        ValueError
            If out_of_range_event is not one of the acceptable values
        ValueError
            If out_of_range_event is 'modulo' and some of the bounds
            were not specified
        ValueError
            If type_ is not one of the acceptable values
        """
        if self.out_of_range_event not in ('fail', 'modulo', 'coerce'):
            raise ValueError(f'Invalid {self.out_of_range_event=}. '
                             'Must be "fail", "coerce", or "modulo"')
        if self.type_ not in (int, float):
            raise ValueError('type_ must be int or float')
        is_modulo = self.out_of_range_event == 'modulo'
        if is_modulo and any(limit is None for limit in self.range_):
            raise ValueError('out_of_range_event "modulo" needs both limits')
        if is_modulo and not all(self.accept_limits):
            raise ValueError('out_of_range_event "modulo" needs both limits '
                             'to be acceptable.')

        # Use object.__setattr__ for frozen data-class
        object.__setattr__(self, '_msg_', self._make_message())

    @property
    def coerce(self):
        """Return whether the out-of-bounds reaction is 'coerce'."""
        return self.out_of_range_event == 'coerce'

    @property
    def fail(self):
        """Return whether the out-of-bounds reaction is 'fail'."""
        return self.out_of_range_event == 'fail'

    @property
    def max(self):
        """Return the lower bound."""
        return self.range_[1]

    @property
    def min(self):
        """Return the lower bound."""
        return self.range_[0]

    @property
    def closed_range(self):
        """A range within which values should fall in, including extremes."""
        increment = 1 if self.type_ is int else 1e-4  # float
        _min, _max = self.range_
        if _min is not None and not self.accept_limits[0]:
            _min += increment
        if _max is not None and not self.accept_limits[1]:
            _max -= increment
        return _min, _max

    @property
    def unlimited(self):
        """Return whether there is any bound."""
        return all(limit is None for limit in self.range_)

    def format_out_of_range(self, value):
        """Return a string for an out-of-range value."""
        return f'Value {value} is {self._msg_}'

    def is_in_range(self, val):
        """Return whether a value is within the acceptable bounds."""
        _min, _max = self.closed_range
        return (_min is None or val >= _min,
                _max is None or val <= _max,)

    def make_in_range(self, value):
        """Return a version of value that is within range."""
        in_range = self.is_in_range(value)
        if self.unlimited or all(in_range):
            return value
        if self.fail and not all(in_range):
            raise RuntimeError

        if self.coerce:
            _, value, _ = sorted((value, *self.closed_range))
            return value

        assert self.out_of_range_event == 'modulo'
        _min, _max = self.range_
        return (value - _min) % (_max - _min) + _min

    def _make_message(self):
        """Return a message for value-out-of-bounds."""
        if self.unlimited:
            return ''

        accept_min, accept_max = self.accept_limits
        range_chars = ({True: '[', False: ']'}, {True: ']', False: '['})
        if all(limit is not None for limit in self.range_):
            return (f'not in range {range_chars[0][accept_min]}{self.min}, '
                    f'{self.max}{range_chars[1][accept_max]}')
        if self.min is not None and accept_min:
            return f'less than {self.min}'
        if self.min is not None:
            return f'less than or equal to {self.min}'

        assert self.max is not None
        if accept_max:
            return f'larger than {self.max}'
        return f'larger than or equal to {self.max}'


_POSITIVE_INT = NumericBounds(type_=int, range_=(1, None))
_POSITIVE_FLOAT = NumericBounds(type_=float, range_=(0, None))


@dataclass(frozen=True)
class Assignment:
    """Class to store the flags and values of a parameter.

    Attributes
    ----------
    values_str : str
        The right-hand side of the parameter assignment. Can also be
        passed as a tuple of strings upon instantiation. In this case
        elements will be joined.
    flags_str : str
        The left-hand side of the parameter assignment, excluding the
        parameter itself. Can also be passed as a tuple of strings upon
        instantiation. In this case elements will be joined.
    parameter : str
        The PARAMETER this assignment corresponds to.
    flags : tuple
        The flags of the parameter as a tuple of strings.
        This is automatically generated from the string
        arguments given at instantiation.
    values : tuple
        The values of the parameter as a tuple of strings.
        This is automatically generated from the string
        arguments given at instantiation.
    """
    values_str: str
    parameter: str
    flags_str: str = ''
    flags: tuple = field(init=False)
    values: tuple = field(init=False)

    def __post_init__(self):
        """Split out left- and right-hand sides into flags and values."""
        flags = self._unpack_assignment_side(self.flags_str)                    # TODO: are there any instances where we can have more than one flag per line?? If not we should raise always in that case!!
        values = self._unpack_assignment_side(self.values_str)

        object.__setattr__(self, 'flags', flags)
        object.__setattr__(self, 'values', values)

        # Make sure values_str and flags_str are actually strings:
        # we also accept Sequence of strings at __init__
        if not isinstance(self.values_str, str):
            object.__setattr__(self, 'values_str', ' '.join(self.values_str))
        if not isinstance(self.flags_str, str):
            object.__setattr__(self, 'flags_str', ' '.join(self.flags_str))
        if not self.parameter:
            raise ValueError('parameter must be a non-empty string')

    @property
    def flag(self):
        """Return the leftmost flag as a string."""
        return self.flags[0] if self.flags else ''

    @property
    def value(self):
        """Return the leftmost value as a string."""
        return self.values[0] if self.values else ''

    @property
    def other_flags(self):
        """Return all the flags except for the first one."""
        return self.flags[1:]

    @property
    def other_values(self):
        """Return all the values except for the first one (as strings)."""
        return self.values[1:]

    def _unpack_assignment_side(self, side):
        """Split the side left or right of the equal into bits."""
        if isinstance(side, str):
            return tuple(side.split())
        if isinstance(side, Sequence):
            return tuple(side)
        raise TypeError