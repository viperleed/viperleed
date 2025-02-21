"""Module utils of viperleed.calc.files.parameters.

Initial version by @fkraushofer on 2020-08-18, major rewrite by
@amimre and @michele-riva in June 2023. This module used to be
part of parameters.py. Refactored in October 2023.

Contains functions and classes used in multiple submodules of
the viperleed.calc.files.parameters package.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-10-16'
__license__ = 'GPLv3+'

from collections.abc import Sequence

from viperleed.calc.lib.dataclass_utils import frozen
from viperleed.calc.lib.dataclass_utils import non_init_field
from viperleed.calc.lib.dataclass_utils import set_frozen_attr


# TODO: some of these classes are probably also useful for other
# files, possibly with little modification. If they are, they
# could go higher up the hierarchy into a file.utils.py or similar


@frozen
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
    _msg_: str = non_init_field()

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
            raise ValueError(
                f'Invalid out_of_range_event={self.out_of_range_event}. '
                'Must be "fail", "coerce", or "modulo"'
                )
        if self.type_ not in (int, float):
            raise ValueError('type_ must be int or float')
        is_modulo = self.out_of_range_event == 'modulo'
        if is_modulo and any(limit is None for limit in self.range_):
            raise ValueError('out_of_range_event "modulo" needs both limits')
        if is_modulo and not all(self.accept_limits):
            raise ValueError('out_of_range_event "modulo" needs both limits '
                             'to be acceptable.')

        set_frozen_attr(self, '_msg_', self._make_message())

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

        if self.out_of_range_event == 'modulo':
            _min, _max = self.range_
            return (value - _min) % (_max - _min) + _min

        assert self.coerce
        _min, _max = self.closed_range
        if _min is not None:
            value = max(value, _min)
        if _max is not None:
            value = min(value, _max)
        return value

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


POSITIVE_INT = NumericBounds(type_=int, range_=(1, None))
POSITIVE_FLOAT = NumericBounds(type_=float, range_=(0, None))


@frozen
class Assignment:
    """Class to store the flags and values of a parameter.

    Attributes
    ----------
    values_str : str
        The right-hand side of the parameter assignment. Can also be
        passed as a tuple of strings upon instantiation. In this case
        elements will be joined.
    parameter : str
        The PARAMETER this assignment corresponds to.
    raw_line : str
        The full line, as read from the PARAMETERS file, possibly
        excluding comments.
    flags_str : str
        The left-hand side of the parameter assignment, excluding the
        parameter itself. Can also be passed as a tuple of strings upon
        instantiation. In this case elements will be joined.
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
    raw_line: str
    flags_str: str = ''
    flags: tuple = non_init_field()
    values: tuple = non_init_field()

    def __post_init__(self):
        """Split out left- and right-hand sides into flags and values."""
        try:
            value = self.parameter.strip()
        except AttributeError as exc:  # Not a string
            raise TypeError('parameter must be a string') from exc
        if not value:
            raise ValueError('parameter must contain printable characters')
        set_frozen_attr(self, 'parameter', value)

        flags = self._unpack_assignment_side(self.flags_str)
        values = self._unpack_assignment_side(self.values_str)

        set_frozen_attr(self, 'flags', flags)
        set_frozen_attr(self, 'values', values)

        # Make sure values_str and flags_str are actually strings:
        # we also accept Sequence of strings at __init__
        if not isinstance(self.values_str, str):
            set_frozen_attr(self, 'values_str', ' '.join(self.values_str))
        if not isinstance(self.flags_str, str):
            set_frozen_attr(self, 'flags_str', ' '.join(self.flags_str))

    @property
    def flag(self):
        """Return the leftmost flag as a string."""
        # pylint: disable-next=unsubscriptable-object  # Can't infer
        return self.flags[0] if self.flags else ''

    @property
    def value(self):
        """Return the leftmost value as a string."""
        # pylint: disable-next=unsubscriptable-object  # Can't infer
        return self.values[0] if self.values else ''

    @property
    def other_flags(self):
        """Return all the flags except for the first one."""
        # pylint: disable-next=unsubscriptable-object  # Can't infer
        return self.flags[1:]

    @property
    def other_values(self):
        """Return all the values except for the first one (as strings)."""
        # pylint: disable-next=unsubscriptable-object  # Can't infer
        return self.values[1:]

    @staticmethod
    def _unpack_assignment_side(side):
        """Split the side left or right of the equal into bits."""
        if isinstance(side, str):
            return tuple(side.split())
        if isinstance(side, Sequence):
            return tuple(side)
        raise TypeError
