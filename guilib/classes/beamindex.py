"""Module beamindex of viperleed.guilib.classes.

======================================
  ViPErLEED Graphical User Interface
======================================

Created: 2024-03-22
Author: Michele Riva (@michele-riva)

Defines the BeamIndex class, representing Miller indices of a LEED beam.
Used to be part of the guilib.base module.
"""

import re

from quicktions import Fraction  # Faster than standard Fraction

from viperleed.guilib.base import format_floats
from viperleed.guilib.base import integer_part_length


class BeamIndex(tuple):
    """A 2-tuple of the Miller indices of a LEED beam.

    Each index is a Fraction.

    Instances are cached for performance reasons.
    """

    separators = ',|'
    __cache = {}

    def __new__(cls, *indices, denominator=1, from_numerators=False):
        """
        Parameters
        ----------
        indices : str, iterable of str, or iterable of numbers
            indices of the beam. Can be passed as a single argument or as two
            arguments.
            When a single argument, it should be either a string of the form
            'idx1, idx' or 'idx1 | idx2' (spaces don't count), or a 2-element
            iterable with indices.
        denominator : int (default=1)
            this is used for speeding up instantiation when indices are given
            as numbers rather than strings, and it is not used at all for string
            inputs, nor for those indices that are given as Fraction.
            It should be the largest common denominator between the
            indices. This is mandatory when passing true floating-point indices.
        from_numerators : bool (default=False)
            Use True when passing only the numerators as indices. The
            denominator is taken from the denominator optional parameter. When
            passing numerators only, it is most efficient to give the indices
            as ints rather than floats
        """
        indices = cls.__process_indices(indices)

        # We will hash the object only if it does not contain
        # '-1' as this is a special value for hashing. In fact,
        # hash(-1) = -2 Hashing stuff that contains -1 would
        # thus produce a significant number of collisions.
        # Moreover hash(obj) never returns -1.
        input_hash = -1
        for_hash = (*indices, denominator)
        if -1 not in for_hash:
            input_hash = hash(for_hash)
        if input_hash in cls.__cache:
            return cls.__cache[input_hash]

        instance = (cls.__index_to_fraction(
                        indices[0], denominator=denominator,
                        from_numerator=from_numerators
                        ),
                    cls.__index_to_fraction(
                        indices[1], denominator=denominator,
                        from_numerator=from_numerators
                        ))
        instance = super().__new__(cls, instance)
        if input_hash != -1:
            cls.__cache[input_hash] = instance
        return instance

    @staticmethod
    def clear_cache():
        """Fully clear cache. Use only if memory usage becomes an issue."""
        BeamIndex.__cache = {}

    @staticmethod
    def __process_indices(indices):
        """
        Checks and processes the indices as needed

        Returns
        -------
        indices : tuple, 2 elements

        Raises
        ------
        TypeError
        ValueError
        """
        n_indices = len(indices)

        if n_indices == 1:
            indices = indices[0]
            if isinstance(indices, str):
                indices = BeamIndex.__indices_from_string(indices)
            try:
                n_indices = len(indices)
            except TypeError as err:
                raise TypeError('BeamIndex: when one argument given, '
                                'it should be a string or a 2-element '
                                'array-like.') from err

        if n_indices != 2:
            raise ValueError('BeamIndex: too many/few indices. '
                             'Exactly 2 indices should be given. '
                             f'Found {n_indices} instead.')
        return tuple(indices)

    @staticmethod
    def __indices_from_string(str_indices):
        """
        """
        for separator in BeamIndex.separators:
            indices = str_indices.split(separator)
            if len(indices) == 2:
                # found an acceptable separator
                return (Fraction(indices[0]), Fraction(indices[1]))
        raise ValueError('BeamIndex: too many/few indices in '
                         f'{str_indices!r}, or incorrect '
                         "separator (acceptable: ',', '|').")

    @staticmethod
    def __index_to_fraction(index, denominator=1, from_numerator=False):
        # if isinstance(index, Fraction):
            # return index
        try:
            index.limit_denominator  # attribute of Fraction only
        except AttributeError:
            pass
        else:
            return index

        if from_numerator:
            try:
                int_numerator = int(index)
            except TypeError as err:
                raise TypeError('BeamIndex: when using from_numerators, '
                                'the indices should be integers.') from err
        else:
            # The index passed is fractional, get the numerator
            float_numerator = index*denominator
            int_numerator = round(float_numerator)
            if abs(int_numerator - float_numerator) > 1e-6:
                raise ValueError(f'Fractional index {index} is not consistent '
                                 f'with denominator {denominator}.')
        return Fraction(int_numerator, denominator)

    def __str__(self):
        return ', '.join(str(index) for index in self)

    def __repr__(self):
        return f'BeamIndex({self})'  # Delegate to __str__

    def __mul__(self, factor):
        """Override tuple.__mul__().

        Make it such that self*number = (self[0]*number, self[1]*number)
        """
        cls = type(self)
        if isinstance(factor, (int, Fraction)):
            return cls(self[0]*factor, self[1]*factor)
        raise TypeError('unsupported operand type(s) for *: '
                        f'{cls.__name__!r} and {type(factor).__name__!r}')

    def __rmul__(self, factor):
        return self.__mul__(factor)

    def __imul__(self, factor):
        return self * factor

    def __format__(self, format_spec):
        """
        Customized formatting of BeamIndex. format_spec is a standard format
        specifier in the form:

        [[fill]align][sign][#][0][minimumwidth][.precision][type]

        The following formats apply:
        - type == 's':
            returns '(num/den|num/den)' where the width of each of the h,k
            fields is dictated by the minimumwidth specifier in format_spec. In
            this case, minimumwidth should be of the form '(w_num,w_den)', so
            that the indices can be aligned at the slash. If this is omitted,
            the width of both fields is equal, and equal to the longest among
            the two.
            If h or k are integers, their '/den' is omitted, and replaced with
            white spaces if the other is fractional.
            Negative signs are printed, positive ones are replaced with spaces
            if the other index is negative
        - type == 'f':
            returns f'({float(h)},{float(k)})'. If .precision is not given,
            uses 5 digits. If minimumwidth is given, it is treated as the
            minimum width of the integer part only. The two indices are aligned
            on the decimal point.
        - all others return f'({str(self)}:{format_spec}})'
        """
        if not format_spec:
            return str(self)
        if format_spec[-1] not in 'sf':
            # basic format
            return f'({str(self):{format_spec}})'
        if format_spec.endswith('f'):
            return f'({format_floats(format_spec, *self)})'

        # Case 's':
        num_min_len, den_min_len = self.get_format_lengths('s')

        # now search the specifier for something like '(\d+,\d+)' to be
        # interpreted as the minimum widths of the two fields, to update
        # the minimum lengths of the fields
        m = re.search(r'((?P<num>\d+),(?P<den>\d+))', format_spec)
        if m is not None:
            num_min_len = max(num_min_len, int(m['num']))
            den_min_len = max(den_min_len, int(m['den']))
            format_spec = format_spec.replace(f'({m["num"]},{m["den"]})', '')

        raws = ['', '']
        for i, hk in enumerate(self):
            # numerator is right-justified in its field
            raws[i] = f'{hk.numerator:>{num_min_len}}'

            # denominator a bit more complicated, as it depends on whether
            # it is == 1 or not
            if hk.denominator == 1:
                n_white = den_min_len
                n_white += 1 if den_min_len else 0  # slash if needed
                raws[i] += ' ' * n_white
            else:
                raws[i] += f'/{hk.denominator:<{den_min_len}}'
        return f"{f'({raws[0]}|{raws[1]})':{format_spec}}"

    @property
    def numerators(self):
        return tuple(index.numerator for index in self)

    def get_format_lengths(self, str_or_float):
        """
        Returns the minimum number of characters that can represent the
        BeamIndex as a fraction or as a float. The two cases behave differently:
        - str_or_float = 's':
            returns a 2-tuple with the minimum lengths of numerator and
            denominator, such that both h and k are aligned at the slash (if
            fractional) or at the lowest significant digit of the numerator
        - str_or_float = 'f':
            returns a 1-tuple with the minimum length of the string
            representation of the integer parts that can represent both
        """
        if str_or_float not in ('s', 'f'):
            raise ValueError("Invalid format specifier. Should be 's' or 'f'.")
        if str_or_float == 'f':
            return (integer_part_length(*self),)
        # case 's'
        num_min_len = max(len(str(hk.numerator)) for hk in self)
        dens = [hk.denominator for hk in self]
        if all(den == 1 for den in dens):
            den_min_len = 0
        else:
            den_min_len = max(len(str(den)) for den in dens)
        return num_min_len, den_min_len
