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

from viperleed.guilib.helpers import format_floats


_SEPARATORS = ',|'


class BeamIndex(tuple):
    """A 2-tuple of the Miller indices (as Fraction) of a LEED beam.

    Instances are cached for performance reasons.
    """

    _cache = {}

    def __new__(cls, *indices, denominator=1, from_numerators=False):
        """Return a new BeamIndex instance.

        Parameters
        ----------
        *indices : str, or Sequence of numbers
            indices of the beam. Can be passed as a single argument
            or as two arguments. When a single argument, it should be
            either a string of the form 'idx1, idx' or 'idx1 | idx2'
            (spaces don't count), or a 2-element iterable with indices.
        denominator : int, default=1
            This is used for speeding up instantiation when indices are
            given as numbers rather than strings, and it is not used at
            all for string inputs, nor for those indices that are given
            as Fraction. It should be the largest common denominator
            between the indices. This is mandatory when passing true
            floating-point indices.
        from_numerators : bool, default=False
            Use True when passing only the numerators as indices.
            `denominator` is used as the denominator. This argument
            is unused for those indices that are Fraction (even when
            generated from strings). If given, only the integer part
            of `indices` is used.

        Raises
        ------
        TypeError
            If `indices` --- or its first and only element --- is
            not a sequence.
        TypeError
            If `from_numerator` is True-thy but one of the indices
            is not int-able.
        ValueError
            If `indices` has a single string items, but it does not
            consist of two indices separated by one of the acceptable
            separators.
        ValueError
            If `indices` --- or its first and only element --- does
            not evaluate to exactly two items.
        ValueError
            If any of the indices is a floating-point number that
            cannot be expressed as <numerator> / `denominator`, with
            <numerator> integer.
        """
        indices = cls._process_indices(indices)

        # We will hash the object only if it does not contain '-1' as
        # this is a special value for hashing. In fact, hash(-1) == -2.
        # Hashing stuff that contains -1 would thus produce a large
        # number of collisions. Moreover hash(obj) never returns -1.
        for_hash = (*indices, denominator)
        input_hash = hash(for_hash) if -1 not in for_hash else -1
        if input_hash in cls._cache:
            return cls._cache[input_hash]

        values = (cls._index_to_fraction(i, denominator, from_numerators)
                  for i in indices)
        instance = super().__new__(cls, values)
        if input_hash != -1:
            cls._cache[input_hash] = instance
        return instance

    @property
    def numerators(self):
        """Return a tuple of the numerators in this BeamIndex."""
        return tuple(index.numerator for index in self)

    @classmethod
    def clear_cache(cls):
        """Fully clear cache. Use only if memory usage becomes an issue."""
        cls._cache = {}

    def __add__(self, value):
        """Disallow making BeamIndex longer."""
        cls = type(self)
        raise TypeError(f'unsupported operand type(s) for +: {cls.__name__!r} '
                        f'and {type(value).__name__!r}')

    def __format__(self, format_spec):
        """Return a formatted string version of this BeamIndex.

        Parameters
        ----------
        format_spec : str
            A standard format specifier in the form:

            [[fill]align][sign][#][0][minimumwidth][.precision][type]

            The following formats apply:
            - type == 's':
                returns '(num/den|num/den)' where the width of each
                of the h,k fields is dictated by the `minimumwidth`
                field of format_spec. `minimumwidth` should have
                form '(w_num,w_den)', so that the indices can be
                aligned at the slash. If this is omitted, the width
                of both fields is equal, and equal to the longest
                among the two.
                If h or k are integers, their '/den' is omitted, and
                replaced with white spaces if the other is fractional.
                Negative signs are printed, positive ones are replaced
                with spaces if the other index is negative.
            - type == 'f':
                returns f'({float(h)},{float(k)})'. Five decimal places
                are used as the default `precision`. If `minimumwidth`
                is given, it is treated as the minimum width of the
                integer part only. The two indices are aligned on the
                decimal point.
            - all others return f'({str(self)}:{format_spec}})'

        Return
        ------
        str

        Raises
        ------
        TypeError
            If `format_spec` is incorrectly formed.
        """
        if not format_spec:
            return str(self)
        if format_spec[-1] not in 'sf':  # Basic format
            return f'({str(self):{format_spec}})'
        if format_spec.endswith('f'):
            return f'({format_floats(format_spec, *self)})'

        # Case 's':
        num_min_len, den_min_len = self.get_format_widths()

        # Search the specifier for something like '(\d+,\d+)' to
        # be interpreted as the minimum widths of the numerator
        # and denominator fields
        found = re.search(r'(?P<minwidths>\((?P<num>\d+),(?P<den>\d+)\))',
                          format_spec)
        if found:
            num_min_len = max(num_min_len, int(found['num']))
            den_min_len = max(den_min_len, int(found['den']))
            format_spec = format_spec.replace(found['minwidths'], '')

        indices = '|'.join(self._fmt_one_index(i, num_min_len, den_min_len)
                           for i in self)
        return f'{indices:{format_spec}}'

    def __imul__(self, value):
        """Return a scaled version of this BeamIndex."""
        return self * value  # Delegate to __mul__

    def __mul__(self, value):
        """Return a scaled version of this BeamIndex.

        Parameters
        ----------
        value : int or Fraction
            The scaling factor to be applied identically to both of
            the Miller indices in this BeamIndex.

        Returns
        -------
        scaled : BeamIndex
            A new BeamIndex with Miller indices multiplied by `value`.

        Raises
        ------
        TypeError
            If value is not int or Fraction.
        """
        cls = type(self)
        if isinstance(value, (int, Fraction)):
            return cls(self[0]*value, self[1]*value)
        raise TypeError('unsupported operand type(s) for *: '
                        f'{cls.__name__!r} and {type(value).__name__!r}')

    def __repr__(self):
        """Return a representation string for this BeamIndex."""
        return f'BeamIndex({self})'  # Delegate to __str__

    def __rmul__(self, value):
        """Return a scaled version of this BeamIndex."""
        return self * value  # Delegate to __mul__

    def __str__(self):
        """Return a string version of this BeamIndex."""
        return ', '.join(str(index) for index in self)

    def get_format_widths(self):
        """Return the minimum number of characters to format this BeamIndex.

        Returns
        -------
        num_min_len, den_min_len : int
            The minimum widths of numerator and denominator, such that
            both h and k Miller indices are aligned at the slash (if
            fractional), or at the lowest significant digit of the
            numerator.
        """
        num_min_len = max(len(str(hk.numerator)) for hk in self)
        dens = [hk.denominator for hk in self]
        den_min_len = (0 if all(den == 1 for den in dens)
                       else max(len(str(den)) for den in dens))
        return num_min_len, den_min_len

    @staticmethod
    def _fmt_one_index(index, num_min_len, den_min_len):
        """Return a "num/den" string for index with defined field lengths."""
        # Numerator is right-justified in its field
        numerator = f'{index.numerator:>{num_min_len}}'

        # Denominator is a bit more complicated, as
        # it depends on whether it is == 1 or not
        if index.denominator != 1:
            return numerator + f'/{index.denominator:<{den_min_len}}'

        n_white = den_min_len
        n_white += 1 if den_min_len else 0  # slash if needed
        return numerator + ' ' * n_white

    @classmethod
    def _index_to_fraction(cls, index, denominator=1, from_numerator=False):
        """Return a Fraction version of index. See also __new__."""
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
                raise TypeError(f'{cls.__name__}: when using from_numerators, '
                                'the indices should be integers.') from err
        else:
            # The index passed is fractional, get the numerator
            float_numerator = index*denominator
            int_numerator = round(float_numerator)
            if abs(int_numerator - float_numerator) > 1e-6:
                raise ValueError(f'Fractional index {index} is not consistent '
                                 f'with denominator {denominator}.')
        return Fraction(int_numerator, denominator)

    @classmethod
    def _indices_from_string(cls, str_indices):
        """Return two Fractions from `str_indices`."""
        for separator in _SEPARATORS:
            indices = str_indices.split(separator)
            if len(indices) == 2:
                # found an acceptable separator
                return (Fraction(indices[0]), Fraction(indices[1]))
        raise ValueError(
            f'{cls.__name__}: too many/few indices in '
            f'{str_indices!r}, or incorrect separator (acceptable: '
            + 'or '.join(repr(s) for s in _SEPARATORS) + ').'
            )

    @classmethod
    def _process_indices(cls, indices):
        """Return a 2-items tuple from indices after checking them.

        Parameters
        ----------
        indices : Sequence
            One or two items. If a single item, it must be a string
            or a 2-item sequence. When a single string, it is expected
            to contain both Miller indices.

        Returns
        -------
        processed_indices : tuple
            Two items, corresponding to the two Miller indices in
            `indices`. The only processing of `indices` concerns
            the case of a single-string value.

        Raises
        ------
        TypeError
            If `indices` --- or its first and only element --- is not a
            sequence (i.e., it has no __len__).
        ValueError
            If `indices` has a single string items, but it does not
            consist of two indices separated by one of the acceptable
            separators.
        ValueError
            If `indices` --- or its first and only element --- does not
            evaluate to exactly two items.
        """
        n_indices = len(indices)
        if n_indices == 1:
            indices = indices[0]
            if isinstance(indices, str):
                indices = cls._indices_from_string(indices)
            try:
                n_indices = len(indices)
            except TypeError as err:
                raise TypeError(f'{cls.__name__}: when one argument given, '
                                'it should be a string or a 2-element '
                                'sequence.') from err
        if n_indices != 2:
            raise ValueError(f'{cls.__name__}: too many/few indices. '
                             'Exactly 2 indices should be given. '
                             f'Found {n_indices} instead.')
        return tuple(indices)
