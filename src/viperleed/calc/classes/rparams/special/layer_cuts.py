"""Module layer_cuts of viperleed.calc.classes.rparams.special.

Defines the LayerCutTokenType, LayerCutToken, and LayerCuts classes.
They are convenience classes for handling user input of slab cut
positions used to generate layers of a Slab. The code is a rewrite
of the functionality, originally by @fkraushofer, previously split
partly in ParameterInterpreter, partly in Slab.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-10-21'
__license__ = 'GPLv3+'

from collections.abc import Sequence
from enum import IntEnum, auto
import itertools
from numbers import Real
import re
from typing import Any

from viperleed.calc.lib.base import pairwise
from viperleed.calc.lib.dataclass_utils import frozen
from viperleed.calc.lib.dataclass_utils import set_frozen_attr

from .base import SpecialParameter


def threewise(iterable):
    """Yield triplets of items from iterable."""
    orig, shift_one, shift_two = itertools.tee(iterable, 3)
    next(shift_one, None)
    next(shift_two, None)
    next(shift_two, None)
    yield from zip(orig, shift_one, shift_two)


class LayerCutTokenType(IntEnum):
    """Possible types of layer cut tokens."""

    NUMERIC = auto()
    ORDERING = auto()  # < or >
    AUTO_DC = auto()
    AUTO_DZ = auto()
    INVALID = auto()


_AUTO_CUT_RE = re.compile(r'(?P<type>dz|dc)\((?P<cutoff>[\d.]+)\)')


@frozen
class LayerCutToken:
    """A single token of a LayerCuts object."""

    type_: LayerCutTokenType
    value: Any

    # Two attributes that can only be set once, and only for auto-cut
    lower: 'LayerCutToken' = None
    upper: 'LayerCutToken' = None

    def __eq__(self, other):
        """Return whether this LayerCutToken is equal to other."""
        if isinstance(other, str):
            other = LayerCutToken.from_string(other)
        elif isinstance(other, Real):
            other = LayerCutToken.make_numeric(other)
        if not isinstance(other, LayerCutToken):
            return NotImplemented
        equal = (self.type_ is other.type_ and self.value == other.value
                 and self.lower == other.lower and self.upper == other.upper)
        return equal or NotImplemented

    def __format__(self, format_spec):
        """Return a formatted version of this LayerCutToken.

        Parameters
        ----------
        format_spec : str
            The formatting specification for this token. Currently
            format_spec is used only if token.is_numeric. It is
            then passed unchanged to self.value.__format__, i.e.,
            to float.__format__.

        Returns
        -------
        formatted : str
            Formatted version of this token.
        """
        if not self.is_numeric:
            return str(self)
        return format(self.value, format_spec)

    def __str__(self):
        """Return a string version of this token."""
        if self.is_ordering:
            return self.value
        if self.is_numeric:
            return str(self.value)
        if self.is_auto_cut:
            dc_dz = self.type_.name.replace('AUTO_', '').lower()
            return f'{dc_dz}({self.value})'
        return 'INVALID'

    def __post_init__(self):
        """Check that the initialization parameters are acceptable."""
        if ((self.is_numeric or self.is_auto_cut)
                and not isinstance(self.value, Real)):
            raise TypeError('Numeric tokens must have numeric values')
        if self.is_numeric and not 0 <= self.value <= 1:
            raise ValueError(f'Numeric token value {self.value} '
                             'is not in range [0, 1]')
        if self.is_ordering and self.value not in '<>':
            raise ValueError('Ordering tokens can only have < or > as value')
        if not self.is_auto_cut and (self.lower, self.upper) != (None, None):
            raise ValueError('Only dc/dz can have bounds')
        if isinstance(self.lower, Real):
            set_frozen_attr(self, 'lower', self.make_numeric(self.lower))
        if isinstance(self.upper, Real):
            set_frozen_attr(self, 'upper', self.make_numeric(self.upper))
        if self.is_auto_cut:
            self._check_valid_bound_type(self.lower, self.upper)

    def __repr__(self):
        """Return a representation string for this LayerCutToken."""
        repr_ = f'LayerCutToken({self.type_.name}, {self.value}'
        if self.is_auto_cut and (self.lower, self.upper) != (None, None):
            repr_ += f', {self.lower}, {self.upper}'
        return repr_ + ')'

    @property
    def is_auto_cut(self):
        """Return whether this token is of type 'dc'/'dz'."""
        return self.type_ in (LayerCutTokenType.AUTO_DC,
                              LayerCutTokenType.AUTO_DZ)

    @property
    def is_auto_cut_dc(self):
        """Return whether this token is of type 'dc'."""
        return self.type_ is LayerCutTokenType.AUTO_DC

    @property
    def is_ordering(self):
        """Return whether this token is '>' or '<'."""
        return self.type_ is LayerCutTokenType.ORDERING

    @property
    def is_larger_than(self):
        """Return whether this token is a '>'."""
        return self.is_ordering and self.value == '>'

    @property
    def is_numeric(self):
        """Return whether this token is a number."""
        return self.type_ is LayerCutTokenType.NUMERIC

    @classmethod
    def from_string(cls, string):
        """Return a LayerCutToken from a `string`."""
        string = re.sub(r'\s', '', string).lower()
        _invalid = cls(LayerCutTokenType.INVALID, string)
        if len(string) == 1 and string in '<>':
            return cls(LayerCutTokenType.ORDERING, string)
        try:
            value = float(string)
        except ValueError:
            pass
        else:
            try:
                return cls(LayerCutTokenType.NUMERIC, value)
            except ValueError:
                return _invalid
        auto_match = _AUTO_CUT_RE.match(string)
        if not auto_match:
            return _invalid
        try:
            value = float(auto_match['cutoff'])
        except ValueError:
            return _invalid
        type_name = f'AUTO_{auto_match["type"].upper()}'
        return cls(LayerCutTokenType[type_name], value)

    @classmethod
    def make_numeric(cls, value):
        """Return a LayerCutToken of type NUMERIC from a numeric `value`."""
        return cls(LayerCutTokenType.NUMERIC, float(value))

    def get_bounds(self, other_lower_bounds=()):
        """Return the numeric value of lower and upper bounds.

        This method assumes that this token is used internally in
        a LayerCuts object. This means that the scope of this token
        is bound only with '<', and not with '>'.

        Parameters
        ----------
        other_lower_bounds : Sequence, optional
            Other bounds to be considered for the lower side.
            Items are numbers. Default is no extra bounds.

        Raises
        ------
        AttributeError
            If this method is called on a non-auto-cut LayerCutToken
        RuntimeError
            If this token has its lower bound larger than its higher
            bound.
        ValueError
            If `other_lower_bounds` contains values larger than the
            upper bound of this token.
        """
        if not self.is_auto_cut:
            raise AttributeError(f'{self.type_.name} has no bounds')

        self_lower = self.lower.value if self.lower is not None else 0.0
        self_upper = self.upper.value if self.upper is not None else 1.0
        if self_lower >= self_upper:
            raise RuntimeError(f'{self} has swaped upper and lower bounds')
        if any(b > self_upper for b in other_lower_bounds):
            raise ValueError(
                'All other_lower_bounds must be at most as large '
                f'as the upper bound ({self_upper}) of this token'
                )
        return max((self_lower, *other_lower_bounds)), self_upper

    def with_bounds(self, lower, upper):
        """Return a copy of self with `lower` and `upper` bounds."""
        if not self.is_auto_cut:
            raise TypeError('Only auto-cut can have bounds')
        if any(attr is not None for attr in (self.lower, self.upper)):
            raise RuntimeError(f'{self} already has bounds')
        if isinstance(lower, Real):
            lower = self.make_numeric(lower)
        if isinstance(upper, Real):
            upper = self.make_numeric(upper)
        self._check_valid_bound_type(lower, upper)
        cls = type(self)
        return cls(self.type_, self.value, lower, upper)

    @staticmethod
    def _check_valid_bound_type(*bounds):
        """Ensure bounds are NUMERIC LayerCutToken or None."""
        try:
            all_numbers = all((b is None or b.is_numeric) for b in bounds)
        except AttributeError:  # Not a LayerCutToken
            raise TypeError('Bounds must be LayerCutToken objects') from None
        if not all_numbers:
            raise TypeError('Only LayerCutTokenType.NUMERIC-type '
                            'tokens can be used as bounds')


class LayerCuts(SpecialParameter, param='LAYER_CUTS'):
    """A container of LayerCutToken objects."""

    def __init__(self, *tokens):
        """Initialize instance from some LayerCutToken objects."""
        self._tokens = None
        self._tokens_with_ordering = None   # Also, sorted as original
        self._had_larger_than = False  # Internally we change all to <

        self._check_tokens(tokens)
        ori_tokens, tokens = self._clean_up_tokens(tokens)
        self._tokens = tokens
        self._tokens_with_ordering = ori_tokens

    def __bool__(self):
        """Return whether there is any cut in this LayerCuts."""
        return bool(self._tokens)

    def __format__(self, format_spec):
        """Return a formatted version of this LayerCuts object.

        Parameters
        ----------
        format_spec : str
            Format specification. This is passed unchanged to
            each of the tokens in this LayerCuts. See also
            help(LayerCutToken.__format__).

        Returns
        -------
        formatted : str
            Formatted version of this LayerCuts.
        """
        return ' '.join(format(t, format_spec)
                        for t in self._tokens_with_ordering)

    def __iter__(self):
        """Return an iterator of tokens in this LayerCuts object."""
        return iter(self._tokens)

    def __len__(self):
        """Return the number of cut specifications (excluding '<'/'>')."""
        return len(self._tokens)

    def __str__(self):
        """Return a string version of this LayerCuts object."""
        return ' '.join(str(t) for t in self._tokens_with_ordering)

    def __repr__(self):
        """Return a string representation of this LayerCuts."""
        return f'LayerCuts({str(self)!r})'

    @classmethod
    def as_layer_cuts(cls, other):
        """Return a LayerCuts from other, raise otherwise."""
        if isinstance(other, LayerCuts):
            return other
        if isinstance(other, str):
            return cls.from_string(other)     # May raise ValueError
        if isinstance(other, Sequence):
            cuts = cls()
            cuts.update_from_sequence(other)  # May Value/TypeError
            return cuts
        raise TypeError('Cannot produce a LayerCuts from '
                        f'{type(other).__name__!r}. Only '
                        'str and Sequence allowed')

    # Override parent's method
    from_value = as_layer_cuts

    @classmethod
    def from_string(cls, string):
        """Return LayerCuts from a `string`."""
        if all(c in string for c in '<>'):
            raise ValueError('Cannot parse list with both "<" and ">"')
        tokens = cls._tokenize_string(string)
        instance = cls(*tokens)
        return instance

    def update_from_sequence(self, sequence):
        """Clear this LayerCuts, and fill it with values from `sequence`.

        Parameters
        ----------
        sequence : Sequence
            Items must be LayerCutToken objects, real-valued numbers or
            tokenized strings. Each item is used to create the correct
            LayerCutToken.
        """
        new_tokens = []
        for item in sequence:
            if isinstance(item, LayerCutToken):
                token = item
            elif isinstance(item, Real):
                token = LayerCutToken.make_numeric(item)
            elif isinstance(item, str):
                token = LayerCutToken.from_string(item)
            else:
                raise TypeError('Expected Real or str, '
                                f'found {type(item).__name__}')
            new_tokens.append(token)
        self._check_tokens(new_tokens)
        ori_tokens, tokens = self._clean_up_tokens(new_tokens)
        self._tokens = tokens
        self._tokens_with_ordering = ori_tokens

    @staticmethod
    def _check_consistency_of_auto_cuts(tokens):
        """Raise if auto-cut tokens are inconsistent.

        Parameters
        ----------
        tokens : Sequence
            The tokens to be checked.

        Raises
        ------
        ValueError
            If auto-cut tokens do not have an ORDERING operator next
            to them (except when there's only one auto-cut item).
        """
        has_auto_cut = any(t.is_auto_cut for t in tokens)
        if not has_auto_cut or len(tokens) == 1:
            return
        err = 'dc/dz always requires < or > unless it is the only cut'
        for pair in pairwise(tokens):
            if not any(tok.is_auto_cut for tok in pair):
                continue
            left, right = pair
            if left.is_auto_cut and right.is_ordering:
                continue
            if left.is_ordering and right.is_auto_cut:
                continue
            raise ValueError(err)

    @staticmethod
    def _check_consistency_of_bounds(tokens):
        """Complain if any of the bounds of `tokens` is inconsistent.

        Parameters
        ----------
        tokens : Sequence
            The LayerCutToken objects to be checked. They are assumed
            to be already sorted so that only '<' bounds are used.

        Raises
        ------
        ValueError
            If any of the auto-cuts has its lower bound larger than
            its upper one
        ValueError
            If any fractional cuts falls within an 'auto-cut' range
        """
        auto_cuts = [t for t in  tokens if t.is_auto_cut]
        if not auto_cuts:
            return

        # First check that all the dc/dz have lower <= upper. While
        # going through, also collect the bounds for later checks.
        bounds = []
        for token in auto_cuts:
            try:
                bounds.append(token.get_bounds())
            except RuntimeError:
                raise ValueError(f'{token} has swapped lower '
                                 'and upper bounds') from None

        # Now check that none of the numeric cuts falls within
        # any of the auto-cut bounds
        for token in tokens:
            if not token.is_numeric:
                continue
            try:
                invalid = next(i for i, (lower, upper) in enumerate(bounds)
                               if lower < token.value < upper)
            except StopIteration:
                continue
            raise ValueError(f'Numeric cut {token.value} falls within '
                             f'the bounds of {auto_cuts[invalid]}')

    @staticmethod
    def _check_consistency_of_order_operators(tokens):
        """Raise if ordering operators are inconsistent.

        This method filters out weird input like the following:
        - '< dc [...]' and '[...] dc <'
        - 'N N < < dc'
        - 'N < N' and 'dc < dc'

        Parameters
        ----------
        tokens : Sequence
            Elements are LayerCutToken objects.

        Raises
        ------
        ValueError
            If an ORDERING token appears as the first or last
            item in tokens.
        ValueError
            If an ORDERING token does not separate a NUMERIC
            from an AUTO token.
        """
        has_ordering = any(t.is_ordering for t in tokens)
        if not has_ordering:
            return

        # Complain if '<' is first or last. Notice that this also
        # implies that, if this test passes, there are at least
        # three items.
        if any(tokens[i].is_ordering for i in (0, -1)):
            raise ValueError('Cannot have < or > as first or last items')

        # Check that left and right of a '<' or '>' there's
        # a number on one side and a dc/dz on the other
        for left, middle, right in threewise(tokens):
            if not middle.is_ordering:
                continue
            if left.is_numeric and right.is_auto_cut:
                continue
            if right.is_numeric and left.is_auto_cut:
                continue
            invalid = f'{left.value} {middle.value} {right.value}'
            raise ValueError(f'Invalid token triplet {invalid!r}. > / < '
                             'must separate a number and a dc / dz')

    @classmethod
    def _check_tokens(cls, tokens):
        """Run some basic consistency checks on `tokens`."""
        if any(not isinstance(t, LayerCutToken) for t in tokens):
            raise TypeError('Can only contain LayerCutToken items')

        invalid = [repr(t.value) for t in tokens
                   if t.type_ is LayerCutTokenType.INVALID]
        if invalid:
            raise ValueError('Could not parse ' + ', '.join(invalid))
        # Make sure that ordering operators are OK, i.e., they do not
        # appear first and last, and they separate a number and dc/dz
        cls._check_consistency_of_order_operators(tokens)
        cls._check_consistency_of_auto_cuts(tokens)

    def _clean_up_tokens(self, tokens):
        """Return tokens and the filtered copy used internally."""
        # We assume that tokens have already gone through basic
        # checks, i.e., they passed _check_tokens. Store a copy
        # to return before messing with tokens.
        ori_tokens = list(tokens)

        # Turn the tokens around if there is a '>'. This allows
        # us not to store the ordering operators at all (except
        # for the purpose of returning a string).
        self._had_larger_than = any(t.is_larger_than for t in tokens)
        if self._had_larger_than:
            tokens = reversed(tokens)

        # Now get rid of all the '<'/'>'
        tokens = [t for t in tokens if not t.is_ordering]

        # Assign 'bounds' to the dc/dz items. This also filters
        # out weird input like 'dz dc < N'. It is important to
        # realize that tokens is modified in-place and its items
        # are replaced by new ones. This makes the ori_tokens
        # above and the ones in tokens different. Should we ever
        # need to do anything other than str on the ori_tokens,
        # we need to deal with this appropriately here.
        self._set_bounds_of_auto_tokens(tokens)

        # And check (i) lower <= upper, (ii) no numbers fall in bounds
        self._check_consistency_of_bounds(tokens)
        return ori_tokens, tokens

    @staticmethod
    def _set_bounds_of_auto_tokens(tokens):
        """Store bounds information in the dc/dz items of `tokens`."""
        if len(tokens) <= 1:
            return
        shift_left, shift_right = (None, *tokens[:-1]), (*tokens[1:], None)
        triplets = zip(shift_left, tokens, shift_right)
        for i, (lower, middle, upper) in enumerate(triplets):
            if not middle.is_auto_cut:
                continue
            tokens[i] = middle.with_bounds(lower, upper)

    @staticmethod
    def _tokenize_string(string):
        """Return a list of LayerCutToken from a `string`."""
        # Make sure ">"/"<" are split off from the rest
        string = re.sub(r'([<>])', r' \1 ', string.strip()).lower()
        return [LayerCutToken.from_string(s) for s in string.split()]
