"""Module target of viperleed.calc.files.new_displacements.tokens."""

__authors__ = ('Alexandra Mia Imre (@alexmiame)',)
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2024-10-14'
__license__ = 'GPLv3+'

import re

from viperleed.calc.lib.string_utils import read_int_range

from .base import DisplacementsFileToken, TokenParserError

LAYER_REGEX = r'L\s*\(\s*(?P<ranges>[\d\-:\s]+)\s*\)'
MALFORMED_LAYER_CHARS = ('L', '(', ')')


class TargetingError(TokenParserError):
    """Base class for errors in the targeting module."""


class TargetToken(DisplacementsFileToken):
    """Class to handle the <target> token in the DISPLACEMENTS file.

    A target is a site label (with optional '*' wildcards), followed
    by optional atom numbers and L(...) layer specifications, in any
    order (e.g., 'O 1 3 L(2)'). The target selects the atoms matching
    the label that have one of the numbers or are in one of the
    layers (their union). Without numbers and layers, the label alone
    selects the atoms.

    Attributes
    ----------
    nums : tuple of int or None
        The atom numbers given, None if there are none.
    layers : tuple of int or None
        The layers given (1-based) in all L(...), None if there are
        none.
    """

    def __init__(self, target_str):
        self.target_str = target_str.strip()
        self.nums = None
        self.layers = None
        self._regex = None
        self._parse_target()

    @property
    def regex(self):
        """Return the compiled regex pattern for matching site labels."""
        return self._regex

    def _parse_target(self):
        """Parse the site label and optional nums and layers."""
        parts = self.target_str.split(maxsplit=1)
        if not parts:
            raise TargetingError('Target string is empty')
        site_str = parts[0]

        # do not allow numeric or layer-only targets
        if site_str[0].isdigit() or 'L(' in site_str or '[' in site_str:
            msg = f'Target must start with a non-numeric label, got: "{site_str}"'
            raise TargetingError(msg)

        self._regex = _generate_label_match_regex(site_str)

        if len(parts) == 1:
            # only site is specified, no nums or layers
            return

        self._parse_nums_and_layers(parts[1])

    def _parse_nums_and_layers(self, specs):
        """Parse atom numbers and L(...) layers following the label.

        Numbers and layers may come in any order, with any number of
        L(...) groups. The target is their union: the atoms with the
        given numbers, plus those in the given layers.
        """
        layer_specs = [m['ranges'] for m in re.finditer(LAYER_REGEX, specs)]
        nums_spec = re.sub(LAYER_REGEX, ' ', specs).strip()
        if any(c in nums_spec for c in MALFORMED_LAYER_CHARS):
            msg = f'Invalid layer specification: "{specs}".'
            raise TargetingError(msg)
        if layer_specs:
            self.layers = _read_ints(' '.join(layer_specs), specs)
        if nums_spec:
            self.nums = _read_ints(nums_spec, specs)

    def __str__(self):
        """Return a string representation of the TargetToken object."""
        return f'TargetToken(target_str={self.target_str})'

    def __eq__(self, other):
        """Compare self to other TargetToken object.

        Note this check may not be reliable. Comparison of TargetToken objects
        should be made using a selection of Atoms.
        """
        if not isinstance(other, TargetToken):
            return NotImplemented
        return (
            other.nums == self.nums
            and self.layers == other.layers
            and self.regex == other.regex
        )


def _read_ints(ranges, specs):
    """Return a tuple of the integers in `ranges`, part of `specs`."""
    try:
        return tuple(read_int_range(ranges))
    except ValueError as err:
        msg = f'Invalid target specification: "{specs}".'
        raise TargetingError(msg) from err


def _generate_label_match_regex(label):
    """Generate a regex pattern to match variations of the given label.

    The label can contain wildcards, where '*' matches any number of characters,
    including none.
    """
    # Escape any special characters in the label, except for '*'
    escaped_label = re.escape(label).replace(r'\*', r'\w*')

    pattern = rf'^{escaped_label}'

    # Compile the final regex pattern
    return re.compile(pattern)
