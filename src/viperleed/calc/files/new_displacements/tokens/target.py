"""Module for the <target> token in the DISPLACEMENTS file."""

__authors__ = ("Alexander M. Imre (@amimre)",)
__copyright__ = "Copyright (c) 2019-2025 ViPErLEED developers"
__created__ = "2024-10-14"
__license__ = "GPLv3+"

import re

from .base import DisplacementsFileToken, TokenParserError


class TargetingError(TokenParserError):
    """Base class for errors in the targeting module."""


class TargetToken(DisplacementsFileToken):
    """Class to handle the <target> token in the DISPLACEMENTS file."""

    def __init__(self, target_str):
        self.target_str = target_str.strip()
        self.nums = None
        self.layers = None
        self._parse_target()

    def _parse_target(self):
        """Parse the site, and optional nums or layers from the target string."""
        parts = self.target_str.split(maxsplit=1)
        if not parts:
            raise TargetingError('Target string is empty')
        site_str = parts[0]

        # do not allow numeric or layer-only targets
        if site_str[0].isdigit() or 'L(' in site_str or '[' in site_str:
            msg = (
                f'Target must start with a non-numeric label, got: "{site_str}"'
            )
            raise TargetingError(msg)

        self.regex = _generate_label_match_regex(site_str)

        if len(parts) == 1:
            # only site is specified, no nums or layers
            return

        # Check if we have a layer specification
        layer_match = re.match(r'L\((\d+)(-(\d+))?\)', parts[1])
        if layer_match:
            start_layer = int(layer_match.group(1))
            end_layer = (
                int(layer_match.group(3))
                if layer_match.group(3)
                else start_layer
            )
            self.layers = list(range(start_layer, end_layer + 1))
            # make sure there is nothing else after the layer
            if len(parts[1].split()) > 1:
                msg = f'Invalid target layer specification: "{parts[1]}".'
                raise TargetingError(msg)
            return

        # Check for a range like "1-4"
        range_match = re.match(r'(\d+)-(\d+)', parts[1])
        if range_match:
            start_num = int(range_match.group(1))
            end_num = int(range_match.group(2))
            self.nums = list(range(start_num, end_num + 1))
            return

        # Check for a list of numbers
        try:
            self.nums = list(map(int, parts[1].split()))
        except ValueError as err:
            msg = (
                f'Invalid target specification: "{parts[1]}".'
            )
            raise TargetingError(msg) from err



    def __str__(self):
        """Return a string representation of the TargetToken object."""
        return f'TargetToken(target_str={self.target_str})'

    def __eq__(self, other):
        """Compare self to other TargetToken object.

        Note this check may not be reliable. Comparison of TargetToken objects
        should be made using a selection of Atoms.
        """
        if not isinstance(other, TargetToken):
            return False
        return (
            (other.nums == self.nums if self.nums is not None else True)
            and (other.layers == self.layers if self.layers is not None else True)
            and (other.regex == self.regex)
        )


def _generate_label_match_regex(label):
    """Generate a regex pattern to match variations of the given label.

    The label can contain wildcards, where '*' matches any number of characters,
    including none.
    """
    # Escape any special characters in the label, except for '*'
    escaped_label = re.escape(label).replace(r'\*', r'\w*')

    # Append `\w*` at the end to match strings starting with the pattern
    pattern = rf'^{escaped_label}\w*'

    # Compile the final regex pattern
    return re.compile(pattern)
