"""Module targeting."""

__authors__ = ('Alexander M. Imre (@amimre)',)
__created__ = '2024-10-14'

import re

import numpy as np

from .base import DisplacementsFileToken, TokenParserError

class TargetingError(TokenParserError):
    """Base class for errors in the targeting module."""


class Targets:
    """Class to handle multiple <target> tokens."""

    def __init__(self, target_str):
        self.target_str = target_str
        self.subtargets = []
        self._parse_target(target_str)

    def _parse_target(self, target_str):
        """Parse multiple subtargets separated by commas."""
        subtarget_strs = target_str.split(',')
        self.subtargets = [TargetToken(sub.strip()) for sub in subtarget_strs]

    def select(self, atom_basis):
        """Take the 'or' of all subtargets, combining masks."""
        combined_mask = np.full(len(atom_basis), False)
        for subtarget in self.subtargets:
            combined_mask = combined_mask | subtarget.select(atom_basis)
        return combined_mask

    def __repr__(self):
        """Return the string representation of the target."""
        return f'Target({self.target_str})'


class TargetToken(DisplacementsFileToken):
    """Class to handle the <target> token in the DISPLACEMENTS file."""

    def __init__(self, target_str):
        self.target_str = target_str
        self.nums = None
        self.layers = None
        self._parse_target()

    def _parse_target(self):
        """Parse the site, and optional nums or layers from the target string."""
        parts = self.target_str.split()
        if not parts:
            raise TargetingError('Target string is empty')
        site_str = parts[0]
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
        else:
            # Check for a range like "1-4"
            range_match = re.match(r'(\d+)-(\d+)', parts[1])
            if range_match:
                start_num = int(range_match.group(1))
                end_num = int(range_match.group(2))
                self.nums = list(range(start_num, end_num + 1))
            else:
                # It's a list of numbers
                self.nums = list(map(int, parts[1:]))



    def __repr__(self):
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
            and (other.nums == self.layers if self.layers is not None else True)
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
