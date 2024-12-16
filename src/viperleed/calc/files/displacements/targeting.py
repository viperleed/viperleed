"""Module targeting."""

__authors__ = ('Alexander M. Imre (@amimre)',)
__created__ = '2024-10-14'

import re

import numpy as np

class TargetingError(ValueError):
    """Base class for errors in the targeting module."""


def generate_label_match_regex(label):
    """Generate a regex pattern to match variations of the given label,
    with '*' acting as a wildcard for word characters, and matching prefixes.
    """
    # Escape any special characters in the label, except for '*'
    escaped_label = re.escape(label).replace(r'\*', r'\w*')

    # Append `\w*` at the end to match strings starting with the pattern
    pattern = rf'^{escaped_label}\w*'

    # Compile the final regex pattern
    return re.compile(pattern)


class BSSubtarget:
    def __init__(self, target_str):
        self.target_str = target_str
        self.nums = None
        self.layers = None
        self._parse_target()

    def _parse_target(self):
        """Parse the site, and optional nums or layers from the target string."""
        parts = self.target_str.split()
        if not parts:
            raise ValueError('Subtarget string is empty')
        site_str = parts[0]
        self.regex = generate_label_match_regex(site_str)

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

    def select(self, atom_basis):
        """Select base scatterers that match the subtarget specification."""
        mask = np.full(len(atom_basis), True)

        # mask based on the site
        matches = np.array(
            [self.regex.match(bs.site) is not None for bs in atom_basis]
        )
        mask = mask & matches

        # mask based on the labels
        label_mask = mask.copy()

        # If nums are specified, apply the selection based on nums
        if self.nums is not None:
            # check range for nums
            if any(num < 1 or num > len(atom_basis) for num in self.nums):
                msg = f'Invalid atom number for subtarget: {self.target_str}'
                raise TargetingError(msg)
            num_mask = np.array([bs.num in self.nums for bs in atom_basis])
            # check if any of the given nums have the wrong label
            wrong_label = np.logical_and(num_mask, ~label_mask)
            if np.any(wrong_label):
                msg = (
                    'Atom numbers do not match label for subtarget: '
                    f'{self.target_str}'
                )
                raise TargetingError(msg)
            mask = mask & num_mask

        # If layers are specified, apply the selection based on layers
        if self.layers is not None:
            mask = mask & np.array(
                # TODO: layer counting from 1; can we unify this somewhere?
                [bs.layer+1 in self.layers for bs in atom_basis]
            )

        if mask.sum() == 0:
            msg = f'No atoms selected for subtarget: {self.target_str}'
            raise TargetingError(msg)

        return mask

    def __eq__(self, other):
        # Technically, different strings could refer to the same targets due to
        # implicit symmetry, but we'll ignore that for now
        if not isinstance(other, BSSubtarget):
            return False
        return not self.target_str != other.target_str


class BSTarget:
    def __init__(self, target_str):
        self.target_str = target_str
        self.subtargets = []
        self._parse_target(target_str)

    def _parse_target(self, target_str):
        """Parse multiple subtargets separated by commas."""
        subtarget_strs = target_str.split(',')
        self.subtargets = [BSSubtarget(sub.strip()) for sub in subtarget_strs]

    def select(self, atom_basis):
        """Take the 'or' of all subtargets, combining masks."""
        combined_mask = np.full(len(atom_basis), False)
        for subtarget in self.subtargets:
            combined_mask = combined_mask | subtarget.select(atom_basis)
        return combined_mask

    def __eq__(self, other):
        if not isinstance(other, BSTarget):
            return False
        if len(self.subtargets) != len(other.subtargets):
            return False
        for sub1, sub2 in zip(self.subtargets, other.subtargets):
            if sub1 != sub2:
                return False
        return True

    def __repr__(self):
        return f'BSTarget({self.target_str})'
