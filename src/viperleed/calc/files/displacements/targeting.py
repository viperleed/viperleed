"""Module targeting."""
__authors__ = ("Alexander M. Imre (@amimre)",)
__created__ = "2024-10-14"

import re
import numpy as np


def generate_label_match_regex(label):
    """Generate a regex pattern to match variations of the given label."""
    # Escape any special characters in the label
    escaped_label = re.escape(label)

    # Replace '*' in the label with a regex pattern that matches any characters
    pattern = escaped_label.replace(r"\*", r"\w*")

    # Compile the final regex pattern
    return re.compile(rf"^{pattern}$")


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
            raise ValueError("Subtarget string is empty")
        site_str = parts[0]
        self.regex = generate_label_match_regex(site_str)

        if len(parts) == 1:
            # only site is specified, no nums or layers
            return

        # Check if we have a layer specification
        layer_match = re.match(r"L\((\d+)(-(\d+))?\)", parts[1])
        if layer_match:
            start_layer = int(layer_match.group(1))
            end_layer = (
                int(layer_match.group(3))
                if layer_match.group(3)
                else start_layer
            )
            self.layers = list(range(start_layer, end_layer + 1))
        else:
            # It's a list of numbers
            self.nums = list(map(int, parts[1:]))

    def select(self, base_scatterers):
        """Selects base scatterers that match the subtarget specification."""
        mask = np.full(len(base_scatterers), True)

        # mask based on the site
        matches = np.array([self.regex.match(bs.site) is not None 
                            for bs in base_scatterers])
        mask = mask & matches

        # mask based on the labels
        label_mask = mask.copy()

        # If nums are specified, apply the selection based on nums
        if self.nums is not None:
            # check range for nums
            if any(num < 1 or num > len(base_scatterers) for num in self.nums):
                raise ValueError(
                    "Invalid atom number for subtarget: " f"{self.target_str}"
                )
            num_mask = np.array([bs.num in self.nums for bs in base_scatterers])
            # check if any of the given nums have the wrong label
            wrong_label = np.logical_and(num_mask, ~label_mask)
            if np.any(wrong_label):
                raise ValueError(
                    "Atom numbers do not match label for subtarget: "
                    f"{self.target_str}"
                )
            mask = mask & num_mask

        # If layers are specified, apply the selection based on layers
        if self.layers is not None:
            mask = mask & np.array(
                [bs.layer in self.layers for bs in base_scatterers]
            )

        if mask.sum() == 0:
            raise ValueError(
                f"No atoms selected for subtarget: {self.target_str}"
            )

        return mask

    def __eq__(self, other):
        # Technically, different strings could refer to the same targets due to
        # implicit symmetry, but we'll ignore that for now
        if not isinstance(other, BSSubtarget):
            return False
        if self.target_str != other.target_str:
            return False


class BSTarget:
    def __init__(self, target_str):
        self.target_str = target_str
        self.subtargets = []
        self._parse_target(target_str)

    def _parse_target(self, target_str):
        """Parse multiple subtargets separated by commas."""
        subtarget_strs = target_str.split(",")
        self.subtargets = [BSSubtarget(sub.strip()) for sub in subtarget_strs]

    def select(self, base_scatterers):
        """Take the 'or' of all subtargets, combining masks."""
        combined_mask = np.full(len(base_scatterers), False)
        for subtarget in self.subtargets:
            combined_mask = combined_mask | subtarget.select(base_scatterers)
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
        return f"BSTarget({self.target_str})"
