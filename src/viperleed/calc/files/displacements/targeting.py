__authors__ = ("Alexander M. Imre (@amimre)",)
__created__ = "2024-10-14"

import re
import numpy as np


class ASESubtarget:
    def __init__(self, target_str):
        self.target_str = target_str
        self.element = None
        self.exact_label = None
        self.partial_label = None
        self.nums = None
        self.layers = None
        self._parse_target()

    def _parse_target(self):
        """Parse the site, and optional nums or layers from the target string."""
        parts = self.target_str.split()
        if not parts:
            raise ValueError("Subtarget string is empty")
        site_str = parts[0]

        if site_str.endswith("*"):
            # It's an inexact site label
            self.partial_label = site_str[:-1]
        elif "_" in site_str:
            self.exact_label = site_str
        else:
            self.element = site_str

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

        # Match element
        if self.element is not None:
            mask = mask & np.array(
                [bs.element == self.element for bs in base_scatterers]
            )
        # Match exact label
        if self.exact_label is not None:
            mask = mask & np.array(
                [bs.site == self.exact_label for bs in base_scatterers]
            )
        # Match partial label
        if self.partial_label is not None:
            mask = mask & np.array(
                [
                    bs.atom.site.starts_with(self.partial_label)
                    for bs in base_scatterers
                ]
            )

        # mask based on the labels
        label_mask = mask.copy()

        # If nums are specified, apply the selection based on nums
        if self.nums is not None:
            num_mask = np.array([bs.num in self.nums for bs in base_scatterers])
            # check if any of the given nums have the wrong label
            if np.any(np.logical_xor(mask, num_mask)):
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


class ASETarget:
    def __init__(self, target_str):
        self.subtargets = []
        self._parse_target(target_str)

    def _parse_target(self, target_str):
        """Parse multiple subtargets separated by commas."""
        subtarget_strs = target_str.split(",")
        self.subtargets = [ASESubtarget(sub.strip()) for sub in subtarget_strs]

    def select(self, base_scatterers):
        """Take the 'or' of all subtargets, combining masks."""
        combined_mask = np.full(len(base_scatterers), False)
        for subtarget in self.subtargets:
            combined_mask = combined_mask | subtarget.select(base_scatterers)
        return combined_mask
