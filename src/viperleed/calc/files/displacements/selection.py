"""Target selection module."""

__authors__ = ('Alexander M. Imre (@amimre)',)
__created__ = '2025-05-12'

import numpy as np

from .tokens.target import TargetToken


class TargetSelectionError(ValueError):
    """Error raised if target selection goes wrong."""


class TargetSelection:
    """Class to apply the targeting by TargetTokens to a basis."""

    def __init__(self, targets):
        if not all(isinstance(target, TargetToken) for target in targets):
            raise TargetSelectionError(
                'All supplied targets must be of TargetToken type.')
        self.targets = targets

    def select(self, atom_basis):
        """Take the 'or' of all targets, combining masks."""
        combined_mask = np.full(len(atom_basis), fill_value=False)
        for target in self.targets:
            combined_mask = (combined_mask |
                             create_target_selection_mask(atom_basis, target))
        return combined_mask


def create_target_selection_mask(atom_basis, target_token):
    """Select base scatterers that match the target specification."""
    mask = np.full(len(atom_basis), fill_value=True)

    # mask based on the site
    matches = np.array(
        [target_token.regex.match(bs.site) is not None for bs in atom_basis]
    )
    mask = mask & matches

    # mask based on the labels
    label_mask = mask.copy()

    # If nums are specified, apply the selection based on nums
    if target_token.nums is not None:
        # check range for nums
        if any(num < 1 or num > len(atom_basis) for num in target_token.nums):
            msg = f'Invalid atom number for target: {target_token.target_str}'
            raise TargetSelectionError(msg)
        num_mask = np.array([bs.num in target_token.nums for bs in atom_basis])
        # check if any of the given nums have the wrong label
        wrong_label = np.logical_and(num_mask, ~label_mask)
        if np.any(wrong_label):
            msg = (
                'Atom numbers do not match label for target: '
                f'{target_token.target_str}'
            )
            raise TargetSelectionError(msg)
        mask = mask & num_mask

    # If layers are specified, apply the selection based on layers
    if target_token.layers is not None:
        mask = mask & np.array(
            # TODO: layer counting from 1; can we unify this somewhere?
            [bs.layer + 1 in target_token.layers for bs in atom_basis]
        )

    if mask.sum() == 0:
        msg = f'No atoms selected for TargetToken: {target_token}'
        raise TargetSelectionError(msg)

    return mask
