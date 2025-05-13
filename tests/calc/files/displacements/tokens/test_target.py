# test_target_token.py
import re

import pytest

from viperleed_jax.files.displacements.tokens.target import (
    TargetingError,
    TargetToken,
    _generate_label_match_regex,
)


# Minimal fake scatterer for tests
class FakeScatterer:
    def __init__(self, site, num, layer):
        self.site = site
        self.num = num
        self.layer = layer


# Build a small atom_basis:
#  A1@layer0, A2@layer1, A3@layer2, B1@layer0, C5@layer1
atom_basis = [
    FakeScatterer('A', 1, 0),
    FakeScatterer('A', 2, 1),
    FakeScatterer('A_extra', 3, 2),
    FakeScatterer('B', 4, 0),
    FakeScatterer('C', 5, 1),
]


def mask_to_indices(mask):
    """Helper to list the 1-based nums of selected atoms for easier asserts."""
    return [atom_basis[i].num for i, m in enumerate(mask) if m]


def test_simple_label_matching():
    # Exact label
    t = Targets('A')
    expected = [True, True, True, False, False]
    mask = t.select(atom_basis)
    assert np.all(mask == expected)
    assert mask_to_indices(mask) == [1, 2, 3]
    # Wildcard postfix
    t2 = Targets('A_*')
    expected2 = [False, False, True, False, False]
    mask2 = t2.select(atom_basis)
    assert np.all(mask2 == expected2)
    assert mask_to_indices(mask2) == [3]
    # Bare '*' matches all
    t3 = Targets('*')
    expected3 = [True, True, True, True, True]
    mask3 = t3.select(atom_basis)
    assert np.all(mask3 == expected3)
    assert mask_to_indices(mask3) == [1, 2, 3, 4, 5]


def test_numeric_subtarget_list_and_range():
    # List of atom numbers
    t = Targets('A 1 3')
    mask = t.select(atom_basis)
    assert mask_to_indices(mask) == [1, 3]
    # Range notation
    t2 = Targets('A 1-3')
    mask2 = t2.select(atom_basis)
    assert mask_to_indices(mask2) == [1, 2, 3]


def test_layer_selection():
    # L(1) means layer==0 internally (layer+1==1)
    t = Targets('A L(1)')
    mask = t.select(atom_basis)
    # Only A1 is on layer0
    assert mask_to_indices(mask) == [1]
    # Multi-layer range
    t2 = Targets('A L(1-2)')
    mask2 = t2.select(atom_basis)
    assert mask_to_indices(mask2) == [1, 2]


def test_comma_separated_subtargets_union():
    # Combine two subtargets
    t = Targets('A 1, B')
    mask = t.select(atom_basis)
    # A1 plus all B
    assert sorted(mask_to_indices(mask)) == [1, 4]


def test_invalid_atom_number_raises():
    # 10 is out of bounds
    with pytest.raises(TargetingError):
        Targets('A 10').select(atom_basis)
    # Range reversed gives empty nums → error
    with pytest.raises(TargetingError):
        Targets('A 3-1').select(atom_basis)


def test_no_matches_raises():
    # Label that matches nothing
    with pytest.raises(TargetingError):
        Targets('Z').select(atom_basis)
    # Valid label but wrong num for that label
    with pytest.raises(TargetingError):
        Targets('B 2').select(atom_basis)


def test_empty_and_bad_inputs():
    # Empty subtarget string
    with pytest.raises(ValueError):
        TargetToken('').select(atom_basis)
    # Bad layer syntax
    with pytest.raises(ValueError):
        TargetToken('A L()')
    # Non-integer num token
    with pytest.raises(ValueError):
        TargetToken('A one').select(atom_basis)


def test_repr():
    # repr(Targets)
    t = Targets('A 1')
    assert repr(t).startswith('Target(')
