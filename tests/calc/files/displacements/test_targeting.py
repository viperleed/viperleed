import numpy as np
import pytest
from pytest_cases import case, parametrize_with_cases

from viperleed_jax.atom_basis import AtomBasis
from viperleed_jax.files.displacements.targeting import BSTarget

# Test cases for valid target with the expected mask
CU_111_VALID_TARGETS = {
    'Cu_surf': [True, False, False, False, False],
    'Cu_def 2': [False, True, False, False, False],
    'Cu_def 2 4': [False, True, False, True, False],
    'Cu_def L(2-4)': [False, False, True, True, True],
}

# Test cases for invalid targeting
CU_111_INVALID_TARGETS = [
    'Cu_nothere',
    'Cu_def 0',
    'Cu_surf L(2)',
]


class ValidTargets:
    @pytest.mark.parametrize(
        'target_str_mask',
        CU_111_VALID_TARGETS.items(),
        ids=CU_111_VALID_TARGETS.keys(),
    )
    @case(tags='cu_111')
    def case_cu_111_exact(
        self, cu_111_fixed_l_max_state_after_init, target_str_mask
    ):
        slab, _ = cu_111_fixed_l_max_state_after_init
        atom_basis = AtomBasis(slab)
        target_str, expected_mask = target_str_mask
        return atom_basis, target_str, expected_mask


class InvalidTargets:
    @pytest.mark.parametrize(
        'target_str', CU_111_INVALID_TARGETS, ids=CU_111_INVALID_TARGETS
    )
    @case(tags='cu_111')
    def case_cu_111_exact(
        self, cu_111_fixed_l_max_state_after_init, target_str
    ):
        slab, _ = cu_111_fixed_l_max_state_after_init
        atom_basis = AtomBasis(slab)
        return atom_basis, target_str


@parametrize_with_cases(
    'atom_basis,target_str,expected_mask', cases=ValidTargets
)
def test_target_lexing(atom_basis, target_str, expected_mask):
    target = BSTarget(target_str)


@parametrize_with_cases(
    'atom_basis,target_str,expected_mask', cases=ValidTargets
)
def test_target_selection(atom_basis, target_str, expected_mask):
    target = BSTarget(target_str)
    mask = target.select(atom_basis)
    assert np.all(mask == expected_mask)


@parametrize_with_cases('atom_basis,target_str', cases=InvalidTargets)
def test_invalid_target(atom_basis, target_str):
    target = BSTarget(target_str)
    with pytest.raises(ValueError):
        target.select(atom_basis)
