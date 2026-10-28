"""Which tests belong to the slow and TensErLEED tiers.

Data only; the hooks that use it are in slow_marks.py.

Attributes
----------
TENSORLEED_FIXTURES : frozenset
    Names of fixtures that run, or need the sources of, TensErLEED.
SLOW_THRESHOLD : float
    Duration (seconds, call phase) from which a test counts as slow.
SLOW_CASES : dict
    Maps tests (node-id prefix) to the cases that are marked 'slow'.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-10-23'
__license__ = 'GPLv3+'

SLOW_THRESHOLD = 0.5

TENSORLEED_FIXTURES = frozenset((
    'tensorleed_path',
    'run_from_ase_refcalc',
    ))

# Cases that take about SLOW_THRESHOLD or more (call phase, serial
# run). Each test (node-id prefix) is mapped to the case-id parts,
# i.e., the '-'-separated fixture and slab names, that mark its slow
# cases. A case is listed if its call took about SLOW_THRESHOLD or
# more. A whole (test, slab) group is listed if at least 40 % of its
# cases did, or if one took 1.4 times SLOW_THRESHOLD or more, as the
# cost of a first use of a slab moves between the cases of a group.
# Only these cases are marked 'slow'.
_SYMMETRY_TESTS = 'tests/calc/symmetry/test_symmetry.py::'
_SLAB_TESTS = 'tests/calc/classes/test_slab.py::'
_POSCAR_WRITE = 'tests/calc/files/test_poscar.py::TestPOSCARWrite::'
_FIND_SYMMETRY = (
    'tests/calc/utilities/poscar/test_find_symmetry.py::TestFindSymmetry::'
    )
_MODIFY_VACUUM = (
    'tests/calc/utilities/poscar/test_modify_vacuum.py::TestModifyVacuum::'
    )
_WOODS_TESTS = 'tests/gui/leedsim/classes/woods/test_woods.py::'
_SUPERCELL_TESTS = _SLAB_TESTS + 'TestSuperAndSubCell::'
_PLANE_GROUP = _SYMMETRY_TESTS + 'TestPlaneGroupFinding::'
_CONSTRAINTS = _SYMMETRY_TESTS + 'TestSymmetryConstraints::'
_REDUCTION = _SYMMETRY_TESTS + 'TestSymmetryReduction::'
_SYMMETRIZATION = _SYMMETRY_TESTS + 'TestSlabSymmetrization::'
SLOW_CASES = {
    _SUPERCELL_TESTS + 'test_subcell_reverses_supercell': (
        '2x1-infoless_poscar-Al2O3_NiAl(111)_cHole_20061025',
        '2x2-infoless_poscar-Al2O3_NiAl(111)_cHole_20061025',
        '2x2-infoless_poscar-TiO2_supercell',
        'off-diagonal-infoless_poscar-Al2O3_NiAl(111)_cHole_20061025',
        ),
    _SUPERCELL_TESTS + 'test_supercell_does_not_change_original_slab': (
        '2x1-infoless_poscar-TiO2_supercell',
        '2x2-infoless_poscar-Al2O3_NiAl(111)_cHole_20061025',
        '2x2-infoless_poscar-Fe3O4_(001)_cod1010369',
        'off-diagonal-infoless_poscar-In2O3_(111)',
        ),
    _SUPERCELL_TESTS + 'test_supercell_valid': (
        '2x1-infoless_poscar-Al2O3_NiAl(111)_cHole_20061025',
        '2x1-infoless_poscar-TiO2_supercell',
        '2x2-infoless_poscar-Al2O3_NiAl(111)_cHole_20061025',
        '2x2-infoless_poscar-TiO2_supercell',
        'identity-infoless_poscar-Al2O3_NiAl(111)_cHole_20061025',
        'off-diagonal-infoless_poscar-Al2O3_NiAl(111)_cHole_20061025',
        ),
    _SLAB_TESTS + 'TestEquivalence::test_equivalent_copy': (
        'infoless_poscar-Al2O3_NiAl(111)_cHole_20061025',
        ),
    _SLAB_TESTS + 'TestEquivalence::test_slab_translated_not_equivalent': (
        'infoless_poscar-Al2O3_NiAl(111)_cHole_20061025',
        ),
    _SLAB_TESTS + 'TestExtraBulk::test_double_thickness_twice': (
        'poscar_tio2_small',
        ),
    _SLAB_TESTS + 'TestRevertUnitCell::test_revert_in_two_steps': (
        'infoless_poscar-STO(110)-4x1-n_undo=4',
        ),
    _SLAB_TESTS + 'TestUnitCellReduction::test_ab_cell_already_minimal': (
        'poscar-Al2O3_NiAl(111)_cHole_20061025',
        ),
    _POSCAR_WRITE + 'test_group_conserved': (
        'infoless_poscar-36C_cm', 'infoless_poscar-In2O3_(111)',
        'infoless_poscar-TiO2_supercell',
        'poscar-Al2O3_NiAl(111)_cHole_20061025', 'poscar-TiO2_supercell',
        'poscar_36carbon_atoms_cm',
        ),
    _POSCAR_WRITE + 'test_nr_atoms_conserved': (
        'infoless_poscar-Al2O3_NiAl(111)_cHole_20061025',
        ),
    _FIND_SYMMETRY + 'test_pipe_poscar': (
        'poscar-Al2O3_NiAl(111)_cHole_20061025', 'poscar-TiO2_supercell',
        ),
    _FIND_SYMMETRY + 'test_plane_group_to_outfile': (
        'poscar-Al2O3_NiAl(111)_cHole_20061025', 'poscar-TiO2_supercell',
        ),
    _FIND_SYMMETRY + 'test_plane_group_to_stdout': (
        'poscar-Al2O3_NiAl(111)_cHole_20061025', 'poscar-TiO2_supercell',
        ),
    _MODIFY_VACUUM + 'test_accept_small_gap': (
        'infoless_poscar-Al2O3_NiAl(111)_cHole_20061025',
        ),
    _MODIFY_VACUUM + 'test_gap_absolute': (
        'infoless_poscar-Al2O3_NiAl(111)_cHole_20061025',
        ),
    _MODIFY_VACUUM + 'test_gap_relative': (
        'infoless_poscar-Al2O3_NiAl(111)_cHole_20061025',
        ),
    _MODIFY_VACUUM + 'test_not_enough_vacuum_error': (
        'infoless_poscar-TiO2_supercell',
        ),
    _MODIFY_VACUUM + 'test_shifted_slab_down': (
        'infoless_poscar-Al2O3_NiAl(111)_cHole_20061025',
        'infoless_poscar-TiO2_supercell',
        ),
    _MODIFY_VACUUM + 'test_zero_gap_raises': (
        'infoless_poscar-TiO2_supercell',
        ),
    _WOODS_TESTS + 'TestRaises::test_invalid_init_args': (
        'gamma too complex',
        ),
    _PLANE_GROUP + 'test_any_plane_group_found': (
        'infoless_poscar-36C_cm', 'infoless_poscar-36C_p6m',
        'infoless_poscar-Al2O3_NiAl(111)_cHole_20061025',
        'infoless_poscar-In2O3_(111)', 'infoless_poscar-TiO2_supercell',
        'infoless_poscar-graphene', 'poscar-Al2O3_NiAl(111)_cHole_20061025',
        'poscar-TiO2_supercell', 'poscar_36carbon_atoms_cm',
        'poscar_36carbon_atoms_p6m', 'poscar_fe3o4_001_cod',
        'poscar_sb_si_111',
        ),
    _PLANE_GROUP + 'test_correct_plane_group': (
        'hex_cm_10', 'infoless_poscar-36C_cm',
        'infoless_poscar-Al2O3_NiAl(111)_cHole_20061025',
        'infoless_poscar-TiO2_supercell',
        'poscar-Al2O3_NiAl(111)_cHole_20061025', 'poscar-TiO2_supercell',
        'poscar_36carbon_atoms_cm', 'poscar_36carbon_atoms_p6m',
        'poscar_ru0001_rt3_te', 'poscar_sb_si_111', 'poscar_sto110_4x1',
        ),
    _SYMMETRIZATION + 'test_rattled_group': (
        'hex_cmm_11', 'infoless_poscar-36C_cm', 'infoless_poscar-36C_p6m',
        'infoless_poscar-Al2O3_NiAl(111)_cHole_20061025',
        'infoless_poscar-Cu2O_111', 'infoless_poscar-Fe3O4_SCV',
        'infoless_poscar-In2O3_(111)', 'infoless_poscar-Sb_Si(111)_rect',
        'infoless_poscar-TiO2_supercell', 'infoless_poscar-diamond',
        'infoless_poscar-graphene', 'poscar-Al2O3_NiAl(111)_cHole_20061025',
        'poscar-Fe3O4_SCV', 'poscar-TiO2_supercell',
        'poscar_36carbon_atoms_cm', 'poscar_36carbon_atoms_p6m',
        'poscar_fe3o4_001_cod', 'poscar_sb_si_111', 'square_pg_10',
        ),
    _CONSTRAINTS + 'test_correct_atom_linking': (
        'bulk_repeat_poscar-Cu2O_111', 'infoless_poscar-36C_cm',
        'infoless_poscar-36C_p6m',
        'infoless_poscar-Al2O3_NiAl(111)_cHole_20061025',
        'infoless_poscar-Cu2O(111)_1x1_surplus_oxygen',
        'infoless_poscar-Fe3O4_(001)_cod1010369', 'infoless_poscar-Fe3O4_SCV',
        'infoless_poscar-In2O3_(111)', 'infoless_poscar-STO(110)-4x1',
        'infoless_poscar-Sb_Si(111)_rect', 'infoless_poscar-TiO2_supercell',
        'infoless_poscar-diamond', 'infoless_poscar-graphene', 'p31m', 'p6m',
        'poscar-Al2O3_NiAl(111)_cHole_20061025',
        'poscar-Cu2O(111)_1x1_surplus_oxygen', 'poscar-Fe3O4_SCV',
        'poscar-TiO2_supercell', 'poscar_36carbon_atoms_cm',
        'poscar_36carbon_atoms_p6m', 'poscar_diamond', 'poscar_fe3o4_001_cod',
        'poscar_lsmo_001_rt2', 'poscar_sb_si_111', 'poscar_sto110_4x1',
        'square_cmm',
        ),
    _CONSTRAINTS + 'test_correct_free_directions': (
        'hex_cm_12', 'infoless_poscar-36C_cm',
        'infoless_poscar-Al2O3_NiAl(111)_cHole_20061025',
        'infoless_poscar-TiO2_supercell', 'p31m',
        'poscar-Al2O3_NiAl(111)_cHole_20061025', 'poscar-TiO2_supercell',
        'poscar_36carbon_atoms_cm', 'poscar_diamond', 'poscar_fe3o4_001_cod',
        'poscar_lsmo_001_rt2', 'square_pm_01',
        ),
    _CONSTRAINTS + 'test_correct_group_after_constraint': (
        'hex_cm_12', 'infoless_poscar-36C_cm',
        'infoless_poscar-Al2O3_NiAl(111)_cHole_20061025',
        'infoless_poscar-Fe3O4_(001)_cod1010369',
        'infoless_poscar-Ir(100)-(2x1)-O', 'infoless_poscar-TiO2_supercell',
        'p6m', 'pmg-ucell=square',
        'poscar-Al2O3_NiAl(111)_cHole_20061025',
        'poscar-Fe3O4_SCV', 'poscar-TiO2_supercell',
        'poscar_36carbon_atoms_cm', 'poscar_36carbon_atoms_p6m',
        'poscar_diamond', 'poscar_fe3o4_001_cod', 'poscar_mgo',
        'poscar_ptrh_100_1x3_o', 'poscar_sb_si_111', 'poscar_sto110_4x1',
        'square_pm_01',
        ),
    _CONSTRAINTS + 'test_preserve_displaced_symmetry': (
        'bulk_repeat_poscar-Cu2O_111', 'fe3o4_bulk', 'hex_cm_01',
        'infoless_poscar-36C_cm', 'infoless_poscar-36C_p6m',
        'infoless_poscar-Al2O3_NiAl(111)_cHole_20061025',
        'infoless_poscar-In2O3_(111)', 'infoless_poscar-LSMO_001_A_site',
        'infoless_poscar-TiO2_supercell', 'infoless_poscar-graphene', 'p6',
        'poscar-Al2O3_NiAl(111)_cHole_20061025',
        'poscar-Cu2O(111)_1x1_surplus_oxygen', 'poscar-Fe3O4_SCV',
        'poscar-TiO2_supercell', 'poscar_36carbon_atoms_cm',
        'poscar_36carbon_atoms_p6m', 'poscar_sto110_4x1',
        ),
    _REDUCTION + 'test_invalid_group': (
        'infoless_poscar-36C_cm', 'infoless_poscar-36C_p6m',
        'infoless_poscar-Al2O3_NiAl(111)_cHole_20061025',
        'infoless_poscar-TiO2_supercell',
        'poscar-Al2O3_NiAl(111)_cHole_20061025', 'poscar-TiO2_supercell',
        'poscar_36carbon_atoms_cm', 'poscar_lsmo_001_rt2',
        ),
    _REDUCTION + 'test_invalid_subgroup': (
        'infoless_poscar-Al2O3_NiAl(111)_cHole_20061025',
        'infoless_poscar-TiO2_supercell',
        'poscar-Al2O3_NiAl(111)_cHole_20061025', 'poscar_36carbon_atoms_cm',
        ),
    _REDUCTION + 'test_missing_direction': (
        'hex_cmm_01', 'hex_cmm_10', 'infoless_poscar-36C_cm',
        'infoless_poscar-36C_p6m',
        'infoless_poscar-Al2O3_NiAl(111)_cHole_20061025',
        'infoless_poscar-Fe3O4_SCV', 'infoless_poscar-STO(110)-4x1',
        'infoless_poscar-TiO2_supercell', 'infoless_poscar-graphene',
        'poscar-Al2O3_NiAl(111)_cHole_20061025', 'poscar-TiO2_supercell',
        'poscar_36carbon_atoms_cm', 'poscar_36carbon_atoms_p6m',
        'poscar_diamond', 'poscar_sb_si_111', 'poscar_sto110_4x1',
        ),
    }
