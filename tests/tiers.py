"""Which tests belong to the slow and TensErLEED tiers.

Data only; the hooks that use it are in slow_marks.py.

Attributes
----------
TENSORLEED_FIXTURES : frozenset
    Names of fixtures that run, or need the sources of, TensErLEED.
SLOW_THRESHOLD : float
    Duration (seconds, call phase) from which a test counts as slow.
SLOW_BAND : float
    Relative margin around SLOW_THRESHOLD in which a test keeps its
    marking. A test that is not marked 'slow' is listed from
    (1 + SLOW_BAND) * SLOW_THRESHOLD, a marked one is dropped below
    (1 - SLOW_BAND) * SLOW_THRESHOLD.
SLOW_CASES : dict
    Maps tests (node-id prefix) to the cases that are marked 'slow'.
RELEASE_CASES : dict
    Maps fixture names to the case-id parts for which the tests that
    use the fixture are marked 'release'. These are TensErLEED cases
    that take long and add little over others, and only run with the
    --release option. An entry is reported as stale if tests use its
    fixture, but none has the case part. A fixture that no collected
    test uses is not checked, as this cannot be judged when running
    a subset of the suite.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-10-23'
__license__ = 'GPLv3+'

SLOW_THRESHOLD = 0.5
SLOW_BAND = 0.2

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
# Cases and groups are counted if they have a call phase, also if they
# skip or xfail in the test body. The rule has a band of SLOW_BAND
# (20 %) around SLOW_THRESHOLD, because chunks of 0.5-0.8 s move
# between unrelated cases whenever the code changes: a case that is
# not marked is listed from 1.2 times SLOW_THRESHOLD, a listed one is
# dropped below 0.8 times; in between it keeps its marking.
# Only these cases are marked 'slow'.
_SYMMETRY_TESTS = 'tests/calc/symmetry/test_symmetry.py::'
_SLAB_TESTS = 'tests/calc/classes/test_slab.py::'
_POSCAR_WRITE = 'tests/calc/files/test_poscar.py::TestPOSCARWrite::'
_DELETE = 'tests/calc/utilities/poscar/test_delete.py::'
_FIND_SYMMETRY = (
    'tests/calc/utilities/poscar/test_find_symmetry.py::TestFindSymmetry::'
    )
_MODIFY_VACUUM = (
    'tests/calc/utilities/poscar/test_modify_vacuum.py::TestModifyVacuum::'
    )
_SORT_BY_Z = 'tests/calc/utilities/poscar/test_sort_by_z.py::'
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
        '2x2-infoless_poscar-Al2O3_NiAl(111)_cHole_20061025',
        ),
    _SUPERCELL_TESTS + 'test_supercell_valid': (
        'infoless_poscar-Al2O3_NiAl(111)_cHole_20061025',
        'infoless_poscar-TiO2_supercell',
        ),
    _SLAB_TESTS + 'TestEquivalence::test_slab_translated_not_equivalent': (
        'infoless_poscar-Al2O3_NiAl(111)_cHole_20061025',
        ),
    _POSCAR_WRITE + 'test_group_conserved': (
        'infoless_poscar-36C_cm',
        'infoless_poscar-Al2O3_NiAl(111)_cHole_20061025',
        'infoless_poscar-TiO2_supercell', 'poscar-TiO2_supercell',
        'poscar_36carbon_atoms_cm',
        ),
    _DELETE + 'TestDeleteBetween::test_delete_between_cli': (
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
    _MODIFY_VACUUM + 'test_gap_relative': (
        'infoless_poscar-Al2O3_NiAl(111)_cHole_20061025',
        'infoless_poscar-TiO2_supercell-vacuum_gap_size=-1.0',
        ),
    _MODIFY_VACUUM + 'test_not_enough_vacuum_error': (
        'infoless_poscar-Al2O3_NiAl(111)_cHole_20061025',
        ),
    _MODIFY_VACUUM + 'test_shifted_slab_down': (
        'infoless_poscar-Al2O3_NiAl(111)_cHole_20061025',
        'infoless_poscar-In2O3_(111)-vacuum_gap_size=15.0',
        ),
    _SORT_BY_Z + 'test_overwrite_file': (
        'infoless_poscar-Al2O3_NiAl(111)_cHole_20061025',
        ),
    _WOODS_TESTS + 'TestRaises::test_invalid_init_args': (
        'gamma too complex',
        ),
    _PLANE_GROUP + 'test_any_plane_group_found': (
        'infoless_poscar-36C_cm', 'infoless_poscar-36C_p6m',
        'infoless_poscar-Al2O3_NiAl(111)_cHole_20061025',
        'infoless_poscar-Cu2O_111', 'infoless_poscar-In2O3_(111)',
        'infoless_poscar-TiO2_supercell', 'infoless_poscar-graphene',
        'poscar-Al2O3_NiAl(111)_cHole_20061025', 'poscar-TiO2_supercell',
        'poscar_36carbon_atoms_cm', 'poscar_36carbon_atoms_p6m',
        'poscar_sb_si_111',
        ),
    _PLANE_GROUP + 'test_correct_plane_group': (
        'infoless_poscar-36C_cm',
        'infoless_poscar-Al2O3_NiAl(111)_cHole_20061025',
        'poscar-Al2O3_NiAl(111)_cHole_20061025', 'poscar-TiO2_supercell',
        'poscar_36carbon_atoms_cm', 'poscar_36carbon_atoms_p6m',
        'poscar_sb_si_111',
        ),
    _SYMMETRIZATION + 'test_rattled_group': (
        'infoless_poscar-36C_cm', 'infoless_poscar-36C_p6m',
        'infoless_poscar-Al2O3_NiAl(111)_cHole_20061025',
        'infoless_poscar-STO(110)-4x1', 'infoless_poscar-Sb_Si(111)_rect',
        'infoless_poscar-TiO2_supercell', 'infoless_poscar-graphene',
        'poscar-Al2O3_NiAl(111)_cHole_20061025', 'poscar-TiO2_supercell',
        'poscar_36carbon_atoms_cm', 'poscar_36carbon_atoms_p6m',
        'poscar_fe2o3_012_1x1', 'poscar_sb_si_111',
        ),
    _CONSTRAINTS + 'test_correct_atom_linking': (
        'infoless_poscar-36C_cm', 'infoless_poscar-36C_p6m',
        'infoless_poscar-Al2O3_NiAl(111)_cHole_20061025',
        'infoless_poscar-Fe3O4_(001)_cod1010369',
        'infoless_poscar-In2O3_(111)', 'infoless_poscar-STO(110)-4x1',
        'infoless_poscar-Sb_Si(111)_rect', 'infoless_poscar-TiO2_supercell',
        'infoless_poscar-diamond', 'infoless_poscar-graphene',
        'poscar-Al2O3_NiAl(111)_cHole_20061025', 'poscar-Fe3O4_SCV',
        'poscar-TiO2_supercell', 'poscar_36carbon_atoms_cm',
        'poscar_36carbon_atoms_p6m', 'poscar_diamond', 'poscar_fe3o4_001_cod',
        'poscar_mgo', 'poscar_sb_si_111', 'rcmm',
        ),
    _CONSTRAINTS + 'test_correct_free_directions': (
        'infoless_poscar-36C_cm',
        'infoless_poscar-Al2O3_NiAl(111)_cHole_20061025',
        'infoless_poscar-TiO2_supercell',
        'poscar-Al2O3_NiAl(111)_cHole_20061025', 'poscar-TiO2_supercell',
        'poscar_36carbon_atoms_cm', 'poscar_fe3o4_001_cod',
        ),
    _CONSTRAINTS + 'test_correct_group_after_constraint': (
        'double_bulk-fe3o4', 'hex_cm_11', 'infoless_poscar-36C_cm',
        'infoless_poscar-Al2O3_NiAl(111)_cHole_20061025',
        'infoless_poscar-Fe3O4_(001)_cod1010369',
        'infoless_poscar-TiO2_supercell',
        'poscar-Al2O3_NiAl(111)_cHole_20061025',
        'poscar-Cu2O(111)_1x1_surplus_oxygen', 'poscar-Fe3O4_SCV',
        'poscar-TiO2_supercell', 'poscar_36carbon_atoms_cm',
        'poscar_36carbon_atoms_p6m', 'poscar_diamond', 'poscar_fe3o4_001_cod',
        'poscar_sb_si_111', 'poscar_sto110_4x1',
        ),
    _CONSTRAINTS + 'test_preserve_displaced_symmetry': (
        'displacement=([0.2, 0, 0],)-poscar_lsmo_001_rt2', 'hex_cm_01',
        'infoless_poscar-36C_cm', 'infoless_poscar-36C_p6m',
        'infoless_poscar-Al2O3_NiAl(111)_cHole_20061025',
        'infoless_poscar-Cu2O_111', 'infoless_poscar-LSMO_001_A_site',
        'infoless_poscar-TiO2_supercell', 'infoless_poscar-graphene',
        'poscar-Al2O3_NiAl(111)_cHole_20061025', 'poscar-TiO2_supercell',
        'poscar_36carbon_atoms_cm', 'poscar_36carbon_atoms_p6m',
        ),
    _REDUCTION + 'test_invalid_group': (
        'infoless_poscar-36C_cm', 'infoless_poscar-graphene',
        'poscar_36carbon_atoms_cm',
        ),
    _REDUCTION + 'test_invalid_subgroup': (
        'infoless_poscar-Al2O3_NiAl(111)_cHole_20061025',
        'poscar_36carbon_atoms_cm',
        ),
    _REDUCTION + 'test_missing_direction': (
        'infoless_poscar-36C_cm', 'infoless_poscar-36C_p6m',
        'infoless_poscar-Al2O3_NiAl(111)_cHole_20061025',
        'infoless_poscar-Cu2O(111)_1x1_surplus_oxygen',
        'infoless_poscar-In2O3_(111)', 'infoless_poscar-TiO2_supercell',
        'infoless_poscar-graphene', 'p3m1',
        'poscar-Al2O3_NiAl(111)_cHole_20061025', 'poscar-TiO2_supercell',
        'poscar_36carbon_atoms_cm', 'poscar_36carbon_atoms_p6m',
        'poscar_lsmo_001_rt2',
        ),
    }

# Cases that only run with --release (and are also 'tensorleed'): each
# fixture is mapped to case-id parts, in the same format as for
# SLOW_CASES. A test is a 'release' case if it uses the fixture and its
# case id has one of the parts. This makes any consumer of the
# run_phaseshift fixture for the Al2O3 slab a release case: its
# PHASESHIFTS calculation is by far the longest one of the tier.
_AL2O3_CHOLE = 'poscar-Al2O3_NiAl(111)_cHole_20061025'
RELEASE_CASES = {
    'run_phaseshift': (_AL2O3_CHOLE,),
    }
