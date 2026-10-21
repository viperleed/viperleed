"""Test configuration for viperleed.tests.

Defines fixtures and fixture factories used in multiple tests.

Random-input tests (those requesting `draws`) check N_DRAWS draws
seeded from the test id, so that outcomes are reproducible and known
failures can be listed (xfail). The opt-in option --rng-seed=<int> or
--rng-seed=random derives the seeds from a different base, to explore
other draws locally; it is not meant for CI. With the option, the
xfail marks of random-input tests are removed, and `draws.fresh` is
True: tests that expect a known failure (e.g., via may_fail or inline
pytest.xfail) must then just assert, so that the real outcome is
reported. With `random`, the base is picked once by the controller
(and passed on to pytest-xdist workers) and is printed in the report
header and in the terminal summary, to reproduce a run.

Fixtures
--------
check_log_records (factory)
    Raise unless caplog records are exactly as expected.
data_path
    Path to the top-level folder containing test data.
draws
    N_DRAWS seeded random generators for one random-input test.
    Its `fresh` attribute tells whether --rng-seed was given.
mock_import_module (factory)
    Replace importlib.import_module to return selected fake modules.
    Retain normal behavior for all the others.
first_case
    The first of the current pytest-cases.
re_match (factory)
    Return a match object from a pattern and a string.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    'Alexandra Mia Imre (@alexmiame)',
)
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2023-02-28'
__license__ = 'GPLv3+'

import collections
from importlib import import_module
import random
import re
import zlib

import numpy as np
import pytest
from pytest_cases import fixture

from . import optional_dependencies
from .helpers import TEST_DATA

pytest_plugins = ('pytester',)  # For testing the options themselves

# Fixture of optional_dependencies that needs no wrapping
fixture_missing_dependency = optional_dependencies.fixture_missing_dependency

# Number of seeded draws a random-input test checks, all within one
# test (one id). This is the only place where the number is defined.
N_DRAWS = 5

Draw = collections.namedtuple('Draw', ('index', 'seed', 'rng'))

# Fixtures that run, or need the sources of, TensErLEED. Any test
# that has one of them in its fixture closure is marked 'tensorleed'
TENSORLEED_FIXTURES = frozenset((
    'tensorleed_path',
    'run_from_ase_refcalc',
    ))

# Cases that take about half a second or more (call phase, serial
# run). Each test (node-id prefix) is mapped to the case-id parts,
# i.e., the '-'-separated fixture and slab names, that mark its slow
# cases. A case is listed if its call took about 0.5 s or more. A
# whole (test, slab) group is listed if at least 40 % of its cases
# did, or if one took 0.7 s or more, as the cost of a first use of a
# slab moves between the cases of a group. Only these cases are
# marked 'slow'.
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

RNG_SEED_OPTION = '--rng-seed'
WORKER_SEED_KEY = 'rng_seed'  # In pytest-xdist's workerinput

class Draws(tuple):
    """A tuple of Draw objects, knowing whether the base seed is fresh."""

    fresh = False


def pytest_addoption(parser):
    """Add the options of the test suite."""
    optional_dependencies.pytest_addoption(parser)
    parser.addoption(
        RNG_SEED_OPTION,
        dest='rng_seed',
        metavar='INT|random',
        default=None,
        help=('Derive the seeds of random-input tests from this base '
              'integer, or from a random one if "random" (printed at '
              'the end of the run, to reproduce it). Known failures '
              '(xfail) of random-input tests then report their real '
              'outcome. For local exploration, not for CI.'),
        )


def pytest_configure(config):
    """Resolve --rng-seed to an integer, once, in the controller.

    A pytest-xdist worker does not resolve it: it takes the base that
    the controller resolved and handed over (pytest_configure_node).
    """
    workerinput = getattr(config, 'workerinput', None)
    if workerinput is not None:
        config.option.rng_seed = workerinput[WORKER_SEED_KEY]
        return
    value = config.option.rng_seed
    if value is None:
        return  # No option
    if value == 'random':
        value = random.SystemRandom().randrange(2**32)
    try:
        config.option.rng_seed = int(value)
    except ValueError:
        raise pytest.UsageError(
            f'{RNG_SEED_OPTION} must be an integer or "random", '
            f'not {value!r}'
            ) from None


@pytest.hookimpl(optionalhook=True)  # Only exists with pytest-xdist
def pytest_configure_node(node):
    """Hand the resolved base seed over to a pytest-xdist worker."""
    node.workerinput[WORKER_SEED_KEY] = node.config.option.rng_seed


def _rng_seed_message(config):
    """Return a line with the base seed, or None if not used."""
    base = config.option.rng_seed
    if base is None:
        return None
    return f'rng-seed base: {base} (reproduce with --rng-seed={base})'


def pytest_report_header(config):
    """Report the base seed, if used."""
    return _rng_seed_message(config)


def pytest_terminal_summary(terminalreporter, config):
    """Report the base seed, if used, also with -q."""
    message = _rng_seed_message(config)
    if message:
        terminalreporter.write_line(message)


def _is_slow(item):
    """Return whether `item` is one of the SLOW_CASES."""
    test, bracket, params = item.nodeid.partition('[')
    if not bracket:
        return False
    # Pad, to match whole '-'-separated parts of the case id
    params = f'-{params[:-1]}-'
    return any(f'-{part}-' in params for part in SLOW_CASES.get(test, ()))


def pytest_collection_modifyitems(items):
    """Mark TensErLEED and slow tests; no xfail with fresh seeds."""
    for item in items:
        if TENSORLEED_FIXTURES.intersection(item.fixturenames):
            item.add_marker(pytest.mark.tensorleed)
        if _is_slow(item):
            item.add_marker(pytest.mark.slow)
        if (item.config.option.rng_seed is not None
                and 'draws' in item.fixturenames):
            item.own_markers[:] = [
                mark for mark in item.own_markers if mark.name != 'xfail'
                ]


@fixture
def check_log_records(caplog):
    """Raise unless log records are exactly as expected."""
    def _check(expected_records):
        logged = tuple(r.getMessage() for r in caplog.records)
        assert len(logged) == len(expected_records)
        for log, expect in zip(logged, expected_records):
            if isinstance(expect, str):
                assert log == expect
            else:
                assert expect.fullmatch(log)
    return _check


@fixture(scope='session')
def re_match():  # This is actually a fixture factory
    """Return a re.match object from a pattern and a string."""
    def _match(pattern, string):
        return re.match(pattern, string)
    return _match


@fixture(scope='session')
def data_path():
    """Return the Path to the top-level folder containing test data."""
    return TEST_DATA


@fixture(name='first_case')
def first_case(current_cases):
    """Return the first of the current cases."""
    def _find_case(cases_dict):
        for value in cases_dict.values():
            if isinstance(value, dict):
                try:
                    return _find_case(value)
                except ValueError:
                    pass
            try:
                value.id
            except AttributeError:
                pass
            else:
                return value
        raise ValueError('No case found')
    return _find_case(current_cases)


@fixture
def draws(request):
    """Return N_DRAWS Draw tuples, seeded from the test id."""
    # Seed from the test id (without the file path), so that the
    # random values of a test do not depend on which other tests
    # run, in which order, or in which process: outcomes are then
    # the same in every run. Not hash(): it depends on PYTHONHASHSEED.
    _, test_id = request.node.nodeid.split('::', 1)
    base = request.config.option.rng_seed
    prefix = '' if base is None else f'{base}:'
    seeds = (
        zlib.crc32(f'{prefix}{test_id}#{index}'.encode())
        for index in range(N_DRAWS)
        )
    draws = Draws(
        Draw(index, seed, np.random.default_rng(seed))
        for index, seed in enumerate(seeds)
        )
    draws.fresh = base is not None
    return draws


@fixture
def mock_import_module(mocker):
    """Make importlib.import_module return fake modules."""
    def _patch(**fake_modules):
        def _fake_import(module_name, *args, **kwargs):
            try:
                return fake_modules[module_name]
            except KeyError:
                return import_module(module_name, *args, **kwargs)
        mocker.patch('importlib.import_module', side_effect=_fake_import)
        mocker.patch.dict('sys.modules', fake_modules)
    return _patch
