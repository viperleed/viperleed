"""Tests for the input files of the search.

The PARAM, search.steu, and restrict.f files written by
viperleed.calc.files.iosearch.generateSearchInput are compared to
stored references, for each version group of TensErLEED (see
conftest.py). The search.steu is the one written before any search
has run. No TensErLEED tree, compiler, tensors, or deltas are needed:
the delta files only exist as names (empty stubs) that the atoms
register as their current deltas, and the execution is never started.

Everything that depends on the machine or on random numbers is pinned:
the number of cores, the population, the starting configuration (a
deterministic one rather than the default 'crandom'), the system name
and the timestamp. The seed of the random-number generator written to
search.steu is a constant (0, i.e., use the system time) and needs no
patching.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-11'
__license__ = 'GPLv3+'

from itertools import product

import pytest

from viperleed.calc.classes.beam import Beam
from viperleed.calc.files import iosearch
from viperleed.calc.lib.context import execute_in_dir

from ...helpers import TEST_DATA
from .conftest import DISPLACEMENTS
from .conftest import GROUPS

SEARCH_FILES = 'PARAM', 'search.steu', 'restrict.f'
SEARCH_POPULATION = 4
SEARCH_START = 'centered'
HUGE_MAX_GEN = 10**9  # Above what the input file format allows
MAX_GEN_WRITTEN = {'1.6.0-1.6.1': 999999, '1.7.3-': 9999999}
EXP_ENERGIES = tuple(range(50, 351))  # Synthetic experimental beams
EXP_BEAMS = ((1, 0), (1, 1), (2, 0))
_AG100_DISPLACEMENTS = TEST_DATA / 'Ag(100)/displacements'
_PTRH_INPUTS = TEST_DATA / 'stored_inputs/inputs/ptrh_100_1x3_o'

# Variant name -> (case as in conftest, DISPLACEMENTS or None). A
# DISPLACEMENTS of None means that the file of the case is used.
VARIANTS = {
    'ag100_z': ('ag100', _AG100_DISPLACEMENTS / 'DISPLACEMENTS_z'),
    'ag100_vib': ('ag100', _AG100_DISPLACEMENTS / 'DISPLACEMENTS_vib'),
    'ag100_z+vib': ('ag100', _AG100_DISPLACEMENTS / 'DISPLACEMENTS_z+vib'),
    'ag100_mixed': ('ag100_mixed', None),
    'telluride_stacking': ('telluride_stacking', None),
    'ptrh_100_1x3_o': ('ptrh_100_1x3_o', _PTRH_INPUTS / DISPLACEMENTS),
    }

PARAMS = [
    pytest.param((group, variant, name), id=f'{group}-{variant}-{name}')
    for group, variant, name in product(GROUPS, VARIANTS, SEARCH_FILES)
    ]


def _synthetic_exp_beams():
    """Return experimental beams, for cases without an EXPBEAMS."""
    beams = []
    for hk in EXP_BEAMS:
        beam = Beam(hk)
        beam.intens = dict.fromkeys(EXP_ENERGIES, 1.0)
        beams.append(beam)
    return beams


def _write_search_inputs(prepared, max_gen=None):
    """Write the search inputs of PreparedDeltas, return the folder."""
    rpars = prepared.rpars
    if max_gen:
        rpars.SEARCH_MAX_GEN = max_gen
    parts = [(prepared.slab, prepared.folder)]
    parts.extend((d.slab, d.workdir) for d in rpars.domainParams)
    for slab, folder in parts:
        slab.deltas_initialized = True
        for atom in slab:  # Name-only stubs of the delta files
            for name in atom.current_deltas:
                (folder / name).touch()
    rpars.SEARCH_POPULATION = SEARCH_POPULATION
    rpars.SEARCH_START = SEARCH_START
    if not rpars.expbeams:
        rpars.expbeams = _synthetic_exp_beams()
    with execute_in_dir(prepared.folder):
        iosearch.generateSearchInput(prepared.slab, rpars)
    return prepared.folder


@pytest.fixture(scope='session', name='search_folders')
def factory_search_folders(prepared_deltas):
    """Return a function that writes the search inputs (cached)."""
    cache = {}

    def _write(group, variant):
        key = group, variant
        if key not in cache:
            case, displacements = VARIANTS[variant]
            prepared = prepared_deltas(case, group, displacements)
            cache[key] = _write_search_inputs(prepared)
        return cache[key]
    return _write


@pytest.fixture(name='written_and_reference', params=PARAMS)
def fixture_written_and_reference(request, search_folders, data_path):
    """Return a written file and its reference."""
    group, variant, name = request.param
    reference = data_path / 'stored_inputs/search' / variant / group / name
    return search_folders(group, variant) / name, reference


def test_search_input(written_and_reference, check_reference):
    """Check a file written for a search."""
    written, reference = written_and_reference
    check_reference(written, reference)


def test_max_generations_limited(group, prepared_deltas):
    """Check that too many generations are cut to what fits."""
    case, displacements = VARIANTS['ag100_z']
    prepared = prepared_deltas(case, group, displacements)
    folder = _write_search_inputs(prepared, max_gen=HUGE_MAX_GEN)
    lines = (folder / 'search.steu').read_text().splitlines()
    generations = next(line for line in lines if 'desired number' in line)
    assert int(generations.split()[0]) == MAX_GEN_WRITTEN[group]
