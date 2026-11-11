"""Tests for the input files of the search.

The PARAM, search.steu, restrict.f and searchpars.info files written by
viperleed.calc.files.iosearch.generateSearchInput are compared to
stored references, for each version group of TensErLEED (see
conftest.py). The search.steu is the one written before any search
has run. No TensErLEED tree, compiler, tensors, or deltas are needed:
the delta files only exist as names (empty stubs) that the atoms
register as their current deltas, and the execution is never started.

Everything that depends on the machine or on random numbers is pinned:
the number of cores, the system name and the timestamp. Most variants
also pin the population and use a deterministic starting
configuration. One variant (CRANDOM_VARIANT) instead keeps the
default 'crandom' start and the default population (computed from the
number of search parameters and cores), with the random numbers drawn
by the 'crandom' start seeded. The seed of the random-number generator
written to search.steu is a constant (0, i.e., use the system time)
and needs no patching.
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
from viperleed.calc.classes.rparams import rparams
from viperleed.calc.files import iosearch
from viperleed.calc.lib.context import execute_in_dir

from ...helpers import TEST_DATA
from ...helpers import private_random
from .conftest import DISPLACEMENTS
from .conftest import GROUPS
from .conftest import reference_group

SEARCH_FILES = 'PARAM', 'search.steu', 'restrict.f', 'searchpars.info'
SEARCH_POPULATION = 4
SEARCH_START = 'centered'
CRANDOM_VARIANT = 'ag100_z_crandom'
CRANDOM_SEED = 1
# The default population is rounded up to a multiple of the cores. A
# number of cores that does not divide it makes the rounding visible.
CRANDOM_N_CORES = 3
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
    CRANDOM_VARIANT: ('ag100', _AG100_DISPLACEMENTS / 'DISPLACEMENTS_z'),
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


def _write_search_inputs(prepared, max_gen=None, crandom=False):
    """Write the search inputs of PreparedDeltas, return the folder.

    Parameters
    ----------
    prepared : PreparedDeltas
        The prepared delta inputs.
    max_gen : int, optional
        Value for SEARCH_MAX_GEN. Default is None, i.e., leave it.
    crandom : bool, optional
        Whether to keep the default population and 'crandom' start,
        with seeded random numbers. Default is False, i.e., pin both.

    Returns
    -------
    folder : Path
        Where the files were written.
    """
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
    if crandom:
        for rp in (rpars, *(d.rpars for d in rpars.domainParams)):
            rp.N_CORES = CRANDOM_N_CORES
    else:
        rpars.SEARCH_POPULATION = SEARCH_POPULATION
        rpars.SEARCH_START = SEARCH_START
    if not rpars.expbeams:
        rpars.expbeams = _synthetic_exp_beams()
    with execute_in_dir(prepared.folder), private_random(rparams,
                                                         CRANDOM_SEED):
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
            cache[key] = _write_search_inputs(
                prepared,
                crandom=variant == CRANDOM_VARIANT,
                )
        return cache[key]
    return _write


@pytest.fixture(name='written_and_reference', params=PARAMS)
def fixture_written_and_reference(request, search_folders, data_path):
    """Return a written file and its reference."""
    group, variant, name = request.param
    reference = (data_path / 'stored_inputs/search' / variant
                 / reference_group(group) / name)
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
    expected = MAX_GEN_WRITTEN[reference_group(group)]
    assert int(generations.split()[0]) == expected
