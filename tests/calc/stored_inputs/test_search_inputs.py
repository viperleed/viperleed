"""Tests for the input files of the search.

The PARAM, search.steu, restrict.f and searchpars.info files written by
viperleed.calc.files.iosearch.generateSearchInput are compared to
stored references, for each version group of TensErLEED (see
conftest.py). The search.steu is the one written before any search
has run. No TensErLEED tree, compiler, tensors, or deltas are needed:
the delta files only exist as names (empty stubs) that the atoms
register as their current deltas, and the execution is never started.

Everything that depends on the machine or on random numbers is pinned:
the number of cores, the system name and the timestamp. Most cases
also pin the population and use a deterministic starting
configuration. One case (CRANDOM_CASE) instead keeps the
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

import contextlib
import functools

from pytest_cases import case
from pytest_cases import fixture
from pytest_cases import get_current_cases
from pytest_cases import parametrize
from pytest_cases import parametrize_with_cases

from viperleed.calc.classes.beam import Beam
from viperleed.calc.classes.rparams import rparams
from viperleed.calc.files import iosearch
from viperleed.calc.lib.context import execute_in_dir

from ...helpers import private_random
from . import cases_stored_inputs
from .cases_stored_inputs import SEARCH
from .conftest import reference_group

SEARCH_FILES = 'PARAM', 'search.steu', 'restrict.f', 'searchpars.info'
SEARCH_POPULATION = 4
SEARCH_START = 'centered'
CRANDOM_CASE = 'ag100_z_crandom'
CRANDOM_SEED = 1
# The default population is rounded up to a multiple of the cores. A
# number of cores that does not divide it makes the rounding visible.
CRANDOM_N_CORES = 3
MAX_GEN_LINE = 'desired number'  # In the comment of the value
HUGE_MAX_GEN = 10**9  # Above what the input file format allows
MAX_GEN_WRITTEN = {'1.6.0-1.6.1': 999999, '1.7.3-': 9999999}
EXP_ENERGIES = tuple(range(50, 351))  # Synthetic experimental beams
EXP_BEAMS = ((1, 0), (1, 1), (2, 0))


@case(id=CRANDOM_CASE, tags=SEARCH)
def case_ag100_z_crandom():
    """Return the inputs of ag100_z, for the CRANDOM_CASE."""
    return cases_stored_inputs.case_ag100_z()


def _synthetic_exp_beams():
    """Return experimental beams, for cases without an EXPBEAMS."""
    beams = []
    for indices in EXP_BEAMS:
        beam = Beam(indices)
        beam.intens = dict.fromkeys(EXP_ENERGIES, 1.0)
        beams.append(beam)
    return beams


def _write_search_inputs(prepared, *, crandom=False):
    """Write the search inputs of PreparedDeltas, return the folder.

    Parameters
    ----------
    prepared : PreparedDeltas
        The prepared delta inputs.
    crandom : bool, optional
        Whether to keep the default population and 'crandom' start,
        with seeded random numbers. Default is False, i.e., pin both.

    Returns
    -------
    folder : Path
        Where the files were written.
    """
    rpars = prepared.rpars
    domains = [(prepared.slab, prepared.folder)]
    domains.extend((d.slab, d.workdir) for d in rpars.domainParams)
    for slab, folder in domains:
        slab.deltas_initialized = True
        for atom in slab:  # Name-only stubs of the delta files
            for name in atom.current_deltas:
                (folder / name).touch()
    if crandom:
        for domain_rpars in (rpars, *(d.rpars for d in rpars.domainParams)):
            domain_rpars.N_CORES = CRANDOM_N_CORES
    else:
        rpars.SEARCH_POPULATION = SEARCH_POPULATION
        rpars.SEARCH_START = SEARCH_START
    if not rpars.expbeams:
        rpars.expbeams = _synthetic_exp_beams()
    seeded = (
        private_random(rparams, CRANDOM_SEED)
        if crandom else
        contextlib.nullcontext()
        )
    with execute_in_dir(prepared.folder), seeded:
        iosearch.generateSearchInput(prepared.slab, rpars)
    return prepared.folder


@fixture(scope='session', name='search_folders')
def factory_search_folders(prepared_deltas):
    """Return a function that writes the search inputs (cached)."""
    @functools.lru_cache(maxsize=None)
    def _write(group, inputs, crandom):
        prepared = prepared_deltas(inputs, group)
        return _write_search_inputs(prepared, crandom=crandom)
    return _write


@fixture(name='written_and_reference')
@parametrize(name=SEARCH_FILES)
@parametrize_with_cases('inputs',
                        cases=(cases_stored_inputs, case_ag100_z_crandom),
                        has_tag=SEARCH)
# pylint: disable-next=too-many-arguments,too-many-positional-arguments
def fixture_written_and_reference(group, inputs, name,
                                  search_folders, data_path, request):
    """Return a written file and its reference."""
    case_id = get_current_cases(request)[request.fixturename]['inputs'].id
    reference = (data_path / 'stored_inputs/search' / case_id
                 / reference_group(group) / name)
    folder = search_folders(group, inputs, case_id == CRANDOM_CASE)
    return folder / name, reference


def test_search_input(written_and_reference, check_reference):
    """Check a file written for a search."""
    written, reference = written_and_reference
    check_reference(written, reference)


def test_max_generations_limited(group, prepared_deltas):
    """Check that too many generations are cut to what fits."""
    prepared = prepared_deltas(cases_stored_inputs.case_ag100_z(), group)
    prepared.rpars.SEARCH_MAX_GEN = HUGE_MAX_GEN
    folder = _write_search_inputs(prepared)
    lines = (folder / 'search.steu').read_text(encoding='utf-8').splitlines()
    generations = next(line for line in lines if MAX_GEN_LINE in line)
    expected = MAX_GEN_WRITTEN[reference_group(group)]
    assert int(generations.split()[0]) == expected
