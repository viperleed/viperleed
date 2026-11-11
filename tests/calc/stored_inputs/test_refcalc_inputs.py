"""Tests for the input files of a reference calculation.

The files written by viperleed.calc.sections.refcalc and
viperleed.calc.files.iorefcalc are compared to stored references, for
each version group of TensErLEED (see conftest.py). Neither a
TensErLEED tree nor a compiler is needed: the execution is suppressed
after the input files are written, and PARAM is created directly.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-11'
__license__ = 'GPLv3+'

import pytest

from viperleed.calc.files import iorefcalc
from viperleed.calc.lib.context import execute_in_dir
from viperleed.calc.sections import refcalc

from ...helpers import TEST_DATA
from .conftest import GROUPS
from .conftest import PRESET_CASES
from .conftest import STORED_CASES
from .conftest import copy_case_state
from .conftest import reference_group
from .conftest import reference_name

# Case (as in conftest) -> parts. The case 'ag100_mixed' is not here:
# its VIBROCC offsets and DISPLACEMENTS only matter after a search, so
# its reference calculation inputs are identical to those of 'ag100'.
CASES = {
    'ag100': ('ag100',),
    'ag100_el_rename': ('ag100_el_rename',),
    'telluride_stacking': ('STA', 'FAU'),
    'ptrh_100_1x3_o': ('ptrh_100_1x3_o',),
    }
assert set(CASES) <= set(STORED_CASES) | set(PRESET_CASES)

PARAM_FILE = 'refcalc-PARAM'
FIN_FILE = 'refcalc-FIN'
PHASESHIFTS = 'PHASESHIFTS'
# PHASESHIFTS is an input that refcalc-FIN contains verbatim. It is
# stored as a placeholder, rather than copied to each reference.
PHASESHIFTS_PLACEHOLDER = (
    f'<{PHASESHIFTS}: verbatim copy of the work-dir {PHASESHIFTS}>\n'
    )
# TensErLEED v1.7.3 and later do not use an AUXBEAMS file
_FILES = ('AUXLATGEO', 'AUXNONSTRUCT', 'AUXGEO', 'BEAMLIST', FIN_FILE,
          PARAM_FILE)
FILES_BY_GROUP = {
    '1.6.0-1.6.1': (*_FILES, 'AUXBEAMS'),
    '1.7.3-': _FILES,
    }
assert set(FILES_BY_GROUP) == set(map(reference_group, GROUPS))


PARAMS = [
    pytest.param((group, case, part, name),
                 id=f'{group}-{reference_name(case, parts, part)}-{name}')
    for group in GROUPS
    for names in (FILES_BY_GROUP[reference_group(group)],)
    for case, parts in CASES.items()
    for part in parts
    for name in names
    ]


@pytest.fixture(scope='session', name='refcalc_folder')
def factory_refcalc_folder(initialized_case, tmp_path_factory):
    """Return a function that writes the refcalc inputs (cached)."""
    cache = {}

    def _write(group, case, part):
        key = group, case, part
        if key in cache:
            return cache[key]
        slab, rpars, source = initialized_case(case, group).parts[part]
        folder = tmp_path_factory.mktemp('refcalc')
        slab, rpars = copy_case_state(source, folder, slab, rpars)
        rpars.SUPPRESS_EXECUTION = True
        with execute_in_dir(folder):
            refcalc.refcalc(slab, rpars)
            lmax = rpars.LMAX.max
            param = iorefcalc.writePARAM(slab, rpars, lmax=lmax)
            (folder / PARAM_FILE).write_text(
                f'### PARAM file for LMAX = {lmax} ###\n\n{param}\n\n',
                encoding='utf-8',
                )
        cache[key] = folder
        return folder
    return _write


@pytest.fixture(name='written_and_reference', params=PARAMS)
def fixture_written_and_reference(request, refcalc_folder, data_path):
    """Return a written file and its reference.

    Preparing the inputs (once per case) happens in this fixture, so
    that the time of the test itself is that of the comparison.
    """
    group, case, part, name = request.param
    reference = (data_path / 'stored_inputs/refcalc'
                 / reference_name(case, CASES[case], part)
                 / reference_group(group) / name)
    return refcalc_folder(group, case, part) / name, reference


def _with_phaseshifts_placeholder(fin_file, tmp_path):
    """Return a copy of `fin_file` with PHASESHIFTS as a placeholder."""
    phaseshifts = (fin_file.parent / PHASESHIFTS).read_text(encoding='utf-8')
    if not phaseshifts.endswith('\n'):  # As in collectFIN
        phaseshifts += '\n'
    fin = fin_file.read_text(encoding='utf-8')
    assert fin.count(phaseshifts) == 1
    result = tmp_path / FIN_FILE
    result.write_text(fin.replace(phaseshifts, PHASESHIFTS_PLACEHOLDER),
                      encoding='utf-8')
    return result


def test_refcalc_input(written_and_reference, check_reference, tmp_path):
    """Check a file written for a reference calculation."""
    written, reference = written_and_reference
    if written.name == FIN_FILE:
        written = _with_phaseshifts_placeholder(written, tmp_path)
    check_reference(written, reference)


# The PHASESHIFTS of the preset cases are written by the case set-up
# itself, and have no stored input to compare with.
_STORED_PHASESHIFTS = [
    pytest.param((group, case, part), id=f'{group}-{case}-{part}')
    for group in GROUPS
    for case, parts in CASES.items()
    if case in STORED_CASES
    for part in parts
    ]


@pytest.fixture(name='phaseshifts_and_stored', params=_STORED_PHASESHIFTS)
def fixture_phaseshifts_and_stored(request, refcalc_folder):
    """Return the work-dir PHASESHIFTS and the stored input."""
    group, case, part = request.param
    source = TEST_DATA / STORED_CASES[case].source
    if STORED_CASES[case].domains:
        source /= part
    written = refcalc_folder(group, case, part) / PHASESHIFTS
    return written, source / PHASESHIFTS


def test_phaseshifts_is_stored_input(phaseshifts_and_stored):
    """Check that the PHASESHIFTS in refcalc-FIN is the stored input."""
    written, stored = phaseshifts_and_stored
    assert written.read_bytes() == stored.read_bytes()
