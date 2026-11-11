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

import copy
import shutil

import pytest

from viperleed.calc.files import iorefcalc
from viperleed.calc.lib.context import execute_in_dir
from viperleed.calc.sections import refcalc

from .conftest import GROUPS
from .conftest import PRESET_CASES
from .conftest import STORED_CASES

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
# TensErLEED v1.7.3 and later do not use an AUXBEAMS file
_FILES = ('AUXLATGEO', 'AUXNONSTRUCT', 'AUXGEO', 'BEAMLIST', FIN_FILE,
          PARAM_FILE)
FILES_BY_GROUP = {
    '1.6.0-1.6.1': (*_FILES, 'AUXBEAMS'),
    '1.7.3-': _FILES,
    }
assert set(FILES_BY_GROUP) == set(GROUPS)


def _reference_name(case, part):
    """Return the name of the reference folder for a case and part."""
    return case if len(CASES[case]) == 1 else f'{case}_{part}'


PARAMS = [
    pytest.param((group, case, part, name),
                 id=f'{group}-{_reference_name(case, part)}-{name}')
    for group, names in FILES_BY_GROUP.items()
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
        shutil.copytree(source, folder, dirs_exist_ok=True,
                        ignore=shutil.ignore_patterns('*.log'))
        slab, rpars = copy.deepcopy((slab, rpars))
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
                 / _reference_name(case, part) / group / name)
    return refcalc_folder(group, case, part) / name, reference


def test_refcalc_input(written_and_reference, check_reference):
    """Check a file written for a reference calculation."""
    check_reference(*written_and_reference)
