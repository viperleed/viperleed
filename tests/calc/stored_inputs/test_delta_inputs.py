"""Tests for the input files of delta-amplitude calculations.

The delta-input file and the PARAM files written by
viperleed.calc.sections.deltas and viperleed.calc.files.iodeltas are
compared to stored references, for each version group of TensErLEED
(see conftest.py). No TensErLEED tree, compiler, tensors or deltas are
needed: the Tensors folder only has to exist (its contents are never
read, as the original-state files are those of the initialization), no
delta files exist yet (so all the requested ones are 'to calculate'),
and the execution is suppressed once the inputs are written.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-11'
__license__ = 'GPLv3+'

from itertools import product

import pytest

from ...helpers import TEST_DATA
from .conftest import DISPLACEMENTS
from .conftest import GROUPS
from .conftest import reference_group

DELTA_INPUT = 'delta-input'
_AG100_DISPLACEMENTS = TEST_DATA / 'Ag(100)/displacements'
_PTRH_INPUTS = TEST_DATA / 'stored_inputs/inputs/ptrh_100_1x3_o'

# Variant name -> (case as in conftest, parts, DISPLACEMENTS or None,
# number of PARAM files, i.e., different array dimensions, of each
# part). A DISPLACEMENTS of
# None means that the file of the case (in its main folder) is used.
VARIANTS = {
    'ag100_z': ('ag100', ('ag100',),
                _AG100_DISPLACEMENTS / 'DISPLACEMENTS_z', 1),
    'ag100_vib': ('ag100', ('ag100',),
                  _AG100_DISPLACEMENTS / 'DISPLACEMENTS_vib', 1),
    'ag100_z+vib': ('ag100', ('ag100',),
                    _AG100_DISPLACEMENTS / 'DISPLACEMENTS_z+vib', 1),
    'ag100_mixed': ('ag100_mixed', ('ag100_mixed',), None, 3),
    'telluride_stacking': ('telluride_stacking', ('STA', 'FAU'), None, 2),
    'ptrh_100_1x3_o': ('ptrh_100_1x3_o', ('ptrh_100_1x3_o',),
                       _PTRH_INPUTS / DISPLACEMENTS, 3),
    }


def _reference_name(variant, part):
    """Return the name of the reference folder of a variant and part."""
    parts = VARIANTS[variant][1]
    return variant if len(parts) == 1 else f'{variant}_{part}'


def _param_name(index):
    """Return the name of the stored PARAM file of a compile task."""
    return f'Delta_Compile_{index}-PARAM'


PARAMS = [
    pytest.param((group, variant, part, name),
                 id=f'{group}-{_reference_name(variant, part)}-{name}')
    for group in GROUPS
    for variant, (_, parts, _, n_param) in VARIANTS.items()
    for part in parts
    for name in (DELTA_INPUT, *(_param_name(i) for i in range(n_param)))
    ]


def _param_files(compile_tasks):
    """Write the PARAM of each compile task, return the paths."""
    paths = []
    for task in compile_tasks:
        path = task.exec_dir.parent / f'{task.foldername}-PARAM'
        path.write_text(task.param, encoding='utf-8')
        paths.append(path)
    return paths


@pytest.fixture(scope='session', name='delta_folders')
def factory_delta_folders(prepared_deltas):
    """Return a function that writes the delta inputs (cached).

    The function returns a (PreparedDeltas, {PARAM names}) tuple.
    """
    cache = {}

    def _write(group, variant):
        key = group, variant
        if key not in cache:
            case, _, displacements, _ = VARIANTS[variant]
            prepared = prepared_deltas(case, group, displacements)
            written = _param_files(prepared.tasks[0])
            cache[key] = prepared, {p.name for p in written}
        return cache[key]
    return _write


@pytest.fixture(name='written_and_reference', params=PARAMS)
def fixture_written_and_reference(request, delta_folders, data_path):
    """Return a written file and its reference."""
    group, variant, part, name = request.param
    reference = (data_path / 'stored_inputs/deltas'
                 / _reference_name(variant, part) / reference_group(group)
                 / name)
    prepared, _ = delta_folders(group, variant)
    return prepared.folders[part] / name, reference


def test_delta_input(written_and_reference, check_reference):
    """Check a file written for a delta-amplitudes calculation."""
    written, reference = written_and_reference
    check_reference(written, reference)


@pytest.fixture(name='written_param_names',
                params=list(product(GROUPS, VARIANTS)),
                ids='-'.join)
def fixture_written_param_names(request, delta_folders):
    """Return the names of the PARAM files written for a variant."""
    group, variant = request.param
    _, written = delta_folders(group, variant)
    return variant, written


@pytest.fixture(name='run_tasks', params=list(product(GROUPS, VARIANTS)),
                ids='-'.join)
def fixture_run_tasks(request, delta_folders):
    """Return the delta-run tasks of a variant."""
    prepared, _ = delta_folders(*request.param)
    return prepared.tasks[1]


def _read_with_newline(path):
    """Return the contents of a file, with a terminating newline."""
    contents = path.read_text(encoding='utf-8')
    return contents if contents.endswith('\n') else contents + '\n'


def test_piped_delta_input(run_tasks):
    """Check the input piped to the delta executables.

    It is the abbreviated input that is written to delta-input, with
    the AUXBEAMS and PHASESHIFTS tags replaced by the files that are in
    the work folder of the task.
    """
    assert run_tasks
    for task in run_tasks:
        auxbeams = _read_with_newline(task.io_dir / 'AUXBEAMS')
        phaseshifts = _read_with_newline(task.io_dir / 'PHASESHIFTS')
        expected = task.din_short.replace('[AUXBEAMS]\n', auxbeams)
        expected = expected.replace('[PHASESHIFTS]\n', phaseshifts)
        assert task.din == expected


def test_no_other_param_files(written_param_names):
    """Check that no PARAM file is written that has no reference."""
    variant, written = written_param_names
    *_, n_param = VARIANTS[variant]
    assert written == {_param_name(i) for i in range(n_param)}
