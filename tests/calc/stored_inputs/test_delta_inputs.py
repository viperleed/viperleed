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
from typing import NamedTuple

import pytest

from .conftest import GROUPS
from .conftest import VARIANTS
from .conftest import reference_group
from .conftest import reference_name

DELTA_INPUT = 'delta-input'


class Layout(NamedTuple):
    """The files written for a variant.

    Attributes
    ----------
    parts : tuple
        The names of the parts of the case, as in CaseState.
    n_param : int
        The number of PARAM files, i.e., of different array dimensions,
        of each part.
    """

    parts: tuple
    n_param: int


LAYOUTS = {
    'ag100_z': Layout(('ag100',), 1),
    'ag100_vib': Layout(('ag100',), 1),
    'ag100_z+vib': Layout(('ag100',), 1),
    'ag100_mixed': Layout(('ag100_mixed',), 3),
    'telluride_stacking': Layout(('STA', 'FAU'), 2),
    'ptrh_100_1x3_o': Layout(('ptrh_100_1x3_o',), 3),
    }
assert LAYOUTS.keys() == VARIANTS.keys()


def _param_name(index):
    """Return the name of the stored PARAM file of a compile task."""
    return f'Delta_Compile_{index}-PARAM'


PARAMS = [
    pytest.param((group, variant, part, name),
                 id=f'{group}-{reference_name(variant, parts, part)}-{name}')
    for group in GROUPS
    for variant, (parts, n_param) in LAYOUTS.items()
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
            case, displacements = VARIANTS[variant]
            prepared = prepared_deltas(case, group, displacements)
            written = _param_files(prepared.tasks[0])
            cache[key] = prepared, {p.name for p in written}
        return cache[key]
    return _write


@pytest.fixture(name='written_and_reference', params=PARAMS)
def fixture_written_and_reference(request, delta_folders, data_path):
    """Return a written file and its reference."""
    group, variant, part, name = request.param
    parts = LAYOUTS[variant].parts
    reference = (data_path / 'stored_inputs/deltas'
                 / reference_name(variant, parts, part)
                 / reference_group(group) / name)
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
    n_param = LAYOUTS[variant].n_param
    assert written == {_param_name(i) for i in range(n_param)}
