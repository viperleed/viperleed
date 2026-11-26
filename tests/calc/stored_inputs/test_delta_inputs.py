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

import functools

from pytest_cases import fixture
from pytest_cases import parametrize
from pytest_cases import parametrize_with_cases

from . import cases_stored_inputs
from .cases_stored_inputs import DELTAS
from .conftest import reference_group
from .conftest import reference_params

DELTA_INPUT = 'delta-input'


def _param_name(index):
    """Return the name of the stored PARAM file of a compile task."""
    return f'Delta_Compile_{index}-PARAM'


def _file_names(_, inputs):
    """Return the names of the reference files of each domain."""
    return (DELTA_INPUT,
            *(_param_name(i) for i in range(inputs.n_param)))


PARAMS = reference_params(DELTAS, _file_names)


def _param_files(compile_tasks):
    """Write the PARAM of each compile task, return the paths."""
    paths = []
    for task in compile_tasks:
        path = task.exec_dir.parent / f'{task.foldername}-PARAM'
        path.write_text(task.param, encoding='utf-8')
        paths.append(path)
    return paths


@fixture(scope='session', name='delta_folders')
def factory_delta_folders(prepared_deltas):
    """Return a function that writes the delta inputs (cached).

    The function returns a (PreparedDeltas, {PARAM names}) tuple.
    """
    @functools.lru_cache(maxsize=None)
    def _write(group, inputs):
        prepared = prepared_deltas(inputs, group)
        written = _param_files(prepared.compile_tasks)
        return prepared, {p.name for p in written}
    return _write


@fixture(name='written_and_reference')
@parametrize('group,inputs,domain,ref_name,name', PARAMS)
# pylint: disable-next=too-many-arguments,too-many-positional-arguments
def fixture_written_and_reference(group, inputs, domain, ref_name, name,
                                  delta_folders, data_path):
    """Return a written file and its reference."""
    reference = (data_path / 'stored_inputs/deltas' / ref_name
                 / reference_group(group) / name)
    prepared, _ = delta_folders(group, inputs)
    return prepared.folders[domain] / name, reference


def test_delta_input(written_and_reference, check_reference):
    """Check a file written for a delta-amplitudes calculation."""
    written, reference = written_and_reference
    check_reference(written, reference)


@fixture(name='written_deltas')
@parametrize_with_cases('inputs', cases=cases_stored_inputs, has_tag=DELTAS)
def fixture_written_deltas(group, inputs, delta_folders):
    """Return the StoredInputs of a case and its delta_folders."""
    prepared, written = delta_folders(group, inputs)
    return inputs, prepared, written


def _read_with_newline(path):
    """Return the contents of a file, with a terminating newline."""
    contents = path.read_text(encoding='utf-8')
    return contents if contents.endswith('\n') else contents + '\n'


def test_piped_delta_input(written_deltas):
    """Check the input piped to the delta executables.

    It is the abbreviated input that is written to delta-input, with
    the AUXBEAMS and PHASESHIFTS tags replaced by the files that are in
    the work folder of the task.
    """
    _, prepared, _ = written_deltas
    assert prepared.run_tasks
    for task in prepared.run_tasks:
        auxbeams = _read_with_newline(task.io_dir / 'AUXBEAMS')
        phaseshifts = _read_with_newline(task.io_dir / 'PHASESHIFTS')
        expected = task.din_short.replace('[AUXBEAMS]\n', auxbeams)
        expected = expected.replace('[PHASESHIFTS]\n', phaseshifts)
        assert task.din == expected


def test_no_other_param_files(written_deltas):
    """Check that no PARAM file is written that has no reference."""
    inputs, _, written = written_deltas
    assert written == {_param_name(i) for i in range(inputs.n_param)}
