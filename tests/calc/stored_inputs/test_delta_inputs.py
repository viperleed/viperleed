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

import copy
from itertools import product
import shutil

import pytest

from viperleed.calc.constants import DEFAULT_TENSORS
from viperleed.calc.files.displacements import readDISPLACEMENTS
from viperleed.calc.lib.context import execute_in_dir
from viperleed.calc.sections import deltas

from ...helpers import TEST_DATA
from .conftest import GROUPS

DELTA_INPUT = 'delta-input'
TENSOR_INDEX = 1  # As after a reference calculation
DISPLACEMENTS = 'DISPLACEMENTS'
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
def factory_delta_folders(initialized_case, tmp_path_factory):
    """Return a function that writes the delta inputs (cached).

    The function returns a {part: folder} dictionary.
    """
    cache = {}

    def _write(group, variant):
        key = group, variant
        if key in cache:
            return cache[key]
        case, _, displacements, _ = VARIANTS[variant]
        state = initialized_case(case, group)
        rpars, source = state.main
        folder = tmp_path_factory.mktemp('deltas')
        shutil.copytree(source, folder, dirs_exist_ok=True,
                        ignore=shutil.ignore_patterns('*.log'))
        slab, rpars = copy.deepcopy((next(iter(state.parts.values()))[0],
                                     rpars))
        folders = dict.fromkeys(state.parts, folder)
        for domain in rpars.domainParams:
            folders[domain.name] = folder / domain.workdir.relative_to(source)
            domain.workdir = folders[domain.name]
        if displacements:
            shutil.copy2(displacements, folder / DISPLACEMENTS)
        for rp in (rpars, *(d.rpars for d in rpars.domainParams)):
            rp.TENSOR_INDEX = TENSOR_INDEX
            if 1 not in rp.runHistory:  # Original state is the initial
                rp.runHistory.append(1)
        with execute_in_dir(folder):
            readDISPLACEMENTS(rpars)
            for part_folder in set(folders.values()):
                tensors = f'{DEFAULT_TENSORS}_{TENSOR_INDEX:03d}'
                (part_folder / DEFAULT_TENSORS / tensors).mkdir(parents=True)
            if rpars.domainParams:
                tasks = deltas._prepare_deltas_for_domains(rpars)
            else:
                tasks = deltas._prepare_deltas_for_one_domain(slab, rpars)
        written = _param_files(tasks[0])
        cache[key] = folders, {p.name for p in written}
        return cache[key]
    return _write


@pytest.fixture(name='written_and_reference', params=PARAMS)
def fixture_written_and_reference(request, delta_folders, data_path):
    """Return a written file and its reference."""
    group, variant, part, name = request.param
    reference = (data_path / 'stored_inputs/deltas'
                 / _reference_name(variant, part) / group / name)
    folders, _ = delta_folders(group, variant)
    return folders[part] / name, reference


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


def test_no_other_param_files(written_param_names):
    """Check that no PARAM file is written that has no reference."""
    variant, written = written_param_names
    *_, n_param = VARIANTS[variant]
    assert written == {_param_name(i) for i in range(n_param)}
