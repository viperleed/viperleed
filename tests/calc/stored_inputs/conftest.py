"""Test configuration for tests.calc.stored_inputs.

Defines the cases of the stored-input tests, and a session-scoped
factory fixture that brings each case to the state it has right after
an INITIALIZATION, without a TensErLEED tree, a compiler, or any
execution of TensErLEED code. The inputs written by later sections
(reference calculation, deltas, search) can then be produced by calling
the writers on the slabs and Rparams of the case.

Fixtures
--------
group
    Run of the tests, i.e., a version group (or a version of it)
    (parametrized).
initialized_case (factory)
    Return the CaseState of a case in a version group.
prepared_deltas (factory)
    Prepare the delta-amplitude inputs of a case in a copy of its
    working folder.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-11'
__license__ = 'GPLv3+'

import copy
from dataclasses import dataclass
from dataclasses import field
from pathlib import Path
import shutil
from types import SimpleNamespace
import time
from typing import NamedTuple

import pytest
from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.calc.classes.rparams.special.energy_range import TheoEnergies
from viperleed.calc.constants import DEFAULT_TENSORS
from viperleed.calc.constants import DEFAULT_WORK
from viperleed.calc.files.displacements import readDISPLACEMENTS
from viperleed.calc.files import beams as iobeams
from viperleed.calc.files import vibrocc
from viperleed.calc.lib import time_utils
from viperleed.calc.lib.context import execute_in_dir
from viperleed.calc.lib.fs_utils import copytree_exists_ok
from viperleed.calc.lib.time_utils import DateTimeFormat
from viperleed.calc.lib.version import Version
from viperleed.calc.run import run_calc
from viperleed.calc.sections import deltas
from viperleed.calc.sections import initialization
from viperleed.calc.sections import run_sections

from ...helpers import TEST_DATA
from ..poscar_slabs import CasePOSCARSlabs

SYSTEM_NAME = 'stored-input'
# The clock is pinned to this time, so that the timestamps written to
# the inputs (format DateTimeFormat.FILE_SUFFIX) are '260101-000000'.
PINNED_TIME = time.struct_time((2026, 1, 1, 0, 0, 0, 3, 1, -1))
N_CORES = 1
THEO_ENERGIES = 50, 350, 3  # As in the PARAMETERS of the stored cases
TENSOR_INDEX = 1  # As after a reference calculation
DISPLACEMENTS = 'DISPLACEMENTS'
AG100_DISPLACEMENTS = TEST_DATA / 'Ag(100)/displacements'
PTRH_INPUTS = TEST_DATA / 'stored_inputs/inputs/ptrh_100_1x3_o'

# Runs of the tests: name of the run and the TensErLEED version used.
# The reference files are grouped by version group, i.e., by the range
# of TensErLEED versions that share them. Each group runs at the highest
# version it contains, and, if it contains more than one TensErLEED
# version with different code paths, also at the lowest one. The latter
# are named '<group>@<version>' and use the references of <group>.
GROUPS = {
    '1.6.0-1.6.1': '1.6.1',
    '1.7.3-': '2.1.0',
    '1.7.3-@1.7.4': '1.7.4',
    }


def reference_group(group):
    """Return the name of the reference folder of a run."""
    return group.split('@')[0]


def reference_name(name, parts, part):
    """Return the name of the reference folder of a part of a case.

    Parameters
    ----------
    name : str
        The name of the case or variant.
    parts : Sequence
        The names of all parts of the case.
    part : str
        The part to return the reference folder of.

    Returns
    -------
    str
        `name` for single-part cases, `name`_`part` otherwise.
    """
    return name if len(parts) == 1 else f'{name}_{part}'


@dataclass
class StoredCase:
    """A case with stored inputs.

    Attributes
    ----------
    source : str
        A folder in the test data. All files in it are copied to the
        folders of the calculation.
    extra : dict
        {destination: source} pairs of files that are added on top.
    domains : tuple
        For domain cases, the sub-folders of `source` to copy.
    """

    source: str
    extra: dict = field(default_factory=dict)
    domains: tuple = ()


STORED_CASES = {
    'ag100': StoredCase('Ag(100)/initialization'),
    'ag100_mixed': StoredCase(
        'Ag(100)/initialization',
        extra={'VIBROCC': 'Ag(100)/mergeDisp/VIBROCC',
               'DISPLACEMENTS': 'Ag(100)/mergeDisp/DISPLACEMENTS_mixed'},
        ),
    'ag100_el_rename': StoredCase('Ag(100)_el_rename/initialization'),
    'telluride_stacking': StoredCase('domains/telluride_stacking/deltas',
                                     domains=('STA', 'FAU')),
    }


class Variant(NamedTuple):
    """The case and DISPLACEMENTS of a variant of the tests.

    Attributes
    ----------
    case : str
        The name of the case, as in STORED_CASES or PRESET_CASES.
    displacements : Path or None
        The DISPLACEMENTS file to use. None means that the file of the
        case (in its main folder) is used.
    """

    case: str
    displacements: object


# Variant name -> Variant. Shared by the tests of the delta and of the
# search inputs.
VARIANTS = {
    'ag100_z': Variant('ag100', AG100_DISPLACEMENTS / 'DISPLACEMENTS_z'),
    'ag100_vib': Variant('ag100', AG100_DISPLACEMENTS / 'DISPLACEMENTS_vib'),
    'ag100_z+vib': Variant('ag100',
                           AG100_DISPLACEMENTS / 'DISPLACEMENTS_z+vib'),
    'ag100_mixed': Variant('ag100_mixed', None),
    'telluride_stacking': Variant('telluride_stacking', None),
    'ptrh_100_1x3_o': Variant('ptrh_100_1x3_o', PTRH_INPUTS / DISPLACEMENTS),
    }
_NOT_COPIED = {'Tensors'}  # Tensors are never read: no stored data
PRESET_CASES = {  # POSCAR_SLABS case methods, with synthetic PHASESHIFTS
    'ptrh_100_1x3_o': 'case_poscar_ptrh_100_1x3_o',
    }
# Synthetic PHASESHIFTS of the preset cases: header (Rundgren
# parameters of the muffin-tin potential, copied from a stored file),
# energies (Hartree) that cover the theoretical energy range, number
# of L values.
_PS_HEADER = '{:3d}   -10.85    0.00  -67.83   13.72      synthetic'
_PS_ENERGIES = (0.1, 2.0, 5.0, 10.0, 15.0)
_PS_N_L = 18


@dataclass
class PreparedDeltas:
    """Result of preparing the delta inputs of a case.

    Attributes
    ----------
    folder : Path
        The working folder of the (main) calculation.
    folders : dict
        Maps the name of each part (see CaseState.parts) to the folder
        in which its delta inputs were written.
    slab : Slab
        The slab of the (main) calculation, with the delta files
        registered in the `current_deltas` of its atoms.
    rpars : Rparams
        The (main) run parameters.
    compile_tasks : Sequence
        The tasks that compile the delta executables.
    run_tasks : Sequence
        The tasks that run the delta executables.
    """

    folder: object
    folders: dict
    slab: object
    rpars: object
    compile_tasks: object
    run_tasks: object


@dataclass
class CaseState:
    """Slabs and Rparams of a case after an initialization.

    Attributes
    ----------
    parts : dict
        Maps a name to a (slab, rpars, folder) tuple. Single-surface
        cases have one entry, named after the case. Domain cases have
        one entry per domain, named after the domain. `folder` is the
        directory in which the INITIALIZATION ran, i.e., where the
        input files for later sections can be written.
    main : tuple
        The (rpars, folder) of the main calculation. For single-surface
        cases, this is the same as the Rparams and folder in `parts`.
    """

    parts: dict = field(default_factory=dict)
    main: tuple = ()


@fixture(name='group')
@parametrize(run=tuple(GROUPS), ids=[f'group={run}' for run in GROUPS])
def fixture_group(run):
    """Return the name of a run of the tests (see GROUPS).

    A run is a version group or one more version of it, whose
    references are those of reference_group(run).
    """
    return run


@fixture(scope='session', name='initialized_case')
def factory_initialized_case(tmp_path_factory):
    """Return a function that initializes a case (cached)."""
    cache = {}

    def _initialize(case, group):
        tl_version = GROUPS[group]
        key = case, tl_version
        if key not in cache:
            tree = tmp_path_factory.mktemp('tree')
            # An empty folder is enough to find a TensErLEED version
            (tree / f'TensErLEED-v{tl_version}').mkdir()
            work = tmp_path_factory.mktemp(case)
            if case in PRESET_CASES:
                cache[key] = _init_preset_case(case, work, tree, tl_version)
            else:
                cache[key] = _init_stored_case(case, work, tree, tl_version)
        return cache[key]
    return _initialize


def copy_case_state(source, folder, slab, rpars):
    """Copy a working folder, and return copies of a slab and Rparams.

    Parameters
    ----------
    source : Path
        The working folder to copy. Log files are not copied.
    folder : Path
        The (existing) destination folder.
    slab : Slab
        The slab to copy.
    rpars : Rparams
        The run parameters to copy.

    Returns
    -------
    slab, rpars : tuple
        Deep copies of the arguments.
    """
    copytree_exists_ok(source, folder,
                       ignore=shutil.ignore_patterns('*.log'))
    return copy.deepcopy((slab, rpars))


def _prepare_folder(folder, spec):
    """Copy the stored inputs of a case to `folder`."""
    source = TEST_DATA / spec.source
    for item in source.iterdir():
        if item.name in _NOT_COPIED:
            continue
        if item.is_dir() and item.name not in spec.domains:
            continue
        if item.is_dir():
            copytree_exists_ok(
                item,
                folder / item.name,
                ignore=shutil.ignore_patterns(*_NOT_COPIED),
                )
        else:
            shutil.copy2(item, folder)
    for destination, extra_source in spec.extra.items():
        shutil.copy2(TEST_DATA / extra_source, folder / destination)


def pin_clock(patcher):
    """Make viperleed.calc.lib.time_utils see a fixed time."""
    patcher.setattr(
        time_utils,
        'time',
        SimpleNamespace(strftime=time.strftime,
                        localtime=lambda *_: PINNED_TIME,
                        gmtime=lambda *_: PINNED_TIME),
        )


def _init_stored_case(case, tmp, tree, tl_version):
    """Run an INITIALIZATION on the stored inputs of `case`."""
    spec = STORED_CASES[case]
    work = tmp / DEFAULT_WORK
    work.mkdir()
    _prepare_folder(tmp, spec)
    _prepare_folder(work, spec)
    captured = {}
    with pytest.MonkeyPatch.context() as patcher:
        pin_clock(patcher)
        _record_and_keep_files(patcher, captured)
        with execute_in_dir(work):
            failed, _ = run_calc(
                source=tree,
                preset_params={'RUN': [0],
                               'TL_VERSION': Version(tl_version),
                               'N_CORES': N_CORES},
                system_name=SYSTEM_NAME,
                home=tmp,
                )
    assert not failed
    rpars, slab = captured['args']
    parts = {case: (slab, rpars, work)}
    if spec.domains:
        parts = {dom.name: (dom.slab, dom.rpars, dom.workdir)
                 for dom in rpars.domainParams}
    for _, part_rpars, _ in parts.values():
        part_rpars.timestamp = rpars.timestamp
    return CaseState(parts, (rpars, work))


def _record_and_keep_files(patcher, captured):
    """Make run_calc hand over its slab/Rparams, and not clean up."""
    original = run_sections.section_loop

    def _section_loop(rpars, slab):
        captured['args'] = rpars, slab
        return original(rpars, slab)

    patcher.setattr('viperleed.calc.run.section_loop', _section_loop)
    patcher.setattr(run_sections, 'cleanup', lambda *_: None)


def _write_synthetic_phaseshifts(slab, rpars, filename):
    """Write a deterministic PHASESHIFTS file that fits `slab`."""
    n_blocks = sum(
        len(rpars.ELEMENT_MIX.get(el, [el]))
        * len([s for s in slab.sitelist if s.el == el])
        for el in slab.elements
        )
    lines = [_PS_HEADER.format(n_blocks)]
    for i_energy, energy in enumerate(_PS_ENERGIES):
        lines.append(f'{energy:7.4f}')
        for i_block in range(n_blocks):
            values = [round((0.3*(i_block + 1) + 0.05*i_energy)/(n_l + 1), 4)
                      for n_l in range(_PS_N_L)]
            for start in range(0, _PS_N_L, 10):
                chunk = values[start:start + 10]
                lines.append(''.join(f'{v:7.4f}' for v in chunk))
    Path(filename).write_text('\n'.join(lines) + '\n', encoding='utf-8')


def _init_preset_case(case, tmp, tree, tl_version):
    """Run an INITIALIZATION on a slab from the POSCAR presets."""
    slab, rpars, *_ = getattr(CasePOSCARSlabs(), PRESET_CASES[case])()
    work = tmp / DEFAULT_WORK
    work.mkdir()
    rpars.paths.tensorleed = tree
    rpars.TL_VERSION = Version(tl_version)
    rpars.systemName = SYSTEM_NAME
    with pytest.MonkeyPatch.context() as patcher:
        pin_clock(patcher)
        rpars.timestamp = DateTimeFormat.FILE_SUFFIX.now()
    rpars.N_CORES = N_CORES
    rpars.THEO_ENERGIES = TheoEnergies(*THEO_ENERGIES)
    with execute_in_dir(work):
        _write_synthetic_phaseshifts(slab, rpars, 'PHASESHIFTS')
        for name in ('IVBEAMS', 'VIBROCC'):
            shutil.copy2(TEST_DATA / 'stored_inputs/inputs' / case / name,
                         name)
        slab.full_update(rpars)  # chemical elements of ELEMENT_MIX
        vibrocc.readVIBROCC(rpars, slab)
        rpars.fileLoaded['VIBROCC'] = True
        slab.full_update(rpars)
        rpars.ivbeams = iobeams.readIVBEAMS()
        rpars.ivbeams_sorted = False
        rpars.fileLoaded['IVBEAMS'] = True
        rpars.updateDerivedParams()
        initialization.initialization(slab, rpars)
    return CaseState({case: (slab, rpars, work)}, (rpars, work))


@fixture(scope='session', name='prepared_deltas')
def factory_prepared_deltas(initialized_case, tmp_path_factory):
    """Return a function that prepares the delta inputs (not cached).

    The function takes the case, the version group, and the path to
    a DISPLACEMENTS file to use (None to use the one of the case), and
    returns a PreparedDeltas. Slabs and Rparams are copies, so that
    the cached state of the case is not modified. The Tensors folder
    only has to exist: its contents are never read, as the
    original-state files are those of the initialization, and no delta
    files exist yet.
    """
    def _prepare(case, group, displacements):
        state = initialized_case(case, group)
        rpars, source = state.main
        folder = tmp_path_factory.mktemp('deltas')
        main_slab, *_ = next(iter(state.parts.values()))
        slab, rpars = copy_case_state(source, folder, main_slab, rpars)
        folders = dict.fromkeys(state.parts, folder)
        for domain in rpars.domainParams:
            folders[domain.name] = folder / domain.workdir.relative_to(source)
            domain.workdir = folders[domain.name]
        if displacements:
            shutil.copy2(displacements, folder / DISPLACEMENTS)
        for part_rpars in (rpars, *(d.rpars for d in rpars.domainParams)):
            part_rpars.TENSOR_INDEX = TENSOR_INDEX
            if 1 not in part_rpars.runHistory:  # Original state: initial
                part_rpars.runHistory.append(1)
        with execute_in_dir(folder):
            readDISPLACEMENTS(rpars)
            for part_folder in set(folders.values()):
                tensors = f'{DEFAULT_TENSORS}_{TENSOR_INDEX:03d}'
                (part_folder / DEFAULT_TENSORS / tensors).mkdir(parents=True)
            # There is no public entry point to the step under test
            # pylint: disable=protected-access
            tasks = (
                deltas._prepare_deltas_for_domains(rpars)
                if rpars.domainParams else
                deltas._prepare_deltas_for_one_domain(slab, rpars)
                )
            # pylint: enable=protected-access
        return PreparedDeltas(folder, folders, slab, rpars, *tasks)
    return _prepare
