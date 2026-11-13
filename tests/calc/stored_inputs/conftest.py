"""Test configuration for tests.calc.stored_inputs.

Defines the version groups of the stored-input tests, and a
session-scoped factory fixture that brings each case (see
cases_stored_inputs.py) to the state it has right after an
INITIALIZATION, without a TensErLEED tree, a compiler, or any
execution of TensErLEED code. The inputs written by later sections
(reference calculation, deltas, search) can then be produced by calling
the writers on the slabs and Rparams of the case.

Fixtures
--------
pinned_clock (autouse)
    Pin the time seen by viperleed.calc.lib.time_utils, so that the
    timestamps written to the inputs are reproducible.
group
    Run of the tests, i.e., a version group (or a version of it)
    (parametrized).
initialized_case (factory)
    Return the InitializedCase of the StoredInputs of a case in a
    version group.
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

import pytest
from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.calc.classes.rparams.special.energy_range import TheoEnergies
from viperleed.calc.constants import DEFAULT_TENSORS
from viperleed.calc.constants import DEFAULT_WORK
from viperleed.calc.files.displacements import readDISPLACEMENTS
from viperleed.calc.files import beams as iobeams
from viperleed.calc.files import vibrocc
from viperleed.calc.lib import checksums
from viperleed.calc.lib import time_utils
from viperleed.calc.lib.context import execute_in_dir
from viperleed.calc.lib.fs_utils import copytree_exists_ok
from viperleed.calc.lib.time_utils import DateTimeFormat
from viperleed.calc.lib.version import Version
from viperleed.calc.run import run_calc
from viperleed.calc.sections import deltas
from viperleed.calc.sections import initialization
from viperleed.calc.sections import run_sections

from ..poscar_slabs import CasePOSCARSlabs
from .cases_stored_inputs import DISPLACEMENTS
from .cases_stored_inputs import PRESET_FILES
from .cases_stored_inputs import tagged_inputs

SYSTEM_NAME = 'stored-input'
# The clock is pinned to this time, so that the timestamps written to
# the inputs (format DateTimeFormat.FILE_SUFFIX) are '260101-000000'.
PINNED_TIME = time.struct_time((2026, 1, 1, 0, 0, 0, 3, 1, -1))
N_CORES = 1
THEO_ENERGIES = 50, 350, 3  # As in the PARAMETERS of the stored cases
TENSOR_INDEX = 1  # As after a reference calculation

# Runs of the tests: name of the run and the TensErLEED version used.
# The reference files are grouped by version group, i.e., by the range
# of TensErLEED versions that share them. Each group runs at the highest
# version it contains, and, if it contains more than one TensErLEED
# version with different code paths, also at the lowest one. The latter
# are named '<group>@<version>' and use the references of <group>.
# The newest group runs at the newest version known to the checksums:
# a newer version that changes the inputs fails these tests, to flag
# that the version groups need review.
GROUPS = {
    '1.6.0-1.6.1': '1.6.1',
    '1.7.3-': max(checksums.KNOWN_TL_VERSIONS, key=Version),
    '1.7.3-@1.7.4': '1.7.4',
    }


def reference_group(group):
    """Return the name of the reference folder of a run."""
    return group.split('@')[0]


def reference_name(name, domain):
    """Return the name of the reference folder of a domain of a case.

    Parameters
    ----------
    name : str
        The id of the case.
    domain : str or None
        The domain to return the reference folder of, None for
        single-surface cases (see InitializedCase.domains).

    Returns
    -------
    str
        `name` for single-surface cases, `name`_`domain` otherwise.
    """
    return name if domain is None else f'{name}_{domain}'


def group_case_domains(tag):
    """Yield every version group, case and domain of the tagged cases.

    Parameters
    ----------
    tag : str
        The tag of the cases (see cases_stored_inputs.py).

    Yields
    ------
    group : str
        A version group (see GROUPS).
    case_id : str
        The id of the case.
    inputs : StoredInputs
        The inputs of the case.
    domain : str or None
        A domain of the case, None for single-surface cases.
    """
    for group in GROUPS:
        for case_id, inputs in tagged_inputs(tag):
            for domain in inputs.domain_names:
                yield group, case_id, inputs, domain


def reference_params(tag, file_names):
    """Return the parameters of a test per reference file.

    Parameters
    ----------
    tag : str
        The tag of the cases to test (see cases_stored_inputs.py).
    file_names : callable
        Called with a group and StoredInputs, returns the names of the
        reference files of each domain.

    Returns
    -------
    list
        One pytest.param of (group, inputs, domain, reference name,
        file name) per group, case, domain and file name, with id
        '<group>-<reference name>-<file name>'. The reference name is
        that of reference_name.
    """
    return [
        pytest.param(group,
                     inputs,
                     domain,
                     ref_name,
                     name,
                     id=f'{group}-{ref_name}-{name}')
        for group, case_id, inputs, domain in group_case_domains(tag)
        for ref_name in (reference_name(case_id, domain),)
        for name in file_names(group, inputs)
        ]


_NOT_COPIED = {'Tensors'}  # Tensors are never read: no stored data
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
        Maps the name of each domain (see InitializedCase.domains) to
        the folder in which its delta inputs were written.
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
class InitializedCase:
    """Slabs and Rparams of a case after an initialization.

    Attributes
    ----------
    domains : dict
        Maps a name to a (slab, rpars, folder) tuple. Single-surface
        cases have one entry, with name None. Domain cases have one
        entry per domain, named after the domain. `folder` is the
        directory in which the INITIALIZATION ran, i.e., where the
        input files for later sections can be written.
    main : tuple
        The (rpars, folder) of the main calculation. For single-surface
        cases, this is the same as the Rparams and folder in `domains`.
    """

    domains: dict = field(default_factory=dict)
    main: tuple = ()


@fixture(autouse=True, name='pinned_clock')
def fixture_pinned_clock(mocker):
    """Make viperleed.calc.lib.time_utils see the time PINNED_TIME.

    The session-scoped factories run inside tests, hence they, too,
    see the pinned clock.
    """
    mocker.patch.object(
        time_utils,
        'time',
        SimpleNamespace(strftime=time.strftime,
                        localtime=lambda *_: PINNED_TIME,
                        gmtime=lambda *_: PINNED_TIME),
        )


@fixture(name='group')
@parametrize(run=tuple(GROUPS), ids=[f'group={run}' for run in GROUPS])
def fixture_group(run):
    """Return the name of a run of the tests (see GROUPS).

    A run is a version group or one more version of it, whose
    references are those of reference_group(run).
    """
    return run


@fixture(scope='session', name='initialized_case')
def factory_initialized_case(tmp_path_factory, data_path):
    """Return a function that initializes a case (cached)."""
    cache = {}

    def _initialize(inputs, group):
        tl_version = GROUPS[group]
        key = inputs.initialization, tl_version
        if key not in cache:
            tree = tmp_path_factory.mktemp('tree')
            # An empty folder is enough to find a TensErLEED version
            (tree / f'TensErLEED-v{tl_version}').mkdir()
            work = tmp_path_factory.mktemp('case')
            init = _init_preset_case if inputs.preset else _init_stored_case
            cache[key] = init(inputs, work, tree, tl_version, data_path)
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
    copytree_exists_ok(source,
                       folder,
                       ignore=shutil.ignore_patterns('*.log'))
    return copy.deepcopy((slab, rpars))


def _prepare_folder(folder, inputs, data_path):
    """Copy the stored inputs of a case to `folder`."""
    source = data_path / inputs.source
    for item in source.iterdir():
        if item.name in _NOT_COPIED:
            continue
        if item.is_dir() and item.name not in inputs.domains:
            continue
        if item.is_dir():
            copytree_exists_ok(
                item,
                folder / item.name,
                ignore=shutil.ignore_patterns(*_NOT_COPIED),
                )
        else:
            shutil.copy2(item, folder)
    for destination, extra_source in inputs.extra:
        shutil.copy2(data_path / extra_source, folder / destination)


def _init_stored_case(inputs, tmp, tree, tl_version, data_path):
    """Run an INITIALIZATION on StoredInputs from a folder."""
    work = tmp / DEFAULT_WORK
    work.mkdir()
    _prepare_folder(tmp, inputs, data_path)
    _prepare_folder(work, inputs, data_path)
    captured = {}
    with pytest.MonkeyPatch.context() as patcher:
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
    domains = {None: (slab, rpars, work)}
    if inputs.domains:
        domains = {dom.name: (dom.slab, dom.rpars, dom.workdir)
                   for dom in rpars.domainParams}
    for _, domain_rpars, _ in domains.values():
        domain_rpars.timestamp = rpars.timestamp
    return InitializedCase(domains, (rpars, work))


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


def _init_preset_case(inputs, tmp, tree, tl_version, data_path):
    """Run an INITIALIZATION on a slab from the POSCAR presets."""
    slab, rpars, *_ = getattr(CasePOSCARSlabs(), inputs.preset)()
    work = tmp / DEFAULT_WORK
    work.mkdir()
    rpars.paths.tensorleed = tree
    rpars.TL_VERSION = Version(tl_version)
    rpars.systemName = SYSTEM_NAME
    rpars.timestamp = DateTimeFormat.FILE_SUFFIX.now()
    rpars.N_CORES = N_CORES
    rpars.THEO_ENERGIES = TheoEnergies(*THEO_ENERGIES)
    with execute_in_dir(work):
        _write_synthetic_phaseshifts(slab, rpars, 'PHASESHIFTS')
        for name in PRESET_FILES:
            shutil.copy2(data_path / inputs.source / name, name)
        slab.full_update(rpars)  # chemical elements of ELEMENT_MIX
        vibrocc.readVIBROCC(rpars, slab)
        rpars.fileLoaded['VIBROCC'] = True
        slab.full_update(rpars)
        rpars.ivbeams = iobeams.readIVBEAMS()
        rpars.ivbeams_sorted = False
        rpars.fileLoaded['IVBEAMS'] = True
        rpars.updateDerivedParams()
        initialization.initialization(slab, rpars)
    return InitializedCase({None: (slab, rpars, work)}, (rpars, work))


@fixture(scope='session', name='prepared_deltas')
def factory_prepared_deltas(initialized_case, tmp_path_factory, data_path):
    """Return a function that prepares the delta inputs (not cached).

    The function takes the StoredInputs of a case and the version
    group, and returns a PreparedDeltas. Slabs and Rparams are copies,
    so that the cached state of the case is not modified. The Tensors
    folder only has to exist: its contents are never read, as the
    original-state files are those of the initialization, and no delta
    files exist yet.
    """
    def _prepare(inputs, group):
        state = initialized_case(inputs, group)
        rpars, source = state.main
        folder = tmp_path_factory.mktemp('deltas')
        main_slab, *_ = next(iter(state.domains.values()))
        slab, rpars = copy_case_state(source, folder, main_slab, rpars)
        folders = dict.fromkeys(state.domains, folder)
        for domain in rpars.domainParams:
            folders[domain.name] = folder / domain.workdir.relative_to(source)
            domain.workdir = folders[domain.name]
        if inputs.displacements:
            shutil.copy2(data_path / inputs.displacements,
                         folder / DISPLACEMENTS)
        for domain_rpars in (rpars, *(d.rpars for d in rpars.domainParams)):
            domain_rpars.TENSOR_INDEX = TENSOR_INDEX
            if 1 not in domain_rpars.runHistory:  # Original state: initial
                domain_rpars.runHistory.append(1)
        with execute_in_dir(folder):
            readDISPLACEMENTS(rpars)
            for domain_folder in set(folders.values()):
                tensors = f'{DEFAULT_TENSORS}_{TENSOR_INDEX:03d}'
                (domain_folder / DEFAULT_TENSORS / tensors).mkdir(parents=True)
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
