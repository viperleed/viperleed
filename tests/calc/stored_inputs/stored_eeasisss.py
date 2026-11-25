"""Run psgen with eeasisss replaced by its stored raw outputs.

The PHASESHIFTS calculation of some slabs takes long with the real
eeasisss. For these slabs, the raw outputs of eeasisss (PS.r.* files)
are stored in the test data, in an archive per slab, together with
the eeasisss input from which they were produced and the PHASESHIFTS
that viperleed.calc.psgen reads from them. A mocked run of psgen
'produces' the stored raw outputs instead of running eeasisss, and
checks that the eeasisss input written is the stored one (so that the
stored outputs are not stale) and that the PHASESHIFTS read are the
stored ones. No TensErLEED tree is needed. The energy range (ENERGIES)
is narrower than the default one to keep the stored data small.

Regenerate the stored data of all slabs of MOCKED_CASES (or of some,
by giving their case ids) with the real eeasisss, from the root of
the repository:

    VIPERLEED_TENSORLEED=<tree> python -m \\
        tests.calc.stored_inputs.stored_eeasisss [case_id ...]

and review the changes (e.g., with git) before committing them.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-25'
__license__ = 'GPLv3+'

from collections import namedtuple
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
from types import SimpleNamespace
from unittest import mock
import zipfile

from pytest_cases import case

from viperleed.calc import psgen
from viperleed.calc.files.phaseshifts import writePHASESHIFTS
from viperleed.calc.files.tenserleed import get_tensorleed_path
from viperleed.calc.lib.context import execute_in_dir
from viperleed.calc.lib.periodic_table import PERIODIC_TABLE

from ...helpers import PSGEN_SEED
from ...helpers import TEST_DATA
from ...helpers import private_random
from ...harness.reference_files import fail_with_diff
from ..poscar_slabs import CasePOSCARSlabs as POSCARSlabs
from ..poscar_slabs import POSCARS_WITH_LITTLE_SYMMETRY_INFO

STORED_PATH = TEST_DATA / 'stored_inputs/eeasisss_output'
ENERGIES = 20, 30, 10  # THEO_ENERGIES of the stored runs
SYSTEM_NAME = 'stored-input'
TIMESTAMP = '260101-000000'  # As in the other stored inputs
INPUT_NAME = 'eeasisss-input'
REFERENCE = 'PHASESHIFTS'
RAW_OUTPUT_ARCHIVE = 'PS.r.zip'
RAW_OUTPUT_GLOB = 'PS.r.*'
TREE_NAME = 'TensErLEED'
EXECUTABLE = 'eeasisss'  # Without WINDOWS_SUFFIX
EXCHANGE_CORRELATION = 'seSernelius'
ATOM_DENSITIES = 'atom_density_files'
WORK_NAME = 'work'
WINDOWS = 'nt'  # Value of os.name on Windows, as psgen checks it
WINDOWS_SUFFIX = '.exe'  # Of EXECUTABLE on Windows, as psgen adds it
POSIX_SEP = '/'  # Path separator in the stored inputs
_ENCODING = 'utf-8'
_ZIP_DATE = (1980, 1, 1, 0, 0, 0)  # Earliest allowed: reproducible

# Case ids (as in the run_phaseshift fixture) of the slabs whose
# eeasisss outputs are stored. Each id is also the name of the folder
# of its stored data. The ids of the mocked cases are these with
# MOCKED_SUFFIX appended, so that they differ in every '-'-separated
# part (see tests/harness/tiers.py) from the ids of the real runs.
MOCKED_SUFFIX = '_mocked'
_AL2O3 = 'poscar-Al2O3_NiAl(111)_cHole_20061025'
_FE3O4 = 'poscar_fe3o4_001_cod'
_SB_SI = 'poscar_sb_si_111'
_STO = 'poscar_sto110_4x1'
_TIO2 = 'poscar_tio2_supercell'
MOCKED_CASES = (_AL2O3, _FE3O4, _SB_SI, _STO, _TIO2)

MockedSlab = namedtuple(  # noqa: PYI024  # As the other test tuples
    'MockedSlab',
    ('slab', 'rpars', 'info', 'stored'),
    )
MockedSlab.__doc__ = 'A slab whose eeasisss outputs are in `stored`.'


def _al2o3_case():
    """Return a slab, an Rparams and a TestInfo for Al2O3/NiAl(111)."""
    poscar_name = 'POSCAR_' + _AL2O3.split('-', 1)[1]
    info = next(info for info in POSCARS_WITH_LITTLE_SYMMETRY_INFO
                if info.poscar.name == poscar_name)
    return POSCARSlabs().case_poscar(info)


# How to make the slab, Rparams and TestInfo of each case
_MAKE_CASE = {
    _AL2O3: _al2o3_case,
    _FE3O4: POSCARSlabs().case_poscar_fe3o4_001_cod,
    _SB_SI: POSCARSlabs().case_poscar_sb_si_111,
    _STO: POSCARSlabs().case_poscar_sto110_4x1,
    _TIO2: POSCARSlabs().case_poscar_tio2_supercell,
    }


def _mocked(case_id):
    """Return a MockedSlab for one of MOCKED_CASES."""
    return MockedSlab(*_MAKE_CASE[case_id](), STORED_PATH / case_id)


def _mocked_case(case_id):
    """Return a case function giving the MockedSlab of `case_id`."""
    @case(id=case_id + MOCKED_SUFFIX)
    def _case():
        return _mocked(case_id)
    _case.__doc__ = f'Return {case_id} with its stored eeasisss outputs.'
    return _case


# One case function per slab of MOCKED_CASES, for the cases argument
# of pytest_cases. A case function parametrized over MOCKED_CASES
# would prepend its own name to the case ids.
CASES_MOCKED_EEASISSS = tuple(_mocked_case(case_id)
                              for case_id in MOCKED_CASES)


def _executable_name():
    """Return the name of the eeasisss executable that psgen runs."""
    if os.name == WINDOWS:
        return EXECUTABLE + WINDOWS_SUFFIX
    return EXECUTABLE


def make_stub_tree(parent):
    """Return a fake TensErLEED directory with stub eeasisss/data.

    Parameters
    ----------
    parent : Path
        The directory in which the TREE_NAME directory is created.

    Returns
    -------
    tree : Path
        The fake TensErLEED tree.
    """
    tree = parent / TREE_NAME
    tree.mkdir()
    (tree / _executable_name()).touch()
    (tree / EXCHANGE_CORRELATION).touch()
    for element in PERIODIC_TABLE:
        element_dir = tree / ATOM_DENSITIES / element
        element_dir.mkdir(parents=True)
        (element_dir / f'chgden{element}').touch()
    return tree


def copy_real_tree(parent, real_tree):
    """Return a TREE_NAME folder in `parent` with files of `real_tree`.

    psgen writes the paths of the files of its tree into the eeasisss
    input relative to the folder in which it runs. A copy in the place
    of the stub tree of make_stub_tree makes psgen write the same input
    while running the real eeasisss. The files are copied rather than
    linked, as links need privileges on Windows.

    Parameters
    ----------
    parent : Path
        The directory in which the TREE_NAME directory is created.
    real_tree : Path
        The TensErLEED tree with a compiled eeasisss.

    Returns
    -------
    tree : Path
        The copy of the files of `real_tree` that psgen uses.
    """
    tree = parent / TREE_NAME
    tree.mkdir()
    for name in (_executable_name(), EXCHANGE_CORRELATION):
        shutil.copy2(real_tree / name, tree / name)
    shutil.copytree(real_tree / ATOM_DENSITIES, tree / ATOM_DENSITIES)
    return tree


def use_forward_slashes(written):
    """Replace OS-native path separators in `written` with '/'.

    psgen writes the paths of the charge-density files with the
    separators of the OS (backslashes on Windows, where eeasisss
    accepts them). The stored inputs use forward slashes, so that
    they hold on every OS.

    Parameters
    ----------
    written : Path
        The file to edit in place.
    """
    text = written.read_text(encoding=_ENCODING)
    written.write_text(text.replace(os.sep, POSIX_SEP), encoding=_ENCODING)


def run_psgen(slab, rpars, tree, tmp_path, energies=ENERGIES):
    """Run psgen in a WORK_NAME folder of `tmp_path`, with `tree`.

    The eeasisss input is written with forward slashes (see
    use_forward_slashes), and the PHASESHIFTS read into the folder.

    Parameters
    ----------
    slab : Slab
        The slab whose phase shifts are calculated.
    rpars : Rparams
        Its parameters. Edited in place.
    tree : Path
        The TensErLEED tree that psgen uses.
    tmp_path : Path
        An empty directory in which the WORK_NAME folder is created.
    energies : Sequence, optional
        The THEO_ENERGIES. Default is ENERGIES, those of the stored
        eeasisss outputs of MOCKED_CASES.

    Returns
    -------
    work : Path
        The folder in which psgen ran.
    firstline : str
        The first line of the PHASESHIFTS.
    phaseshifts : list
        The phase shifts read from the eeasisss outputs.
    """
    rpars.paths.tensorleed = tree
    rpars.systemName = SYSTEM_NAME
    rpars.timestamp = TIMESTAMP
    rpars.initTheoEnergies()
    rpars.THEO_ENERGIES = rpars.THEO_ENERGIES.from_value(energies)
    work = tmp_path / WORK_NAME
    work.mkdir()
    with execute_in_dir(work), private_random(psgen, PSGEN_SEED):
        firstline, phaseshifts = psgen.runPhaseshiftGen_old(slab, rpars)
    use_forward_slashes(work / INPUT_NAME)
    writePHASESHIFTS(firstline, phaseshifts, work / REFERENCE)
    return work, firstline, phaseshifts


def produce_stored_outputs(stored):
    """Return a replacement of subprocess.run 'producing' `stored`.

    Parameters
    ----------
    stored : Path
        The folder of the stored raw outputs of eeasisss: either a
        RAW_OUTPUT_ARCHIVE, or loose files matching RAW_OUTPUT_GLOB.

    Returns
    -------
    run : callable
        A function with the call signature of subprocess.run, that
        writes the raw outputs to the current directory.
    """
    archive_path = stored / RAW_OUTPUT_ARCHIVE
    def _run(args, **_):
        if archive_path.is_file():
            with zipfile.ZipFile(archive_path) as archive:
                archive.extractall()
        else:
            for raw_output in stored.glob(RAW_OUTPUT_GLOB):
                shutil.copy2(raw_output, raw_output.name)
        return subprocess.CompletedProcess(args, 0, 'stdout', '')
    return _run


def _check_stored(written, stored):
    """Fail if file `written` differs from file `stored`.

    The diff and failure are those of the check_reference fixture
    (tests/harness/reference_files.py). The fixture itself is not
    used, as run_mocked also runs in session-scoped fixtures, and
    the stored data is regenerated with the real eeasisss, not with
    --rewrite-references.
    """
    written_text = written.read_text(encoding=_ENCODING)
    stored_text = stored.read_text(encoding=_ENCODING)
    if written_text == stored_text:
        return
    fail_with_diff(f'{written.name} differs from the stored one. If the '
                   'change is intended, regenerate the stored eeasisss '
                   f'outputs (see {__name__}).',
                   stored_text,
                   written_text,
                   fromfile=f'stored: {stored}',
                   tofile=f'written: {written}')


def run_mocked(mocked, tmp_path):
    """Run psgen for a MockedSlab, with its stored eeasisss outputs.

    Parameters
    ----------
    mocked : MockedSlab
        The slab, its Rparams, and the folder of its stored data.
    tmp_path : Path
        An empty directory in which psgen runs.

    Returns
    -------
    work : Path
        The folder in which psgen ran.
    firstline : str
        The first line of the PHASESHIFTS.
    phaseshifts : list
        The phase shifts read from the stored eeasisss outputs.

    Raises
    ------
    pytest.fail.Exception
        If the eeasisss input written or the PHASESHIFTS read differ
        from the stored ones.
    """
    tree = make_stub_tree(tmp_path)
    with mock.patch.object(psgen.subprocess, 'run',
                           side_effect=produce_stored_outputs(mocked.stored)
                           ) as mock_run:
        work, *results = run_psgen(mocked.slab, mocked.rpars,
                                   tree, tmp_path)
    mock_run.assert_called_once()
    _check_stored(work / INPUT_NAME, mocked.stored / INPUT_NAME)
    _check_stored(work / REFERENCE, mocked.stored / REFERENCE)
    return (work, *results)


def _store_raw_outputs(store):
    """Return a replacement of subprocess.run that stores PS.r.* files.

    The real eeasisss is run, then its PS.r.* outputs are written,
    sorted by name and with a fixed date, to an archive in `store`.
    """
    def _run(*args, **kwargs):
        # Runs eeasisss with the arguments of psgen
        completed = subprocess.run(  # noqa: S603
            *args,
            check=False,
            **kwargs,
            )
        if completed.returncode:
            message = f'eeasisss failed:\n{completed.stderr}'
            raise RuntimeError(message)
        archive_path = store / RAW_OUTPUT_ARCHIVE
        with zipfile.ZipFile(archive_path, 'w') as archive:
            for raw_output in sorted(Path().glob(RAW_OUTPUT_GLOB)):
                info = zipfile.ZipInfo(raw_output.name, _ZIP_DATE)
                info.compress_type = zipfile.ZIP_DEFLATED
                archive.writestr(info, raw_output.read_bytes())
        return completed
    return _run


def regenerate(case_id, real_tree):
    """Rewrite the stored data of a case with the real eeasisss.

    Parameters
    ----------
    case_id : str
        One of MOCKED_CASES.
    real_tree : Path
        The TensErLEED tree with a compiled eeasisss.
    """
    mocked = _mocked(case_id)
    mocked.stored.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir)
        tree = copy_real_tree(tmp_path, real_tree)
        # Replace the subprocess module seen by psgen, not its run
        # function, which _store_raw_outputs calls
        run_real = _store_raw_outputs(mocked.stored)
        with mock.patch.object(psgen, 'subprocess',
                               SimpleNamespace(run=run_real)):
            work, *_ = run_psgen(mocked.slab, mocked.rpars,
                                 tree, tmp_path)
        for name in (INPUT_NAME, REFERENCE):
            shutil.copyfile(work / name, mocked.stored / name)


def main(case_ids):
    """Regenerate the stored data of `case_ids` (default: all)."""
    real_tree = get_tensorleed_path()
    for case_id in case_ids or MOCKED_CASES:
        print(f'Regenerating {case_id}...',  # noqa: T201  # Script
              flush=True)
        regenerate(case_id, real_tree)


if __name__ == '__main__':
    main(sys.argv[1:])
