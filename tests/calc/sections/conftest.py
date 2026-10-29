"""Test configuration for tests.calc.sections.

Defines fixtures useful for testing the successful execution of
different segments of viperleed.calc.

Fixtures
--------
delta_files_ag100
    INITIALIZATION + DELTAS run, for Ag(100) only
make_section_tempdir (factory)
    Temporary directory for a system name and a viperleed.calc section.
init_files
    INITIALIZATION run
refcalc_files
    INITIALIZATION + REFCALC run
search_files_ag100
    INITIALIZATION + SEARCH run, for Ag(100) only
"""

__authors__ = (
    'Alexandra Mia Imre (@alexmiame)',
    'Michele Riva (@michele-riva)',
)
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2023-07-26'
__license__ = 'GPLv3+'

from contextlib import contextmanager
from pathlib import Path
import re
import shutil
import subprocess
from zipfile import ZipFile

import pytest
from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.calc.classes.rparams.defaults import DEFAULTS
from viperleed.calc.constants import DEFAULT_WORK
from viperleed.calc.constants import LOG_PREFIX
from viperleed.calc.files import parameters
from viperleed.calc.files import tenserleed
from viperleed.calc.lib import leedbase
from viperleed.calc.lib.context import execute_in_dir
from viperleed.calc.lib.fs_utils import copytree_exists_ok
from viperleed.calc.lib.version import Version
from viperleed.calc.run import run_calc

from ...helpers import TEST_DATA


ALWAYS_REQUIRED_FILES = ('PARAMETERS', 'EXPBEAMS.csv', 'POSCAR')

TENSERLEED_TEST_VERSIONS = (Version('2.0.0'), Version('1.6.1'))

INIT_SURFACES = (
    'Ag(100)',
    'Ag(100)_el_rename',
    )
INIT_DOMAINS = (   # In the domains subfolder of TEST_DATA
    'silver_and_bismuth',
    )
DELTA_DOMAINS = (  # Also in the domain subfolder of TEST_DATA
    'telluride_stacking',
    )
REFCALC_SURFACES = ('Ag(100)',)
# Seconds between reads of SD.TL during a SEARCH. The default (60 s) is
# longer than what the search of Ag(100) needs to reach SEARCH_MAX_GEN
# on a fast, idle machine. Then, no evaluation happens, and the
# Search-progress/Search-report PDFs are not written. See also #582.
SEARCH_EVAL_TIME = 1
CHAIN_TL_VERSION = Version('1.6.1')  # Older than 1.7.0
CHAIN_DISPLACEMENTS = 'DISPLACEMENTS_z'
AG_100_DISPLACEMENTS = {  # For DELTAS and SEARCH
    'DISPLACEMENTS_z': 'Deltas_z.zip',
    'DISPLACEMENTS_vib': 'Deltas_vib.zip',
    'DISPLACEMENTS_z+vib': 'Deltas_z+vib.zip'
    }

with_tl_versions = parametrize(
    tl_version=TENSERLEED_TEST_VERSIONS,
    ids=[str(v) for v in TENSERLEED_TEST_VERSIONS],
    )


@fixture(scope='session', name='make_section_tempdir')
def fixture_factory_make_section_tempdir(tmp_path_factory):
    """Return a temporary directory for a surface and a calc section."""
    def _make(surface, section, *other_specifiers):
        surface = surface.replace('/', '_')
        tmp_dir_name = f'{surface}_{section}' + '_'.join(other_specifiers)
        return tmp_path_factory.mktemp(basename=tmp_dir_name, numbered=True)
    return _make


def fail_if_tl_version_missing(source, version):
    """Fail the current test/fixture if `source` has no TensErLEED `version`.

    The tests that call this are in the 'tensorleed' tier, where a
    missing TensErLEED version is an error rather than a reason to skip.

    Parameters
    ----------
    source : path-like
        Folder with TensErLEED sources (e.g., $VIPERLEED_TENSORLEED).
    version : Version or str
        The TensErLEED version that is needed.

    Raises
    ------
    Failed
        If `version` is not among the TensErLEED sources in `source`.
    Exception
        Any error from looking up the TensErLEED sources in `source`
        (e.g., an entry without a version).
    """
    version = Version(str(version))
    available = [s.version
                 for s in tenserleed.get_tenserleed_sources(source)]
    if version not in available:
        pytest.fail(f'TensErLEED v{version} is not available in {source}. '
                    'Set VIPERLEED_TENSORLEED to a folder that contains it.')


class BaseCalcFilesSetup:
    """Utility class for collecting files and running viperleed.calc tests."""

    def __init__(self, surface_dir, tmp_test_path,
                 required_files=(), copy_dirs=()):
        """Initialize instance.

        Parameters
        ----------
        surface_dir : str
            Name of the directory in the tests/_test_data tree
            from which input files should be collected.
        tmp_test_path : Path
            Path to the (temporary) directory in which
            viperleed.calc should be executed.
        required_files : Iterable of str, optional
            Name of files that should be present. Files that
            are common to all sections do not need to be listed.
            Default is an empty tuple.
        copy_dirs : Iterable of str, optional
            Subfolders of the input path that should be copied
            before execution. Default is an empty tuple.

        Returns
        -------
        None.
        """
        self.surface_name = surface_dir
        self.required_files = set(ALWAYS_REQUIRED_FILES)
        self.required_files.update(required_files)
        # pylint: disable-next=magic-value-comparison
        if 'domains' in surface_dir:
            self.required_files.remove('POSCAR')
        self.test_path = tmp_test_path

        base_dir = self.inputs_path
        self.input_files_paths = [base_dir / pth for pth in copy_dirs]
        for input_dir in self.input_files_paths:
            copytree_exists_ok(input_dir, self.test_path)
            copytree_exists_ok(input_dir, self.work_path)

        self.failed = -1
        self.records = None
        self.work_files_after_run = []

    @property
    def inputs_path(self):
        """Return the path containing input files."""
        return TEST_DATA / self.surface_name

    @property
    def work_path(self):
        """Return the path to the work directory."""
        return self.test_path / DEFAULT_WORK

    def run_calc_from_setup(self, source, preset_params):
        """Move to work folder, execute, collect outcome, go back home."""
        fail_if_tl_version_missing(source, preset_params['TL_VERSION'])
        with execute_in_dir(self.work_path):
            self.failed, self.records = run_calc(
                source=source,
                preset_params=preset_params,
                home=self.test_path,
                )
        self.work_files_after_run = [f.name for f in self.work_path.glob('*')]

    def expected_file_exists(self, expected_file):
        """Return whether the work directory contains `expected_file`."""
        return (self.work_path / expected_file).exists()

    def section_ran(self, section):
        """Return whether `section` ran, and nothing failed."""
        return (
            not self.failed
            and self.records is not None
            and bool(self.records.get_last_state_for_section(section))
            )

    def out_suffixed_files(self):
        """Return the files in work with an old-style _OUT suffix."""
        # Skip R_OUT files
        return [f for f in self.work_path.rglob('*_OUT*')
                # pylint: disable-next=magic-value-comparison
                if 'R=' not in f.name]

    def final_r_factor(self, section):
        """Return the overall R factor of `section` in the calc log.

        Parameters
        ----------
        section : str
            Name of the section, as in 'Final R (<section>)'.

        Returns
        -------
        r_factor : float or None
            The R factor that was logged at the end of the run, or
            None if there is none (e.g., no log file, section did not
            run).
        """
        pattern = re.compile(rf'Final R \({section}\): (\d+\.\d+)')
        for log in sorted(self.work_path.glob(f'{LOG_PREFIX}-*.log')):
            match = pattern.search(log.read_text())
            if match:
                return float(match.group(1))
        return None

    def copy_displacements(self, displacements_path):
        """Copy displacements_path to work as a DISPLACEMENTS file."""
        shutil.copy(displacements_path, self.work_path / 'DISPLACEMENTS')

    def copy_deltas(self, deltas_path):
        """Copy delta-amplitude files to test and work directories."""
        shutil.copy(deltas_path, self.test_path / 'Deltas' / 'Deltas_001.zip')
        shutil.copy(deltas_path, self.work_path / 'Deltas' / 'Deltas_001.zip')
        shutil.copy(deltas_path, self.work_path / 'Deltas_001.zip')
        with ZipFile(deltas_path, 'r') as archive:
            archive.extractall(self.work_path)

    def read_manifest(self):
        """Collect the contents of the main manifest file."""
        return (self.work_path/'manifest').read_text().splitlines()


class DomainsCalcFilesSetup(BaseCalcFilesSetup):
    """Utility for executing multi-domain calc runs."""

    def __init__(self, surface_dir, *args, **kwargs):
        """Initialize instance."""
        super().__init__(f'domains/{surface_dir}', *args, **kwargs)
        self.src_folders = ()    # Paths (as str) given in PARAMETERS
        self.work_domains = {}   # {src_path_as_str: work_path}

    def find_most_recent_tensors(self):
        """Collect the most recent tensor number for each domain."""
        self._collect_domain_info()
        src_tensors = {src: leedbase.getMaxTensorIndex(src)
                       for src in self.src_folders}
        work_tensors = {src: leedbase.getMaxTensorIndex(work)
                        for src, work in self.work_domains.items()}
        return {src: max(src_tensors[src], work_tensors[src])
                for src in self.src_folders}

    def run_calc_from_setup(self, *args, **kwargs):
        """Move to work folder, execute, collect outcome, go back home."""
        super().run_calc_from_setup(*args, **kwargs)
        self._collect_domain_info()

    def _collect_domain_info(self):
        """Collect information about the domains in the calculation."""
        src_parameters = [file
                          for src_path in self.input_files_paths
                          for file in src_path.glob('PARAMETERS')]
        assert len(src_parameters) == 1
        rpars = parameters.read(src_parameters[0])
        domains = rpars.readParams['DOMAIN']
        src_to_folders_created = {
            domain.values_str: {f'Domain_{domain.flags_str}',
                                Path(domain.values_str).name}
            for domain in domains
            }
        self.src_folders = tuple(src_to_folders_created)
        self.work_domains = {
            src: self.work_path/n
            for src, paths in src_to_folders_created.items()
            for n in paths
            if (self.work_path/n).is_dir()
            }


@fixture(scope='session')
@parametrize(surface=INIT_SURFACES, ids=INIT_SURFACES)
@with_tl_versions
def init_files(surface, tl_version, make_section_tempdir, tensorleed_path):
    """Collect files and run an initialization."""
    files = BaseCalcFilesSetup(
        surface_dir=surface,
        tmp_test_path=make_section_tempdir(surface, 'init'),
        required_files=['PHASESHIFTS',],
        copy_dirs=['initialization'],
        )
    files.run_calc_from_setup(
        source=tensorleed_path,
        preset_params={'RUN': [0,],  # only initialization
                       'TL_VERSION': tl_version,}
        )
    return files


@fixture(scope='session')
@parametrize(domains=INIT_DOMAINS)
@with_tl_versions
def init_domains(domains, tl_version, make_section_tempdir, tensorleed_path):
    """Collect input files for a DOMAINS calculation and run initialization."""
    setup = DomainsCalcFilesSetup(
        surface_dir=domains,
        tmp_test_path=make_section_tempdir(domains, 'init'),
        required_files=['PHASESHIFTS',],
        copy_dirs=['initialization'],
        )
    setup.run_calc_from_setup(
        source=tensorleed_path,
        preset_params={'RUN': [0,],  # only initialization
                       'TL_VERSION': tl_version,}
        )
    return setup



_NON_INIT_TL_VERSION = tenserleed.CURRENT_TL_VERSION  # i.e., most recent       # TODO: to prevent regressions like #101, it's probably better to run this stuff also for other versions!


@fixture(scope='session')
@parametrize(surface=REFCALC_SURFACES, ids=REFCALC_SURFACES)
def refcalc_files(surface, make_section_tempdir, tensorleed_path):
    """Collect files and execute a reference calculation."""
    files = BaseCalcFilesSetup(
        surface_dir=surface,
        tmp_test_path=make_section_tempdir(surface, 'refcalc'),
        required_files=['PHASESHIFTS',],
        copy_dirs=['initialization']
        )
    files.run_calc_from_setup(
        source=tensorleed_path,
        preset_params={'RUN': [0, 1],                                           # TODO: calc should probably automatically inject INIT!
                       'TL_VERSION': _NON_INIT_TL_VERSION,}
        )
    return files


@contextmanager
def _short_search_eval_time():
    """Temporarily make SEARCH read SD.TL every SEARCH_EVAL_TIME s."""
    with pytest.MonkeyPatch.context() as patcher:
        # SEARCH_EVAL_TIME is not a PARAMETERS keyword
        patcher.setitem(DEFAULTS, 'SEARCH_EVAL_TIME', SEARCH_EVAL_TIME)
        yield


@fixture(scope='session')
@parametrize(displacements=AG_100_DISPLACEMENTS, ids=AG_100_DISPLACEMENTS)
def delta_files_ag100(displacements, make_section_tempdir, tensorleed_path):
    """Collect files, and run a delta-amplitude calculation for Ag(100)."""
    surface = 'Ag(100)'

    # correct DISPLACEMENTS
    files = BaseCalcFilesSetup(
        surface_dir=surface,
        tmp_test_path=make_section_tempdir(surface, 'deltas', displacements),
        required_files=['PHASESHIFTS',],
        copy_dirs=['initialization', 'deltas']
        )
    disp_source = files.inputs_path / 'displacements' / displacements
    files.copy_displacements(displacements_path=disp_source)
    files.run_calc_from_setup(
        source=tensorleed_path,
        preset_params={'RUN': [0, 2],  # init and deltas
                       'TL_VERSION': _NON_INIT_TL_VERSION,}
        )
    return files


@fixture(scope='session')
@parametrize(domains=DELTA_DOMAINS)
@with_tl_versions
def delta_domains(domains, tl_version, make_section_tempdir, tensorleed_path):
    """Collect files, and run a multi-domain delta-amplitude calculation."""
    setup = DomainsCalcFilesSetup(
        surface_dir=domains,
        tmp_test_path=make_section_tempdir(domains, 'deltas'),
        copy_dirs=['deltas'],
        )
    setup.run_calc_from_setup(
        source=tensorleed_path,
        preset_params={'RUN': [0, 2],  # init and deltas
                       'TL_VERSION': tl_version,}
        )
    return setup


@fixture(scope='session')
@parametrize('displacements,deltas',
             AG_100_DISPLACEMENTS.items(),
             ids=AG_100_DISPLACEMENTS)
def search_files_ag100(displacements, deltas,
                       make_section_tempdir,
                       tensorleed_path):
    """Collect input files and run a structure optimization."""
    surface = 'Ag(100)'
    files = BaseCalcFilesSetup(
        surface_dir=surface,
        tmp_test_path=make_section_tempdir(surface, 'search', displacements),
        required_files=[],
        copy_dirs=['initialization', 'deltas', 'search']
        )
    disp_source = files.inputs_path / 'displacements' / displacements
    deltas_source = files.inputs_path / 'search' / 'Deltas' / deltas
    files.copy_displacements(disp_source)
    files.copy_deltas(deltas_source)
    with _short_search_eval_time():
        files.run_calc_from_setup(
            source=tensorleed_path,
            preset_params={'RUN': [0, 3],  # init and search
                           'TL_VERSION': _NON_INIT_TL_VERSION}
            )
    return files


@fixture(scope='session', name='tensorleed_chain_path')
def fixture_tensorleed_chain_path(tensorleed_path, tmp_path_factory):
    """Return a TensErLEED tree whose CHAIN_TL_VERSION can be linked.

    The precompiled randomizer objects (MPIrandom_.o, random_.o) in
    the lib folder of TensErLEED versions up to 1.7.3 are not
    position-independent. Linking them fails with the default (PIE)
    options of recent toolchains. The returned tree is a copy of
    `tensorleed_path` in which they are recompiled for CHAIN_TL_VERSION.
    Everything else is a link to the original tree.

    Returns
    -------
    tree : Path
        The path to the new TensErLEED tree.
    """
    fail_if_tl_version_missing(tensorleed_path, CHAIN_TL_VERSION)
    chain_source, = (
        s.path for s in tenserleed.get_tenserleed_sources(tensorleed_path)
        if s.version == CHAIN_TL_VERSION
        )
    tree = tmp_path_factory.mktemp('tensorleed_chain')
    chain_source_name = None
    for entry in Path(tensorleed_path).iterdir():
        if entry.resolve() == chain_source:
            chain_source_name = entry.name
            shutil.copytree(entry, tree / entry.name)
        else:
            (tree / entry.name).symlink_to(entry)
    lib = tree / chain_source_name / 'lib'
    source_c = lib / 'random_.c'
    compilers = [('gcc', 'random_.o')]
    if shutil.which('mpicc'):  # Otherwise, calc uses a serial search
        compilers.append(('mpicc', 'MPIrandom_.o'))
    for compiler, obj in compilers:
        command = [compiler, '-fPIC', '-c', str(source_c),
                   '-o', str(lib / obj)]
        try:
            # S603: the command is built from fixed names and paths
            subprocess.run(  # noqa: S603
                command,
                check=True,
                capture_output=True,
                text=True,
                )
        except FileNotFoundError:
            pytest.fail(f'Cannot rebuild {obj} for TensErLEED '
                        f'v{CHAIN_TL_VERSION}: {compiler} not found.')
        except subprocess.CalledProcessError as exc:
            command_line = ' '.join(command)
            pytest.fail(f'Failed to rebuild {obj} for TensErLEED '
                        f'v{CHAIN_TL_VERSION}: {command_line}\n'
                        f'{exc.stderr}')
    return tree


@fixture(scope='session', name='chain_files_ag100_v161')
def fixture_chain_files_ag100_v161(make_section_tempdir,
                                   tensorleed_chain_path):
    """Run INITIALIZATION to SEARCH for Ag(100) with TL 1.6.1.

    Nothing stored is used apart from the inputs of the INITIALIZATION
    and the displacements: the tensors for the delta-amplitude
    calculation come from the reference calculation, and the delta
    amplitudes for the search from the DELTAS section.
    """
    surface = 'Ag(100)'
    files = BaseCalcFilesSetup(
        surface_dir=surface,
        tmp_test_path=make_section_tempdir(surface, 'chain', 'v161'),
        required_files=['PHASESHIFTS'],
        copy_dirs=['initialization']
        )
    disp_source = (files.inputs_path / 'displacements'
                   / CHAIN_DISPLACEMENTS)
    files.copy_displacements(disp_source)
    with _short_search_eval_time():
        files.run_calc_from_setup(
            source=tensorleed_chain_path,
            preset_params={'RUN': [0, 1, 2, 3],
                           'TL_VERSION': CHAIN_TL_VERSION}
            )
    return files
