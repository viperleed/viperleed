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
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-07-26'
__license__ = 'GPLv3+'

import shutil
from zipfile import ZipFile

from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.calc.constants import DEFAULT_WORK
from viperleed.calc.files import tenserleed
from viperleed.calc.lib.base import copytree_exists_ok
from viperleed.calc.lib.context import execute_in_dir
from viperleed.calc.lib.version import Version
from viperleed.calc.run import run_calc

from ...helpers import TEST_DATA


ALWAYS_REQUIRED_FILES = ('PARAMETERS', 'EXPBEAMS.csv', 'POSCAR')

TENSERLEED_TEST_VERSIONS = (Version('2.0.0'), Version('1.6.1'))

INIT_SURFACES = (
    'Ag(100)',
    'Ag(100)_el_rename',
    )
INIT_DOMAINS = (  # In the domains subfolder of TEST_DATA
    'silver_and_bismuth',
    )
REFCALC_SURFACES = ('Ag(100)',)
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
    setup = BaseCalcFilesSetup(
        surface_dir=f'domains/{domains}',
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
    files.run_calc_from_setup(
        source=tensorleed_path,
        preset_params={'RUN': [0, 3],  # init and search
                       'TL_VERSION': _NON_INIT_TL_VERSION,}
        )
    return files
