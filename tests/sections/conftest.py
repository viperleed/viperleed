"""Test configuration for tests/sections.

Created on 2023-07-26

@author: Alexander M. Imre (@amimre)
@author: Michele Riva (@michele-riva)

Defines fixtures useful for testing the successful execution of
different segments of tleedm.

Fixtures
--------
delta_files_ag100
    INITIALIZATION + DELTAS run, for Ag(100) only
make_section_tempdir (factory)
    Temporary directory for a system name and a TLEEDM section.
init_files
    INITIALIZATION run
refcalc_files
    INITIALIZATION + REFCALC run
search_files_ag100
    INITIALIZATION + SEARCH run, for Ag(100) only
"""

import os
from pathlib import Path
import shutil
import sys
from zipfile import ZipFile

import pytest
import pytest_cases

VPR_PATH = str(Path(__file__).resolve().parents[3])
if VPR_PATH not in sys.path:
    sys.path.append(VPR_PATH)

# pylint: disable=wrong-import-position
# Will be fixed in installable version
from viperleed.tleedm import run_tleedm

from ..helpers import TEST_DATA
# pylint: enable=wrong-import-position


SOURCE_STR = str(Path(VPR_PATH) / 'viperleed')
ALWAYS_REQUIRED_FILES = ('PARAMETERS', 'EXPBEAMS.csv', 'POSCAR')

TENSERLEED_TEST_VERSIONS = (1.61, 0.0)  # (0 == Newest)

INIT_SURFACES = ('Ag(100)', 'Ag(100)_el_rename')
REFCALC_SURFACES = ('Ag(100)',)
AG_100_DISPLACEMENTS = {  # For DELTAS and SEARCH
    'DISPLACEMENTS_z': 'Deltas_z.zip',
    'DISPLACEMENTS_vib': 'Deltas_vib.zip',
    'DISPLACEMENTS_z+vib': 'Deltas_z+vib.zip'
    }


@pytest_cases.fixture(scope='session', name='make_section_tempdir')
def fixture_factory_make_section_tempdir(tmp_path_factory):
    """Return a temporary directory for a surface and a TLEEDM section."""
    def _make(surface, section, *other_specifiers):
        tmp_dir_name = f'{surface}_{section}' + '_'.join(other_specifiers)
        return tmp_path_factory.mktemp(basename=tmp_dir_name, numbered=True)
    return _make


class BaseTleedmFilesSetup:
    """Utility class for collecting files an running TLEEDM tests."""

    def __init__(self, surface_dir, tmp_test_path,
                 required_files=(), copy_dirs=()):
        """Initialize instance.

        Parameters
        ----------
        surface_dir : str
            Name of the directory in the tests/fixtures tree
            from which input files should be collected.
        tmp_test_path : Path
            Path to the (temporary) directory in which TLEEDM
            should be executed.
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
        self.test_path = tmp_test_path

        base_dir = self.inputs_path
        self.input_files_paths = [base_dir / pth for pth in copy_dirs]
        for input_dir in self.input_files_paths:
            shutil.copytree(input_dir, self.test_path, dirs_exist_ok=True)
            shutil.copytree(input_dir, self.work_path, dirs_exist_ok=True)

        self.home = Path().resolve()
        self.failed = -1
        self.work_files_after_run = []

    @property
    def inputs_path(self):
        """Return the path containing input files."""
        return TEST_DATA / self.surface_name

    @property
    def work_path(self):
        """Return the path to the work directory."""
        return self.test_path / 'work'

    def run_tleedm_from_setup(self, source, preset_params):
        """Move to work folder, execute, collect outcome, go back home."""
        os.chdir(self.work_path)
        try:
            self.failed = run_tleedm(source=source,
                                     preset_params=preset_params)
        finally:
            os.chdir(self.home)
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


@pytest_cases.fixture(scope='session')
@pytest.mark.parametrize('surface', INIT_SURFACES, ids=INIT_SURFACES)
@pytest.mark.parametrize('tl_version', TENSERLEED_TEST_VERSIONS,
                         ids=(str(v) for v in TENSERLEED_TEST_VERSIONS))
def init_files(surface, tl_version, make_section_tempdir):
    """Collect files and run an initialization."""
    files = BaseTleedmFilesSetup(
        surface_dir=surface,
        tmp_test_path=make_section_tempdir(surface, 'init'),
        required_files=['PHASESHIFTS',],
        copy_dirs=['initialization']
        )
    files.run_tleedm_from_setup(
        source=SOURCE_STR,
        preset_params={'RUN': [0,],  # only initialization
                       'TL_VERSION': tl_version,}
        )
    return files


_NON_INIT_TL_VERSION = 0.0  # i.e., most recent


@pytest_cases.fixture(scope='session')
@pytest.mark.parametrize('surface', REFCALC_SURFACES, ids=REFCALC_SURFACES)
def refcalc_files(surface, make_section_tempdir):
    """Collect files and execute a reference calculation."""
    files = BaseTleedmFilesSetup(
        surface_dir=surface,
        tmp_test_path=make_section_tempdir(surface, 'refcalc'),
        required_files=['PHASESHIFTS',],
        copy_dirs=['initialization']
        )
    files.run_tleedm_from_setup(
        source=SOURCE_STR,
        preset_params={'RUN': [0, 1],                                           # TODO: tleedm should probably automatically inject INIT!
                       'TL_VERSION': _NON_INIT_TL_VERSION,}
        )
    return files


@pytest_cases.fixture(scope='session')
@pytest.mark.parametrize('displacements', AG_100_DISPLACEMENTS,
                         ids=AG_100_DISPLACEMENTS)
def delta_files_ag100(displacements, make_section_tempdir):
    """Collect files, and run a delta-amplitude calculation for Ag(100)."""
    surface = 'Ag(100)'

    # correct DISPLACEMENTS
    files = BaseTleedmFilesSetup(
        surface_dir=surface,
        tmp_test_path=make_section_tempdir(surface, 'deltas', displacements),
        required_files=['PHASESHIFTS',],
        copy_dirs=['initialization', 'deltas']
        )
    disp_source = files.inputs_path / 'displacements' / displacements
    files.copy_displacements(displacements_path=disp_source)
    files.run_tleedm_from_setup(
        source=SOURCE_STR,
        preset_params={'RUN': [0, 2],  # init and deltas
                       'TL_VERSION': _NON_INIT_TL_VERSION,}
        )
    return files


@pytest_cases.fixture(scope='session')
@pytest.mark.parametrize('displacements, deltas', AG_100_DISPLACEMENTS.items(),
                         ids=AG_100_DISPLACEMENTS)
def search_files_ag100(displacements, deltas, make_section_tempdir):
    """Collect input files and run a structure optimization."""
    surface = 'Ag(100)'
    files = BaseTleedmFilesSetup(
        surface_dir=surface,
        tmp_test_path=make_section_tempdir(surface, 'search', displacements),
        required_files=[],
        copy_dirs=['initialization', 'deltas', 'search']
        )
    disp_source = files.inputs_path / 'displacements' / displacements
    deltas_source = files.inputs_path / 'search' / 'Deltas' / deltas
    files.copy_displacements(disp_source)
    files.copy_deltas(deltas_source)
    files.run_tleedm_from_setup(
        source=SOURCE_STR,
        preset_params={'RUN': [0, 3],  # init and search
                       'TL_VERSION': _NON_INIT_TL_VERSION,}
        )
    return files
