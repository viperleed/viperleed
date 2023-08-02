"""Module helpers of viperleed.tests.

Created on 2023-02-28

@author: Michele Riva
@author: Alexander M. Imre

Contains some useful general definitions that can be used when creating
or running tests.
"""


# Think about a decorator for injecting fixtures.
# Some ideas at
# https://github.com/pytest-dev/pytest/issues/2424
# https://github.com/pytest-dev/pytest/issues/6322
# https://github.com/nteract/testbook/issues/4

from pathlib import Path
import os
import sys
import shutil
from zipfile import ZipFile

import pytest
import numpy as np

vpr_path = str(Path(__file__).parent.parent.parent)
if os.path.abspath(vpr_path) not in sys.path:
    sys.path.append(os.path.abspath(vpr_path))

from viperleed.calc import run_tleedm
from viperleed.lib import symmetry
from viperleed.lib.classes.atom import Atom
from viperleed.lib.classes.rparams import Rparams
from viperleed.lib.classes.slab import Slab
from viperleed.lib.files import parameters, poscar
from viperleed.lib.files.vibrocc import readVIBROCC
from viperleed.lib.files.displacements import readDISPLACEMENTS, readDISPLACEMENTS_block


_FIXTURES_PATH = Path('tests/fixtures/')
_POSCARs_PATH = _FIXTURES_PATH / 'POSCARs'

_EXAMPLE_POSCAR_EXPECTATIONS = [("POSCAR_Ag(100)", 6, 'p4m', 0),
                                ("POSCAR_STO(100)-4x1", 136, 'pm', 0),
                                ("POSCAR_TiO2", 540, 'pmm', -1),
                                ("POSCAR_diamond", 96, 'pm', 89),
                                ("POSCAR_36C_p6m", 36, 'p6m', 0),
                                ("POSCAR_36C_cm", 36,'cm', 0),]
                               #("POSCAR_Fe3O4_SCV", 83, 'cmm', 50)]            #TODO: Phaseshift generation fails. Why? @Fkraushofer (worked in fkpCurie:Florian_OldLocalTests/Fe3O4-001-SCV/history/t000.r013_211220-133452)

_EXAMPLE_POSCARs = [file.name for file in _POSCARs_PATH.glob('POSCAR*')]

SOURCE_STR = str(Path(vpr_path) / "viperleed")
TENSORLEED_PATH = Path(vpr_path) / "viperleed" / "tensorleed"
ALWAYS_REQUIRED_FILES = ('PARAMETERS', 'EXPBEAMS.csv', 'POSCAR')
INPUTS_ORIGIN = Path(__file__).parent / "fixtures"
POSCAR_PATHS = INPUTS_ORIGIN / "POSCARs"

TENSERLEED_TEST_VERSIONS = ('1.72', '1.73', '1.74')

AG_100_DISPLACEMENTS_NAMES = ['DISPLACEMENTS_z', 'DISPLACEMENTS_vib', 'DISPLACEMENTS_z+vib']
AG_100_DELTAS_NAMES = ['Deltas_z.zip', 'Deltas_vib.zip', 'Deltas_z+vib.zip']


@pytest.fixture
def poscars_path():
    """Return the path to a POSCAR file."""
    return POSCAR_PATHS

@pytest.fixture
def inputs_path():
    """Return the path to an input file."""
    return INPUTS_ORIGIN

@pytest.fixture
def tensorleed_path():
    """Return the path to the tensorleed executable."""
    return TENSORLEED_PATH  


@pytest.fixture(params=[('Ag(100)')], ids=['Ag(100)',])
def init_files(request, tmp_path_factory, scope="function"):
    surface_name = request.param
    tmp_dir_name = f'{surface_name}_init'
    tmp_path = tmp_path_factory.mktemp(basename=surface_name, numbered=True)
    run = [0,] # only initialization
    files = BaseTleedmFilesSetup(surface_dir=surface_name,
                                tmp_test_path=tmp_path,
                                required_files=["PHASESHIFTS",],
                                copy_dirs=["initialization"])
    files.run_tleedm_from_setup(source=SOURCE_STR,
                                preset_params={
                                    "RUN":run,
                                    "TL_VERSION":1.73,
                                })
    return files


class BaseTleedmFilesSetup():
    def __init__(self, surface_dir, tmp_test_path, required_files=(), copy_dirs=()):
        self.surface_name = surface_dir
        self.required_files = set(ALWAYS_REQUIRED_FILES)
        self.required_files.update(required_files)
        self.copy_dirs = copy_dirs
        self.test_path = tmp_test_path

        self.work_path = self.test_path / "work"
        self.inputs_path = INPUTS_ORIGIN / self.surface_name
        self.input_files_paths = []

        for pth in self.copy_dirs:
            cur_input_dir = self.inputs_path / pth
            self.input_files_paths.append(cur_input_dir)
            shutil.copytree(cur_input_dir, self.test_path, dirs_exist_ok=True)
            shutil.copytree(cur_input_dir, self.work_path, dirs_exist_ok=True)

        self.home = os.getcwd()
        self.exit_code = None

    def run_tleedm_from_setup(self, source, preset_params):
        os.chdir(self.work_path)
        exit_code = run_tleedm(source=source,
                                preset_params=preset_params)
        self.exit_code = exit_code
        self.work_files_after_run = [file.name for file in Path(self.work_path).glob('*')]
        os.chdir(self.home)

    def expected_file_exists(self, expected_file):
        expected_path = Path(self.work_path) / expected_file
        return expected_path.exists()

    def copy_displacements(self, displacements_path):
        shutil.copy(displacements_path, self.work_path / 'DISPLACEMENTS')

    def copy_deltas(self, deltas_path):
        shutil.copy(deltas_path, self.test_path / "Deltas" / "Deltas_001.zip")
        shutil.copy(deltas_path, self.work_path / "Deltas" / "Deltas_001.zip")
        shutil.copy(deltas_path, self.work_path / "Deltas_001.zip")
        ZipFile(deltas_path, 'r').extractall(self.work_path)



class TestSetup:
    def test_work_path_exists(self, init_files):
        """Check that work_path was created properly."""
        assert init_files.work_path.is_dir()

    def test_files_copied_correctly(self, init_files):
        """Check if all files were copied correctly."""
        input_files = [file.name for file in init_files.test_path.glob('*')]
        source_copied = all(file in input_files for file in init_files.required_files)
        input_files = [file.name for file in Path(init_files.work_path).glob('*')]
        work_copied = all(file in input_files for file in init_files.required_files)
        assert source_copied and work_copied

@pytest.fixture()
def ag100_parameters_example():
    # read Ag(100) POSCAR and PARAMETERS files
    slab = poscar.readPOSCAR(_FIXTURES_PATH / 'Ag(100)' / 'initialization' / 'POSCAR')
    rpars = parameters.readPARAMETERS(_FIXTURES_PATH / 'Ag(100)' / 'initialization' / 'PARAMETERS')
    # interpret PARAMETERS file
    interpreter = parameters.ParameterInterpreter(rpars)
    interpreter.interpret(slab)
    symmetry.findSymmetry(slab, rpars)
    symmetry.findBulkSymmetry(slab, rpars)
    return (rpars, slab)


@pytest.fixture()
def ag100_rename_ax(request, tmp_path):
    dir = INPUTS_ORIGIN / "Ag(100)_el_rename"
    run = [0,] # only initialization
    files = BaseTleedmFilesSetup(surface_dir=dir,
                                tmp_test_path=tmp_path,
                                required_files=["PHASESHIFTS",],
                                copy_dirs=["initialization"])
    files.run_tleedm_from_setup(source=SOURCE_STR,
                                preset_params={
                                    "RUN":run,
                                    "TL_VERSION":1.73,
                                })
    return files


@pytest.fixture(scope="function", params=_EXAMPLE_POSCARs)
def example_poscars(request):
    file_path = _POSCARs_PATH / request.param
    slab = poscar.readPOSCAR(str(file_path))
    return slab

@pytest.fixture(scope="function", params=_EXAMPLE_POSCAR_EXPECTATIONS)
def slab_and_expectations(request):
    filename, expected_n_atoms, expected_pg, offset_at = request.param
    file_path = _POSCARs_PATH / filename
    pos_slab = poscar.readPOSCAR(str(file_path))
    return (pos_slab, expected_n_atoms, expected_pg, offset_at)

@pytest.fixture(scope="function")
def slab_pg_rp(slab_and_expectations):
    slab, *_ = slab_and_expectations
    rp = Rparams()
    slab.fullUpdate(rp)
    pg = symmetry.findSymmetry(slab, rp, output=False)
    symmetry.enforceSymmetry(slab, rp)
    return slab, pg, rp


@pytest.fixture(scope='function')
def manual_slab_3_atoms():
    slab = Slab()
    slab.ucell = np.diag([3., 4., 5.])
    positions = (np.array([-0.25, 0, 0]),
                 np.array([0.00, 0, 0]),
                 np.array([0.25, 0, 0]))
    slab.atlist = [Atom('C', pos, i+1, slab)
                   for i, pos in enumerate(positions)]
    param = Rparams()
    slab.fullUpdate(param)
    return slab


@pytest.fixture()
def manual_slab_1_atom_trigonal():
    slab = Slab()
    slab.ucell = np.array([[ 1, 0, 0],
                           [-2, 3, 0],
                           [ 1, 2, 3]],dtype=float)
    slab.atlist = [Atom('C', np.array([0.2, 0.7, 0.1]), 1, slab),]  # "random" position
    param = Rparams()
    slab.fullUpdate(param)
    return slab


@pytest.fixture()
def ag100_slab_param(poscars_path):
    slab = poscar.readPOSCAR(poscars_path /"POSCAR_Ag(100)")
    param = Rparams()
    param.N_BULK_LAYERS = 1
    slab.fullUpdate(param)
    return slab, param


@pytest.fixture()
def ag100_slab_with_displacements_and_offsets(ag100_slab_param, inputs_path):
    slab, param = ag100_slab_param
    vibrocc_path = inputs_path / "Ag(100)" / "mergeDisp" / "VIBROCC"
    displacements_path = inputs_path / "Ag(100)" / "mergeDisp" / "DISPLACEMENTS_mixed"
    readVIBROCC(param, slab, str(vibrocc_path))
    readDISPLACEMENTS(param, str(displacements_path))
    readDISPLACEMENTS_block(param, slab, param.disp_blocks[param.search_index])
    return slab, param


@pytest.fixture(params=[('Ag(100)')], ids=['Ag(100)',])
def refcalc_files(request, tmp_path_factory, scope="session"):
    surface_name = request.param
    tmp_dir_name = f'{surface_name}_refcalc'
    tmp_path = tmp_path_factory.mktemp(basename=tmp_dir_name, numbered=True)
    run = [0, 1] # initialization and refcalc
    files = BaseTleedmFilesSetup(surface_dir=surface_name,
                                tmp_test_path=tmp_path,
                                required_files=["PHASESHIFTS",],
                                copy_dirs=["initialization"])
    files.run_tleedm_from_setup(source=SOURCE_STR,
                                preset_params={
                                    "RUN":run,
                                    "TL_VERSION":1.74,
                                })
    return files

@pytest.fixture(params=AG_100_DISPLACEMENTS_NAMES, ids=AG_100_DISPLACEMENTS_NAMES)
def delta_files_ag100(request, tmp_path_factory, scope="session"):
    displacements_name = request.param
    surface_name = 'Ag(100)'
    tmp_dir_name = tmp_dir_name = f'{surface_name}_deltas_{displacements_name}'
    tmp_path = tmp_path_factory.mktemp(basename=tmp_dir_name, numbered=True)
    run = [0, 2] # init and deltas
    required_files = ["PHASESHIFTS",]
    copy_dirs=["initialization", "deltas"]
    # correct DISPLACEMENTS
    files = BaseTleedmFilesSetup(surface_dir=surface_name,
                                tmp_test_path=tmp_path,
                                required_files=required_files,
                                copy_dirs=copy_dirs)
    disp_source = files.inputs_path / "displacements" / displacements_name
    files.copy_displacements(displacements_path=disp_source)
    files.run_tleedm_from_setup(source=SOURCE_STR,
                                preset_params={
                                    "RUN":run,
                                    "TL_VERSION":1.74,
                                })
    return files


@pytest.fixture(params=list(zip(AG_100_DISPLACEMENTS_NAMES, AG_100_DELTAS_NAMES)),
                ids=AG_100_DISPLACEMENTS_NAMES)
def search_files_ag100(request, tmp_path_factory, scope="session"):
    surface_name = 'Ag(100)'
    displacements_name, deltas_name = request.param
    tmp_dir_name = tmp_dir_name = f'{surface_name}_search_{displacements_name}'
    tmp_path = tmp_path_factory.mktemp(basename=tmp_dir_name, numbered=True)
    run = [0, 3] # init and search
    required_files = []
    copy_dirs=["initialization", "deltas", "search"]
    files = BaseTleedmFilesSetup(surface_dir=surface_name,
                                tmp_test_path=tmp_path,
                                required_files=required_files,
                                copy_dirs=copy_dirs)
    disp_source = files.inputs_path / "displacements" / displacements_name
    deltas_source = files.inputs_path / "search" / "Deltas" / deltas_name
    files.copy_displacements(disp_source)
    files.copy_deltas(deltas_source)
    files.run_tleedm_from_setup(source=SOURCE_STR,
                                preset_params={
                                    "RUN":run,
                                    "TL_VERSION":1.74,
                                })
    return files