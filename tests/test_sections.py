import pytest
import shutil, tempfile
import sys
import os
from pathlib import Path
from zipfile import ZipFile
from copy import deepcopy
import numpy as np

vpr_path = str(Path(__file__).parent.parent.parent)
if os.path.abspath(vpr_path) not in sys.path:
    sys.path.append(os.path.abspath(vpr_path))


import viperleed
import viperleed.tleedmlib
from viperleed.tleedm import run_tleedm

from tleedmlib.files.poscar import readPOSCAR
from viperleed import tleedmlib as tl
from tleedmlib.symmetry import findSymmetry

SOURCE_STR = str(Path(vpr_path) / "viperleed")
ALWAYS_REQUIRED_FILES = ('PARAMETERS', 'EXPBEAMS.csv', 'POSCAR')
INPUTS_ORIGIN = Path(__file__).parent / "fixtures"

TENSERLEED_TEST_VERSIONS = ('1.71', '1.72', '1.73', '1.74')

AG_100_DISPLACEMENTS_NAMES = ['DISPLACEMENTS_z', 'DISPLACEMENTS_vib', 'DISPLACEMENTS_z+vib']
AG_100_DELTAS_NAMES = ['Deltas_z.zip', 'Deltas_vib.zip', 'Deltas_z+vib.zip']


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


@pytest.fixture(params=[('Ag(100)')], ids=['Ag(100)',])
def init_files(request, tmp_path_factory, scope="session"):
    surface_name = request.param
    tmp_dir_name = f'{surface_name}_inti'
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
                                    "TL_VERSION":1.73,
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
                                    "TL_VERSION":1.73,
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
                                    "TL_VERSION":1.73,
                                })
    return files


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


class TestInitialization(TestSetup):
    def test_exit_code_0(self, init_files):
        """Test if initialization gives exit code 0."""
        assert init_files.exit_code == 0

    @pytest.mark.parametrize('expected_file', (('IVBEAMS'), ('BEAMLIST'), ('VIBROCC')))
    def test_init_files_present(self, init_files, expected_file):
        """Checks if files are present after initialization"""
        assert init_files.expected_file_exists(expected_file)


class TestRefCalc(TestInitialization):
    @pytest.mark.parametrize('expected_file', (('THEOBEAMS.csv',)))
    def test_refcalc_files_present(self, refcalc_files, expected_file):
        assert refcalc_files.expected_file_exists(expected_file)


class TestDeltasAg100(TestSetup):
    def test_delta_input_written(self, delta_files_ag100):
        assert delta_files_ag100.expected_file_exists("delta-input")


    def test_exit_code_0(self, delta_files_ag100):
        assert delta_files_ag100.exit_code == 0


    def test_deltas_zip_created(self, delta_files_ag100):
        assert delta_files_ag100.expected_file_exists(Path("Deltas") / "Deltas_001.zip")


class TestSearchAg100(TestSetup):
    def test_exit_code_0(self, search_files_ag100):
        assert search_files_ag100.exit_code == 0
        
    @pytest.mark.parametrize('expected_file', ('search.steu',))
    def test_search_input_exist(self, search_files_ag100, expected_file):
        assert search_files_ag100.expected_file_exists(expected_file)

    @pytest.mark.parametrize('expected_file', ('SD.TL', 'control.chem'))
    def test_search_raw_files_exist(self, search_files_ag100, expected_file):
        assert search_files_ag100.expected_file_exists(expected_file)

    @pytest.mark.parametrize('expected_file', ('Search-report.pdf', 'Search-progress.pdf'))
    def test_search_pdf_files_exist(self, search_files_ag100, expected_file):
        assert search_files_ag100.expected_file_exists(expected_file)



@pytest.fixture(scope="class", params=[('POSCAR_STO(100)-4x1', 136, 'pm'),
                                       ("POSCAR_TiO2", 540, 'pmm'),
                                       ("POSCAR_diamond", 96, 'pm'),
                                       ("POSCAR_p6m_36C", 36, 'p6m')])
def slab_and_expectations(request, scope="session"):
    filename, expected_n_atoms, expected_pg = request.param
    file_path = Path(__file__).parent / "fixtures" / "POSCARs" / filename
    pos_slab = readPOSCAR(str(file_path))
    return (pos_slab, expected_n_atoms, expected_pg)

@pytest.fixture()
def slab_pg(slab_and_expectations, scope="session"):
    slab, *_ = slab_and_expectations
    rp = tl.Rparams()
    slab.fullUpdate(rp)
    pg = findSymmetry(slab, rp, output=False)
    return pg


class TestPOSCARRead:
    def test_read_in_atoms(self, slab_and_expectations):
        slab, *_ = slab_and_expectations
        assert len(slab.atlist) > 0

    def test_n_atom_correct(self, slab_and_expectations):
        slab, expected_n_atoms, _ = slab_and_expectations
        assert len(slab.atlist) == expected_n_atoms


class TestPOSCARSymmetry(TestPOSCARRead):
    def test_pg_found(self, slab_pg):
        assert slab_pg != 'unknown'

    def test_pg_correct(self, slab_and_expectations, slab_pg):
        _, _, expected_pg = slab_and_expectations
        assert slab_pg == expected_pg


# Test Slab with p6m symmetry by Michele and Alex

_CARBON_SLABS = (('POSCAR_p6m_36C','p6m'),
                 ('POSCAR_36C_slanted_cm','cm'))
_CARBON_SLAB_NAMES = [name for name, _ in _CARBON_SLABS]

@pytest.fixture(params=_CARBON_SLABS, ids=_CARBON_SLAB_NAMES, scope="function")
def carbon_slab(request):
    poscar_name, expected_group = request.param
    poscar_path = Path(__file__).parent / "fixtures" / "POSCARs" / poscar_name
    slab = readPOSCAR(str(poscar_path))
    slab.expected_group = expected_group
    return slab


@pytest.fixture(scope="function")
def carbon_setup(carbon_slab):
    param = tl.classes.rparams.Rparams()
    param.BULK_REPEAT = 3.0
    tl.files.parameters.interpretPARAMETERS(param, carbon_slab)
    param.updateDerivedParams()
    carbon_slab.fullUpdate(param)
    return carbon_slab, param

def test_recognize_carbon_symmetry(carbon_setup):
    slab, param = carbon_setup
    assert tl.symmetry.findSymmetry(slab, param) == slab.expected_group


@pytest.mark.parametrize("displacement", [(4, (np.array([0.2, 0, 0]),)),
                                          (4, (np.array([0, 0.2, 0]),)),
                                          (4, (np.array([0, 0, 0.2]),)),
                                          ])
def test_preserve_carbon_symmetry_with_displacement(displacement, carbon_setup):
    carbon_slab, param = carbon_setup
    tl.symmetry.findSymmetry(carbon_slab, param)
    tl.symmetry.enforceSymmetry(carbon_slab, param)
    sl_copy = deepcopy(carbon_slab)

    # manually assign displacements
    sl_copy.atlist[0].assignDisp(*displacement)

    for at in sl_copy.atlist:
        disp = at.disp_geo_offset['all'][0]
        at.cartpos += disp
    sl_copy.getFractionalCoordinates()

    assert tl.symmetry.findSymmetry(sl_copy, deepcopy(param)) == sl_copy.expected_group
