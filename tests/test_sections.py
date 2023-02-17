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
from tleedmlib.symmetry import findSymmetry, enforceSymmetry
from tleedmlib.psgen import runPhaseshiftGen_old

SOURCE_STR = str(Path(vpr_path) / "viperleed")
TENSORLEED_PATH = Path(vpr_path) / "viperleed" / "tensorleed"
ALWAYS_REQUIRED_FILES = ('PARAMETERS', 'EXPBEAMS.csv', 'POSCAR')
INPUTS_ORIGIN = Path(__file__).parent / "fixtures"

TENSERLEED_TEST_VERSIONS = ('1.71', '1.72', '1.73', '1.74')

AG_100_DISPLACEMENTS_NAMES = ['DISPLACEMENTS_z', 'DISPLACEMENTS_vib', 'DISPLACEMENTS_z+vib']
AG_100_DELTAS_NAMES = ['Deltas_z.zip', 'Deltas_vib.zip', 'Deltas_z+vib.zip']

_EXAMPLE_POSCARs = [("POSCAR_Ag(100)", 6, 'p4m', 0),
                    ("POSCAR_STO(100)-4x1", 136, 'pm', 0),
                    ("POSCAR_TiO2", 540, 'pmm', -1),
                    ("POSCAR_diamond", 96, 'pm', 89),
                    ("POSCAR_36C_p6m", 36, 'p6m', 0),
                    ("POSCAR_36C_cm", 36,'cm', 0)]

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


@pytest.fixture(scope="function", params=_EXAMPLE_POSCARs)
def slab_and_expectations(request):
    filename, expected_n_atoms, expected_pg, offset_at = request.param
    file_path = Path(__file__).parent / "fixtures" / "POSCARs" / filename
    pos_slab = readPOSCAR(str(file_path))
    return (pos_slab, expected_n_atoms, expected_pg, offset_at)

@pytest.fixture(scope="function")
def slab_pg_rp(slab_and_expectations):
    slab, *_ = slab_and_expectations
    rp = tl.Rparams()
    slab.fullUpdate(rp)
    pg = findSymmetry(slab, rp, output=False)
    enforceSymmetry(slab, rp)
    return slab, pg, rp


class TestPOSCARRead:
    def test_read_in_atoms(self, slab_and_expectations):
        slab, *_ = slab_and_expectations
        assert len(slab.atlist) > 0

    def test_n_atom_correct(self, slab_and_expectations):
        slab, expected_n_atoms, *_ = slab_and_expectations
        assert len(slab.atlist) == expected_n_atoms


class TestPOSCARSymmetry(TestPOSCARRead):
    def test_any_pg_found(self, slab_pg_rp):
        _, slab_pg, _ = slab_pg_rp
        assert slab_pg != 'unknown'

    def test_pg_correct(self, slab_and_expectations, slab_pg_rp):
        _, _, expected_pg, _ = slab_and_expectations
        _, slab_pg, _ = slab_pg_rp
        assert slab_pg == expected_pg

    @pytest.mark.parametrize("displacement", [(4, (np.array([0.2, 0, 0]),)),
                                            (4, (np.array([0, 0.2, 0]),)),
                                            (4, (np.array([0, 0, 0.2]),)),
                                            ])
    def test_preserve_symmetry_with_displacement(self, displacement, slab_and_expectations, slab_pg_rp):
        slab, _, expected_pg, offset_at = slab_and_expectations
        _, _, rp = slab_pg_rp
        sl_copy = deepcopy(slab)
        
        # manually assign displacements
        sl_copy.atlist[offset_at].assignDisp(*displacement)

        for at in sl_copy.atlist:
            disp = at.disp_geo_offset['all'][0]
            at.cartpos += disp
        sl_copy.getFractionalCoordinates()

        assert findSymmetry(sl_copy, rp) == expected_pg


@pytest.fixture(scope='function')
def manual_slab_3_atoms():
    slab = tl.Slab()
    slab.ucell = np.diag([3., 4., 5.])
    positions = (np.array([-0.25, 0, 0]),
                 np.array([0.00, 0, 0]),
                 np.array([0.25, 0, 0]))
    slab.atlist = [tl.Atom('C', pos, i+1, slab)
                   for i, pos in enumerate(positions)]
    param = tl.Rparams()
    slab.fullUpdate(param)
    return slab


@pytest.fixture()
def manual_slab_1_atom_trigonal():
    slab = tl.Slab()
    slab.ucell = np.array([[ 1, 0, 0],
                           [-2, 3, 0],
                           [ 1, 2, 3]],dtype=float)
    slab.atlist = [tl.Atom('C', np.array([0.2, 0.7, 0.1]), 1, slab),]  # "random" position
    param = tl.Rparams()
    slab.fullUpdate(param)
    return slab

class TestSlabTransforms:
    def test_mirror(self, manual_slab_3_atoms):
        slab = manual_slab_3_atoms
        mirrored_slab = deepcopy(slab)
        symplane = tl.classes.slab.SymPlane(pos=(0,0),
                                            dr=np.array([0,1]),
                                            abt=slab.ucell.T[:2,:2])
        mirrored_slab.mirror(symplane)
        mirrored_slab.collapseCartesianCoordinates()
        assert all(at.isSameXY(mir_at.cartpos[:2])
                for at, mir_at in
                zip(slab.atlist, reversed(mirrored_slab.atlist)))

    def test_180_rotation(self, manual_slab_3_atoms):
        slab = manual_slab_3_atoms
        rotated_slab = deepcopy(slab)
        rotated_slab.rotateAtoms((0,0), order=2)
        rotated_slab.collapseCartesianCoordinates()
        assert all(at.isSameXY(mir_at.cartpos[:2])
                for at, mir_at in
                zip(slab.atlist, reversed(rotated_slab.atlist)))


@pytest.fixture()
def run_phaseshift(slab_pg_rp, tmp_path_factory):
    slab, _,  param = slab_pg_rp
    param.workdir = tmp_path_factory.mktemp(basename="phaseshifts", numbered=True)
    # run EEASISSS
    firstline, phaseshift = runPhaseshiftGen_old(slab,
                                                 param,
                                                 psgensource = TENSORLEED_PATH/'EEASiSSS.x',
                                                 excosource=TENSORLEED_PATH/'seSernelius',
                                                 atdenssource=TENSORLEED_PATH/'atom_density_files')
    return param, slab, firstline, phaseshift


class TestPhaseshifts:
    def test_phaseshifts_firstline_not_empty(self, run_phaseshift):
        _, _, firstline, _ = run_phaseshift
        assert firstline

    def test_phaseshifts_firstline_len(self, run_phaseshift):
        _, _, firstline, _ = run_phaseshift
        potential_param = firstline.split()
        assert len(potential_param) >= 4


    def test_phaseshift_log_exists(self, run_phaseshift):
        param, _, _, _ = run_phaseshift
        assert len(list(param.workdir.glob('phaseshift*.log'))) > 0


    def test_write_phaseshifts(self, run_phaseshift):
        from tleedmlib.files.phaseshifts import writePHASESHIFTS
        param, _, firstline, phaseshift = run_phaseshift
        writePHASESHIFTS(firstline, phaseshift, file_path=param.workdir/'PHASESHIFTS')
        assert len(list(param.workdir.glob('PHASESHIFTS'))) > 0


    def test_phaseshifts_not_empty(self, run_phaseshift):
        _, _, _, phaseshift = run_phaseshift
        assert len(phaseshift) > 0


# Slab Matrix operations

def test_rotation_on_trigonal_slab(manual_slab_1_atom_trigonal):
    rot_15 = np.array([[ 0.96592583, -0.25881905,  0.        ],
                       [ 0.25881905,  0.96592583,  0.        ],
                       [ 0.        ,  0.        ,  1.        ]])
    expected_cell = np.array([[ 0.44828774, -2.1906707 ,  1.        ],
                              [ 0.77645714,  2.89777748,  2.        ],
                              [ 0.        ,  0.        ,  3.        ]])
    expected_atom_cartpos = [0.63317754, 1.5903101]
    slab = manual_slab_1_atom_trigonal
    slab.apply_matrix_transformation(rot_15)
    assert np.allclose(slab.ucell.T, expected_cell)
    assert np.allclose(slab.atlist[0].cartpos[:2], expected_atom_cartpos)