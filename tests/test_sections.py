import unittest   # The test framework
from parameterized import parameterized, parameterized_class
import shutil, tempfile
import sys
import os
from pathlib import Path
from zipfile import ZipFile

vpr_path = str(Path(__file__).parent.parent.parent)
if os.path.abspath(vpr_path) not in sys.path:
    sys.path.append(os.path.abspath(vpr_path))

import viperleed
import viperleed.tleedmlib
from viperleed.tleedm import run_tleedm

ALWAYS_REQUIRED_FILES = ('PARAMETERS', 'EXPBEAMS.csv', 'POSCAR')


class TestTleedmFromFiles(unittest.TestCase):
    @classmethod
    def setUpClass(cls, surface_name, required_files, section_dirs):
        cls.surface_name = surface_name
        cls.required_files = set(ALWAYS_REQUIRED_FILES)
        cls.section_dirs = section_dirs
        # add additional files
        cls.required_files.update(required_files)

        # Create a temporary directory to run in
        cls.test_dir = Path(tempfile.mkdtemp())
        cls.work_path = str(cls.test_dir / "work")
        os.makedirs(cls.work_path, exist_ok=True)
        
        cls.fixtures_dir = Path(__file__).parent / "fixtures" / cls.surface_name
        cls.input_files_dirs = []

        for dir in section_dirs:
            cur_input_dir = cls.fixtures_dir / dir
            cls.input_files_dirs.append(cur_input_dir)
            shutil.copytree(cur_input_dir, cls.test_dir, dirs_exist_ok=True)
            shutil.copytree(cur_input_dir, cls.work_path, dirs_exist_ok=True)

        # get current work directory
        cls.home = os.getcwd()

    @classmethod
    def tearDownClass(cls):
        # Remove the directory after the test
        #shutil.rmtree(cls.test_dir)
        os.chdir(cls.home)


    def expected_file_exists(self, expected_file):
        expected_path = Path(self.work_path) / expected_file
        return expected_path.exists()


class TestInitializationAg(TestTleedmFromFiles):
    @classmethod
    def setUpClass(cls):
        super().setUpClass(surface_name="Ag(100)",
                           required_files=('PHASESHIFTS',),
                           section_dirs=("initialization",))


class TestCaseInitializationAgRun(TestInitializationAg):
    def test_work_path_exists(self):
        """Check that work_path was created properly."""
        assert Path(self.work_path).is_dir()

    def test_files_copied_correctly(self):
        """Check if all files were copied correctly."""
        input_files = [file.name for file in self.test_dir.glob('*')]
        source_copied = all(file in input_files for file in self.required_files)
        input_files = [file.name for file in Path(self.work_path).glob('*')]
        work_copied = all(file in input_files for file in self.required_files)
        assert source_copied and work_copied


class TestCaseInitializationAgFileChecks(TestInitializationAg):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        os.chdir(cls.work_path)
        
        # run tleedm
        cls.exit_code = run_tleedm(source=str(Path(vpr_path) / "viperleed"))
        cls.files_after_run = [file.name for file in Path(cls.work_path).glob('*')]

    def test_initialization_runs(self):
        """Test if initialization gives exit code 0."""
        assert self.exit_code == 0

    def test_files_present(self):
        """Checks if files are present after initialization"""
        files_to_test = ('IVBEAMS', 'BEAMLIST', 'VIBROCC')
        for file in files_to_test:
            with self.subTest(file=file):
                self.assertIn(file, self.files_after_run)


class AgRefCalc(TestInitializationAg):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        os.chdir(cls.work_path)
        # run tleedm
        cls.exit_code = run_tleedm(source=str(Path(vpr_path) / "viperleed"),
                   preset_params={"RUN":[0, 1],
                                  "TL_VERSION":1.73})

        cls.files_after_run = [file.name for file in Path(cls.work_path).glob('*')]

    def test_THEOBEAMS_csv_present(self):
        """Check if THEOBEAMS.csv is present
        """
        assert 'THEOBEAMS.csv' in self.files_after_run


    def test_exit_code_0(self):
        """Checks if tleedm exits with exit code 0.
        """
        assert self.exit_code == 0


@parameterized_class(('displacements_name'), [('DISPLACEMENTS_z',), ('DISPLACEMENTS_vib',), ('DISPLACEMENTS_z+vib',)])
class AgDeltas(TestTleedmFromFiles):
    @classmethod
    def setUpClass(cls):
        super().setUpClass(surface_name="Ag(100)",
                           required_files=('PHASESHIFTS','Tensors/Tensors001.zip'),
                           section_dirs=("initialization", "deltas"))
        os.chdir(cls.work_path)

        # copy DISPLACEMENTS
        disp_f = cls.fixtures_dir / "displacements" / cls.displacements_name
        shutil.copy(disp_f, Path(cls.work_path) / 'DISPLACEMENTS')

        # run Delta calculation
        cls.exit_code = run_tleedm(source=str(Path(vpr_path) / "viperleed"),
                   preset_params={"RUN":[0, 2],
                                  "TL_VERSION":1.73})

    def test_delta_input_written(self):
        assert self.expected_file_exists("delta-input")


    def test_exit_code_0(self):
        assert self.exit_code == 0


    def test_deltas_zip_created(self):
        assert self.expected_file_exists(Path("Deltas") / "Deltas_001.zip")


@parameterized_class(('displacements_name', 'deltas_name'),
                     [('DISPLACEMENTS_z', 'Deltas_z.zip'),
                      ('DISPLACEMENTS_vib', 'Deltas_vib.zip'),
                      ('DISPLACEMENTS_z+vib', 'Deltas_z+vib.zip')])
class AgSearch(TestTleedmFromFiles):
    @classmethod
    def setUpClass(cls):
        super().setUpClass(surface_name="Ag(100)",
                           required_files=('PHASESHIFTS'),
                           section_dirs=("initialization", "deltas", "search"))
        os.chdir(cls.work_path)

        # copy DISPLACEMENTS
        disp_f = cls.fixtures_dir / "displacements" / cls.displacements_name
        shutil.copy(disp_f, Path(cls.work_path) / 'DISPLACEMENTS')


        # copy Deltas
        deltas_f = cls.fixtures_dir / "search" / "Deltas" / cls.deltas_name
        shutil.copy(deltas_f, Path(cls.test_dir) / "Deltas" / "Deltas_001.zip")
        #ZipFile(deltas_f, 'r').extractall(cls.work_path)
        pass

        # run Delta calculation
        cls.exit_code = run_tleedm(source=str(Path(vpr_path) / "viperleed"),
                   preset_params={"RUN":[0, 3],
                                  "TL_VERSION":1.73})


    @parameterized.expand(('SD.TL', 'control.chem'))
    def test_search_file_exists(self, req_file):
        print(self.work_path)
        self.assertTrue(self.expected_file_exists(req_file))


# POSCAR and symmetry tests
class TestPOSCAR(unittest.TestCase):
    @classmethod
    def setUpClass(cls, filename):
        # import reading and wrting POSCARS
        from tleedmlib.files.poscar import readPOSCAR, writePOSCAR

        file_path = Path(__file__).parent / "fixtures" / "POSCARs" / filename
        cls.slab = readPOSCAR(str(file_path))

    def test_read_in_atoms(self):
        assert len(self.slab.atlist) > 0

    def test_n_atom_correct(self):
        assert len(self.slab.atlist) == self.expected_n_atoms

# for tests with Symmetry recognition
class TestPOSCARSymmetry(TestPOSCAR):
    @classmethod
    def setUpClass(cls, filename):
        super().setUpClass(filename)
        cls.filename = filename

        import tleedmlib as tl
        from tleedmlib.symmetry import findSymmetry, findBulkSymmetry
        # create dummy Rparams object
        cls.rp = tl.Rparams()
        cls.slab.fullUpdate(cls.rp)
        cls.pg = findSymmetry(cls.slab, cls.rp, output=False)

    def test_pg_found(self):
        assert self.pg != 'unknown'

    def test_pg_correct(self):
        self.assertEqual(self.pg, self.expected_pg,
                         f"POSCAR file {self.filename}: "
                         f"expected planegroup {self.expected_pg}, found {self.pg}.")

class read_Ag100(TestPOSCARSymmetry):
    @classmethod
    def setUpClass(cls):
        cls.filename = "POSCAR_Ag(100)"
        cls.expected_n_atoms = 6
        cls.expected_pg = 'p4m'
        super().setUpClass(cls.filename)

    def test_all_atoms_are_Ag(self):
        atoms = self.slab.atlist
        atom_elems = [atom.el for atom in atoms]
        atom_is_Ag = [el == 'Ag' for el in atom_elems]
        assert all(atom_is_Ag)

class read_STO_4x1(TestPOSCARSymmetry):
    @classmethod
    def setUpClass(cls):
        cls.filename = "POSCAR_STO(100)-4x1"
        cls.expected_n_atoms = 136
        cls.expected_pg = 'pm'
        super().setUpClass(cls.filename)

# try huge unit cell
class read_TiO2(TestPOSCARSymmetry):
    @classmethod
    def setUpClass(cls):
        cls.filename = "POSCAR_TiO2"
        cls.expected_n_atoms = 540
        cls.expected_pg = 'pmm'
        super().setUpClass(cls.filename)

class read_diamond(TestPOSCARSymmetry):
    @classmethod
    def setUpClass(cls):
        cls.filename = "POSCAR_diamond"
        cls.expected_n_atoms = 96
        cls.expected_pg = 'pm'
        super().setUpClass(cls.filename)

class read_graphene(TestPOSCARSymmetry):
    @classmethod
    def setUpClass(cls):
        cls.filename = "POSCAR_graphene"
        cls.expected_n_atoms = 36
        cls.expected_pg = 'pmm'
        super().setUpClass(cls.filename)


if __name__ == '__main__':
    unittest.main()
