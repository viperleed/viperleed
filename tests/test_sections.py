import unittest   # The test framework
import shutil, tempfile
import sys
import os
from pathlib import Path

vpr_path = str(Path(__file__).parent.parent.parent)
if os.path.abspath(vpr_path) not in sys.path:
    sys.path.append(os.path.abspath(vpr_path))

import viperleed
import viperleed.tleedmlib
from viperleed.tleedm import run_tleedm

ALWAYS_REQUIRED_FILES = ('PARAMETERS', 'EXPBEAMS.csv', 'POSCAR')


class TesTTleedmFromFiles(unittest.TestCase):
    @classmethod
    def setUpClass(cls, surface_name, required_files):
        cls.surface_name = surface_name
        cls.required_files = set(ALWAYS_REQUIRED_FILES)
        # add additional files
        cls.required_files.update(required_files)

        cls.input_files_dir = Path(__file__).parent / "fixtures" / cls.surface_name / "initialization"
        
        # Create a temporary directory to run in
        cls.test_dir = Path(tempfile.mkdtemp())
        cls.work_path = str(cls.test_dir / "work")
        os.makedirs(cls.work_path, exist_ok=True)
        
        shutil.copytree(cls.input_files_dir, cls.test_dir, dirs_exist_ok=True)
        shutil.copytree(cls.input_files_dir, cls.work_path, dirs_exist_ok=True)
        
        # get current work directory
        cls.home = os.getcwd()

    @classmethod
    def tearDownClass(cls):
        # Remove the directory after the test
        shutil.rmtree(cls.test_dir)
        os.chdir(cls.home)

class TestInitializationAg(TesTTleedmFromFiles):
    @classmethod
    def setUpClass(cls):
        super().setUpClass(surface_name="Ag(100)", required_files=('PHASESHIFTS',))


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
        run_tleedm(source=str(Path(vpr_path) / "viperleed"),
                   preset_params={"RUN":[0, 1],
                                  "TL_VERSION":1.73})

        cls.files_after_run = [file.name for file in Path(cls.work_path).glob('*')]

    def test_THEOBEAMS_csv_present(self):
        print(self.files_after_run)
        assert 'THEOBEAMS.csv' in self.files_after_run



if __name__ == '__main__':
    unittest.main()
