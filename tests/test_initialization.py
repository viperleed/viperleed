"""Test section initialization.

Created on 2023-07-19

@author: Alexander M. Imre
"""

from pathlib import Path
import os
import sys

import pytest

vpr_path = str(Path(__file__).parent.parent.parent)
if os.path.abspath(vpr_path) not in sys.path:
    sys.path.append(os.path.abspath(vpr_path))


from viperleed.tests.helpers import (ag100_parameters_example,
                                     BaseTleedmFilesSetup,
                                     SOURCE_STR,
                                     TestSetup,
                                     init_files,  # fixtures must be imported explicitly!
                                     INPUTS_ORIGIN
                                     )

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


class TestInitialization(TestSetup):
    def test_exit_code_0(self, init_files):
        """Test if initialization gives exit code 0."""
        assert init_files.exit_code == 0

    @pytest.mark.parametrize('expected_file', (('IVBEAMS'), ('BEAMLIST'), ('VIBROCC')))
    def test_init_files_present(self, init_files, expected_file):
        """Checks if files are present after initialization"""
        assert init_files.expected_file_exists(expected_file)

    def test_PARAMETERS_was_updated(self, init_files):
        """Checks if PARAMETERS file was updated."""
        assert init_files.expected_file_exists('PARAMETERS')
        with open(init_files.work_path / 'PARAMETERS', 'r') as param_file:
            param_content = param_file.read()
        assert "line commented out automatically" in param_content

    @pytest.mark.parametrize('expected_file', (('IVBEAMS'), ('BEAMLIST'), ('VIBROCC')))
    def test_init_with_element_rename(self, ag100_rename_ax, expected_file):
        assert ag100_rename_ax.expected_file_exists(expected_file)