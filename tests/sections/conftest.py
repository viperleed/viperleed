"""sections/conftest.py

Created on 2023-07-26

@author: Alexander M. Imre
"""

import pytest
import sys
import os
from pathlib import Path
from copy import deepcopy
import numpy as np

vpr_path = str(Path(__file__).parent.parent.parent.parent)
if os.path.abspath(vpr_path) not in sys.path:
    sys.path.append(os.path.abspath(vpr_path))



from viperleed.tleedmlib.files.poscar import readPOSCAR
from viperleed.tleedmlib.classes.rparams import Rparams


TENSERLEED_TEST_VERSIONS = ('1.71', '1.72', '1.73', '1.74')

AG_100_DISPLACEMENTS_NAMES = ['DISPLACEMENTS_z', 'DISPLACEMENTS_vib', 'DISPLACEMENTS_z+vib']
AG_100_DELTAS_NAMES = ['Deltas_z.zip', 'Deltas_vib.zip', 'Deltas_z+vib.zip']

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



@pytest.fixture()
def fe3o4_bulk_slab(poscar_path):
    file_name = "POSCAR_Fe3O4_(001)_cod1010369"
    file_path = poscar_path(file_name)
    slab = readPOSCAR(str(file_path))
    param = Rparams()
    param.LAYER_CUTS = [0.1, 0.2, '<', 'dz(1.0)']
    param.N_BULK_LAYERS = 2
    param.SYMMETRY_EPS =0.3
    param.SYMMETRY_EPS_Z = 0.3
    param.BULK_REPEAT = np.array([-0.0, -4.19199991, 4.19199991])
    slab.fullUpdate(param)
    bulk_slab = slab.makeBulkSlab(param)
    return slab, bulk_slab, param


@pytest.fixture()
def fe3o4_thick_bulk_slab(fe3o4_bulk_slab):
    slab, thin_bulk, param = fe3o4_bulk_slab
    thick_bulk = thin_bulk.doubleBulkSlab()
    return slab, thick_bulk, param
