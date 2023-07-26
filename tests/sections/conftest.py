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


from viperleed.tleedmlib.files.displacements import readDISPLACEMENTS, readDISPLACEMENTS_block
from viperleed.tleedmlib.files.poscar import readPOSCAR
from viperleed.tleedmlib.files.vibrocc import readVIBROCC
from viperleed.tleedmlib.symmetry import findSymmetry, enforceSymmetry
from viperleed.tleedmlib.psgen import runPhaseshiftGen_old
from viperleed.tleedmlib.classes.atom import Atom
from viperleed.tleedmlib.classes.rparams import Rparams
from viperleed.tleedmlib.classes.slab import Slab


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

@pytest.fixture(scope='function')
def manual_slab_3_atoms():
    slab = Slab()
    slab.ucell = np.diag([3., 4., 5.])
    positions = (np.array([-0.25, 0, 0]),
                 np.array([0.00, 0, 0]),
                 np.array([0.25, 0, 0]))
    slab.atlist = [Atom('C', pos, i+1, slab)
                   for i, pos in enumerate(positions)]
    param = tl.Rparams()
    slab.fullUpdate(param)
    return slab


@pytest.fixture()
def manual_slab_1_atom_trigonal():
    slab = Slab()
    slab.ucell = np.array([[ 1, 0, 0],
                           [-2, 3, 0],
                           [ 1, 2, 3]],dtype=float)
    slab.atlist = [Atom('C', np.array([0.2, 0.7, 0.1]), 1, slab),]  # "random" position
    param = tl.Rparams()
    slab.fullUpdate(param)
    return slab

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

@pytest.fixture()
def fe3o4_bulk_slab():
    file_name = "POSCAR_Fe3O4_(001)_cod1010369"
    file_path = POSCAR_PATHS / file_name
    slab = readPOSCAR(str(file_path))
    param = tl.Rparams()
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

@pytest.fixture(scope="function")
def atom_with_disp_and_offset():
    slab = readPOSCAR(POSCAR_PATHS / "POSCAR_STO(100)-4x1")
    atom = slab.atlist[0]
    el = atom.el
    atom.disp_geo[el] = [-0.2, 0.0, 0.2]
    atom.disp_vib[el] = [-0.1, 0.0, 0.1]
    atom.disp_occ[el] = [0.7, 0.8, 0.9, 1.0]
    return atom

@pytest.fixture()
def ag100_slab_param():
    slab = readPOSCAR(POSCAR_PATHS / "POSCAR_Ag(100)")
    param = Rparams()
    param.N_BULK_LAYERS = 1
    slab.fullUpdate(param)
    return slab, param

@pytest.fixture()
def ag100_slab_with_displacements_and_offsets(ag100_slab_param):
    slab, param = ag100_slab_param
    vibrocc_path = INPUTS_ORIGIN / "Ag(100)" / "mergeDisp" / "VIBROCC"
    displacements_path = INPUTS_ORIGIN / "Ag(100)" / "mergeDisp" / "DISPLACEMENTS_mixed"
    readVIBROCC(param, slab, str(vibrocc_path))
    readDISPLACEMENTS(param, str(displacements_path))
    readDISPLACEMENTS_block(param, slab, param.disp_blocks[param.search_index])
    return slab, param
