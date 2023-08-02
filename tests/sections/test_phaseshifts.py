"""sections/test_phaseshifts.py

Created on 2023-07-28

@author: Alexander M. Imre
"""

import pytest
import sys
import os
from pathlib import Path
from copy import deepcopy
import numpy as np


from viperleed.lib.files.displacements import readDISPLACEMENTS, readDISPLACEMENTS_block
from viperleed.lib.files.poscar import readPOSCAR
from viperleed.lib.files.vibrocc import readVIBROCC
from viperleed.lib.symmetry import findSymmetry, enforceSymmetry
from viperleed.lib.psgen import runPhaseshiftGen_old
from viperleed.lib.classes.atom import Atom
from viperleed.lib.classes.rparams import Rparams
from viperleed.lib.classes.slab import Slab


@pytest.fixture()
def run_phaseshift(slab_pg_rp, source_path, tmp_path_factory):
    slab, _,  param = slab_pg_rp
    param.workdir = tmp_path_factory.mktemp(basename="phaseshifts",
                                            numbered=True)
    param.source_dir = source_path
    param.THEO_ENERGIES = [10, 500, 3]
    # run EEASISSS
    firstline, phaseshift = runPhaseshiftGen_old(
        slab,
        param,
        psgensource ='EEASiSSS.x',
        excosource='seSernelius',
        atdenssource='atom_density_files'
    )
    return param, slab, firstline, phaseshift



def test_phaseshifts_firstline_not_empty(run_phaseshift):
    _, _, firstline, _ = run_phaseshift
    assert firstline

def test_phaseshifts_firstline_len(run_phaseshift):
    _, _, firstline, _ = run_phaseshift
    potential_param = firstline.split()
    assert len(potential_param) >= 4


def test_phaseshift_log_exists(run_phaseshift):
    param, _, _, _ = run_phaseshift
    assert len(list(param.workdir.glob('phaseshift*.log'))) > 0


def test_write_phaseshifts(run_phaseshift):
    from tleedmlib.files.phaseshifts import writePHASESHIFTS
    param, _, firstline, phaseshift = run_phaseshift
    writePHASESHIFTS(firstline, phaseshift, file_path=param.workdir/'PHASESHIFTS')
    assert len(list(param.workdir.glob('PHASESHIFTS'))) > 0


def test_phaseshifts_not_empty(run_phaseshift):
    _, _, _, phaseshift = run_phaseshift
    assert len(phaseshift) > 0