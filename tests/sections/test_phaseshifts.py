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

from viperleed.calc.files.displacements import (readDISPLACEMENTS,
                                                readDISPLACEMENTS_block)
from viperleed.calc.files.phaseshifts import writePHASESHIFTS
from viperleed.calc.files.vibrocc import readVIBROCC
from viperleed.calc.lib.symmetry import findSymmetry, enforceSymmetry
from viperleed.calc.files.psgen import runPhaseshiftGen_old
from viperleed.calc.classes.atom import Atom
from viperleed.calc.classes.rparams import Rparams
from viperleed.calc.classes.slab import Slab


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


@pytest.mark.xfail(reason="Relative path problem. Fixed on installable branch.")
def test_phaseshifts_firstline_not_empty(run_phaseshift):
    _, _, firstline, _ = run_phaseshift
    assert firstline


@pytest.mark.xfail(reason="Relative path problem. Fixed on installable branch.")
def test_phaseshifts_firstline_len(run_phaseshift):
    _, _, firstline, _ = run_phaseshift
    potential_param = firstline.split()
    assert len(potential_param) >= 4


@pytest.mark.xfail(reason="Relative path problem. Fixed on installable branch.")
def test_phaseshift_log_exists(run_phaseshift):
    param, _, _, _ = run_phaseshift
    assert len(list(param.workdir.glob('phaseshift*.log'))) > 0


@pytest.mark.xfail(reason="Relative path problem. Fixed on installable branch.")
def test_write_phaseshifts(run_phaseshift):
    param, _, firstline, phaseshift = run_phaseshift
    writePHASESHIFTS(firstline, phaseshift, file_path=param.workdir/'PHASESHIFTS')
    assert len(list(param.workdir.glob('PHASESHIFTS'))) > 0


@pytest.mark.xfail(reason="Relative path problem. Fixed on installable branch.")
def test_phaseshifts_not_empty(run_phaseshift):
    _, _, _, phaseshift = run_phaseshift
    assert len(phaseshift) > 0