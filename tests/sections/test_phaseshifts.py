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

vpr_path = str(Path(__file__).parent.parent.parent.parent)
if os.path.abspath(vpr_path) not in sys.path:
    sys.path.append(os.path.abspath(vpr_path))

from viperleed.tleedmlib.files.displacements import (readDISPLACEMENTS,
                                                     readDISPLACEMENTS_block)
from viperleed.tleedmlib.files.vibrocc import readVIBROCC
from viperleed.tleedmlib.symmetry import findSymmetry, enforceSymmetry
from viperleed.tleedmlib.psgen import runPhaseshiftGen_old
from viperleed.tleedmlib.classes.atom import Atom
from viperleed.tleedmlib.classes.rparams import Rparams
from viperleed.tleedmlib.classes.slab import Slab



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
    from tleedmlib.files.phaseshifts import writePHASESHIFTS
    param, _, firstline, phaseshift = run_phaseshift
    writePHASESHIFTS(firstline, phaseshift, file_path=param.workdir/'PHASESHIFTS')
    assert len(list(param.workdir.glob('PHASESHIFTS'))) > 0


@pytest.mark.xfail(reason="Relative path problem. Fixed on installable branch.")
def test_phaseshifts_not_empty(run_phaseshift):
    _, _, _, phaseshift = run_phaseshift
    assert len(phaseshift) > 0