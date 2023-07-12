"""Module helpers of viperleed.tests.

Created on 2023-02-28

@author: Michele Riva

Contains some useful general definitions that can be used when creating
or running tests.
"""


# Think about a decorator for injecting fixtures.
# Some ideas at
# https://github.com/pytest-dev/pytest/issues/2424
# https://github.com/pytest-dev/pytest/issues/6322
# https://github.com/nteract/testbook/issues/4

from pathlib import Path

import pytest

from viperleed.tleedmlib import symmetry
from viperleed.tleedmlib.classes.rparams import Rparams
from viperleed.tleedmlib.files import parameters, poscar


_FIXTURES_PATH = Path('tests/fixtures/')
_POSCARs_PATH = _FIXTURES_PATH / 'POSCARs'

_EXAMPLE_POSCAR_EXPECTATIONS = [("POSCAR_Ag(100)", 6, 'p4m', 0),
                                ("POSCAR_STO(100)-4x1", 136, 'pm', 0),
                                ("POSCAR_TiO2", 540, 'pmm', -1),
                                ("POSCAR_diamond", 96, 'pm', 89),
                                ("POSCAR_36C_p6m", 36, 'p6m', 0),
                                ("POSCAR_36C_cm", 36,'cm', 0),
                                ("POSCAR_Fe3O4_SCV", 83, 'cmm', 50)]            #TODO: Phaseshift generation fails. Why? @Fkraushofer (worked in fkpCurie:Florian_OldLocalTests/Fe3O4-001-SCV/history/t000.r013_211220-133452)

_EXAMPLE_POSCARs = [file.name for file in _POSCARs_PATH.glob('POSCAR*')]


@pytest.fixture()
def ag100_parameters_example():
    # read Ag(100) POSCAR and PARAMETERS files
    slab = poscar.readPOSCAR(_FIXTURES_PATH / 'Ag(100)' / 'initialization' / 'POSCAR')
    rpars = parameters.readPARAMETERS(_FIXTURES_PATH / 'Ag(100)' / 'initialization' / 'PARAMETERS')
    # interpret PARAMETERS file
    interpreter = parameters.ParameterInterpreter(rpars)
    interpreter.interpret(slab)
    symmetry.findSymmetry(slab, rpars)
    symmetry.findBulkSymmetry(slab, rpars)
    return (rpars, slab)

@pytest.fixture(scope="function", params=_EXAMPLE_POSCARs)
def example_poscars(request):
    file_path = _POSCARs_PATH / request.param
    slab = poscar.readPOSCAR(str(file_path))
    return slab

@pytest.fixture(scope="function", params=_EXAMPLE_POSCAR_EXPECTATIONS)
def slab_and_expectations(request):
    filename, expected_n_atoms, expected_pg, offset_at = request.param
    file_path = _POSCARs_PATH / filename
    pos_slab = poscar.readPOSCAR(str(file_path))
    return (pos_slab, expected_n_atoms, expected_pg, offset_at)

@pytest.fixture(scope="function")
def slab_pg_rp(slab_and_expectations):
    slab, *_ = slab_and_expectations
    rp = Rparams()
    slab.fullUpdate(rp)
    pg = symmetry.findSymmetry(slab, rp, output=False)
    symmetry.enforceSymmetry(slab, rp)
    return slab, pg, rp