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
from viperleed.tleedmlib.files import parameters, poscar

_FIXTURES_PATH = Path('tests/fixtures/')

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
