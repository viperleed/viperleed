"""Test configuration for viperleed.tests.

Created on 2023-02-28

@author: Michele Riva (@michele-riva)
@author: Alexander M. Imre (@amimre)

Defines fixtures and fixture factories used in multiple tests.

Fixtures
--------
ag100
    A Ag(100) slab, an Rparams, and a TestInfo.
ag100_with_displacements_and_offsets
    A Ag(100) Slab and an Rparams, after reading a DISPLACEMENTS block.
data_path
    Path to the top-level folder containing test data.
make_poscar (factory)
    Return a Slab from POSCAR, an Rparams and a TestInfo from a TestInfo.
poscars_path
    Path to the data directory containing POSCAR files.
re_match (factory)
    Return a match object from a pattern and a string.
run_phaseshift
    An Rparams, a Slab, and the results of generating phaseshifts
    with them.
tensorleed_path
    Path to the top-level tree with tensor-LEED source code.
"""

import os
from pathlib import Path
import re
import sys

import pytest
import pytest_cases

VPR_PATH = str(Path(__file__).resolve().parents[2])
if VPR_PATH not in sys.path:
    sys.path.append(VPR_PATH)

# pylint: disable=wrong-import-position
# Will be fixed in installable version
from viperleed.tleedmlib.files import displacements, vibrocc
from viperleed.tleedmlib import psgen

from .helpers import TEST_DATA, POSCAR_PATH
from .helpers import CaseTag, exclude_tags
from . import poscar_slabs
# pylint: enable=wrong-import-position


TENSORLEED_PATH = Path(VPR_PATH) / 'viperleed' / 'tensorleed'                   # TODO: this will need to be dynamic!


@pytest.fixture(scope='session')
def re_match():  # This is actually a fixture factory
    """Return a re.match object from a pattern and a string."""
    def _match(pattern, string):
        return re.match(pattern, string)
    return _match


@pytest.fixture(scope='session')
def poscars_path():
    """Return the Path to the directory containing POSCAR files."""
    return POSCAR_PATH


@pytest.fixture(scope='session', name='data_path')
def fixture_data_path():
    """Return the Path to the top-level folder containing test data."""
    return TEST_DATA


@pytest.fixture(scope='session', name='tensorleed_path')
def fixture_tensorleed_path():
    """Return the Path to the top-level tree with tensor-LEED source code."""
    return TENSORLEED_PATH


@pytest.fixture(name='make_poscar', scope='session')
def factory_make_poscar():
    """Return a POSCAR from a TestInfo."""
    def _make(info):
        return poscar_slabs.CasePOSCARSlabs().case_poscar(info)
    return _make


@pytest.fixture(name='ag100')
def fixture_ag100(make_poscar):
    """Return a Ag(100) slab, an Rparams, and a TestInfo."""
    return make_poscar(poscar_slabs.AG_100)


@pytest.fixture
def ag100_with_displacements_and_offsets(ag100, data_path):
    """Return a Slab and Rparams after reading a DISPLACEMENTS block."""
    slab, param, *_ = ag100

    inputs_path = data_path / 'Ag(100)/mergeDisp'
    vibrocc_path = inputs_path / 'VIBROCC'
    displacements_path = inputs_path / 'DISPLACEMENTS_mixed'
    vibrocc.readVIBROCC(param, slab, str(vibrocc_path))
    displacements.readDISPLACEMENTS(param, str(displacements_path))
    displacements.readDISPLACEMENTS_block(param, slab, param.disp_blocks[0])
    return slab, param


# Notice that we need to exclude POSCARs without information as some
# are known to make this fail because of "too small interlayer spacing"
# For Fe3O4 this is because parameter presets are needed. For graphene,
# only one layer is present in the POSCAR.
_PHASESHIFT_SETTINGS = {'cases': poscar_slabs.CasePOSCARSlabs,
                        'filter': exclude_tags(CaseTag.NO_INFO),
                        'scope': 'session'}

@pytest_cases.fixture(scope='session')
@pytest_cases.parametrize_with_cases('args', **_PHASESHIFT_SETTINGS)
def run_phaseshift(args, tensorleed_path, tmp_path_factory):
    """Execute a PHASESHIFTS calculation.

    Parameters
    ----------
    args : pytest_cases.Case
        Slab, Rparams and other (unused) info.
    tensorleed_path : pytest.fixture
        The path to the directory containing the Fortran code.
    tmp_path_factory : pytest.fixture
        To run the PHASESHIFTS calculation in a fresh directory.

    Yields
    ------
    param : Rparams
        The PARAMETERS used during the PHASESHIFTS calculation.
    slab : Slab
        The Slab for which PHASESHIFTS were calculated.
    firstline : str
        The first line of the PHASESHIFTS file, that contains the
        coefficients for the real part of the inner potential.
    phaseshift : list
        The PHASESHIFTS that were generated.
    """
    slab, param, *_ = args
    param.source_dir = tensorleed_path
    param.workdir = tmp_path_factory.mktemp(basename='phaseshifts',
                                            numbered=True)
    param.initTheoEnergies()
    executable = 'EEASiSSS' + ('.exe' if 'nt' in os.name else '')               # TODO: does this cover it or should we use 'win' in sys.platform()?

    # run EEASISSS in the temporary directory
    home = Path()
    os.chdir(param.workdir)
    results = psgen.runPhaseshiftGen_old(slab, param, psgensource=executable)
    yield param, slab, *results
    os.chdir(home)
