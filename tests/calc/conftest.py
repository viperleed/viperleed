"""Test configuration for viperleed.tests.calc.

Defines fixtures and fixture factories used in multiple tests.

Fixtures
--------
ag100
    A Ag(100) slab, an Rparams, and a TestInfo.
ag100_with_displacements_and_offsets
    A Ag(100) Slab and an Rparams, after reading a DISPLACEMENTS block.
make_poscar (factory)
    Return a Slab from POSCAR, an Rparams and a TestInfo from a TestInfo.
mock_path
    A fake path-like object.
poscars_path
    Path to the data directory containing POSCAR files.
run_phaseshift
    An Rparams, a Slab, and the results of generating phaseshifts
    with them.
tensorleed_path
    Path to the top-level tree with tensor-LEED source code.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    'Alexander M. Imre (@amimre)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2023-02-28'
__license__ = 'GPLv3+'

from pathlib import Path

import pytest
import pytest_cases

from viperleed.calc import psgen
from viperleed.calc.files import displacements
from viperleed.calc.files import vibrocc
from viperleed.calc.files.tenserleed import get_tensorleed_path
from viperleed.calc.lib.context import execute_in_dir

from ..helpers import POSCAR_PATH
from ..helpers import exclude_tags
from . import poscar_slabs
from .tags import CaseTag


@pytest.fixture(scope='session')
def poscars_path():
    """Return the Path to the directory containing POSCAR files."""
    return POSCAR_PATH


@pytest.fixture(scope='session', name='tensorleed_path')
def fixture_tensorleed_path():
    """Return the Path to the top-level tree with tensor-LEED source code."""
    return get_tensorleed_path()


@pytest.fixture(name='make_poscar', scope='session')
def factory_make_poscar():
    """Return a POSCAR from a TestInfo."""
    def _make(info):
        return poscar_slabs.CasePOSCARSlabs().case_poscar(info)
    return _make


@pytest.fixture(name='ag100')
def fixture_ag100():
    """Return a Ag(100) slab, an Rparams, and a TestInfo."""
    return poscar_slabs.CasePOSCARSlabs().case_poscar_ag100()


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


@pytest_cases.fixture
def mock_path(mocker):
    """Return a fake pathlib.Path."""
    return mocker.MagicMock(spec=Path)


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
    rpars : Rparams
        The PARAMETERS used during the PHASESHIFTS calculation.
    slab : Slab
        The Slab for which PHASESHIFTS were calculated.
    tmp_path : Path
        Path to the directory in which the calculation was executed.
    firstline : str
        The first line of the PHASESHIFTS file, that contains the
        coefficients for the real part of the inner potential.
    phaseshift : list
        The PHASESHIFTS that were generated.
    """
    slab, rpars, *_ = args
    rpars.paths.tensorleed = tensorleed_path
    tmp_path = tmp_path_factory.mktemp(basename='phaseshifts', numbered=True)
    rpars.initTheoEnergies()
    executable = 'eeasisss'

    # run eeasisss in the temporary directory
    with execute_in_dir(tmp_path):
        results = psgen.runPhaseshiftGen_old(slab, rpars,
                                             psgensource=executable)
        yield (rpars, slab, tmp_path, *results)
