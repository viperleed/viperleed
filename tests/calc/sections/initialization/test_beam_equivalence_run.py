"""Tests for BEAM_EQUIVALENCE handling in a full initialization."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-25'
__license__ = 'GPLv3+'

import pytest
from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.calc.classes.beam import Beam
from viperleed.calc.classes.beam_correspondence import BeamCorrespondence
from viperleed.calc.files.tenserleed import CURRENT_TL_VERSION
from viperleed.calc.lib.context import execute_in_dir
from viperleed.calc.sections.initialization import BeamEquivalenceError
from viperleed.calc.sections.initialization import initialization

pytestmark = pytest.mark.tensorleed


@fixture(name='init_with_file')
def factory_init_with_file(ag100, make_section_tempdir, tensorleed_path):
    """Return a function running initialization with a given file."""
    slab, rpars, *_ = ag100
    rpars.paths.tensorleed = tensorleed_path
    rpars.TL_VERSION = CURRENT_TL_VERSION
    tmp = make_section_tempdir('Ag(100)', 'init_beam_equivalence')

    def _run(assignments, run):
        rpars.RUN = list(run)
        rpars.expbeams = [Beam(hk) for hk in assignments]
        for beam in rpars.expbeams:
            beam.intens = {50.0: 1.0, 100.0: 1.0}
        rpars.fileLoaded['EXPBEAMS'] = True
        rpars.fileLoaded['BEAM_EQUIVALENCE'] = True
        rpars.beam_correspondence = BeamCorrespondence(assignments)
        with execute_in_dir(tmp):
            initialization(slab, rpars)
        return tmp, rpars
    return _run


class TestInitializationWithBeamEquivalence:
    """Tests for the initialization with a loaded BEAM_EQUIVALENCE."""

    _weighted = {(1, 0): {(1, 0): 1, (0, 1): 2}}

    @parametrize('run', ([1], [3], [31], [4]))
    def test_weights_error_at_initialization(self, init_with_file, run):
        """Check the error at initialization if an R factor runs."""
        with pytest.raises(BeamEquivalenceError):
            init_with_file(self._weighted, run)
