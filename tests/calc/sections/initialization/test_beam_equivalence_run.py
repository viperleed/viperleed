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
from viperleed.calc.files.beam_equivalence import read_beam_equivalence
from viperleed.calc.files.tenserleed import CURRENT_TL_VERSION
from viperleed.calc.lib.context import execute_in_dir
from viperleed.calc.sections.calc_section import BEAM_EQUIVALENCE_USED_NAME
from viperleed.calc.sections.initialization import BeamEquivalenceError
from viperleed.calc.sections.initialization import initialization

_WEIGHTED = {(1, 0): {(1, 0): 1, (0, 1): 2}}

pytestmark = pytest.mark.tensorleed


@fixture(name='init_with_file')
def factory_init_with_file(ag100, make_section_tempdir, tensorleed_path):
    """Return a function running initialization with a given file."""
    slab, rpars, *_ = ag100
    rpars.paths.tensorleed = tensorleed_path
    rpars.TL_VERSION = CURRENT_TL_VERSION
    tmp = make_section_tempdir('Ag(100)', 'init_beam_equivalence')

    def _run(assignments, run, expbeams=None, ivbeams=None):
        rpars.RUN = list(run)
        if expbeams is None:
            expbeams = assignments
        rpars.expbeams = [Beam(hk) for hk in expbeams]
        for beam in rpars.expbeams:
            beam.intens = {50.0: 1.0, 100.0: 1.0}
        rpars.fileLoaded['EXPBEAMS'] = True
        if ivbeams is not None:
            rpars.ivbeams = [Beam(hk) for hk in ivbeams]
            rpars.ivbeams_sorted = False
            rpars.fileLoaded['IVBEAMS'] = True
        rpars.fileLoaded['BEAM_EQUIVALENCE'] = True
        rpars.beam_correspondence = BeamCorrespondence(assignments)
        with execute_in_dir(tmp):
            initialization(slab, rpars)
        return tmp, rpars
    return _run


class TestInitializationWithBeamEquivalence:
    """Tests for the initialization with a loaded BEAM_EQUIVALENCE."""

    @parametrize('run', ([1], [3], [31], [4]))
    def test_weights_error_at_initialization(self, init_with_file, run):
        """Check the error at initialization if an R factor runs."""
        with pytest.raises(BeamEquivalenceError):
            init_with_file(_WEIGHTED, run)

    def test_weights_without_r_factor(self, init_with_file):
        """Check that weights are fine if no R factor runs."""
        tmp, rpars = init_with_file(_WEIGHTED, [2])
        used = read_beam_equivalence(tmp / BEAM_EQUIVALENCE_USED_NAME)
        assert used == rpars.beam_correspondence

    def test_ignored_line_not_checked(self, init_with_file):
        """Check that a line of a beam not in EXPBEAMS is not used."""
        assignments = {(1, 0): {(1, 0): 1},
                       (1, 1): {(1, 1): 1},
                       (5, 5): {(2, 2): 1, (3, 1): 2}}
        tmp, _ = init_with_file(assignments, [1],
                                expbeams=[(1, 0), (1, 1)])
        used = read_beam_equivalence(tmp / BEAM_EQUIVALENCE_USED_NAME)
        assert used == BeamCorrespondence({(1, 0): {(1, 0): 1},
                                           (1, 1): {(1, 1): 1}})

    def test_ignored_line_not_in_ivbeams(self, init_with_file):
        """Check no error for an ignored line not in IVBEAMS."""
        assignments = {(1, 0): {(1, 0): 1}, (9, 9): {(2, 0): 1}}
        _, rpars = init_with_file(assignments, [2],
                                  expbeams=[(1, 0)],
                                  ivbeams=[(1, 0)])
        assert rpars.beam_correspondence == BeamCorrespondence(
            {(1, 0): {(1, 0): 1}}
            )

    @parametrize(run=([0], [1]))
    @parametrize(ivbeams=(None, [(1, 0), (1, 1), (2, 0)]))
    def test_all_lines_ignored(self, init_with_file, run, ivbeams):
        """Check the error if no left-hand side is in EXPBEAMS."""
        assignments = {(9, 9): {(1, 0): 1}, (8, 8): {(1, 1): 1}}
        with pytest.raises(BeamEquivalenceError,
                           match='No left-hand side of BEAM_EQUIVALENCE'):
            init_with_file(assignments, run,
                           expbeams=[(1, 0), (1, 1)],
                           ivbeams=ivbeams)
