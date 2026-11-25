"""Tests for BEAM_EQUIVALENCE checks in the initialization section."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-25'
__license__ = 'GPLv3+'

import logging

import pytest
from pytest_cases import fixture

from viperleed.calc.classes.beam_correspondence import BeamCorrespondence
from viperleed.calc.classes.rparams.rparams import Rparams
from viperleed.calc.sections.initialization import BeamEquivalenceError
from viperleed.calc.sections.initialization import (
    _check_beam_equivalence_in_beamlist,
    )

_MODULE = 'viperleed.calc.sections.initialization'
_BEAMLIST = ['  0  0\n', '  1  0\n', '  0  1\n', '  1  1\n']


@fixture(name='make_rpars')
def factory_make_rpars(make_beams):
    """Return a factory of Rparams with a BEAM_EQUIVALENCE loaded."""
    def _make(assignments, expbeams=(), ivbeams=()):
        rpars = Rparams()
        rpars.beamlist = list(_BEAMLIST)
        rpars.fileLoaded['BEAM_EQUIVALENCE'] = True
        rpars.beam_correspondence = BeamCorrespondence(assignments)
        rpars.expbeams = make_beams(*expbeams)
        rpars.fileLoaded['EXPBEAMS'] = bool(expbeams)
        rpars.ivbeams = make_beams(*ivbeams)
        return rpars
    return _make


class TestBeamlistCheck:
    """Tests for _check_beam_equivalence_in_beamlist."""

    def test_all_in_beamlist(self, make_rpars, caplog):
        """Check that no message is emitted for a compatible file."""
        rpars = make_rpars({(1, 0): {(1, 0): 1, (0, 1): 1}})
        with caplog.at_level(logging.DEBUG):
            _check_beam_equivalence_in_beamlist(rpars)
        assert not caplog.records

    def test_beams_not_in_beamlist(self, make_rpars, caplog):
        """Check the error listing all beams missing in BEAMLIST."""
        rpars = make_rpars({(1, 0): {(1, 0): 1, (2, 0): 1},
                            (1, 1): {(3, 0): 1}})
        with pytest.raises(BeamEquivalenceError):
            _check_beam_equivalence_in_beamlist(rpars)
        errors = [r for r in caplog.records if r.levelno == logging.ERROR]
        assert len(errors) == 1
        assert 'not in BEAMLIST: ( 2| 0), ( 3| 0)' in errors[0].message

    def test_average_beams_warning(self, make_rpars, caplog):
        """Check that a given AVERAGE_BEAMS is reported as ignored."""
        rpars = make_rpars({(1, 0): {(1, 0): 1}})
        rpars.readParams['AVERAGE_BEAMS'].append('T')
        with caplog.at_level(logging.WARNING):
            _check_beam_equivalence_in_beamlist(rpars)
        assert len(caplog.records) == 1
        assert 'AVERAGE_BEAMS is ignored' in caplog.text

