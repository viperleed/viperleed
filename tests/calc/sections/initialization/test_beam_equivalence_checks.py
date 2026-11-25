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
from quicktions import Fraction

from viperleed.calc.classes.beam_correspondence import BeamCorrespondence
from viperleed.calc.classes.rparams.rparams import Rparams
from viperleed.calc.sections.initialization import BeamEquivalenceError
from viperleed.calc.sections.initialization import _apply_beam_equivalence
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


class TestApplyIvbeams:
    """Tests for the IVBEAMS part of _apply_beam_equivalence."""

    def test_missing_beams(self, make_rpars, caplog):
        """Check the error listing beams of the file not in IVBEAMS."""
        rpars = make_rpars({(1, 0): {(1, 0): 1, (0, 1): 1, (1, 1): 1}},
                           ivbeams=[(1, 0)])
        with pytest.raises(BeamEquivalenceError):
            _apply_beam_equivalence(rpars)
        errors = [r for r in caplog.records if r.levelno == logging.ERROR]
        assert len(errors) == 1
        assert 'missing in IVBEAMS: ( 0| 1), ( 1| 1)' in errors[0].message

    def test_calculated_beams_on_no_right_hand_side(self, make_rpars, caplog):
        """Check the info log for unassigned calculated beams."""
        rpars = make_rpars({(1, 0): {(1, 0): 1}},
                           ivbeams=[(1, 0), (0, 1), (1, 1)])
        with caplog.at_level(logging.INFO):
            _apply_beam_equivalence(rpars)
        records = [r for r in caplog.records if 'right-hand' in r.message]
        assert len(records) == 1
        assert records[0].levelno == logging.INFO
        assert '( 0| 1), ( 1| 1)' in records[0].message

    def test_all_assigned_no_log(self, make_rpars, caplog):
        """Check that nothing is logged if everything is assigned."""
        rpars = make_rpars({(1, 0): {(1, 0): 1, (0, 1): 1}},
                           ivbeams=[(1, 0), (0, 1)])
        with caplog.at_level(logging.DEBUG):
            _apply_beam_equivalence(rpars)
        assert not caplog.records


class TestApplyExpbeams:
    """Tests for the EXPBEAMS part of _apply_beam_equivalence."""

    def test_left_hand_side_not_in_expbeams(self, make_rpars, caplog):
        """Check that a line with an unknown beam is ignored."""
        rpars = make_rpars({(1, 0): {(1, 0): 1}, (1, 1): {(0, 1): 1}},
                           expbeams=[(1, 0)], ivbeams=[(1, 0), (0, 1)])
        with caplog.at_level(logging.WARNING):
            _apply_beam_equivalence(rpars)
        assert rpars.beam_correspondence == BeamCorrespondence(
            {(1, 0): {(1, 0): 1}}
            )
        assert len(rpars.expbeams) == 1
        warnings = [r.message for r in caplog.records]
        assert len(warnings) == 1
        assert '( 1| 1), which is not in EXPBEAMS' in warnings[0]

    def test_experimental_beam_without_line(self, make_rpars, caplog):
        """Check removal of an experimental beam not in the file."""
        rpars = make_rpars({(1, 0): {(1, 0): 1}},
                           expbeams=[(1, 0), (2, 0), (1, 1)],
                           ivbeams=[(1, 0)])
        with caplog.at_level(logging.WARNING):
            _apply_beam_equivalence(rpars)
        assert [b.hkfrac for b in rpars.expbeams] == [(1, 0)]
        assert rpars.halt == 1
        assert [r.message for r in caplog.records] == [
            ('Experimental beam ( 2| 0) is not in BEAM_EQUIVALENCE, '
             'and will not be used in the R factor.'),
            'Beam ( 2| 0) was removed from set of experimental beams!',
            ('Experimental beam ( 1| 1) is not in BEAM_EQUIVALENCE, '
             'and will not be used in the R factor.'),
            'Beam ( 1| 1) was removed from set of experimental beams!',
            ]

    def test_filtered_correspondence_stored(self, make_rpars):
        """Check that rp keeps the correspondence of the kept lines."""
        half = Fraction(1, 2)
        rpars = make_rpars({(1, 0): {(1, 0): 1, (0, 1): 3},
                            (half, 0): {(half, 0): 1},
                            (3, 0): {(3, 0): 1}},
                           expbeams=[(1, 0), (half, 0)],
                           ivbeams=[(1, 0), (0, 1), (half, 0), (3, 0)])
        _apply_beam_equivalence(rpars)
        assert rpars.beam_correspondence.assignments == {
            (1, 0): {(1, 0): 0.25, (0, 1): 0.75},
            (half, 0): {(half, 0): 1.0},
            }

    def test_no_expbeams_file(self, make_rpars, caplog):
        """Check that the correspondence is kept without EXPBEAMS."""
        rpars = make_rpars({(1, 0): {(1, 0): 1}}, ivbeams=[(1, 0)])
        before = rpars.beam_correspondence
        with caplog.at_level(logging.WARNING):
            _apply_beam_equivalence(rpars)
        assert rpars.beam_correspondence is before
        assert not caplog.records
