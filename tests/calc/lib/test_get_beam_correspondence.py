"""Characterisation tests for leedbase.getBeamCorrespondence."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-24'
__license__ = 'GPLv3+'

import logging

from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.calc import symmetry
from viperleed.calc.classes.beam import Beam
from viperleed.calc.classes.rparams.rparams import Rparams
from viperleed.calc.files import parameters
from viperleed.calc.files import poscar
from viperleed.calc.files.beams import readOUTBEAMS
from viperleed.calc.lib.leedbase import getBeamCorrespondence


_MODULE = 'viperleed.calc.lib.leedbase'

# A fake symmetry-equivalence list, as returned by getSymEqBeams:
# ((h, k), index of the group of equivalent beams)
_SYMEQ = (
    ((0.0, 0.0), 0),
    ((1.0, 0.0), 1), ((0.0, 1.0), 1), ((-1.0, 0.0), 1), ((0.0, -1.0), 1),
    ((1.0, 1.0), 2), ((-1.0, 1.0), 2), ((1.0, -1.0), 2), ((-1.0, -1.0), 2),
    ((2.0, 0.0), 3), ((0.0, 2.0), 3), ((-2.0, 0.0), 3), ((0.0, -2.0), 3),
    ((0.5, 0.0), 4), ((0.0, 0.5), 4),
    )


def _labels(beams):
    """Return the labels of `beams`."""
    return [b.label for b in beams]


def _make_beams(*indices):
    """Return a list of Beam objects with the given (h, k) indices."""
    return [Beam(hk) for hk in indices]


@fixture(name='synthetic')
def fixture_synthetic(mocker):
    """Return a fake slab and an Rparams, with getSymEqBeams mocked."""
    mocker.patch(f'{_MODULE}.getSymEqBeams', return_value=list(_SYMEQ))
    rpars = Rparams()
    slab = mocker.MagicMock()

    def _make(expbeams, ivbeams):
        rpars.expbeams = _make_beams(*expbeams)
        rpars.ivbeams = _make_beams(*ivbeams)
        return slab, rpars
    return _make


class TestSynthetic:
    """Tests for getBeamCorrespondence with a fake symmetry list."""

    def test_matched_directly(self, synthetic):
        """Check experimental beams matching a calculated beam."""
        slab, rpars = synthetic(expbeams=((1, 0), (1, 1)),
                                ivbeams=((1, 1), (1, 0)))
        assert getBeamCorrespondence(slab, rpars) == [1, 0]
        assert len(rpars.expbeams) == 2
        assert rpars.halt == 0

    def test_matched_via_symmetry(self, synthetic):
        """Check an experimental beam matched via an equivalent."""
        slab, rpars = synthetic(expbeams=((1, 1),),
                                ivbeams=((2, 0), (-1, 1)))
        assert getBeamCorrespondence(slab, rpars) == [-1, 0]
        assert _labels(rpars.expbeams) == ['( 1| 1)']
        assert rpars.halt == 0

    def test_no_counterpart(self, synthetic, caplog):
        """Check removal of an experimental beam without counterpart."""
        slab, rpars = synthetic(expbeams=((3, 0), (1, 0), (2, 1)),
                                ivbeams=((1, 0),))
        removed = rpars.expbeams[0], rpars.expbeams[2]
        with caplog.at_level(logging.WARNING, logger=_MODULE):
            assert getBeamCorrespondence(slab, rpars) == [0]
        assert _labels(rpars.expbeams) == ['( 1| 0)']
        assert not any(b in rpars.expbeams for b in removed)
        assert rpars.halt == 1
        assert [r.getMessage() for r in caplog.records] == [
            ('Experimental beam ( 3| 0) does not have a theoretical '
             'counterpart. Consider adding it to IVBEAMS.'),
            ('Experimental beam ( 2| 1) does not have a theoretical '
             'counterpart. Consider adding it to IVBEAMS.'),
            'Beam ( 3| 0) was removed from set of experimental beams!',
            'Beam ( 2| 1) was removed from set of experimental beams!',
            ]

    def test_calculated_beam_via_symmetry(self, synthetic):
        """Check a calculated beam assigned via an equivalent beam."""
        slab, rpars = synthetic(expbeams=((1, 1), (1, 0)),
                                ivbeams=((0, 1), (1, 1), (-1, -1)))
        assert getBeamCorrespondence(slab, rpars) == [1, 0, 0]

    def test_calculated_beam_unassigned(self, synthetic, caplog):
        """Check calculated beams without experimental beam."""
        slab, rpars = synthetic(expbeams=((1, 0),),
                                ivbeams=((1, 0), (2, 0), (1, 1), (3, 1)))
        with caplog.at_level(logging.DEBUG, logger=_MODULE):
            assert getBeamCorrespondence(slab, rpars) == [0, -1, -1, -1]
        assert rpars.halt == 0
        assert ('No experimental beams found for calculated beams: '
                '( 2| 0), ( 1| 1), ( 3| 1)') in caplog.messages
        assert caplog.messages == [
            ('No experimental beams found for calculated beams: '
             '( 2| 0), ( 1| 1), ( 3| 1)'),
            ]

    def test_two_calculated_on_one_experimental(self, synthetic):
        """Check calculated beams averaged for one experimental."""
        slab, rpars = synthetic(
            expbeams=((1, 0), (2, 0)),
            ivbeams=((0, 1), (1, 0), (0, -2), (-1, 0), (2, 0)),
            )
        assert getBeamCorrespondence(slab, rpars) == [0, 0, 1, 0, 1]

    def test_two_experimental_on_one_calculated(self, synthetic):
        """Check that the last experimental beam claiming one wins."""
        slab, rpars = synthetic(expbeams=((1, 0), (0, 1)),
                                ivbeams=((1, 0),))
        assert getBeamCorrespondence(slab, rpars) == [1]
        assert len(rpars.expbeams) == 2

    def test_fractional_and_tolerance(self, synthetic):
        """Check fractional beams and the matching tolerance."""
        slab, rpars = synthetic(expbeams=((0.5, 0), (1, 0)),
                                ivbeams=((0, 0.5), (1 + 2e-6, -3e-6)))
        assert getBeamCorrespondence(slab, rpars) == [0, 1]

    def test_beyond_tolerance(self, synthetic):
        """Check that beams 2e-5 apart do not match."""
        slab, rpars = synthetic(expbeams=((1, 0), (1, 1)),
                                ivbeams=((1 + 2e-5, 0), (1, 1)))
        assert getBeamCorrespondence(slab, rpars) == [-1, 0]
        assert _labels(rpars.expbeams) == ['( 1| 1)']
        assert rpars.halt == 1

    def test_no_symmetry(self, synthetic, mocker):
        """Check matching without symmetry-equivalent beams."""
        slab, rpars = synthetic(expbeams=((1, 0), (1, 1)),
                                ivbeams=((0, 1), (1, 0), (-1, -1)))
        mocker.patch(f'{_MODULE}.getSymEqBeams', return_value=[])
        assert getBeamCorrespondence(slab, rpars) == [-1, 0, -1]
        assert _labels(rpars.expbeams) == ['( 1| 0)']
        assert rpars.halt == 1


_AG100_IVBEAMS = (
    (1, 0), (0, 1), (1, 1), (-1, 1), (2, 0), (0, -2), (1, 2), (3, 0),
    (-3, 0), (5, 0),
    )


@fixture(name='ag100_init')
def fixture_ag100_init(data_path):
    """Return a Ag(100) slab and Rparams after the symmetry search."""
    inputs = data_path / 'Ag(100)' / 'initialization'
    slab = poscar.read(inputs / 'POSCAR')
    rpars = parameters.read(inputs / 'PARAMETERS')
    parameters.interpret(rpars, slab=slab, silent=True)
    slab.full_update(rpars)
    symmetry.findSymmetry(slab, rpars)
    symmetry.enforceSymmetry(slab, rpars)
    slab.make_bulk_slab(rpars)
    bulk = slab.bulkslab
    symmetry.findSymmetry(bulk, rpars, bulk=True, output=False)
    bulk.revert_unit_cell()
    symmetry.findBulkSymmetry(bulk, rpars)
    rpars.expbeams = readOUTBEAMS(str(inputs / 'EXPBEAMS.csv'))
    rpars.ivbeams = _make_beams(*_AG100_IVBEAMS)
    return slab, rpars


_AG100_EXPECTED = {
    # AVERAGE_BEAMS: (correspondence, remaining experimental beams)
    None: ([0, 0, 1, 1, 2, 2, 3, 4, 4, -1],
           ['( 1| 0)', '( 1| 1)', '( 2| 0)', '( 2| 1)', '( 3| 0)']),
    False: ([0, -1, 1, -1, 2, -1, -1, 3, -1, -1],
            ['( 1| 0)', '( 1| 1)', '( 2| 0)', '( 3| 0)']),
    }


@parametrize(average_beams=_AG100_EXPECTED)
def test_ag100(ag100_init, average_beams):
    """Check Ag(100) with and without beam averaging."""
    slab, rpars = ag100_init
    rpars.AVERAGE_BEAMS = average_beams
    expected, remaining = _AG100_EXPECTED[average_beams]
    assert getBeamCorrespondence(slab, rpars) == expected
    assert _labels(rpars.expbeams) == remaining
    assert rpars.halt == 1
