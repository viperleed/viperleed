"""Tests for checkEXPBEAMS of viperleed.calc.files.beams."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-25'
__license__ = 'GPLv3+'

from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.calc.classes.beam import Beam
from viperleed.calc.classes.rparams.rparams import Rparams
from viperleed.calc.files.beams import checkEXPBEAMS

_MODULE = 'viperleed.calc.files.beams'


@fixture(name='rpars')
def fixture_rpars():
    """Return an Rparams with no experimental beams."""
    rpars = Rparams()
    rpars.HALTING = 3
    return rpars


class TestCheckExpbeamsSameHk:
    """Tests for same-(h, k) experimental beams."""

    @parametrize(average_beams=(False, True))
    def test_duplicate_discarded(self, average_beams, rpars, mocker, caplog):
        """Check same-(h, k) beams are reported once and discarded."""
        rpars.AVERAGE_BEAMS = average_beams
        first, second, other = Beam((1, 0)), Beam((1, 0)), Beam((0, 1))
        rpars.expbeams = [first, other, second]
        symeq = [((1.0, 0.0), 0)] if average_beams else []
        mocker.patch(f'{_MODULE}.leedbase.getSymEqBeams', return_value=symeq)
        with caplog.at_level('WARNING'):
            checkEXPBEAMS(mocker.MagicMock(), rpars)
        assert rpars.expbeams == [first, other]
        assert rpars.expbeams[0] is first
        assert 'equivalent to experimental beam' in caplog.text
        assert len(caplog.records) == 1
        assert 'will be discarded' in caplog.text
        assert rpars.halt >= 2

    def test_distinct_beams_kept(self, rpars, mocker, caplog):
        """Check that beams with different (h, k) survive."""
        beams = [Beam((1, 0)), Beam((0, 1))]
        rpars.AVERAGE_BEAMS = False
        rpars.expbeams = list(beams)
        mocker.patch(f'{_MODULE}.leedbase.getSymEqBeams', return_value=[])
        with caplog.at_level('WARNING'):
            checkEXPBEAMS(mocker.MagicMock(), rpars)
        assert rpars.expbeams == beams
        assert not caplog.records
