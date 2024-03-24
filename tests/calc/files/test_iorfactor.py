"""Tests for module viperleed.calc.files.iorfactor."""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-12-11'
__license__ = 'GPLv3+'

import pytest
from pytest_cases import fixture

from viperleed.calc.classes.rparams import Rparams, TheoEnergies
from viperleed.calc.files import beams
from viperleed.calc.files import iorfactor


@fixture(name='ag100_expbeams')
def fixture_ag100_expbeams(data_path):
    """Return a list of experimental beam energies for Ag(100)."""
    _expbeams_path = data_path / 'Ag(100)' / 'initialization' / 'EXPBEAMS.csv'
    return beams.readOUTBEAMS(str(_expbeams_path))


class TestCheckTheoBeamsEnergies:
    """Collection of tests for checking the theoretical beam energies."""

    def test_valid_range(self, ag100_expbeams):
        """Ensure matching theory energies is correctly identified."""
        rpars = Rparams()
        rpars.THEO_ENERGIES = TheoEnergies(300, 500, 0.5)
        iorfactor.check_theobeams_energies(rpars, ag100_expbeams)

    def test_invalid_range(self, ag100_expbeams):
        """Emulate an Rfactor on a energy range larger than in Refcalc."""
        rpars = Rparams()
        rpars.THEO_ENERGIES = TheoEnergies(300, 800, 0.5)
        with pytest.raises(ValueError):
            iorfactor.check_theobeams_energies(rpars, ag100_expbeams)
