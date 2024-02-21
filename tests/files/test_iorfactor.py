"""Tests for module viperleed.tleedmlib.files.iorfactor.

Created on 2023-12-11

@author: Alexander M. Imre (@amimre)
"""

from pathlib import Path
import sys

import pytest
from pytest_cases import fixture

VPR_PATH = str(Path(__file__).resolve().parents[3])
if VPR_PATH not in sys.path:
    sys.path.append(VPR_PATH)

# pylint: disable=wrong-import-position
# Will be fixed in installable version
from viperleed.tleedmlib.classes.rparams import Rparams, TheoEnergies
from viperleed.tleedmlib.files.beams import readOUTBEAMS

from viperleed.tleedmlib.files import iorfactor
# pylint: enable=wrong-import-position


@fixture(name='ag100_expbeams')
def fixture_ag100_expbeams(data_path):
    """Return a list of experimental beam energies for Ag(100)."""
    _expbeams_path = data_path / 'Ag(100)' / 'initialization' / 'EXPBEAMS.csv'
    return readOUTBEAMS(str(_expbeams_path))


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
