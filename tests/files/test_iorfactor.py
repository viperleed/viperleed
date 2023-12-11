"""Tests for module viperleed.tleedmlib.files.iorfactor.

Created on 2023-12-11

@author: Alexander M. Imre (@amimre)
"""

from pathlib import Path
import sys

import pytest
from pytest_cases import parametrize_with_cases, fixture

VPR_PATH = str(Path(__file__).resolve().parents[3])
if VPR_PATH not in sys.path:
    sys.path.append(VPR_PATH)

# pylint: disable=wrong-import-position
# Will be fixed in installable version
from viperleed.tleedmlib.classes.rparams import Rparams, TheoEnergies
from viperleed.tleedmlib.files.beams import readOUTBEAMS

from viperleed.tleedmlib.files import iorfactor


# pylint: enable=wrong-import-position

@pytest.fixture
def ag100_expbeams():
    """Return a list of experimental beam energies for Ag(100)."""
    return readOUTBEAMS(str(Path(VPR_PATH) / 'viperleed' / 'tests' / '_test_data' / 'Ag(100)' / 'initialization' / 'EXPBEAMS.csv'))


class TestCheckTheoBeamsEnergies:
    """Collection of tests for checking the theoretical beam energies."""
    def test_check_theobeams_energies_valid_range(self, ag100_expbeams):
        theobeams = ag100_expbeams
        rpars = Rparams()
        rpars.THEO_ENERGIES = TheoEnergies(300, 500, 0.5)
        iorfactor.check_theobeams_energies(rpars, theobeams)

    def test_check_theobeams_energies_invalid_range(self, ag100_expbeams):
        """Test that emulates a Rfactor calculation on a larger energy range than the Refcalc."""
        theobeams = ag100_expbeams
        rpars = Rparams()
        rpars.THEO_ENERGIES = TheoEnergies(300, 800, 0.5)
        with pytest.raises(ValueError):
            iorfactor.check_theobeams_energies(rpars, theobeams)
