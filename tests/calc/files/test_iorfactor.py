"""Tests for module viperleed.calc.files.iorfactor."""

__authors__ = ('Alexandra Mia Imre (@alexmiame)',)
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2023-12-11'
__license__ = 'GPLv3+'

import numpy as np
import pytest
from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.calc.classes.rparams.rparams import Rparams
from viperleed.calc.classes.rparams.special.r_factor_type import RFactorType
from viperleed.calc.classes.rparams.special.energy_range import TheoEnergies
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


class TestWriteRfactorPdf:
    """Tests for writing R-factor plots from TensErLEED column files."""

    @fixture(name='columns')
    def fixture_columns(self):
        """Return fake theoretical and experimental column data."""
        energies = np.linspace(50, 100, 21)
        theo = np.column_stack([energies, 1 + np.sin(energies / 5)])
        exp = np.column_stack([energies, 1 + np.cos(energies / 5)])
        return [theo], [exp]

    @parametrize(r_factor_type=('pendry', 'smooth'))
    def test_writes_analysis(self, r_factor_type, columns, mocker):
        """Check that the analysis PDF is written for an RFactorType."""
        mocker.patch.object(iorfactor,
                            'read_rfactor_columns',
                            return_value=columns)
        mocker.patch.object(iorfactor, 'plot_iv')
        write_analysis = mocker.patch.object(
            iorfactor,
            '_write_rfactor_analysis_pdf_from_splines',
            )
        which_r = RFactorType(r_factor_type)
        iorfactor.writeRfactorPdf([('(1|0)', 0.1)],
                                  analysisFile='analysis.pdf',
                                  which_r=which_r)
        write_analysis.assert_called_once()
        assert write_analysis.call_args.kwargs['which_r'] == which_r
