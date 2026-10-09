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
from viperleed.calc.classes.rparams.special.energy_range import TheoEnergies
from viperleed.calc.classes.rparams.special.r_factor_type import RFactorType
from viperleed.calc.files import beams
from viperleed.calc.files import iorfactor
from viperleed.calc.lib.version import Version


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
        assert write_analysis.call_args.kwargs['which_r'] is which_r



_WEXPEL_WEIGHTS = {  # r_factor_type, tl_version, (wr_line, warns)
    'pendry': ('pendry', '2.0.0', (' WR=      0.,0.,1.,', False)),
    'r2': ('r2', '2.0.0', (' WR=      1.,0.,0.,', False)),
    'zj': ('zj', '2.0.0', (' WR=      0.,1.,0.,', False)),
    'smooth, old TensErLEED': ('smooth',
                               '2.0.0',
                               (' WR=      0.,0.,1.,', True)),
    'smooth': ('smooth', '2.1.0', (' WR=      0.,0.,0.,', False)),
    }


class TestWriteWEXPELRFactorType:
    """Tests for the R-factor weights written to WEXPEL."""

    @fixture(name='write_wexpel')
    def factory_write_wexpel(self, tmp_path, mocker):
        """Return a function that writes WEXPEL and returns it."""
        energies = mocker.MagicMock(min=50.0, max=100.0)
        iv_shift = mocker.MagicMock(start=-3.0, stop=3.0)
        mocker.patch.object(iorfactor,
                            'prepare_rfactor_energy_ranges',
                            return_value=(None, energies, iv_shift, 0.5))
        mocker.patch.object(iorfactor.leedbase,
                            'getBeamCorrespondence',
                            return_value=[])
        mocker.patch.object(iorfactor, 'writeAUXEXPBEAMS', return_value='')
        wexpel = tmp_path / 'WEXPEL'

        def _write(rpars):
            iorfactor.writeWEXPEL(None, rpars, [], filename=wexpel)
            return wexpel.read_text().splitlines()
        return _write

    @parametrize('r_factor_type,tl_version,expect',
                 _WEXPEL_WEIGHTS.values(),
                 ids=_WEXPEL_WEIGHTS)
    def test_weights(self, r_factor_type, tl_version, expect,
                     write_wexpel, caplog):
        """Check the WR= line written for each R_FACTOR_TYPE."""
        rpars = Rparams()
        rpars.R_FACTOR_TYPE = RFactorType(r_factor_type)
        rpars.TL_VERSION = Version(tl_version)
        wr_line, warns = expect
        assert wr_line in write_wexpel(rpars)
        assert ('Falling back on Pendry' in caplog.text) == warns
