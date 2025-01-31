"""Tests for module viperleed.calc.psgen."""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-07-28'
__license__ = 'GPLv3+'

import numpy as np
import pytest
from pytest_cases import parametrize_with_cases

from viperleed.calc.psgen import adjust_phaseshifts
from viperleed.calc.psgen import runPhaseshiftGen_old

from ..helpers import execute_in_dir
from .poscar_slabs import CasePOSCARSlabs as POSCARSlabs

MIN_ENTRIES_IN_FIRST_LINE = 4


class TestPhaseshiftsGen:
    """Tests for the successful outcome of a PHASESHIFTS calculation."""

    def test_phaseshifts_not_empty(self, run_phaseshift):
        """Assert that the generated PHASESHIFTS contain items."""
        *_, phaseshift = run_phaseshift
        assert phaseshift

    def test_phaseshifts_firstline_not_empty(self, run_phaseshift):
        """Check that the first line contains characters."""
        *_, firstline, _ = run_phaseshift
        assert firstline

    def test_phaseshifts_firstline_len(self, run_phaseshift, subtests):
        """Check that the first line has at least four float coefficients."""
        *_, firstline, _ = run_phaseshift
        _, *potential_param = firstline.split()
        n_floats = 0                                                            # TODO: this calculation is repeated in at least two other places
        for coeff in potential_param:
            try:
                float(coeff)
            except (TypeError, ValueError):
                continue
            n_floats += 1
        with subtests.test('First line: enough items'):
            assert len(potential_param) >= MIN_ENTRIES_IN_FIRST_LINE
        with subtests.test('First line: enough floating-point items'):
            assert n_floats >= MIN_ENTRIES_IN_FIRST_LINE

    def test_phaseshift_log_exists(self, run_phaseshift):
        """Ensure a log file was written to disk."""
        *_, work_path, _, _ = run_phaseshift
        assert any(work_path.glob('phaseshift*.log'))

    @parametrize_with_cases('args', cases=POSCARSlabs.case_poscar_lsmo_001_rt2)
    def test_phaseshift_input_mixed_sites(self, args, tensorleed_path,
                                          tmp_path_factory):
        """Test that phaseshift generation works with mixed sites."""
        tmp_path = tmp_path_factory.mktemp(basename='phaseshifts').resolve()
        slab, rpars, _ = args
        rpars.paths.source = tensorleed_path
        rpars.THEO_ENERGIES = rpars.THEO_ENERGIES.from_value((50,100,5))
        with execute_in_dir(tmp_path):
            firstline, _ = runPhaseshiftGen_old(slab, rpars)
        assert firstline


class TestAdjustPhaseshifts:
    """"Tests for the wrap_phaseshifts function."""

    linear_increasing = np.linspace(0.1, 2*np.pi, 50)
    linear_decreasing = np.linspace(2*np.pi, 0.1, 50)
    sinusoidal = np.sin(np.linspace(0, 2*np.pi, 50))
    mock_phaseshifts = {
    'linear_increasing' : linear_increasing,
    'linear_increasing_offset' : linear_increasing + 10*np.pi,
    'linear_decreasing' : linear_decreasing,
    'linear_through_0_from_pos' : linear_decreasing - 0.5,
    'linear_through_0_from_neg' : linear_increasing - 0.5,
    'linear_through_pi_from_pos' : linear_increasing + np.pi - 0.2,
    'linear_through_pi_from_neg' : np.pi + 0.3 - linear_increasing,
    'linear_with_jumps_from_pos' : (linear_increasing + np.pi - 0.2)  % np.pi,
    'linear_with_jumps_from_neg' : (np.pi + 0.3 - linear_increasing)  % np.pi,
    'sinusoidal' : sinusoidal,
    'sinusoidal_offset' : sinusoidal + 10*np.pi,
    'sinusoidal_start_neg' : sinusoidal - 0.5,
    'sinusoidal_through_pi' : sinusoidal + np.pi - 0.2,
    }

    @pytest.mark.parametrize('phaseshifts', mock_phaseshifts.values(),
                             ids=mock_phaseshifts.keys())
    def test_adjust_phaseshifts_is_constrained(self, phaseshifts):
        """Assert that the adjusted phaseshifts don't diverge."""
        adjusted_ps = adjust_phaseshifts(phaseshifts)
        assert np.max(np.abs(adjusted_ps)) < np.max(np.abs(phaseshifts))*3

    @pytest.mark.parametrize('phaseshifts', mock_phaseshifts.values(),
                                ids=mock_phaseshifts.keys())
    def test_adjust_phaseshifts_is_continuous(self, phaseshifts):
        """Assert that the adjusted phaseshifts are continuous."""
        adjusted_ps = adjust_phaseshifts(phaseshifts)
        assert np.all(np.abs(np.diff(adjusted_ps)) < 0.2)

    @pytest.mark.parametrize('phaseshifts', mock_phaseshifts.values(),
                             ids=mock_phaseshifts.keys())
    def test_adjust_phaseshifts_first_value_is_wrapped(self, phaseshifts):
        """Assert that the first value is wrapped."""
        adjusted_ps = adjust_phaseshifts(phaseshifts)
        assert np.abs(adjusted_ps[0]) < np.pi/2
