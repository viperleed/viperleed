"""Tests for module viperleed.tleedmlib.psgen.

Created on 2023-07-28

@author: Alexander M. Imre (@amimre)
@author: Michele Riva (@michele-riva)
"""


class TestPhaseshiftsGen:
    """Tests for the successful outcome of a PHASESHIFTS calculation."""

    def test_phaseshifts_not_empty(self, run_phaseshift):
        """Assert that the generated PHASESHIFTS contain items."""
        _, _, _, phaseshift = run_phaseshift
        assert phaseshift

    def test_phaseshifts_firstline_not_empty(self, run_phaseshift):
        """Check that the first line contains characters."""
        _, _, firstline, _ = run_phaseshift
        assert firstline

    def test_phaseshifts_firstline_len(self, run_phaseshift, subtests):
        """Check that the first line has at least four float coefficients."""
        _, _, firstline, _ = run_phaseshift
        _, *potential_param = firstline.split()
        n_floats = 0                                                            # TODO: this calculation is repeated in at least two other places
        for coeff in potential_param:
            try:
                float(coeff)
            except (TypeError, ValueError):
                continue
            n_floats += 1
        with subtests.test('First line: enough items'):
            assert len(potential_param) >= 4
        with subtests.test('First line: enough floating-point items'):
            assert n_floats >= 4

    def test_phaseshift_log_exists(self, run_phaseshift):
        """Ensure a log file was written to disk."""
        param, _, _, _ = run_phaseshift
        assert any(param.workdir.glob('phaseshift*.log'))
