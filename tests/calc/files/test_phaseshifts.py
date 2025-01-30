"""Tests for module viperleed.calc.files.phaseshifts."""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-07-28'
__license__ = 'GPLv3+'

from viperleed.calc.files import phaseshifts

from ..test_psgen import MIN_ENTRIES_IN_FIRST_LINE


def test_phaseshifts_firstline_not_empty(run_phaseshift):
    """Ensure a first line of the PHASESHIFTS file has been created."""
    *_, firstline, _ = run_phaseshift
    assert firstline


def test_phaseshifts_firstline_len(run_phaseshift):
    """Ensure the first line of a PHASESHIFTS file has enough items."""
    *_, firstline, _ = run_phaseshift
    potential_param = firstline.split()
    assert len(potential_param) >= MIN_ENTRIES_IN_FIRST_LINE


def test_phaseshift_log_exists(run_phaseshift):
    """Ensure a log file is successfully written to disk."""
    rpars, *_ = run_phaseshift
    assert any(rpars.paths.work.glob('phaseshift*.log'))


def test_write_phaseshifts(run_phaseshift):
    """Ensure a PHASESHIFTS file is successfully written to disk."""
    rpars, _, firstline, phaseshift = run_phaseshift
    phaseshifts.writePHASESHIFTS(firstline, phaseshift,
                                 file_path=rpars.paths.work/'PHASESHIFTS')
    assert any(rpars.paths.work.glob('PHASESHIFTS'))


def test_phaseshifts_not_empty(run_phaseshift):
    """Ensure run_phaseshift has generated some PHASESHIFTS."""
    *_, phaseshift = run_phaseshift
    assert len(phaseshift)


def test_check_consistency_element_order_no_false_positives(run_phaseshift):
    """Check that no example system produces false positives for elements."""
    rpars, slab, _, phaseshift = run_phaseshift
    _check = phaseshifts.__check_consistency_element_order
    inconsistencies = _check(rpars, slab, phaseshift)
    assert not inconsistencies
