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
    assert any(rpars.workdir.glob('phaseshift*.log'))


def test_write_phaseshifts(run_phaseshift):
    """Ensure a PHASESHIFTS file is successfully written to disk."""
    rpars, _, firstline, phaseshift = run_phaseshift
    phaseshifts.writePHASESHIFTS(firstline, phaseshift,
                                 file_path=rpars.workdir/'PHASESHIFTS')
    assert any(rpars.workdir.glob('PHASESHIFTS'))


def test_phaseshifts_not_empty(run_phaseshift):
    """Ensure run_phaseshift has generated some PHASESHIFTS."""
    *_, phaseshift = run_phaseshift
    assert len(phaseshift)


def test_check_consistency_element_order_no_false_positives(run_phaseshift):
    """Test for __check_consitency_element_order().
    Checks that non of the example systems produce false positives."""
    rpars, slab, _, phaseshift = run_phaseshift
    inconsistencies = phaseshifts.__check_consitency_element_order(
        rpars, slab, phaseshift)
    assert not inconsistencies
