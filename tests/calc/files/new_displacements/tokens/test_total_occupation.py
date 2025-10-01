"""Tests for calc.files.new_displacements.tokens.total_occupation."""

__authors__ = ('Alexander M. Imre (@amimre)',)
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-07-18'
__license__ = 'GPLv3+'


import pytest

from viperleed.calc.files.new_displacements.tokens.total_occupation import (
    TotalOccupationToken,
    TotalOccupationTokenParserError,
)

# This token is essentially just a float parser, so not much to test here.


def test_total_occupation_token_valid():
    """Test that a valid total occupation token is parsed correctly."""
    token = TotalOccupationToken('0.5')
    assert token.total_occupation == 0.5


def test_total_occupation_token_invalid():
    """Test that an invalid total occupation token raises an error."""
    with pytest.raises(TotalOccupationTokenParserError):
        TotalOccupationToken('invalid')


def test_total_occupation_token_out_of_range():
    """Test that a total occupation token out of range raises an error."""
    with pytest.raises(TotalOccupationTokenParserError):
        TotalOccupationToken('-0.1')
    with pytest.raises(TotalOccupationTokenParserError):
        TotalOccupationToken('1.1')
