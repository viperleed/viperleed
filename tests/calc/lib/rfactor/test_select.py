"""Tests for selecting R-factor implementations by name."""

import pytest

from viperleed.calc.lib import rfactor


@pytest.mark.parametrize(
    ('name', 'expected'),
    [
        ('pendry', rfactor.r_pendry),
        ('r_p', rfactor.r_pendry),
        ('r_pendry', rfactor.r_pendry),
        ('rp', rfactor.r_pendry),
        ('pendry r-factor', rfactor.r_pendry),
        ('p', rfactor.r_pendry),
        ('r1', rfactor.r_1),
        ('r_1', rfactor.r_1),
        ('r1 factor', rfactor.r_1),
        ('r2', rfactor.r_2),
        ('r_2', rfactor.r_2),
        ('r2 factor', rfactor.r_2),
        ('s', rfactor.r_s),
        ('rs', rfactor.r_s),
        ('r_s', rfactor.r_s),
        ('r_smooth', rfactor.r_s),
        ('smooth', rfactor.r_s),
        ('schmid', rfactor.r_s),
        ('zj', rfactor.R_zj),
        ('zj factor', rfactor.R_zj),
        ('zannazi', rfactor.R_zj),
        ('zannazi jona', rfactor.R_zj),
        ('zannazi-jona', rfactor.R_zj),
    ],
)
def test_select_rfactor_synonyms(name, expected):
    assert rfactor.select_rfactor(name) is expected


def test_select_rfactor_strips_and_lowercases():
    assert rfactor.select_rfactor('  RP  ') is rfactor.r_pendry


def test_select_rfactor_unknown():
    with pytest.raises(ValueError, match="Unknown R-factor name"):
        rfactor.select_rfactor("definitely-not-an-rfactor")
