"""Tests for selecting R-factor implementations by name."""

import pytest

from viperleed.calc.lib import rfactor


@pytest.mark.parametrize(
    ("name", "expected"),
    [
        ("pendry", rfactor.R_pendry),
        ("r_p", rfactor.R_pendry),
        ("r_pendry", rfactor.R_pendry),
        ("rp", rfactor.R_pendry),
        ("pendry r-factor", rfactor.R_pendry),
        ("p", rfactor.R_pendry),
        ("r1", rfactor.R_1),
        ("r_1", rfactor.R_1),
        ("r1 factor", rfactor.R_1),
        ("r2", rfactor.R_2),
        ("r_2", rfactor.R_2),
        ("r2 factor", rfactor.R_2),
        ("s", rfactor.R_s),
        ("rs", rfactor.R_s),
        ("r_s", rfactor.R_s),
        ("r_smooth", rfactor.R_s),
        ("smooth", rfactor.R_s),
        ("schmid", rfactor.R_s),
        ("zj", rfactor.R_zj),
        ("zj factor", rfactor.R_zj),
        ("zannazi", rfactor.R_zj),
        ("zannazi jona", rfactor.R_zj),
        ("zannazi-jona", rfactor.R_zj),
    ],
)
def test_select_rfactor_synonyms(name, expected):
    assert rfactor.select_rfactor(name) is expected


def test_select_rfactor_strips_and_lowercases():
    assert rfactor.select_rfactor("  RP  ") is rfactor.R_pendry


def test_select_rfactor_unknown():
    with pytest.raises(ValueError, match="Unknown R-factor name"):
        rfactor.select_rfactor("definitely-not-an-rfactor")
