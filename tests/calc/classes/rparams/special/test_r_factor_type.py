"""Tests for viperleed.calc.classes.rparams.special.r_factor_type."""

__authors__ = ('Alexander M. Imre (@amimre)',)
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-11-11'
__license__ = 'GPLv3+'


from viperleed.calc.classes.rparams.special.r_factor_type import RFactorType


def test_synonyms():
    assert int(RFactorType('r 2')) == RFactorType.R2
    assert int(RFactorType('Zanazzi Jona')) == RFactorType.ZJ
    assert int(RFactorType('rp')) == RFactorType.PENDRY
    assert int(RFactorType('smooth pendry')) == RFactorType.SMOOTH
