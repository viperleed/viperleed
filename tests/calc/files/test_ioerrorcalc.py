"""Tests for module ioerrorcalc of viperleed.calc.files."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2026-10-08'
__license__ = 'GPLv3+'

from pytest_cases import parametrize

from viperleed.calc.classes.rparams.special.r_factor_type import RFactorType
from viperleed.calc.files import ioerrorcalc


@parametrize(r_factor_type=('pendry', 'r2', 'zj', 'smooth'))
def test_make_errors_figs_any_r_factor_type(r_factor_type):
    """Check that error figures can be made for any R-factor type."""
    figs = ioerrorcalc.make_errors_figs([], RFactorType(r_factor_type))
    assert not figs
