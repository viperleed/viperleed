"""Tests for module ioerrorcalc of viperleed.calc.files."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2026-10-08'
__license__ = 'GPLv3+'

import numpy as np
import pytest
from pytest_cases import parametrize

from viperleed.calc.classes.rparams.special.r_factor_type import RFactorType
from viperleed.calc.files import ioerrorcalc
from viperleed.calc.lib.matplotlib_utils import close_figures


@parametrize('r_factor_type,ylabel', (
    ('pendry', 'Pendry R-factor'),
    ('r2', 'R2 R-factor'),
    ('zj', 'Zanazzi-Jona R-factor'),
    ('smooth', 'Smooth R-factor'),
    ))
def test_make_errors_figs_label(r_factor_type, ylabel, mocker):
    """Check that error figures are labelled with the R-factor type."""
    pyplot = pytest.importorskip('matplotlib.pyplot')
    atom = mocker.MagicMock(site='Ag', num=1)
    error = mocker.MagicMock(mode='geo',
                             lin_disp=np.array([-0.1, 0., 0.1]),
                             rfacs=[0.3, 0.2, 0.3],
                             get_r_min=0.2,
                             var_r=0.05,
                             get_error_estimates=(None, None),
                             atoms=[atom],
                             disp_label='z')
    figs = ioerrorcalc.make_errors_figs([error], RFactorType(r_factor_type))
    try:
        assert figs[0].axes[0].get_ylabel() == ylabel
    finally:
        close_figures(pyplot, *figs)
