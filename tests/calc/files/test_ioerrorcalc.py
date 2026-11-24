"""Tests for module viperleed.calc.files.ioerrorcalc."""

__authors__ = ('Florian Kraushofer (@fkraushofer)',)
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-24'
__license__ = 'GPLv3+'

from types import SimpleNamespace

import numpy as np
import pytest
from pytest_cases import fixture

from viperleed.calc.files import ioerrorcalc
from viperleed.calc.lib.matplotlib_utils import CAN_PLOT

pytestmark = pytest.mark.skipif(not CAN_PLOT, reason='Needs matplotlib')

_LIN_DISP = [0.0, 0.1, 0.2, 0.3, 0.4]
_RFACS = [0.5, 0.3, 0.2, 0.3, 0.5]
_ESTIMATES = (0.1, 0.1)


class _StubError:
    """A minimal stand-in for an R_Error."""

    var_r = 0.15
    get_error_estimates = _ESTIMATES

    def __init__(self, mode):
        """Initialize a stub error curve for `mode`."""
        self.mode = mode
        self.rfacs = _RFACS
        self.lin_disp = _LIN_DISP
        self.disp_label = 'z'
        self.atoms = [SimpleNamespace(num=1, site='A')]

    @property
    def get_r_min(self):
        """Return the minimum R factor."""
        return min(self.rfacs)


@fixture(name='occ_figs')
def fixture_occ_figs():
    """Return the figures for a stub occupation error."""
    return ioerrorcalc.make_errors_figs([_StubError('occ')])


class TestMakeErrorsFigsOcc:
    """Tests for plotting occupation errors."""

    def test_combined_plot(self, occ_figs):
        """Check x data and markers of the combined plot in percent."""
        ax = occ_figs[0].axes[0]
        curve = ax.lines[-1]
        x_data = curve.get_xdata()
        marks = curve.get_markevery()
        assert len(marks) == len(_LIN_DISP)
        assert x_data[marks] == pytest.approx(np.array(_LIN_DISP)*100,
                                              abs=0.05)
        assert len(x_data) > len(marks)
        assert min(x_data) == pytest.approx(0)
        assert max(x_data) == pytest.approx(40, abs=0.05)
        assert ax.get_xlim() == pytest.approx((-2, 42))

    def test_single_plot(self, occ_figs):
        """Check axis limits and error bounds of single plots."""
        ax = occ_figs[1].axes[0]
        assert ax.get_xlim() == pytest.approx((-2, 42))
        curve = ax.lines[-1]
        assert max(curve.get_xdata()) == pytest.approx(40, abs=0.05)
        bounds = sorted(
            seg[0][0] for coll in ax.collections
            for seg in coll.get_segments()
            )
        assert bounds[0] == pytest.approx(10)
        assert bounds[-1] == pytest.approx(30)
