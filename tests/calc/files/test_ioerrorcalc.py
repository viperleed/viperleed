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
from viperleed.calc.lib.matplotlib_utils import close_figures

pytestmark = pytest.mark.skipif(not CAN_PLOT, reason='Needs matplotlib')

_LIN_DISP = [0.0, 0.1, 0.2, 0.3, 0.4]
_RFACS = [0.5, 0.3, 0.2, 0.3, 0.5]
_ESTIMATES = (0.1, 0.1)
_R_VAR_COLOR = 'slategray'  # Color of the R + var(R) lines


class _StubError:
    """A minimal stand-in for an R_Error."""

    error_estimates = _ESTIMATES

    def __init__(self, mode, r_offset=0, var_r=0.15):
        """Initialize a stub error curve for `mode`."""
        self.mode = mode
        self.rfacs = [r + r_offset for r in _RFACS]
        self.var_r = var_r
        self.lin_disp = _LIN_DISP
        self.disp_label = 'z'
        self.atoms = [SimpleNamespace(num=1, site='A')]

    @property
    def r_min(self):
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


def _r_plus_var_levels(axes):
    """Return the y levels of the R + var(R) lines in axes."""
    return [line.get_ydata()[0] for line in axes.lines
            if line.get_color() == _R_VAR_COLOR]


class TestMakeErrorsFigsRPlusVar:
    """Tests for the R + var(R) lines of two error curves."""

    @fixture(name='geo_figs')
    def fixture_geo_figs(self):
        """Yield the figures for two errors with different R minima."""
        errors = (_StubError('geo'),
                  _StubError('geo', r_offset=0.1, var_r=0.1))
        figs = ioerrorcalc.make_errors_figs(errors)
        yield figs
        close_figures(ioerrorcalc.plt, *figs)

    def test_combined_plot(self, geo_figs):
        """Check that the line uses the lowest R minimum and its var."""
        axes = geo_figs[0].axes[0]
        assert _r_plus_var_levels(axes) == pytest.approx([0.35])

    def test_single_plots(self, geo_figs):
        """Check that each plot uses its own R minimum and variance."""
        levels = [_r_plus_var_levels(a) for a in geo_figs[1].axes[:2]]
        assert levels == [pytest.approx([0.35]), pytest.approx([0.4])]


class TestCrossings:
    """Tests for the _crossings helper."""

    # pylint: disable-next=protected-access           # OK in tests
    crossings = staticmethod(ioerrorcalc._crossings)

    def test_crossings(self):
        """Check the points where the curve crosses the line."""
        r_values = [np.nan, 0.5, 0.2, 0.5]
        crossings = self.crossings(range(4), r_values, 0.3)
        assert crossings == [2, 3]

    def test_line_above_curve(self):
        """Check that no crossings are returned if R is always lower."""
        assert not self.crossings(range(3), [0.5, 0.2, 0.5], 1.5)
