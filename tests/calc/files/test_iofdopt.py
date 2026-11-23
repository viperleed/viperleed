"""Tests for module iofdopt of viperleed.calc.files."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-23'
__license__ = 'GPLv3+'

import numpy as np
import pytest
from pytest_cases import fixture

from viperleed.calc.files import iofdopt
from viperleed.calc.lib.matplotlib_utils import CAN_PLOT

pytestmark = pytest.mark.skipif(not CAN_PLOT, reason='Needs matplotlib')


class TestWriteFdOptBeamsPdf:
    """Tests for the write_fd_opt_beams_pdf function."""

    @fixture(name='mock_plot_iv')
    def fixture_mock_plot_iv(self, mocker):
        """Replace reading and plotting of I(V) curves with mocks."""
        mocker.patch.object(iofdopt,
                            'read_rfactor_columns',
                            return_value=('theo', 'exp'))
        return mocker.patch.object(iofdopt, 'plot_iv')

    @fixture(name='rpars')
    def fixture_rpars(self, mocker):
        """Return fake run parameters with one experimental beam."""
        beam = mocker.Mock()
        beam.getLabel.return_value = ('(1|0)', 5)
        return mocker.Mock(PLOT_IV={'overbar': False}, expbeams=[beam])

    def test_colors(self, rpars, mock_plot_iv):
        """Check that curves are colored with viridis, then black."""
        mpl = pytest.importorskip('matplotlib')
        points = np.array(((0.9, 0.3), (1.0, 0.2), (1.1, 0.4)))
        tmpdirs = [f'dir_{i}' for i in range(len(points))]
        iofdopt.write_fd_opt_beams_pdf(rpars, points, 'a', tmpdirs, [0.2])
        _, kwargs = mock_plot_iv.call_args
        colors = kwargs['formatting']['colors']
        viridis = mpl.colormaps['viridis']
        expect = [*viridis(np.linspace(0, 1, len(points))), (0, 0, 0, 1)]
        assert np.allclose(colors, expect)
