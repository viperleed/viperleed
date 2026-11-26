"""Tests for module canvases of viperleed.gui.widgets."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-17'
__license__ = 'GPLv3+'

from pytest_cases import fixture
from pytest_cases import parametrize

from ..known_bugs import MPL_TO_QT_COORDINATES


class TestMPLFigureCanvas:
    """Tests for the MPLFigureCanvas class."""

    @fixture(name='canvas')
    def fixture_canvas(self, import_with_qt, show):
        """Return a shown MPLFigureCanvas."""
        module = import_with_qt('viperleed.gui.widgets.canvases')
        canvas = module.MPLFigureCanvas()
        canvas.resize(200, 150)
        show(canvas)
        return canvas

    @MPL_TO_QT_COORDINATES
    @parametrize(x_y=((0, 0), (30, 40), (199, 149)))
    def test_mpl_to_qt_coordinates(self, x_y, canvas, qtc):
        """Check the inverse of matplotlib's mouseEventCoords."""
        position = qtc.QPoint(*x_y)
        x_mpl, y_mpl = canvas.mouseEventCoords(position)
        assert canvas.mpl_to_qt_coordinates(x_mpl, y_mpl) == position
