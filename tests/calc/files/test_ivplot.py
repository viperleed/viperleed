"""Tests for module viperleed.calc.files.ivplot."""

__authors__ = ('Florian Kraushofer (@fkraushofer)',)
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-18'
__license__ = 'GPLv3+'

import logging

import numpy as np
import pytest
from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.calc.files import ivplot
from viperleed.calc.lib.matplotlib_utils import CAN_PLOT

pytestmark = pytest.mark.skipif(not CAN_PLOT, reason='Needs matplotlib')


@fixture(name='iv_data')
def fixture_iv_data():
    """Return one dataset of two I(V) curves."""
    energies = np.linspace(50, 100, 11)
    curve = np.column_stack((energies, np.sin(energies)))
    return [[curve, curve]]


class TestPlotIVToFile:
    """Tests for writing I(V) curves to a PDF file with plot_iv."""

    def test_success(self, iv_data, tmp_path, caplog):
        """Check that a PDF file is written and None returned."""
        pdf_path = tmp_path / 'iv.pdf'
        with caplog.at_level(logging.ERROR):
            result = ivplot.plot_iv(iv_data, str(pdf_path))
        assert result is None
        assert pdf_path.is_file()
        assert not caplog.records

    def test_close_fails(self, iv_data, tmp_path, caplog, mocker):
        """Check that an error when closing the file is only logged."""
        _patch_close(mocker, RuntimeError('cannot close'))
        with caplog.at_level(logging.ERROR):
            result = ivplot.plot_iv(iv_data, str(tmp_path / 'iv.pdf'))
        assert result is None
        assert 'Error while writing' in caplog.text
        assert 'cannot close' in caplog.text

    @parametrize(exc=(KeyboardInterrupt, SystemExit))
    def test_interrupt_in_close_propagates(self, exc, iv_data,
                                           tmp_path, mocker):
        """Check that an interrupt while closing the file propagates."""
        _patch_close(mocker, exc)
        with pytest.raises(exc):
            ivplot.plot_iv(iv_data, str(tmp_path / 'iv.pdf'))

    @parametrize(exc=(KeyboardInterrupt, SystemExit))
    def test_interrupt_in_write_propagates(self, exc, iv_data,
                                           tmp_path, mocker):
        """Check that an interrupt in writing propagates after close."""
        mocker.patch.object(ivplot.PdfPages, 'savefig', side_effect=exc)
        close = mocker.spy(ivplot.PdfPages, 'close')
        with pytest.raises(exc):
            ivplot.plot_iv(iv_data, str(tmp_path / 'iv.pdf'))
        close.assert_called_once()


def _patch_close(mocker, exc):
    """Make PdfPages.close raise `exc` after closing the file."""
    real_close = ivplot.PdfPages.close

    def _close_then_raise(pdf):
        real_close(pdf)  # Do not leak the open file
        raise exc

    mocker.patch.object(ivplot.PdfPages, 'close', _close_then_raise)
