"""Tests for module matplotlib_utils of viperleed.calc.lib."""


__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-07-10'
__license__ = 'GPLv+3'

import logging
from unittest.mock import Mock

import pytest
from pytest_cases import fixture, parametrize

# Import the module
from viperleed.calc.lib import matplotlib_utils
from viperleed.calc.lib.log_utils import at_level
from viperleed.calc.lib.matplotlib_utils import _LOGGER
from viperleed.calc.lib.matplotlib_utils import HasNoMatplotlibError
from viperleed.calc.lib.matplotlib_utils import close_figures
from viperleed.calc.lib.matplotlib_utils import log_without_matplotlib
from viperleed.calc.lib.matplotlib_utils import prepare_matplotlib_for_calc
from viperleed.calc.lib.matplotlib_utils import raise_without_matplotlib
from viperleed.calc.lib.matplotlib_utils import skip_without_matplotlib
from viperleed.calc.lib.matplotlib_utils import use_calc_style


MOCK_RETURN = 'Function called'


@fixture(name='mock_logger')
def fixture_mock_logger():
    """Return a fake logging.Logger."""
    return Mock(logging.Logger)


@fixture(name='mock_pyplot')
def fixture_mock_pyplot():
    """Return a fake matplotlib.pyplot."""
    return Mock()


@fixture(name='decorated', scope='session')
def fixture_factory_decorated():
    """Return a dummy function decorated by a given decorator."""
    def mock_func():
        return MOCK_RETURN

    def _make(decorator):
        return decorator(mock_func)

    return _make


class TestHasMatplotlibError:
    """Collection of tests for HasNoMatplotlibError."""

    def test_custom_msg(self):
        """Check that a custom message is raised correctly."""
        error_message = 'Custom error message'
        error = HasNoMatplotlibError(error_message)
        assert str(error) == error_message

    def test_default_msg(self):
        """Check that the default message is raised if none is given."""
        error = HasNoMatplotlibError()
        # pylint: disable-next=protected-access
        assert str(error) == error._default_msg


class TestCannotPlot:
    """Collections of tests for missing matplotlib."""

    @fixture(autouse=True)
    def missing_matplotlib(self, monkeypatch):
        """Patch matplotlib_utils to fake the absence of matplotlib."""
        monkeypatch.setattr(matplotlib_utils, 'CAN_PLOT', False)

    @parametrize(decorator=(skip_without_matplotlib,
                            log_without_matplotlib(_LOGGER)))
    def test_returns_none(self, decorated, decorator):
        """Check that decorating returns None without matplotlib."""
        func = decorated(decorator)
        assert func() is None

    @parametrize(level=('debug', 'info', 'warning', 'error'))
    def test_log(self, level, caplog, decorated):
        """Check emits log message."""
        expected_log_msg = 'Necessary modules for plotting not found.'
        func = decorated(log_without_matplotlib(_LOGGER, level=level))
        with at_level(_LOGGER, logging.DEBUG):
            func()
        assert expected_log_msg in caplog.text

    def test_raises(self, decorated):
        """Check that calling a raise-decorated function raises indeed."""
        func = decorated(raise_without_matplotlib)
        with pytest.raises(HasNoMatplotlibError):
            func()


class TestCanPlot:
    """Collection of tests for when matplotlib is available."""

    @fixture(autouse=True)
    def with_matplotlib(self, monkeypatch):
        """Make sure matplotlib is found, even if it's not there."""
        monkeypatch.setattr(matplotlib_utils, 'CAN_PLOT', True)

    @parametrize(decorator=(skip_without_matplotlib,
                            log_without_matplotlib(_LOGGER),
                            raise_without_matplotlib))
    def test_returns_expected(self, decorated, decorator):
        """Check that decorating returns None without matplotlib."""
        func = decorated(decorator)
        assert func() == MOCK_RETURN

    @parametrize(level=('debug', 'info', 'warning', 'error'))
    def test_does_not_log(self, level, decorated, mock_logger):
        """Check that no logging messages are emitted."""
        func = decorated(log_without_matplotlib(mock_logger, level=level))
        func()
        mock_logger.debug.assert_not_called()

    def test_close_figures(self, mock_pyplot):
        """Check that calling close_figures calls pyplot.close."""
        figures = Mock(), Mock()
        close_figures(mock_pyplot, *figures)
        for fig in figures:
            mock_pyplot.close.assert_any_call(fig)

    def test_close_figures_exception_logs(self, mock_pyplot, caplog):
        """Check that exceptions during close_figures are logged."""
        _exc_text = 'Close failed'
        mock_pyplot.close.side_effect = Exception(_exc_text)
        close_figures(mock_pyplot, Mock())
        assert _exc_text in caplog.text

    def test_prepare_matplotlib_for_calc(self, monkeypatch):
        """Check that all the expected functions are called."""
        mock_matplotlib = Mock()
        monkeypatch.setattr(matplotlib_utils, 'matplotlib', mock_matplotlib)
        monkeypatch.setattr(matplotlib_utils,
                            'mpl_style',
                            mock_matplotlib.style)

        prepare_matplotlib_for_calc()
        mock_matplotlib.rcParams.update.assert_called_once_with(
            {'figure.max_open_warning': 0}
            )
        mock_matplotlib.use.assert_called_once_with('Agg')
        mock_matplotlib.style.use.assert_called_once()

    @parametrize(mpl_version=('3.6', '3.7'))
    def test_use_calc_style(self, mpl_version, monkeypatch):
        """Check that the matplotlib style is updated via use_calc_style."""
        mock_mpl_style = Mock()
        monkeypatch.setattr(matplotlib_utils, 'mpl_version', mpl_version)
        monkeypatch.setattr(matplotlib_utils, 'mpl_style', mock_mpl_style)
        use_calc_style()
        mock_mpl_style.use.assert_called_once()
