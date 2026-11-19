"""Tests for plot_iv_from_csv in viperleed.calc.from_ase."""

__authors__ = ('Florian Kraushofer (@fkraushofer)',)
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-19'
__license__ = 'GPLv3+'

import inspect
from pathlib import Path
from unittest.mock import patch

import pytest

from viperleed.calc import from_ase

BEAMS_FILE = (Path(__file__).resolve().parents[1] / '_test_data'
              / 'beam-data' / 'EXPBEAMS-fe2o3-012-1x1.csv')


@pytest.fixture(name='mock_plot')
def fixture_mock_plot():
    """Replace plot_iv with a mock, to avoid producing figures."""
    with patch.object(from_ase, 'plot_iv') as mock:
        yield mock


def test_which_beams_default_is_immutable():
    """Check that which_beams has no mutable default."""
    default = inspect.signature(from_ase.plot_iv_from_csv).parameters[
        'which_beams'
    ].default
    assert default is None or isinstance(default, (str, tuple))


@pytest.mark.parametrize('beam_file', (BEAMS_FILE, str(BEAMS_FILE)))
def test_single_path_is_one_file(beam_file, mock_plot):
    """Check that a lone path is a single file, not a sequence of files."""
    from_ase.plot_iv_from_csv(beam_file)
    all_beam_data = mock_plot.call_args[0][0]
    assert len(all_beam_data) == 1


def test_single_content_string_is_one_file(mock_plot):
    """Check that a lone string with CSV contents is a single file."""
    from_ase.plot_iv_from_csv(BEAMS_FILE.read_text(),
                              beam_file_is_content=True)
    all_beam_data = mock_plot.call_args[0][0]
    assert len(all_beam_data) == 1


@pytest.mark.parametrize('is_content', ((True,), [True, False, True]))
def test_mismatched_lengths(is_content, mock_plot):
    """Check complaints if beam_file_is_content does not fit beam_file."""
    files = (BEAMS_FILE.read_text(), BEAMS_FILE)
    with pytest.raises(ValueError):
        from_ase.plot_iv_from_csv(files, beam_file_is_content=is_content)
    mock_plot.assert_not_called()


@pytest.mark.parametrize('which_beams', (None, [], 'all'))
def test_all_beams(which_beams, mock_plot):
    """Check that None, an empty list, and 'all' select all beams."""
    from_ase.plot_iv_from_csv(BEAMS_FILE, which_beams=which_beams)
    all_beam_data = mock_plot.call_args[0][0]
    assert len(all_beam_data[0]) > 2


def _plotted_labels(mock_plot):
    """Return the labels of the beams passed to the mocked plot_iv."""
    return [beam.label for beam in mock_plot.call_args[0][0][0]]


def test_which_beams_indices(mock_plot):
    """Check that integers select beams by their position."""
    from_ase.plot_iv_from_csv(BEAMS_FILE)
    all_labels = _plotted_labels(mock_plot)
    from_ase.plot_iv_from_csv(BEAMS_FILE, which_beams=[2, 0])
    assert _plotted_labels(mock_plot) == [all_labels[2], all_labels[0]]


def test_which_beams_labels(mock_plot):
    """Check that strings select beams by their label, in file order."""
    from_ase.plot_iv_from_csv(BEAMS_FILE)
    all_labels = _plotted_labels(mock_plot)
    wanted = [all_labels[3], all_labels[1]]
    from_ase.plot_iv_from_csv(BEAMS_FILE, which_beams=wanted)
    assert _plotted_labels(mock_plot) == [all_labels[1], all_labels[3]]


@pytest.mark.parametrize('which_beams', ([0, 'a'], [1.5]))
def test_which_beams_invalid(which_beams, mock_plot):
    """Check complaints for an invalid which_beams."""
    with pytest.raises(ValueError, match='which_beams'):
        from_ase.plot_iv_from_csv(BEAMS_FILE, which_beams=which_beams)
    mock_plot.assert_not_called()
