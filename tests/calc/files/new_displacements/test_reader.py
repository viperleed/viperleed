"""Tests for module viperleed.calc.files.new_displacements.reader."""

__authors__ = ('Alexandra Mia Imre (@alexmiame)',)
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__license__ = 'GPLv3+'


from viperleed.calc.files.new_displacements.lines import (
    GeoDeltaLine,
    LoopMarkerLine,
    ModeHeaderLine,
    SearchHeaderLine,
)
from viperleed.calc.files.new_displacements.reader import DisplacementsReader


def _compare_lines(path, expected_lines):
    with DisplacementsReader(path) as reader:
        parsed_lines = list(reader)
    assert len(parsed_lines) == len(expected_lines)
    for parsed, expected in zip(parsed_lines, expected_lines):
        exp_cls, exp_attrs = expected
        # check correct type
        assert isinstance(parsed, exp_cls)
        # check correct attributes
        for attr, value in exp_attrs.items():
            assert getattr(parsed, attr) == value


def test_displacements_reader(mock_displacements_path_and_lines):
    path, expected_lines = mock_displacements_path_and_lines
    _compare_lines(path, expected_lines)


def test_displacements_reader_simple(
    mock_displacements_cu_111_realistic_path_and_lines,
):
    path, expected_lines = mock_displacements_cu_111_realistic_path_and_lines
    _compare_lines(path, expected_lines)


def test_displacements_reader_loop_markers(tmp_path):
    """Check that loop markers (which contain no '=') are read."""
    path = tmp_path / 'DISPLACEMENTS_loop'
    path.write_text(
        '<loop>\n'
        '== SEARCH a\n'
        '= GEO_DELTA\n'
        'Fe 1 z = 0.0 0.1\n'
        '</loop>\n'
    )
    _compare_lines(
        path,
        [
            (LoopMarkerLine, {'kind': 'start'}),
            (SearchHeaderLine, {'label': 'a'}),
            (ModeHeaderLine, {'section': 'GEO_DELTA'}),
            (GeoDeltaLine, {}),
            (LoopMarkerLine, {'kind': 'end'}),
        ],
    )
