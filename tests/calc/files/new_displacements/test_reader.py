"""Tests for module viperleed.calc.files.new_displacements.reader."""

__authors__ = ('Alexander M. Imre (@amimre)',)
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2024-10-04'
__license__ = 'GPLv3+'


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
