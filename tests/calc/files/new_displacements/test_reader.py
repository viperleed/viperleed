"""Tests for module viperleed.calc.files.new_displacements.reader."""

__authors__ = ('Alexandra Mia Imre (@alexmiame)',)
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2025-10-01'
__license__ = 'GPLv3+'

import pytest
from pytest_cases import parametrize

from viperleed.calc.files.new_displacements.errors import (
    DisplacementsSyntaxError,
)
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


@parametrize(label=('my_sym_delta_run', 'SYM_DELTA'))
def test_displacements_reader_search_label_with_sym_delta(tmp_path, label):
    """Check that search-block labels are not mistaken for SYM_DELTA."""
    path = tmp_path / 'DISPLACEMENTS_label'
    path.write_text(f'== SEARCH {label}\n')
    _compare_lines(path, [(SearchHeaderLine, {'label': label})])


def test_displacements_reader_sym_delta_deprecated(tmp_path):
    """Check that the deprecated SYM_DELTA mode header is rejected."""
    path = tmp_path / 'DISPLACEMENTS_sym_delta'
    path.write_text('== SEARCH a\n= SYM_DELTA\n')
    with pytest.raises(DisplacementsSyntaxError, match='SYMMETRY_FIX'):
        _compare_lines(path, [])


_DROPPED_HEADERS = (
    pytest.param(
        '== S a', '== SEARCH',
        id='search',
    ),
    pytest.param(
        '== SEAR a', '== SEARCH',
        id='search_long',
    ),
    pytest.param(
        '= G', '= GEO_DELTA',
        id='geo',
    ),
    pytest.param(
        '= GEO', '= GEO_DELTA',
        id='geo_long',
    ),
    pytest.param(
        '= V', '= VIB_DELTA',
        id='vib',
    ),
    pytest.param(
        '= O', '= OCC_DELTA',
        id='occ',
    ),
    pytest.param(
        '= C', '= CONSTRAIN',
        id='constrain',
    ),
    pytest.param(
        '= OFF', '= OFFSETS',
        id='offsets',
    ),
    pytest.param(
        '= REGEX on', r'"\*" wildcards',
        id='regex_on',
    ),
    pytest.param(
        '= regex off', r'"\*" wildcards',
        id='regex_off',
    ),
    pytest.param(
        '= R on', r'"\*" wildcards',
        id='regex_r',
    ),
    pytest.param(
        '= REG off', r'"\*" wildcards',
        id='regex_reg',
    ),
    pytest.param(
        '= GEO_D', '= GEO_DELTA',
        id='geo_delta_prefix',
    ),
    pytest.param(
        '== sea a', '== SEARCH',
        id='search_prefix',
    ),
    pytest.param(
        '= constr', '= CONSTRAIN',
        id='constrain_prefix',
    ),
    pytest.param(
        '= Geometry', r'unknown header .*= GEO_DELTA, = VIB',
        id='unknown_word',
    ),
    pytest.param(
        '= GEO_DELTA x', 'unexpected text after header',
        id='text_after_header',
    ),
    pytest.param(
        '==SEARCH a', r'Write "== SEARCH a"',
        id='search_no_space',
    ),
    pytest.param(
        '== SEARCHabc',
        r'space is needed between SEARCH and the label.*"== SEARCH abc"',
        id='search_glued_label',
    ),
)


@parametrize('header, msg', _DROPPED_HEADERS)
def test_displacements_reader_dropped_header(header, msg, tmp_path):
    """Check errors for abbreviated headers and = REGEX."""
    path = tmp_path / 'DISPLACEMENTS'
    path.write_text(f'{header}\nO 1 = 0.1\n')
    raises = pytest.raises(DisplacementsSyntaxError, match=msg)
    with raises as err, DisplacementsReader(path) as reader:
        list(reader)
    assert f'"{header}"' in err.value.message


def test_displacements_reader_regex_in_mode_block(tmp_path):
    """Check that an abbreviated = REGEX is reported inside a block."""
    path = tmp_path / 'DISPLACEMENTS'
    path.write_text('== SEARCH a\n= GEO_DELTA\nO 1 z = 0.1 0.2\n= R on\n')
    raises = pytest.raises(DisplacementsSyntaxError, match=r'"\*" wildcards')
    with raises, DisplacementsReader(path) as reader:
        list(reader)


_DOUBLE_TAGS = (
    pytest.param(
        '<loop><loop>\n== SEARCH a\n', ('start', 'start'),
        id='start',
    ),
    pytest.param(
        '</loop></loop>\n', ('end', 'end'),
        id='end',
    ),
    pytest.param(
        '<loop> <LOOP>\n== SEARCH a\n', ('start', 'start'),
        id='spaced',
    ),
)


@parametrize('text, kinds', _DOUBLE_TAGS)
def test_displacements_reader_several_tags_on_line(text, kinds, tmp_path):
    """Check that several loop tags on one line are read as tags."""
    path = tmp_path / 'DISPLACEMENTS'
    path.write_text(text)
    with DisplacementsReader(path) as reader:
        read = list(reader)
    assert [m.kind for m in read[:2]] == list(kinds)
    assert all(isinstance(m, LoopMarkerLine) for m in read[:2])


def test_displacements_reader_abbreviation_after_loop(tmp_path):
    """Check that an abbreviated header after <loop> is named."""
    path = tmp_path / 'DISPLACEMENTS'
    path.write_text('<loop>\n== S b\n')
    raises = pytest.raises(DisplacementsSyntaxError, match='abbreviated')
    with raises, DisplacementsReader(path) as reader:
        list(reader)


@parametrize('header', ('== SEARCH a', '== search a', '= geo_delta'))
def test_displacements_reader_full_headers_ok(header, tmp_path):
    """Check that full headers, in any case, are still accepted."""
    path = tmp_path / 'DISPLACEMENTS'
    path.write_text(f'{header}\n')
    with DisplacementsReader(path) as reader:
        assert len(list(reader)) == 1
