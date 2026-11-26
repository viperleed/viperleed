"""Tests for module viperleed.calc.files.new_displacements.file."""

__authors__ = ('Alexandra Mia Imre (@alexmiame)',)
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2025-10-01'
__license__ = 'GPLv3+'


import re

import pytest
from pytest_cases import parametrize

from viperleed.calc.files.new_displacements.errors import (
    AlreadyReadError,
    DisplacementsSyntaxError,
    OffsetsNotAtBeginningError,
)
from viperleed.calc.files.new_displacements.file import (
    DisplacementsFile,
    _unparsable_line_message,
)
from viperleed.calc.files.new_displacements.lines import (
    ModeHeaderLine,
    OffsetsHeaderLine,
)
from viperleed.calc.files.new_displacements.segments import (
    ConstraintBlock,
    LoopBlock,
    SearchBlock,
)
from viperleed.calc.files.new_displacements.tokens.direction import (
    FractionalDirectionToken,
)

from ....helpers import (
    REPO_ROOT,
    TEST_DATA,
)

NESTED_LOOPS_PATH = (
    TEST_DATA / 'DISPLACEMENTS/Cu_111/DISPLACEMENTS_nested_loops'
)
LATERAL_1_PATH = (
    TEST_DATA / 'DISPLACEMENTS/Pt_111_10x10_Te/DISPLACEMENTS_lateral_1'
)


def test_read_from_file(displacements_file_path, subtests):
    """Test reading a file and checking its validity."""

    df = DisplacementsFile()
    assert df._has_been_read is False
    df.read(displacements_file_path)
    assert df._has_been_read is True

    with subtests.test('read() called again'), pytest.raises(AlreadyReadError):
        # check that read() cannot be called again
        df.read(displacements_file_path)
PTRH_PATH = TEST_DATA / 'stored_inputs/inputs/ptrh_100_1x3_o/DISPLACEMENTS'


def test_read_constraint_alias_file():
    """Check that a stored file using '= CONSTRAINT' can be read."""
    disp_file = DisplacementsFile()
    disp_file.read(PTRH_PATH)
    assert any(isinstance(b, ConstraintBlock) for b in disp_file.descendants)

ROUGH_PATH = (
    REPO_ROOT / 'doc/_static/example_systems/Cu(111)-Te'
    / 'DISPLACEMENTS_rough_1'
)


def test_read_direction_before_target_file():
    """Check that the old Cu(111)-Te file gives the reordered line."""
    with pytest.raises(DisplacementsSyntaxError) as exc_info:
        DisplacementsFile().read(ROUGH_PATH)
    assert (
        'Write "Cu_def L(2-4) xy = -0.10 0.10 0.01"'
        in exc_info.value.message
    )


_OFFSETS = '== OFFSETS\n    geo Cu_surf z = 0.05\n'
_SEARCH = '== SEARCH a\n    = GEO_DELTA\n    Cu_surf z = -0.1 0.1 0.05\n'


@parametrize(
    'contents, msg',
    (
        (_OFFSETS + _OFFSETS + _SEARCH,
         'There can only be one OFFSETS block at the beginning'),
        (_SEARCH + _OFFSETS,
         'The OFFSETS block must be at the beginning of the file'),
    ),
)
def test_read_offsets_misplaced(contents, msg, tmp_path):
    """Check the error for a duplicated or late OFFSETS block."""
    path = tmp_path / 'DISPLACEMENTS'
    path.write_text(contents)
    with pytest.raises(OffsetsNotAtBeginningError, match=msg):
        DisplacementsFile().read(path)


def test_read_nested_loops():
    """Check that <loop> markers in a file give nested LoopBlocks."""
    df = DisplacementsFile()
    df.read(NESTED_LOOPS_PATH)

    outer, last = df.children
    assert isinstance(outer, LoopBlock)
    assert isinstance(last, SearchBlock)
    assert last.label == 'd'
    first, inner, third = outer.children
    assert isinstance(first, SearchBlock)
    assert isinstance(inner, LoopBlock)
    assert isinstance(third, SearchBlock)
    assert [c.label for c in (first, third)] == ['a', 'c']
    assert [type(c) for c in inner.children] == [SearchBlock]
    assert inner.children[0].label == 'b'


_NESTED_LOOPS_SEQUENCES = (
    # Each loop runs twice: the second pass does not improve R.
    pytest.param((0.5,)*9, 'abbcabbcd', id='constant R'),
    # The inner loop repeats as long as R improves. When the outer
    # loop repeats, the inner one starts afresh and runs twice.
    pytest.param((0.5, 0.4, 0.3) + (0.2,)*7, 'abbbcabbcd', id='improving R'),
    # R improves during the outer loop's second pass: a third pass.
    pytest.param((0.5, 0.5, 0.4, 0.4, 0.4) + (0.3,)*8,
                 'abbcabbcabbcd',
                 id='outer improves'),
)


@parametrize('rfacs,expected', _NESTED_LOOPS_SEQUENCES)
def test_next_nested_loops_sequence(rfacs, expected):
    """Check the search blocks next() returns for nested loops."""
    disp_file = DisplacementsFile()
    disp_file.read(NESTED_LOOPS_PATH)
    labels = ''.join(disp_file.next(current_rfac=rfac).label
                     for rfac in rfacs)
    assert labels == expected
    with pytest.raises(StopIteration):
        disp_file.next(current_rfac=rfacs[-1])


def test_read_fractional_directions():
    """Check the directions read from a file with ab[i j] directions."""
    disp_file = DisplacementsFile()
    disp_file.read(LATERAL_1_PATH)
    search, = disp_file.children
    a_dir = FractionalDirectionToken('ab[1 0]')
    a2b_dir = FractionalDirectionToken('ab[2 1]')
    directions = [line.direction for line in search.geo_delta_lines]
    assert directions == [a_dir, a_dir, a2b_dir, a2b_dir, a_dir, a2b_dir]


def _exhaust(displacements_file):
    """Return all search blocks that next() yields until it stops."""
    blocks = []
    while True:
        try:
            blocks.append(displacements_file.next(current_rfac=0.1))
        except StopIteration:
            return blocks


def test_next(displacements_file_path):
    """Test the next() method of DisplacementsFile."""
    df = DisplacementsFile()
    df.read(displacements_file_path)
    assert all(block is not None for block in _exhaust(df))


def test_reset_iterator(displacements_file_path):
    """Check that reset_iterator restarts a finished file."""
    df = DisplacementsFile()
    df.read(displacements_file_path)
    first_pass = _exhaust(df)
    df.reset_iterator()
    assert _exhaust(df) == first_pass


_MODES = '= GEO_DELTA\n    Cu_surf z = -0.1 0.1 0.05\n'


def _read_text(contents, tmp_path):
    """Return a DisplacementsFile read from a file with `contents`."""
    path = tmp_path / 'DISPLACEMENTS'
    path.write_text(contents)
    disp_file = DisplacementsFile()
    disp_file.read(path)
    return disp_file


def _search_labels(disp_file):
    """Return the labels of all search blocks in `disp_file`."""
    return [
        b.label for b in disp_file.descendants if isinstance(b, SearchBlock)
    ]


def _mode_blocks(disp_file):
    """Return (type, lines) of each mode block, for each search."""
    return [
        [(type(mode), [str(line) for line in mode.lines])
         for mode in search.children]
        for search in disp_file.children
    ]


def test_read_implicit_search_block(tmp_path):
    """Check that mode blocks with no '== SEARCH' make a search."""
    implicit = _read_text(_MODES, tmp_path)
    explicit = _read_text('== SEARCH\n' + _MODES, tmp_path)
    assert len(implicit.children) == 1
    assert isinstance(implicit.children[0], SearchBlock)
    assert _mode_blocks(implicit) == _mode_blocks(explicit)
    assert _search_labels(implicit) == _search_labels(explicit)


def test_read_implicit_then_explicit_search_block(tmp_path):
    """Check that '== SEARCH' after mode lines starts a new block."""
    disp_file = _read_text(_MODES + _SEARCH, tmp_path)
    assert _search_labels(disp_file) == ['search_1', 'a']


def test_read_implicit_search_block_after_offsets(tmp_path):
    """Check that mode blocks right after OFFSETS are a search block."""
    disp_file = _read_text(_OFFSETS + _MODES, tmp_path)
    assert _search_labels(disp_file) == ['search_1']


def test_read_unnamed_search_blocks(tmp_path):
    """Check that unnamed '== SEARCH' blocks are distinct and named."""
    disp_file = _read_text(
        '== SEARCH\n' + _MODES + '== SEARCH\n' + _MODES,
        tmp_path,
    )
    first, second = disp_file.children
    assert first.parent is second.parent is disp_file
    assert [first.label, second.label] == ['search_1', 'search_2']
    assert not set(first.children) & set(second.children)


_BAD_FILES = (
    pytest.param(
        '</loop>\n' + _SEARCH,
        r'Found "</loop>" without a matching "<loop>"',
        id='stray_loop_end',
    ),
    pytest.param(
        '<loop>\n' + _MODES + '</loop>\n',
        r'right after "<loop>"',
        id='loop_without_search',
    ),
    pytest.param(
        '<loop> == SEARCH a\n' + _MODES,
        (
            r'Found "<loop>" together with other text in line "<loop> == '
            r'SEARCH a"'
        ),
        id='loop_on_search_header',
    ),
    pytest.param(
        '== SEARCH a <loop>\n' + _MODES,
        (
            r'Found "<loop>" together with other text in line "== SEARCH '
            r'a <loop>"'
        ),
        id='loop_at_end_of_header',
    ),
    pytest.param(
        _SEARCH + '</loop> == SEARCH b\n',
        r'Found "</loop>" together with other text',
        id='end_tag_with_text',
    ),
    pytest.param(
        '<loop>\n' + _SEARCH,
        r'Found "<loop>" without a matching "</loop>"',
        id='unclosed_loop',
    ),
    pytest.param(
        '<loop>\n' + _SEARCH + '</loop>\n' + _MODES,
        r'Line "= GEO_DELTA" is not inside a search block',
        id='mode_after_loop',
    ),
)


@parametrize('contents, msg', _BAD_FILES)
def test_read_bad_loop_tags(contents, msg, tmp_path):
    """Check the error for misplaced or unbalanced loop tags."""
    with pytest.raises(DisplacementsSyntaxError, match=msg):
        _read_text(contents, tmp_path)


@parametrize(
    'target, number',
    (('O 0 z', '0'), ('O L(0) z', '0'), ('O -2 z', '-2')),
)
def test_read_bad_target_number_in_message(target, number, tmp_path):
    """Check that the error names the bad atom or layer number."""
    contents = f'== SEARCH a\n= GEO_DELTA\n    {target} = -0.1 0.1 0.05\n'
    with pytest.raises(DisplacementsSyntaxError) as exc_info:
        _read_text(contents, tmp_path)
    assert re.search(
        rf'Invalid (atom|layer) number {number}', exc_info.value.message
    )


@parametrize(
    'line',
    (ModeHeaderLine('= GEO_DELTA'), OffsetsHeaderLine('= OFFSETS')),
    ids=('mode', 'offsets'),
)
def test_unparsable_line_message_quotes_line(line):
    """Check that the generic error quotes the line, not its repr."""
    assert _unparsable_line_message(line) == f'Unable to parse line "{line}".'


@parametrize(
    'contents, expected',
    (
        ('== SEARCH\n' + _MODES + _SEARCH.replace('a', 'search_1'),
         ['search_2', 'search_1']),
        (_SEARCH + _SEARCH, ['a', 'a 2']),
        (_SEARCH + _SEARCH + _SEARCH, ['a', 'a 2', 'a 3']),
        (_SEARCH + '== SEARCH\n' + _MODES, ['a', 'search_2']),
        ('== SEARCH\n' + _MODES + _SEARCH, ['search_1', 'a']),
    ),
)
def test_read_unique_search_labels(contents, expected, tmp_path):
    """Check that search labels are unique and skip used ones."""
    assert _search_labels(_read_text(contents, tmp_path)) == expected


def test_search_label_round_trip(tmp_path):
    """Check that a generated label is rendered in the header line."""
    disp_file = _read_text(
        '== SEARCH\n' + _MODES + '== SEARCH\n' + _MODES, tmp_path
    )
    headers = [
        line for line in str(disp_file).splitlines() if '== SEARCH' in line
    ]
    labels = [h.split('== SEARCH', 1)[1].strip() for h in headers]
    assert labels == ['search_1', 'search_2']
    reread = _read_text(
        '\n'.join(f'== SEARCH {label}\n{_MODES}' for label in labels), tmp_path
    )
    assert _search_labels(reread) == labels


def test_str_shows_block_headers(tmp_path):
    """Check that rendering a file names its blocks, not a method."""
    rendered = str(_read_text(_OFFSETS + _SEARCH, tmp_path))
    assert 'bound method' not in rendered
    assert 'GEO_DELTA' in rendered
    assert 'OFFSET' in rendered


def test_read_empty_mode_block(tmp_path):
    """Check that an empty mode block is reported by name."""
    with pytest.raises(
        DisplacementsSyntaxError, match='Found empty GEO_DELTA block'
    ):
        _read_text('== SEARCH a\n= GEO_DELTA\n', tmp_path)
