import pytest

from viperleed_jax.files.displacements.regex import (
    SEARCH_HEADER_PATTERN,
    SECTION_HEADER_PATTERN,
    match_constrain_line,
    match_geo_line,
    match_occ_line,
    match_offsets_line,
    match_vib_line,
)

# Test cases for SECTION_HEADER_PATTERN
TEST_LINES_SECTION = {
    '= GEO_DELTA': 'GEO_DELTA',
    '= VIB_DELTA': 'VIB_DELTA',
    '= OCC_DELTA': 'OCC_DELTA',
    '= CONSTRAIN': 'CONSTRAIN',
    '= OFFSETS': 'OFFSETS',
    '== GEO_DELTA': 'GEO_DELTA',
    'GEO_DELTA': None, # No '=' at the start
    'NOT_A_SECTION': None,
    'GEO_DELTA_EXTRA': None,
}

# Test cases for SEARCH_HEADER_PATTERN
TEST_LINES_SEARCH = {
    '== SEARCH label1': ('label1',),
    '== search another_label': ('another_label',),
    '== SEARCH complicated_label-123': ('complicated_label-123',),
    '== SEARCH With_Spaces In_Label': ('With_Spaces In_Label',),
}

# Test cases for GEO_LINE_PATTERN
TEST_LINES_GEOMETRY = {
    'Ir L(1-6) z = -0.05 0.05 0.01': ('Ir', 'L(1-6)', 'z', -0.05, 0.05, 0.01),
    'Si x = -0.005 0.005 0.0005': ('Si', None, 'x', -0.005, 0.005, 0.0005),
    'Ir L(1-6) xy[1 0] = -0.03 0.03 0.01': (
        'Ir',
        'L(1-6)',
        'xy[1 0]',
        -0.03,
        0.03,
        0.01,
    ),
    'Abc L(1-6) xy = -0.03 0.03 0.01': (
        'Abc',
        'L(1-6)',
        'xy',
        -0.03,
        0.03,
        0.01,
    ),
    'Ir L(1-6) z = -0.03 0.03 0.01': ('Ir', 'L(1-6)', 'z', -0.03, 0.03, 0.01),
    'Au 1 2 4 x = -0.03 0.03 0.01': ('Au', '1 2 4', 'x', -0.03, 0.03, 0.01),
    'B 1 2 4 abc = -0.03 0.03 0.01': ('B', '1 2 4', 'abc', -0.03, 0.03, 0.01),
    'Cd L(1-6) azi(ab[c1 c2]) = -0.03 0.03 0.01': (
        'Cd',
        'L(1-6)',
        'azi(ab[c1 c2])',
        -0.03,
        0.03,
        0.01,
    ),
    'E 5 ab[n1 n2] = -0.005 0': ('E', '5', 'ab[n1 n2]', -0.005, 0.0, None),
}

# Test cases for VIB_DELTA lines
TEST_LINES_VIB = {
    'O 1 = -0.05 0.05 0.02': ('O', '1', -0.05, 0.05, 0.02),
    'Ir_top = -0.05 0.05 0.01': ('Ir_top', None, -0.05, 0.05, 0.01),
    'H 5 = -0.03 0.03': ('H', '5', -0.03, 0.03, None),  # No step
    'C L(1-4) = -0.1 0.1 0.05': ('C', 'L(1-4)', -0.1, 0.1, 0.05),  # With L(1-4)
}

# Test cases for OCC_DELTA lines with optional steps
TEST_LINES_OCC = {
    'O 1 = O 0.8 1.0 0.05': ('O', '1', [('O', 0.8, 1.0, 0.05)]),
    'M_top = Fe 0.4 0.6 0.05, Ni 0.6 0.4 0.05': (
        'M_top',
        None,
        [('Fe', 0.4, 0.6, 0.05), ('Ni', 0.6, 0.4, 0.05)],
    ),
    'M_top = Fe 0.3 0.5, Ni 0.6 0.4, Ti 0.1': (
        'M_top',
        None,
        [
            ('Fe', 0.3, 0.5, None),
            ('Ni', 0.6, 0.4, None),
            ('Ti', 0.1, None, None),
        ],
    ),  # Missing steps
    'Si 1 = Si 0.2 0.5, Ge 0.5': (
        'Si',
        '1',
        [('Si', 0.2, 0.5, None), ('Ge', 0.5, None, None)],
    ),
    'Cu = Cu 0.6 1.0, Zn 0.4 0.0': (
        'Cu',
        None,
        [('Cu', 0.6, 1.0, None), ('Zn', 0.4, 0.0, None)],
    ),
}

# Test cases for CONSTRAIN lines
TEST_LINES_CONSTRAIN = {
    'geo O L(1-2), Ir L(1) = linked': (
        'geo',
        ['O L(1-2)', 'Ir L(1)'],
        'linked',
    ),
    'vib Ir_top = linked': ('vib', ['Ir_top'], 'linked'),
    'vib Ir_top = -0.03': ('vib', ['Ir_top'], -0.03),  # Fix to a specific value
    'geo C L(1), N L(2) = 0.1': (
        'geo',
        ['C L(1)', 'N L(2)'],
        0.1,
    ),  # Geo constraint to 0.1
    'occ Fe L(3), Ni = linked': ('occ', ['Fe L(3)', 'Ni'], 'linked'),
}

# Test cases for OFFSETS lines
TEST_LINES_OFFSETS = {
    'geo O L(1-2) z = 0.05': ('geo', 'O L(1-2) z', 0.05),
    'geo Ir_top x = -0.02': ('geo', 'Ir_top x', -0.02),
    'vib Si 1 = 0.01': ('vib', 'Si 1', 0.01),
    'vib O_top = -0.005': ('vib', 'O_top', -0.005),
    'occ Fe 1 = 0.3': ('occ', 'Fe 1', 0.3),
    'occ M_top = Fe 0.6, Ni 0.4': (
        'occ',
        'M_top',
        'Fe 0.6, Ni 0.4',
    ),  # Complex value (not float)
}


@pytest.mark.parametrize(
    'input, expected',
    TEST_LINES_GEOMETRY.items(),
    ids=TEST_LINES_GEOMETRY.keys(),
)
def test_geo_line_regex(input, expected):
    """Check that the regex for the geometry line works as expected."""
    match = match_geo_line(input)
    assert match is not None
    label, which, dir, start, stop, step = match
    e_label, e_which, e_dir, e_start, e_stop, e_step = expected
    assert label == e_label
    assert which == e_which or e_which is None  # optional
    assert dir == e_dir
    assert start == pytest.approx(e_start)
    assert stop == pytest.approx(e_stop)
    assert step == pytest.approx(e_step) or step is None  # optional


@pytest.mark.parametrize(
    'input, expected', TEST_LINES_SEARCH.items(), ids=TEST_LINES_SEARCH.keys()
)
def test_search_header_regex(input, expected):
    """Check that the regex for the search header works as expected."""
    match = SEARCH_HEADER_PATTERN.match(input)
    assert match is not None
    assert match.group(1) == expected[0]


@pytest.mark.parametrize(
    'input, expected', TEST_LINES_SECTION.items(), ids=TEST_LINES_SECTION.keys()
)
def test_section_header_regex(input, expected):
    """Check that the regex for the section header works as expected."""
    match = SECTION_HEADER_PATTERN.match(input)
    if expected is None:
        assert match is None
        return
    assert match is not None
    assert match.group(1) == expected


@pytest.mark.parametrize(
    'input, expected', TEST_LINES_VIB.items(), ids=TEST_LINES_VIB.keys()
)
def test_vib_line_regex(input, expected):
    """Check that the regex for VIB_DELTA lines works as expected."""
    match = match_vib_line(input)
    assert match is not None
    label, which, start, stop, step = match
    e_label, e_which, e_start, e_stop, e_step = expected
    assert label == e_label
    assert which == e_which or e_which is None  # Optional
    assert start == pytest.approx(e_start)
    assert stop == pytest.approx(e_stop) or stop is None
    assert step == pytest.approx(e_step) or step is None  # Optional


@pytest.mark.parametrize(
    'input, expected', TEST_LINES_OCC.items(), ids=TEST_LINES_OCC.keys()
)
def test_occ_line_regex(input, expected):
    """Check that the regex for the OCC_DELTA line works as expected."""
    match = match_occ_line(input)
    assert match is not None
    label, which, chem_blocks = match
    e_label, e_which, e_chem_blocks = expected
    assert label == e_label
    assert which == e_which or e_which is None
    assert len(chem_blocks) == len(e_chem_blocks)

    for (chem, start, stop, step), (e_chem, e_start, e_stop, e_step) in zip(
        chem_blocks, e_chem_blocks
    ):
        assert chem == e_chem
        assert start == pytest.approx(e_start)
        assert stop == pytest.approx(e_stop) or stop is None
        assert step == pytest.approx(e_step) or step is None


@pytest.mark.parametrize(
    'input, expected',
    TEST_LINES_CONSTRAIN.items(),
    ids=TEST_LINES_CONSTRAIN.keys(),
)
def test_constrain_line_regex(input, expected):
    """Check that the regex for CONSTRAIN lines works as expected."""
    match = match_constrain_line(input)
    assert match is not None
    constraint_type, parameters, value = match
    e_constraint_type, e_parameters, e_value = expected
    assert constraint_type == e_constraint_type
    assert parameters == e_parameters
    assert value == e_value


@pytest.mark.parametrize(
    'input, expected', TEST_LINES_OFFSETS.items(), ids=TEST_LINES_OFFSETS.keys()
)
def test_offsets_line_regex(input, expected):
    """Check that the regex for OFFSETS lines works as expected."""
    match = match_offsets_line(input)
    assert match is not None
    offset_type, parameters, value = match
    e_offset_type, e_parameters, e_value = expected
    assert offset_type == e_offset_type
    assert parameters == e_parameters
    assert value == e_value
