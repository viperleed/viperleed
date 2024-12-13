"""Module displacements/regex."""

__authors__ = ('Alexander M. Imre (@amimre)',)
__created__ = '2024-10-03'

import contextlib
import re

SEARCH_HEADER_PATTERN = re.compile(r'^=+\s+(?i:search)\s+(.*)$')
SECTION_HEADER_PATTERN = re.compile(
    r'^=+\s*(OFFSETS|GEO_DELTA|VIB_DELTA|OCC_DELTA|CONSTRAIN)$'
)

# components
LABEL_PATTERN = (
    r'(?:\*|\*?\w+\*?)'  # Non-capturing group to avoid name conflicts
)
WHICH_PATTERN = r'(?P<which>L\(\d+(-\d+)?\)|\d+(-\d+)?(\s+\d+(-\d+)?)*)?'
DIRECTION_PATTERN = r'(?P<direction>[a-zA-Z]+(?:\[[^\]]+\]|\([^\)]+\))?)'
START_PATTERN = r'(?P<start>-?\d+(\.\d+)?)'
STOP_PATTERN = r'(?P<stop>-?\d+(\.\d+)?)'
STEP_PATTERN = r'(?P<step>-?\d+(\.\d+)?(?:\s+))?'
VALUE_PATTERN = r'(?P<value>-?\d+(\.\d+)?)'

# Use the label pattern to construct the targets pattern
TARGETS_PATTERN = (
    rf'(?P<targets>{LABEL_PATTERN}(?:\s+\d+|\s+\d+-\d+)*'
    rf'(?:\s*,\s*{LABEL_PATTERN}(?:\s+\d+|\s+\d+-\d+)*)*)'
)

CHEM_BLOCKS_PATTERN = (
    r'(?P<chem_blocks>(?P<chem>\w+)\s+'
    + START_PATTERN
    + r'\s+'
    + STOP_PATTERN
    + r'(?:\s+'
    + STEP_PATTERN
    + r'(?:\s*,\s*(?P<additional_blocks>.+))?)?)'
)

# Patterns
OFFSETS_LINE_PATTERN = re.compile(
    rf'^(?P<type>geo|vib|occ)\s+{TARGETS_PATTERN}'
    rf'(?:\s+{DIRECTION_PATTERN})?\s*=\s*{VALUE_PATTERN}$'
)

GEO_LINE_PATTERN = re.compile(
    rf'^(?P<label>{LABEL_PATTERN})(?:\s+{WHICH_PATTERN})?'
    rf'\s+{DIRECTION_PATTERN}\s*=\s*{START_PATTERN}\s+{STOP_PATTERN}'
    rf'(?:\s+{STEP_PATTERN})?$'
)

VIB_LINE_PATTERN = re.compile(
    rf'^(?P<label>{LABEL_PATTERN})(?:\s+{WHICH_PATTERN})?'
    rf'\s*=\s*{START_PATTERN}\s+{STOP_PATTERN}(?:\s+{STEP_PATTERN})?$'
)

OCC_LINE_PATTERN = re.compile(
    rf'^(?P<label>{LABEL_PATTERN})(?:\s+{WHICH_PATTERN})?'
    rf'\s*=\s*{CHEM_BLOCKS_PATTERN}$'
)

CONSTRAIN_LINE_PATTERN = re.compile(
    rf'^(?P<type>geo|vib|occ)\s+{TARGETS_PATTERN}'
    rf'(?:\s+{DIRECTION_PATTERN})?\s*=\s*(?P<value>linked|-?\d+(\.\d+)?)$'
)


def match_geo_line(line):
    """Match and parse a GEO_DELTA line, returning the values."""
    match = GEO_LINE_PATTERN.match(line)
    if match is None:
        return None
    label = match.group('label')
    which = match.group('which')  # optional, can be None
    direction = match.group('dir')
    start = float(match.group('start'))
    stop = (
        float(match.group('stop')) if match.group('stop') is not None else None
    )
    step = (
        float(match.group('step')) if match.group('step') is not None else None
    )
    return label, which, direction, start, stop, step


def match_vib_line(line):
    """Match and parse a VIB_DELTA line, returning the values as floats."""
    match = VIB_LINE_PATTERN.match(line)
    if match is None:
        return None

    label = match.group('label')
    which = match.group('which')
    start = float(match.group('start'))
    stop = (
        float(match.group('stop')) if match.group('stop') is not None else None
    )
    step = (
        float(match.group('step')) if match.group('step') is not None else None
    )

    return label, which, start, stop, step


def match_occ_line(line):
    """Match and parse OCC_DELTA line, returning chemical blocks and ranges."""
    match = OCC_LINE_PATTERN.match(line)
    if match is None:
        return None

    label = match.group('label')
    which = match.group('which')

    # Parse the first chemical block
    chem_blocks = []
    chem = match.group('chem')
    start = float(match.group('start'))
    stop = (
        float(match.group('stop')) if match.group('stop') is not None else None
    )
    step = (
        float(match.group('step')) if match.group('step') is not None else None
    )
    chem_blocks.append((chem, start, stop, step))

    # Handle additional blocks, if present
    additional_blocks = match.group('additional_blocks')
    if additional_blocks:
        for block in additional_blocks.split(','):
            _block = block.strip()
            m = re.match(CHEM_BLOCK_PATTERN, _block)
            if m:
                chem = m.group('chem')
                start = float(m.group('start'))
                stop = (
                    float(m.group('stop'))
                    if m.group('stop') is not None
                    else None
                )
                step = (
                    float(m.group('step'))
                    if m.group('step') is not None
                    else None
                )
                chem_blocks.append((chem, start, stop, step))

    return label, which, chem_blocks


def match_constrain_line(line):
    """Match and parse line, returning type, targets, direction, and value."""
    match = CONSTRAIN_LINE_PATTERN.match(line)
    if match is None:
        return None

    offset_type = match.group(
        'type'
    )  # Type can be 'geo', 'vib', or 'occ'     # TODO: make into Enum
    targets = match.group('targets')  # Multiple comma-separated targets
    direction = match.group('direction')  # Optional complex direction for geo
    value = match.group('value')

    # Convert `value` to float if it's a number; otherwise, keep it as a string
    with contextlib.suppress(ValueError):
        value = float(value)

    return offset_type, targets, direction, value


def match_offsets_line(line):
    """Match and parse line, returning type, targets, direction, and value."""
    match = OFFSETS_LINE_PATTERN.match(line)
    if match is None:
        return None

    offset_type = match.group('type')  # Type can be 'geo', 'vib', or 'occ'
    targets = match.group('targets').strip()  # Multiple comma-separated targets
    direction = match.group('direction')  # Optional complex direction for geo
    value = match.group('value')

    # Convert `value` to float if it's a number; otherwise, keep it as a string
    with contextlib.suppress(ValueError):
        value = float(value)

    return offset_type, targets, direction, value
