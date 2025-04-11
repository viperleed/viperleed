"""Module displacements/regex."""

__authors__ = ('Alexander M. Imre (@amimre)',)
__created__ = '2024-10-03'

import contextlib
import re

from viperleed_jax.perturbation_type import PerturbationType

SEARCH_HEADER_PATTERN = re.compile(r'^=+\s+(?i:search)\s+(.*)$')
SECTION_HEADER_PATTERN = re.compile(
    r'^=+\s*(OFFSETS|GEO_DELTA|VIB_DELTA|OCC_DELTA|CONSTRAIN)$'
)
LOOP_START_PATTERN = re.compile(r'<loop>')
LOOP_END_PATTERN = re.compile(r'<\\loop>|</loop>')

# Common components
LABEL_PATTERN = (
    r'(?P<label>\*|\*?\w+\*?)'  # Match '*' or labels with/without wildcards
)
WHICH_PATTERN = r'(?P<which>L\(\d+(-\d+)?\)|\d+(-\d+)?(\s+\d+(-\d+)?)*)?'

DIRECTION_PATTERN = (
    r'(?P<direction>('
    r'[abcxyz]'  # 1D single-letter shorthand
    r'|ab|xy'    # 2D multi-letter in-plane shorthand
    r'|abc|xyz'  # 3D multi-letter shorthand directions
    r'|\[\s*-?\d+(?:\.\d+)?(?:\s+-?\d+(?:\.\d+)?)*\s*\]'  # [1 0] or [1.0 -2]
    r'|ab\[\s*-?\d+(?:\.\d+)?(?:\s+-?\d+(?:\.\d+)?)*\s*\]'  # ab[1 2]
    r'|xy\[\s*-?\d+(?:\.\d+)?(?:\s+-?\d+(?:\.\d+)?)*\s*\]'  # xy[0 1]
    r'|azi\((?:ab|xy)?\[\s*-?\d+(?:\.\d+)?(?:\s+-?\d+(?:\.\d+)?)*\s*\]\)'  # azi(ab[1 2])
    r'|r\((?:ab|xy)?\[\s*-?\d+(?:\.\d+)?(?:\s+-?\d+(?:\.\d+)?)*\s*\]\)'  # r([1.0 -1])
    r'))'
)


START_PATTERN = r'(?P<start>-?\d+(\.\d+)?)'
STOP_PATTERN = r'(?P<stop>-?\d+(\.\d+)?)'
STEP_PATTERN = r'(?P<step>-?\d+(\.\d+)?)'
VALUE_PATTERN = r'-?\d+(\.\d+)?'

# TARGETS_PATTERN
TARGETS_PATTERN = r'(?P<targets>[^=]+?)'

CHEM_BLOCKS_PATTERN = (
    r'(?P<chem_blocks>(?P<chem>\w+)\s+'
    + START_PATTERN
    + r'\s+'
    + STOP_PATTERN
    + r'(?:\s+'
    + STEP_PATTERN
    + r'(?:\s*,\s*(?P<additional_blocks>.+))?)?)'
)

CHEM_LITERAL_BLOCK_PATTERN = rf'[A-Z][a-z]?(?:\s+{VALUE_PATTERN}){{1,3}}'
CHEM_LITERAL_BLOCKS_PATTERN = (
    rf'{CHEM_LITERAL_BLOCK_PATTERN}(?:\s*,\s*{CHEM_LITERAL_BLOCK_PATTERN})*'
)

OCC_LINE_PATTERN = re.compile(
    rf'^{LABEL_PATTERN}(?:\s+{WHICH_PATTERN})?\s*=\s*'
    rf'(?P<chem_blocks>{CHEM_LITERAL_BLOCKS_PATTERN})$'
)

# Line patterns
OFFSETS_LINE_PATTERN = re.compile(
    rf'^(?P<type>geo|vib|occ)\s+{TARGETS_PATTERN}'
    rf'(?:\s+{DIRECTION_PATTERN})?\s*=\s*'
    rf'(?P<value>{VALUE_PATTERN}|{CHEM_LITERAL_BLOCKS_PATTERN})$'
)

GEO_LINE_PATTERN = re.compile(
    rf'^{LABEL_PATTERN}'         # Match the label
    rf'(?:\s+{WHICH_PATTERN})?'  # Match the optional 'which' group
    rf'\s+{DIRECTION_PATTERN}'   # Match the direction
    rf'\s*=\s*{START_PATTERN}\s+{STOP_PATTERN}'  # Match the ranges
    rf'(?:\s+{STEP_PATTERN})?$'  # Match the optional step size
)

VIB_LINE_PATTERN = re.compile(
    rf'^{LABEL_PATTERN}(?:\s+{WHICH_PATTERN})?'
    rf'\s*=\s*{START_PATTERN}\s+{STOP_PATTERN}(?:\s+{STEP_PATTERN})?$'
)

CONSTRAIN_LINE_PATTERN = re.compile(
    rf'^(?P<type>geo|vib|occ)\s+{TARGETS_PATTERN}'
    rf'(?:\s+{DIRECTION_PATTERN})?\s*=\s*(?P<value>linked|{VALUE_PATTERN})$'
)


def match_geo_line(line):
    """Match and parse a GEO_DELTA line, returning the values."""
    match = GEO_LINE_PATTERN.match(line)
    if match is None:
        return None
    label = match.group('label')
    which = match.group('which')  # optional, can be None
    direction = match.group('direction')
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
    """Match and parse OCC_DELTA line, returning chemical blocks as parametric triples.

    Returns a tuple: (label, which, list of (element, start, stop, step))

    - If only 1 value is given: stop = start, step = None
    - If 2 values: step = None
    - If 3 values: step is set
    """
    match = OCC_LINE_PATTERN.match(line)
    if match is None:
        return None

    label = match.group('label')
    which = match.group('which')

    chem_blocks = []

    # Extract the chemical blocks from the match
    # TODO: refactor to separate function
    # Also TODO: maybe we should introduce a chemical element Enum

    chem_blocks_str = match.group('chem_blocks')
    for block in chem_blocks_str.split(','):
        tokens = block.strip().split()
        if len(tokens) < 2:
            continue  # Not enough data to process

        element = tokens[0]
        try:
            numbers = [float(tok) for tok in tokens[1:]]
        except ValueError:
            continue  # Skip malformed values

        if len(numbers) == 1:
            start = stop = numbers[0]
            step = None
        elif len(numbers) == 2:
            start, stop = numbers
            step = None
        elif len(numbers) == 3:
            start, stop, step = numbers
        else:
            continue  # Invalid number of numeric values

        chem_blocks.append((element, start, stop, step))

    return label, which, chem_blocks


def match_constrain_line(line):
    """Match and parse line, returning type, targets, direction, and value."""
    match = CONSTRAIN_LINE_PATTERN.match(line)
    if match is None:
        return None

    offset_type = PerturbationType.from_string(match.group('type'))
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

    offset_type = PerturbationType.from_string(match.group('type'))
    targets = match.group('targets').strip()  # Multiple comma-separated targets
    direction = match.group('direction')  # Optional complex direction for geo
    value = match.group('value')

    # Convert `value` to float if it's a number; otherwise, keep it as a string
    with contextlib.suppress(ValueError):
        value = float(value)

    return offset_type, targets, direction, value
