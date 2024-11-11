import re

SEARCH_HEADER_PATTERN = re.compile(r"^=+\s+(?i:search)\s+(.*)$")
SECTION_HEADER_PATTERN = re.compile(
    r"^=+\s*(OFFSETS|GEO_DELTA|VIB_DELTA|OCC_DELTA|CONSTRAIN)$"
)

OFFSETS_LINE_PATTERN = re.compile(
    r"^(?P<type>geo|vib|occ)\s+"
    r"(?P<targets>[^\s=,]+(?:\s+\d+)*(?:\s*,\s*[^\s=,]+(?:\s+\d+)*)*)"
    r"(?:\s+(?P<direction>[a-zA-Z]+(?:\[[^\]]+\]|\([^\)]+\))?))?\s*=\s*"
    r"(?P<value>-?\d+(\.\d+)?)$"
)

GEO_LINE_PATTERN = re.compile(
    r"^(?P<label>\*?\w+)"
    r"(?:\s+(?P<which>L\(\d+(-\d+)?\)|\d+(\s+\d+)*))?"
    r"\s+(?P<dir>[a-zA-Z]+(?:\[[^\]]+\]|\([^\)]+\))?)"
    r"\s*=\s*(?P<start>-?\d+(\.\d+)?)"
    r"\s+(?P<stop>-?\d+(\.\d+)?)"
    r"(?:\s+(?P<step>-?\d+(\.\d+)?))?$"
)
VIB_LINE_PATTERN = re.compile(
    r"^(?P<label>\*?\w+)(?:\s+(?P<which>L\(\d+(-\d+)?\)|\d+(\s+\d+)*))?"
    r"\s*=\s*(?P<start>-?\d+(\.\d+)?)"
    r"\s+(?P<stop>-?\d+(\.\d+)?)"
    r"(?:\s+(?P<step>-?\d+(\.\d+)?))?$"
)
OCC_LINE_PATTERN = re.compile(
    r"^(?P<label>\*?\w+)"
    r"(?:\s+(?P<which>L\(\d+(-\d+)?\)|\d+(\s+\d+)*))?"
    r"\s*=\s*(?P<chem_blocks>(?P<chem>\w+)\s+(?P<start>-?\d+(\.\d+)?)"
    r"\s+(?P<stop>-?\d+(\.\d+)?)(?:\s+(?P<step>-?\d+(\.\d+)?))?"
    r"(?:\s*,\s*(?P<additional_blocks>.+))?)$"
)
CONSTRAIN_LINE_PATTERN = re.compile(
    r"^(?P<type>geo|vib|occ)\s+(?P<targets>[^\s=,]+"
    r"(?:\s*[^\s=,]+)*(?:,\s*[^\s=,]+)*)"
    r"(?:\s+(?P<direction>[a-zA-Z]+(?:\[[^\]]+\]|\([^\)]+\))?))?"
    r"\s*=\s*(?P<value>linked|-?\d+(\.\d+)?)$"
)


def match_geo_line(line):
    match = GEO_LINE_PATTERN.match(line)
    if match is None:
        return None
    label = match.group("label")
    which = match.group("which")  # optional, can be None
    dir = match.group("dir")
    start = float(match.group("start"))
    stop = (
        float(match.group("stop")) if match.group("stop") is not None else None
    )
    step = (
        float(match.group("step")) if match.group("step") is not None else None
    )
    return label, which, dir, start, stop, step


def match_vib_line(line):
    """Match and parse a VIB_DELTA line, returning the values as floats."""
    match = VIB_LINE_PATTERN.match(line)
    if match is None:
        return None

    label = match.group("label")
    which = match.group("which")
    start = float(match.group("start"))
    stop = (
        float(match.group("stop")) if match.group("stop") is not None else None
    )
    step = (
        float(match.group("step")) if match.group("step") is not None else None
    )

    return label, which, start, stop, step


def match_occ_line(line):
    """Match and parse an OCC_DELTA line, returning chemical blocks and their ranges."""
    match = OCC_LINE_PATTERN.match(line)
    if match is None:
        return None

    label = match.group("label")
    which = match.group("which")

    # Parse the first chemical block
    chem_blocks = []
    chem = match.group("chem")
    start = float(match.group("start"))
    stop = (
        float(match.group("stop")) if match.group("stop") is not None else None
    )
    step = (
        float(match.group("step")) if match.group("step") is not None else None
    )
    chem_blocks.append((chem, start, stop, step))

    # Handle additional blocks, if present
    additional_blocks = match.group("additional_blocks")
    if additional_blocks:
        for block in additional_blocks.split(","):
            block = block.strip()
            m = re.match(
                r"(?P<chem>\w+)\s+(?P<start>-?\d+(\.\d+)?)(?:\s+(?P<stop>-?\d+(\.\d+)?)(?:\s+(?P<step>-?\d+(\.\d+)?))?)?",
                block,
            )
            if m:
                chem = m.group("chem")
                start = float(m.group("start"))
                stop = (
                    float(m.group("stop"))
                    if m.group("stop") is not None
                    else None
                )
                step = (
                    float(m.group("step"))
                    if m.group("step") is not None
                    else None
                )
                chem_blocks.append((chem, start, stop, step))

    return label, which, chem_blocks


def match_constrain_line(line):
    """Match and parse an OFFSETS line, returning the type, targets, direction, and value."""
    match = CONSTRAIN_LINE_PATTERN.match(line)
    if match is None:
        return None

    offset_type = match.group("type")  # Type can be 'geo', 'vib', or 'occ'     # TODO: make into Enum
    targets = match.group("targets")  # Multiple comma-separated targets
    direction = match.group("direction")  # Optional complex direction for geo
    value = match.group("value")

    # Convert `value` to float if it's a number; otherwise, keep it as a string
    try:
        value = float(value)
    except ValueError:
        pass  # Keep value as a string if it is not a float

    return offset_type, targets, direction, value


def match_offsets_line(line):
    """Match and parse an OFFSETS line, returning the type, targets, direction, and value."""
    match = OFFSETS_LINE_PATTERN.match(line)
    if match is None:
        return None

    offset_type = match.group("type")  # Type can be 'geo', 'vib', or 'occ'
    targets = match.group("targets").strip()  # Multiple comma-separated targets
    direction = match.group("direction")  # Optional complex direction for geo
    value = match.group("value")

    # Convert `value` to float if it's a number; otherwise, keep it as a string
    try:
        value = float(value)
    except ValueError:
        pass  # Keep value as a string if it is not a float

    return offset_type, targets, direction, value
