from enum import Enum
from collections import namedtuple
import re

from .errors import InvalidSyntaxError
from .errors import SymmetryViolationError

from .regex import SEARCH_HEADER_PATTERN
from .regex import SECTION_HEADER_PATTERN
from .regex import match_geo_line
from .regex import match_vib_line
from .regex import match_occ_line
from .regex import match_constrain_line

DisplacementFileSections = Enum('DisplacementFileSections', [
    'GEO_DELTA',
    'VIB_DELTA',
    'OCC_DELTA',
    'CONSTRAIN'
])


LOOP_START_MARKER = 'LOOP_START'
LOOP_END_MARKER = 'LOOP_END'

LoopMarkerLine = namedtuple('LoopMarkerLine', ['type'])
SearchHeaderLine = namedtuple('SearchHeaderLine', ['label'])
SectionLine = namedtuple('SectionLine', ['section', 'line'])
