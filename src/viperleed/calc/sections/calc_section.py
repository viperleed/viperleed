"""Module calc_section of viperleed.calc.sections.

Defines an enumeration of calculation sections.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    'Alexandra Mia Imre (@alexmiame)',
)
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2023-01-23'
__license__ = 'GPLv3+'

from enum import Enum
from itertools import chain
import re

from viperleed.calc.lib.string_utils import read_int_range
from viperleed.calc.lib.string_utils import split_string_range


# Regular expressions for the aliases of each section. They are
# matched against the whole lower-case string first, then against its
# beginning. White space is not allowed, as it separates items in RUN.
# Hyphens, which separate ranges, are allowed only in place of the
# underscore before 'calc(ulation)', 'amplitudes' and 'opt(imization)'
# (as in 'DELTA-AMPLITUDES'). Longer alternatives come first, so that
# matching the beginning of a string finds the longest alias.
_ALIAS_PATTERNS = {
    'INITIALIZATION': r'init(iali[sz](e|ation))?|ini',
    'REFCALC': r'(ref(erence)?|fd|full_?dynamic)([-_]?calc(ulation)?)?',
    'RFACTOR_REFCALC': r'r_?fac(tor)?',
    'RFACTOR_SUPERPOS': r'r_?fac(tor)?_?super(pos)?',
    'DELTAS': r'deltas?([-_]?(calc(ulation)?|amplitudes?))?|del',
    'SEARCH': r'(struct(ure)?_?)?search',
    'SUPERPOS': r'sup(er(pos)?)?',
    'DOMAINS': r'domains?([-_]?calc(ulation)?)?|dom',
    'ERRORCALC': r'err(ors?)?([-_]?calc(ulation)?)?',
    'FD_OPTIMIZATION': (r'opt(imi[sz](e|ation))?'
                        r'|(fd|full_?dynamic)[-_]?opt(imi[sz](e|ation))?'),
    }
_ALIASES = {name: re.compile(pattern)
            for name, pattern in _ALIAS_PATTERNS.items()}
_PREFIX_NAME = re.compile(r'\w+')


class CalcSection(Enum):
    """An enumeration of calculation sections."""
    INITIALIZATION = 0
    REFCALC = 1
    RFACTOR_REFCALC = 11
    RFACTOR_SUPERPOS = 12
    DELTAS = 2
    SEARCH = 3
    SUPERPOS = 31
    DOMAINS = 4
    ERRORCALC = 5
    FD_OPTIMIZATION = 6

    @classmethod
    def _missing_(cls, value):  # Used only with __call__ in EnumMeta!
        """Attempt returning an element from its aliases."""
        if not isinstance(value, str):
            return None
        value = value.lower()
        # See if value is an integer
        try:
            value_int = int(value)
        except ValueError:
            pass
        else:
            try:
                return cls(value_int)
            except ValueError:
                pass

        # Try with the capitalized version of the string
        try:
            return cls[value.upper()]
        except KeyError:
            pass

        # Try with the aliases instead
        for name, alias in _ALIASES.items():
            if alias.fullmatch(value):
                return cls[name]
        return cls._from_alias_prefix(value)

    @classmethod
    def _from_alias_prefix(cls, value):
        """Return the section of the longest alias that starts `value`.

        This keeps names like 'superposition' working, which earlier
        versions accepted as they start with an alias. Only names made
        of letters, digits and underscores are accepted here: hyphens,
        colons and white space separate ranges and items in RUN.
        """
        if not _PREFIX_NAME.fullmatch(value):
            return None
        longest, best_name = 0, None
        for name, alias in _ALIASES.items():
            match = alias.match(value)
            if match and match.end() > longest:
                longest, best_name = match.end(), name
        return cls[best_name] if best_name else None

    @classmethod
    def sequence_from_string(cls, string_):
        """Return a tuple of sections from a string."""
        # See if string is a single section
        try:
            return (cls(string_),)
        except ValueError:
            pass

        # May be a range: convert it to its integer representation;
        try:
            start, stop = split_string_range(string_)
        except ValueError:
            raise ValueError(
                f'{cls.__name__}: {string_!r} is neither a valid '
                'section nor a valid range of sections.'
                ) from None
        try:
            _as_int_range = f'{cls(start).value}-{cls(stop).value}'
        except ValueError:
            raise ValueError(
                f'{cls.__name__}: Could not interpret one of '
                f'{start!r} and/or {stop!r} as a valid section'
                ) from None
        return tuple(cls(v) for v in read_int_range(_as_int_range))

    @property
    def history_tag(self):
        """Return an identifier suitable for history records."""
        _history_sections = (CalcSection.REFCALC,
                             CalcSection.DELTAS,
                             CalcSection.SEARCH)
        return self.name[0] if self in _history_sections else ''

    @property
    def long_name(self):
        """Return a descriptive name for self."""
        return _LONG_NAMES[self]

    def __lt__(self, other):  # Not sure if we need to implement the others too
        """Return whether self comes before other."""
        try:
            other = self.__class__(other)
        except ValueError:
            return NotImplemented
        return _SECTION_ORDERING[self] < _SECTION_ORDERING[other]


_LONG_NAMES = {
    CalcSection.INITIALIZATION: 'INITIALIZATION',
    CalcSection.REFCALC: 'REFERENCE CALCULATION',
    CalcSection.DELTAS: 'DELTA-AMPLITUDES',
    CalcSection.SEARCH: 'SEARCH',
    CalcSection.RFACTOR_REFCALC: 'R-FACTOR CALCULATION',
    CalcSection.RFACTOR_SUPERPOS: 'R-FACTOR CALCULATION',
    CalcSection.SUPERPOS: 'SUPERPOS',
    CalcSection.ERRORCALC: 'ERROR CALCULATION',
    CalcSection.FD_OPTIMIZATION: 'FULL-DYNAMIC OPTIMIZATION',
    }


_SECTION_ORDERING = {  # Default order of execution
    v: i
    for i, v in enumerate((
        CalcSection.INITIALIZATION,
        CalcSection.REFCALC,
        CalcSection.FD_OPTIMIZATION,
        CalcSection.RFACTOR_REFCALC,
        CalcSection.DELTAS,
        CalcSection.SEARCH,
        CalcSection.SUPERPOS,
        CalcSection.RFACTOR_SUPERPOS,
        CalcSection.DOMAINS,
        CalcSection.ERRORCALC
        ))
    }


_REQUIRED_FILES = {  # Required input files per section
    CalcSection.INITIALIZATION: (
        'IVBEAMS',
        'PARAMETERS',
        'POSCAR',
        'VIBROCC',
        ),
    CalcSection.REFCALC: (
        'BEAMLIST',
        'IVBEAMS',
        'PARAMETERS',
        'PHASESHIFTS',
        'POSCAR',
        'VIBROCC',
        ),
    CalcSection.DELTAS: (
        'BEAMLIST',
        'DISPLACEMENTS',
        'IVBEAMS',
        'PARAMETERS',
        'PHASESHIFTS',
        'POSCAR',
        'VIBROCC',
        ),
    CalcSection.SEARCH: (
        'BEAMLIST',
        'DISPLACEMENTS',
        'EXPBEAMS',
        'IVBEAMS',
        'PARAMETERS',
        'PHASESHIFTS',
        'POSCAR',
        'VIBROCC',
        ),
    CalcSection.RFACTOR_REFCALC: (
        'BEAMLIST',
        'EXPBEAMS',
        'IVBEAMS',
        'PARAMETERS',
        ),
    CalcSection.RFACTOR_SUPERPOS: (
        'BEAMLIST',
        'PARAMETERS',
        'IVBEAMS',
        'EXPBEAMS',
        ),
    CalcSection.SUPERPOS: (
        'BEAMLIST',
        'DISPLACEMENTS',
        'IVBEAMS',
        'PARAMETERS',
        'POSCAR',
        'VIBROCC',
        ),
    CalcSection.ERRORCALC: (
        'BEAMLIST',
        'DISPLACEMENTS',
        'EXPBEAMS',
        'IVBEAMS',
        'PARAMETERS',
        'PHASESHIFTS',
        'POSCAR',
        'VIBROCC',
        ),
    CalcSection.FD_OPTIMIZATION: (
        'BEAMLIST',
        'EXPBEAMS',
        'IVBEAMS',
        'PARAMETERS',
        'PHASESHIFTS',
        'POSCAR',
        'VIBROCC',
        ),
    }

# set of all input files
ALL_INPUT_FILES = set(chain.from_iterable(_REQUIRED_FILES.values()))

# allowed names for the file containing the experimental beams
# files will be used in precedence from left to right
EXPBEAMS_NAMES = ('EXPBEAMS.csv', 'EXPBEAMS')
ALL_INPUT_FILES.update(set(EXPBEAMS_NAMES))
