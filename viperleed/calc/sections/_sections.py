"""Module _sections of viperleed.calc.sections.

Defines an enumeration of calculation sections.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    'Alexander M. Imre (@amimre)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-01-23'
__license__ = 'GPLv3+'

from enum import Enum
from itertools import chain

from viperleed.calc.lib.base import split_string_range, readIntRange

_ALIASES = {  # Exact match first, then check .startswith
    'INITIALIZATION': ('ini', 'init', 'initialisation'),
    'REFCALC': ('ref', 'fd'),
    'RFACTOR_REFCALC': ('rfac', 'rfactor'),
    'RFACTOR_SUPERPOS': ('rfacsuper',),
    'DELTAS': ('del', 'delta', 'deltaamplitudes'),
    'SEARCH': (),
    'SUPERPOS': ('sup', 'super'),
    'DOMAINS': ('dom', 'domain'),
    'ERRORCALC': ('err', 'error'),
    'FD_OPTIMIZATION': ('opt', 'optimize', 'fdopt'),
    }


class TLEEDMSection(Enum):
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

        # Try with the aliases instead. First look for an exact match
        for name, aliases in _ALIASES.items():
            if any(alias == value for alias in aliases):
                return cls[name]
        # Alternatively, see if the value begins with an alias
        for name, aliases in _ALIASES.items():
            if any(value.startswith(alias) for alias in aliases):
                return cls[name]
        return None

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
                f"{cls.__name__}: {string_!r} is neither a valid "
                "section nor a valid range of sections."
                ) from None
        try:
            _as_int_range = f"{cls(start).value}-{cls(stop).value}"
        except ValueError:
            raise ValueError(
                f"{cls.__name__}: Could not interpret one of "
                f"{start!r} and/or {stop!r} as a valid section"
                ) from None
        return tuple(cls(v) for v in readIntRange(_as_int_range))

    @property
    def history_tag(self):
        """Return an identifier suitable for history records."""
        _history_sections = (TLEEDMSection.REFCALC,
                             TLEEDMSection.DELTAS,
                             TLEEDMSection.SEARCH)
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
    TLEEDMSection.INITIALIZATION: "INITIALIZATION",
    TLEEDMSection.REFCALC: "REFERENCE CALCULATION",
    TLEEDMSection.DELTAS: "DELTA-AMPLITUDES",
    TLEEDMSection.SEARCH: "SEARCH",
    TLEEDMSection.RFACTOR_REFCALC: "R-FACTOR CALCULATION",
    TLEEDMSection.RFACTOR_SUPERPOS: "R-FACTOR CALCULATION",
    TLEEDMSection.SUPERPOS: "SUPERPOS",
    TLEEDMSection.ERRORCALC: "ERROR CALCULATION",
    TLEEDMSection.FD_OPTIMIZATION: "FULL-DYNAMIC OPTIMIZATION",
    }


_SECTION_ORDERING = {  # Default order of execution
    v: i
    for i, v in enumerate((
        TLEEDMSection.INITIALIZATION,
        TLEEDMSection.REFCALC,
        TLEEDMSection.FD_OPTIMIZATION,
        TLEEDMSection.RFACTOR_REFCALC,
        TLEEDMSection.DELTAS,
        TLEEDMSection.SEARCH,
        TLEEDMSection.SUPERPOS,
        TLEEDMSection.RFACTOR_SUPERPOS,
        TLEEDMSection.DOMAINS,
        TLEEDMSection.ERRORCALC
        ))
    }

_REQUIRED_FILES = {  # Required input files per section
    TLEEDMSection.INITIALIZATION: ["POSCAR", "PARAMETERS", "VIBROCC", "IVBEAMS"],
    TLEEDMSection.REFCALC: ["BEAMLIST", "PHASESHIFTS", "POSCAR", "PARAMETERS",
                         "IVBEAMS", "VIBROCC"],
    TLEEDMSection.DELTAS: ["BEAMLIST", "PHASESHIFTS", "POSCAR", "PARAMETERS",
                         "IVBEAMS", "VIBROCC", "DISPLACEMENTS"],
    TLEEDMSection.SEARCH: ["BEAMLIST", "PHASESHIFTS", "POSCAR", "PARAMETERS",
                         "IVBEAMS", "VIBROCC", "DISPLACEMENTS", "EXPBEAMS"],
    TLEEDMSection.RFACTOR_REFCALC: ["BEAMLIST", "PARAMETERS", "IVBEAMS", "EXPBEAMS"],
    TLEEDMSection.RFACTOR_SUPERPOS: ["BEAMLIST", "PARAMETERS", "IVBEAMS", "EXPBEAMS"],
    TLEEDMSection.SUPERPOS: ["BEAMLIST", "POSCAR", "PARAMETERS", "IVBEAMS",
                          "VIBROCC", "DISPLACEMENTS"],
    TLEEDMSection.ERRORCALC: ["BEAMLIST", "PHASESHIFTS", "POSCAR", "PARAMETERS",
                         "IVBEAMS", "VIBROCC", "DISPLACEMENTS", "EXPBEAMS"],
    TLEEDMSection.FD_OPTIMIZATION: ["BEAMLIST", "PHASESHIFTS", "POSCAR", "PARAMETERS",
                         "IVBEAMS", "VIBROCC", "EXPBEAMS"],
}

# set of all input files
ALL_INPUT_FILES = set(chain.from_iterable(_REQUIRED_FILES.values()))

# allowed names for the file containing the experimental beams
# files will be used in precedence from left to right
EXPBEAMS_NAMES = ("EXPBEAMS.csv", "EXPBEAMS")
ALL_INPUT_FILES.update(set(EXPBEAMS_NAMES))