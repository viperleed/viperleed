"""Module state_recorder of viperleed.calc.classes.

Defines a class to record and retrieve calculation states.
"""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-08-06'
__license__ = 'GPLv3+'

from collections import namedtuple
from collections.abc import MutableSequence
from copy import deepcopy
from typing import Iterator

from viperleed.calc.sections.calc_section import CalcSection

# A state is given by a tuple (slab, rpars).
CalcState = namedtuple('State', 'slab, rpars, section')

class CalcStateSequence(MutableSequence):
    """A sequence of calculation states."""

    def __init__(self):
        """Initialize the sequence with an empty list of states."""
        self._recorded_states = []

    def append(self, state):
        """Add a state to the sequence."""
        self._recorded_states.append(state)

    def insert(self, index, state):
        raise ValueError('CalcStateSequence does not support item insertion.')

    def __getitem__(self, index):
        return self._recorded_states[index]

    def __setitem__(self, index, value):
        raise ValueError('CalcStateSequence does not support item assignment.')

    def __delitem__(self, index):
        raise ValueError('CalcStateSequence does not support item deletion. '
                         'Use pop() instead.')

    def pop(self, index: int = -1):
        return super().pop(index)

    def __len__(self):
        return len(self._recorded_states)

    def __iter__(self) -> Iterator:
        return iter(self._recorded_states)

    def reverse(self) -> None:
        return self._recorded_states.reverse()


class CalcStateRecorder:
    """Class that records and retrieves calculation states."""

    def __init__(self):
        """Initialize the state recorder with an empty list of states."""
        self._recorded_states = CalcStateSequence()

    def record(self, slab, rpars, section):
        """Freeze and record the current state."""
        state = CalcState(slab=deepcopy(slab),
                      rpars=deepcopy(rpars),
                      section=CalcSection(section))
        self._recorded_states.append(state)

    @property
    def last_state(self):
        """Return the last recorded state."""
        return self._recorded_states[-1]

    def get_last_section_state(self, section):
        """Return the last state recorded for a given section."""
        _section = CalcSection(section)
        for state in reversed(self._recorded_states):
            if state.section is _section:
                return state
        raise ValueError(f'No state recorded for section {section.long_name}.')

    def pop_last_state(self):
        """Remove and return the last recorded state, similar to List.pop()."""
        return self._recorded_states.pop()
