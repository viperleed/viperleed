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

    def __getitem__(self, index):
        """Return the state at `index`."""
        return self._recorded_states[index]

    def __iter__(self) -> Iterator:
        """Return an iterator of recorded states."""
        return iter(self._recorded_states)

    def __len__(self):
        """Return the number of recorded states."""
        return len(self._recorded_states)

    def __setitem__(self, index, value):
        """Disallow modifying the state at `index`."""
        raise ValueError('CalcStateSequence does not support item assignment.')

    def __delitem__(self, index):
        """Disallow deleting the state at `index`."""
        raise ValueError('CalcStateSequence does not support item deletion. '
                         'Use pop() instead.')

    def append(self, state):
        """Add a state to the sequence."""
        self._recorded_states.append(state)

    def insert(self, index, state):
        """Disallow inserting states."""
        raise ValueError('CalcStateSequence does not support item insertion.')

    def pop(self, index: int = -1):
        """Remove the last recorded state and return it."""
        return self._recorded_states.pop(index)

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
