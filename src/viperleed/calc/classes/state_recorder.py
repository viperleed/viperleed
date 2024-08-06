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
from copy import deepcopy

from viperleed.calc.sections.calc_section import CalcSection

# A state is given by a tuple (slab, rparams).
State = namedtuple('State', 'slab, rparams, section')

class StateRecorder:
    """Class that records and retrieves calculation states."""

    def __init__(self):
        """Initialize the state recorder with an empty list of states."""
        self.recorded_states = []

    def record(self, slab, rparams, section):
        """Freezes and records the current state."""
        state = State(slab=deepcopy(slab),
                      rparams=deepcopy(rparams),
                      section=CalcSection(section))
        self.recorded_states.append(state)

    def get_last_state(self):
        """Returns the last recorded state."""
        return self.recorded_states[-1]

    def get_last_section_state(self, section):
        """Returns the last state recorded for a given section."""
        _section = CalcSection(section)
        for state in reversed(self.recorded_states):
            if state.rparams.section == _section:
                return state
        raise ValueError(f"No state recorded for section {section.long_name}.")

    def pop_last_state(self, index=None):
        """Removes and returns the last recorded state, similar to List.pop()."""
        return self.recorded_states.pop(index if index is not None else -1)
