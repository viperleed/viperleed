"""Module state_recorder of viperleed.calc.classes.

Defines a class to record and retrieve calculation states.
"""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-08-06'
__license__ = 'GPLv3+'

from abc import ABCMeta
from abc import abstractmethod
from collections import namedtuple
from copy import deepcopy

from viperleed.calc.sections.calc_section import CalcSection


class StateError(Exception):
    """Base class for exceptions concerning states."""


class NoStateError(StateError):
    """No states available."""


CalcState = namedtuple('State', 'slab, rpars, section')


# Do not use Sequence as a base class to avoid risking giving the
# impression that this is an immutable object, as Sequence usually
# is reserved for immutable ones, like tuple or str. However, the
# full set of methods available for a Sequence is implemented
# below.

class _StateSequence(metaclass=ABCMeta):
    """A sequence of states.

    Users can only .record() or .revert() to add/remove states.
    """

    __slots__ = {
        # Prevent adding new attributes at runtime. See
        # also https://stackoverflow.com/questions/472000
        '_recorded_states': 'The internal list used for storage',
        }

    # Tell ABCMeta.__new__ that this class has TPFLAGS_SEQUENCE set.
    __abc_tpflags__ = 1 << 5 # Py_TPFLAGS_SEQUENCE

    def __init__(self):
        """Initialize the sequence with an empty list of states."""
        self._recorded_states = []

    @property
    def last_state(self):
        """Return the last recorded state."""
        if not self:
            raise NoStateError('No states recorded yet.')
        return self[-1]

    def __bool__(self):
        """Return whether there is any state."""
        return bool(self._recorded_states)

    def __contains__(self, state):
        """Return whether `state` is in this sequence."""
        return state in self._recorded_states

    def __getitem__(self, index):
        """Return the state at `index`."""
        return self._recorded_states[index]

    def __iter__(self):
        """Return an iterator of recorded states."""
        return iter(self._recorded_states)

    def __len__(self):
        """Return the number of recorded states."""
        return len(self._recorded_states)

    def __reversed__(self):
        """Return an iterator of recorded states in reverse order."""
        return reversed(self._recorded_states)

    def count(self, state):
        """Return then number of occurrences of state in self."""
        return self._recorded_states.count(state)

    def index(self, state, start=0, stop=None):
        """Return the index of state in this sequence."""
        return self._recorded_states.index(state, start=start, stop=stop)

    def record(self, *args, **kwargs):
        """Freeze and record the current state."""
        state = self._make_new_state(*args, **kwargs)
        self._recorded_states.append(state)

    def revert(self, *args, **kwargs):
        """Revert to the state before the most recent one."""
        if not self:
            raise NoStateError('No states recorded yet.')
        self._revert_last_state(*args, **kwargs)
        self._recorded_states.pop()

    def _check_unexpected_args(self, accepted, *args, **kwargs):
        """Raise if any argument is passed.

        Parameters
        ----------
        accepted : str
            Which arguments are acceptable. Used only for the
            exception message.
        *args : object
            Unexpected positional arguments.
        **kwargs : object
            Unexpected keyword arguments.

        Raises
        ------
        TypeError
            If any `args` or `kwargs` are given.
        """
        if args or kwargs:
            raise TypeError(f'{type(self).__name__!r}: Only '
                            f'{accepted} argument(s) are acceptable.')

    @abstractmethod
    def _make_new_state(self, *args, **kwargs):
        """Return a new state for this _StateSequence."""

    @abstractmethod
    def _revert_last_state(self, *args, **kwargs):
        """Revert to the state before the most recent one."""


class CalcStateRecorder(_StateSequence):
    """Class that records and retrieves calculation states."""

    __slots__ = ()  # Prevent adding new attributes at runtime.

    def get_last_state_for_section(self, section):
        """Return the last state recorded for a given section."""
        section = CalcSection(section)
        try:
            return next(state
                        for state in reversed(self)
                        if state.section is section)
        except StopIteration:
            pass
        raise ValueError(f'No state recorded for section {section.long_name}.')

    # We're just providing a different signature for this concrete
    # implementation. We also check that there are no extra arguments
    # passed, so it's OK to have a different number of arguments. The
    # others will not cause a TypeError from python, but a more
    # understandable one, via _check_unexpected_args.
    # pylint: disable-next=arguments-differ
    def record(self, slab, rpars, section, *args, **kwargs):
        """Freeze and record the current state."""
        self._check_unexpected_args('slab, rpars, and section',
                                    *args,
                                    **kwargs)
        super().record(slab, rpars, section)

    # pylint: disable-next=arguments-differ              # Like .record
    def _make_new_state(self, slab, rpars, section, *args, **kwargs):
        """Return a new CalcState."""
        self._check_unexpected_args('slab, rpars, and section',
                                    *args,
                                    **kwargs)
        return CalcState(slab=deepcopy(slab),
                         rpars=deepcopy(rpars),
                         section=CalcSection(section))

    # pylint: disable-next=arguments-differ              # Like .record
    def _revert_last_state(self, slab, rpars, *args, **kwargs):
        """Revert to the state before the most recent one."""
        self._check_unexpected_args('slab and rpars', *args, **kwargs)
        raise NotImplementedError
