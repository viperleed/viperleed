"""Tests for module viperleed.calc.classes.state_recorder."""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-10-04'
__license__ = 'GPLv3+'

from collections.abc import MutableSequence
from collections.abc import Sequence

import pytest
from pytest_cases import parametrize

from viperleed.calc.sections.calc_section import CalcSection
from viperleed.calc.classes.state_recorder import CalcState
from viperleed.calc.classes.state_recorder import CalcStateRecorder
from viperleed.calc.classes.state_recorder import NoStateError
from viperleed.calc.classes.state_recorder import _StateSequence


@pytest.fixture(name='calc_state')
def fixture_calc_state(ag100):
    """Return a CalcState for a Ag(100) Slab and Rparams."""
    slab, rpars, *_ = ag100
    return CalcState(slab, rpars, CalcSection.INITIALIZATION)


@pytest.fixture(name='state_sequence')
def fixture_state_sequence(calc_state):
    """Return an instance of a concrete subclass of _StateSequence."""
    # We can't directly instantiate a _StateSequence as it's abstract
    class _ConcreteSequence(_StateSequence):
        def _make_new_state(self, *args, **__):
            return args[0] if args else calc_state
        def _revert_last_state(self, *_, **__):
            pass
    return _ConcreteSequence()


@pytest.fixture(name='calc_state_recorder')
def fixture_calc_state_recorder():
    """Return an empty CalcStateRecorder instance."""
    return CalcStateRecorder()


class TestStateSequence:
    """Tests for _StateSequence."""

    def test_is_sequence_like(self, state_sequence):
        """Check that state_sequence has the full Sequence interface."""
        assert not set(dir(Sequence)) - set(dir(state_sequence))
        assert not isinstance(state_sequence, Sequence)

    _mutable_seq_methods = set(dir(MutableSequence)) - set(dir(Sequence))

    @parametrize(method=_mutable_seq_methods)
    def test_not_a_mutable_sequence(self, state_sequence, method):
        """Check that no MutableSequence methods are available."""
        with pytest.raises(AttributeError):
            getattr(state_sequence, method)

    def test_setitem_raises(self, state_sequence):
        """Check that items can't be assigned to."""
        with pytest.raises(TypeError):
            state_sequence[0] = 1

    def test_delitem_raises(self, state_sequence):
        """Check that items can't be deleted."""
        with pytest.raises(TypeError):
            del state_sequence[0]

    def test_last_state(self, state_sequence):
        """Check correct identification of the last state."""
        states = 123, 456, 789
        for state in states:
            state_sequence.record(state)
        assert state_sequence.last_state == states[-1]

    def test_last_state_raises(self, state_sequence):
        """Check complaints when pulling the last state from an empty record."""
        with pytest.raises(NoStateError):
            _ = state_sequence.last_state

    def test_len_and_bool(self, state_sequence):
        """Check expected outcome of __len__ and __bool__."""
        assert not state_sequence
        state_sequence.record(1, 2, 3)
        assert len(state_sequence) == 1
        assert state_sequence

    def test_iter(self, state_sequence):
        """Check correct iteration."""
        states = 'abc', '123'
        for state in states:
            state_sequence.record(state)
        iter_seq = iter(state_sequence)
        assert tuple(iter_seq) == states

    def test_record(self, state_sequence, calc_state):
        """Check that recording a state works correctly."""
        state_sequence.record()
        assert len(state_sequence) == 1
        assert state_sequence[0] == calc_state

    def test_reversed(self, state_sequence, calc_state):
        """Check iteration in reversed order."""
        first_state = 'abc'
        state_sequence.record(first_state)
        state_sequence.record()
        iter_rev = reversed(state_sequence)
        assert next(iter_rev) == calc_state
        assert next(iter_rev) == first_state

    def test_revert(self, state_sequence):
        """Check that reverting removes recorded states."""
        state_sequence.record(1)
        state_sequence.revert()
        assert not state_sequence

    def test_revert_raises(self, state_sequence):
        """Check complaints when attempting to reverting an empty record."""
        with pytest.raises(NoStateError):
            state_sequence.revert()


class TestCalcStateRecorder:
    """Tests for CalcStateRecorder"""

    @staticmethod
    def _check_state(recorded_state, slab, _, section):
        """Check that recorded_state matches expectations."""
        assert CalcSection(recorded_state.section.long_name) is section
        assert recorded_state.section is section
        assert recorded_state.slab.elements == slab.elements
        assert recorded_state.rpars is not None

    def test_last_state(self, calc_state_recorder, calc_state):
        """Check the last state is correctly returned."""
        calc_state_recorder.record(*calc_state)
        last_state = calc_state_recorder.last_state
        self._check_state(last_state, *calc_state)

    def test_get_last_section_state(self, calc_state_recorder, calc_state):
        """Check that the last state stored for a section is correct."""
        *_, section = calc_state
        calc_recorder = calc_state_recorder
        calc_recorder.record(*calc_state)
        last_section_state = calc_recorder.get_last_section_state(section)
        self._check_state(last_section_state, *calc_state)
        assert last_section_state is calc_recorder[-1]

    def test_get_last_section_state_raises(self, calc_state_recorder):
        """Check complaints when accessing state for an unrecorded section."""
        unused_section = CalcSection.REFCALC
        with pytest.raises(ValueError, match='No state recorded for section '
                           f'{unused_section.long_name}.'):
            calc_state_recorder.get_last_section_state(unused_section)

    def test_record(self, calc_state_recorder, calc_state):
        """Check correct recording of a `calc_state`."""
        calc_state_recorder.record(*calc_state)
        assert len(calc_state_recorder) == 1
        recorded_state = calc_state_recorder[0]
        self._check_state(recorded_state, *calc_state)

    @parametrize(n_args=(0, 1, 2, 4, 5, 6))
    def test_record_args(self, calc_state_recorder, n_args):
        """Check complaints with too few/many positional .record arguments."""
        args = (1,) * n_args
        with pytest.raises(TypeError, match='slab|rpars|section'):
            calc_state_recorder.record(*args)

    @parametrize(n_kwargs=(1, 2, 3, 4, 5, 6))
    def test_record_kwargs(self, calc_state_recorder, n_kwargs):
        """Check complaints with too many keyword .record arguments."""
        kwargs = {str(i): i for i in range(n_kwargs)}
        with pytest.raises(TypeError, match='are acceptable'):
            calc_state_recorder.record(1, 2, 3, **kwargs)

    def test_revert(self, calc_state_recorder):
        """Check that reverting a state does what it should."""
        calc_state_recorder.record(1, 2, 3)
        with pytest.raises(NotImplementedError):
            calc_state_recorder.revert(1, 2)

    @parametrize(n_args=(0, 1, 3, 4, 5, 6))
    def test_revert_args(self, calc_state_recorder, n_args):
        """Check complaints with too few/many positional .revert arguments."""
        args = (1,) * n_args
        calc_state_recorder.record(1, 2, 3)
        with pytest.raises(TypeError, match='slab|rpars|section'):
            calc_state_recorder.revert(*args)

    @parametrize(n_kwargs=(1, 2, 3, 4, 5, 6))
    def test_revert_kwargs(self, calc_state_recorder, n_kwargs):
        """Check complaints with too few/many keyword .revert arguments."""
        kwargs = {str(i): i for i in range(n_kwargs)}
        calc_state_recorder.record(1, 2, 3)
        with pytest.raises(TypeError, match='are acceptable'):
            calc_state_recorder.revert(1, 2, **kwargs)
