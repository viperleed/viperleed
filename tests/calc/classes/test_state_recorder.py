"""Tests for module viperleed.calc.classes.state_recorder."""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-10-04'
__license__ = 'GPLv3+'

import pytest
from copy import deepcopy
from collections import namedtuple

# Import the module and classes
from viperleed.calc.classes.state_recorder import CalcStateSequence
from viperleed.calc.classes.state_recorder import CalcStateRecorder, CalcState
from viperleed.calc.sections.calc_section import CalcSection


# Fixtures for test data
@pytest.fixture
def calc_state(ag100):
    slab, rpars, _ = ag100
    return CalcState(slab, rpars, CalcSection('0'))

@pytest.fixture
def calc_state_sequence():
    return CalcStateSequence()

@pytest.fixture
def calc_state_recorder():
    return CalcStateRecorder()

class TestCalcStateSequence:
    """Tests for CalcStateSequence"""
    def test_append(self, calc_state_sequence, calc_state):
        calc_state_sequence.append(calc_state)
        assert len(calc_state_sequence) == 1
        assert calc_state_sequence[0] == calc_state

    def test_insert_raises(self, calc_state_sequence, calc_state):
        with pytest.raises(ValueError, match='CalcStateSequence does not support item insertion.'):
            calc_state_sequence.insert(0, calc_state)

    def test_setitem_raises(self, calc_state_sequence, calc_state):
        calc_state_sequence.append(calc_state)
        with pytest.raises(ValueError, match='CalcStateSequence does not support item assignment.'):
            calc_state_sequence[0] = calc_state

    def test_delitem_raises(self, calc_state_sequence):
        with pytest.raises(ValueError, match='CalcStateSequence does not support item deletion.'):
            del calc_state_sequence[0]

    def test_pop(self, calc_state_sequence, calc_state):
        calc_state_sequence.append(calc_state)
        with pytest.raises(ValueError, match='CalcStateSequence does not support item deletion.'):
            calc_state_sequence.pop()

    def test_len(self, calc_state_sequence, calc_state):
        assert len(calc_state_sequence) == 0
        calc_state_sequence.append(calc_state)
        assert len(calc_state_sequence) == 1

    def test_iter(self, calc_state_sequence, calc_state):
        calc_state_sequence.append(calc_state)
        iter_seq = iter(calc_state_sequence)
        assert next(iter_seq) == calc_state

    def test_reverse(self, calc_state_sequence, calc_state):
        calc_state_sequence.append(calc_state)
        calc_state_sequence.append(calc_state)
        calc_state_sequence.reverse()
        assert calc_state_sequence[0] == calc_state

class TestCalcStateRecorder:
    """Tests for CalcStateRecorder"""
    def test_record(self, calc_state_recorder, calc_state):
        slab, rpars, section = calc_state
        calc_state_recorder.record(slab, rpars, section)
        assert len(calc_state_recorder._recorded_states) == 1
        recorded_state = calc_state_recorder._recorded_states[0]
        assert CalcSection(recorded_state.section.long_name) == section
        assert recorded_state.slab.elements == slab.elements
        assert recorded_state.rpars is not None


    def test_last_state(self, calc_state_recorder, calc_state):
        slab, rpars, section = calc_state
        calc_state_recorder.record(slab, rpars, section)
        last_state = calc_state_recorder.last_state
        assert CalcSection(last_state.section.long_name) == section
        assert last_state.slab.elements == slab.elements
        assert last_state.rpars is not None

    def test_get_last_section_state(self, calc_state_recorder, calc_state):
        slab, rpars, section = calc_state
        calc_state_recorder.record(slab, rpars, section)
        last_section_state = calc_state_recorder.get_last_section_state(section)
        assert CalcSection(last_section_state.section.long_name) == section
        assert last_section_state.slab.elements == slab.elements
        assert last_section_state.rpars is not None

    def test_get_last_section_state_raises(self, calc_state_recorder, calc_state):
        unused_section = CalcSection('1')
        with pytest.raises(ValueError, match='No state recorded for section '
                           f'{unused_section.long_name}.'):
            calc_state_recorder.get_last_section_state(unused_section)

    def test_pop_last_state(self, calc_state_recorder, calc_state):
        slab, rpars, section = calc_state
        calc_state_recorder.record(slab, rpars, section)
        popped_state = calc_state_recorder.pop_last_state()
        assert CalcSection(popped_state.section.long_name) == section
        assert popped_state.slab.elements == slab.elements
        assert popped_state.rpars is not None
        assert len(calc_state_recorder._recorded_states) == 0

