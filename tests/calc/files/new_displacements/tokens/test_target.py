"""Tests for module viperleed.calc.files.new_displacements.tokens.target."""

__authors__ = ('Alexander M. Imre (@amimre)',)
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2024-10-14'
__license__ = 'GPLv3+'

import re

import pytest

from viperleed.calc.files.new_displacements.tokens.target import (
    TargetingError,
    TargetToken,
    _generate_label_match_regex,
)


@pytest.mark.parametrize(
    'raw, expected_regex, expected_nums, expected_layers',
    [
        # only label, no nums or layers
        ('A', re.compile(r'^A'), None, None),
        # list of numbers
        ('A 1 3 5', re.compile(r'^A'), [1, 3, 5], None),
        # numeric range
        ('A 2-4', re.compile(r'^A'), [2, 3, 4], None),
        # numeric range
        ('A 1-3 5 7', re.compile(r'^A'), [1, 2, 3, 5, 7], None),
        # single layer
        ('A L(1)', re.compile(r'^A'), None, [1]),
        # layer range
        ('A L(2-3)', re.compile(r'^A'), None, [2, 3]),
        # layer with spaces
        ('A L( 1-3 )', re.compile(r'^A'), None, [1, 2, 3]),
        # complicated layer numbers
        ('A L(1 4-5 7)', re.compile(r'^A'), None, [1, 4, 5, 7]),
        # wildcard label
        ('* 1-2', re.compile(r'^\w*'), [1, 2], None),
    ],
)
def test_parse_valid_targets(
    raw, expected_regex, expected_nums, expected_layers
):
    tok = TargetToken(raw)
    # regex pattern
    assert isinstance(tok.regex, re.Pattern)
    assert tok.regex.pattern == expected_regex.pattern
    # nums and layers
    assert tok.nums == expected_nums
    assert tok.layers == expected_layers


@pytest.mark.parametrize(
    'raw',
    [
        # combined layer and number
        'A L(1-2) 3',
        # numeric without label
        '1 2 3',
        '2',
        '[1 0 0'
        # layer without label
        'L(1)',
        # wrongly formatted layer
    ],
)
def test_invalid_target_raises(raw):
    with pytest.raises(TargetingError):
        TargetToken(raw)


def test_empty_string_raises():
    with pytest.raises(TargetingError, match='Target string is empty'):
        TargetToken('')


def test_repr():
    token = TargetToken('A')
    repr_str = str(token)
    assert 'TargetToken(target_str=A)' in repr_str


def test_eq_same_and_different():
    a = TargetToken('A 1-2')
    b = TargetToken('A 1 2')
    c = TargetToken('A 1-2')
    # a and c should be equal (same raw string)
    assert a == c
    # a and b: same nums but b parsed as nums=[1,2] as well => they compare regex equal & nums equal
    assert a == b
    # different label
    d = TargetToken('B 1-2')
    assert a != d
    # comparing to different type
    assert a != 'A 1-2'


def test_eq_layer_specifier():
    # regression test for bug
    a = TargetToken('Fe L(1-3)')
    b = TargetToken('Fe L(1-3)')
    assert a == b


def test_selection_regex_matching():
    # test that regex actually matches sites correctly
    tok = TargetToken('A*')
    # should match A, A1, ABCXYZ, but not B
    regex = tok.regex
    assert regex.match('A')
    assert regex.match('A1')
    assert regex.match('ABCXYZ')
    assert not regex.match('B')


def test_invalid_number_format_raises():
    with pytest.raises(TargetingError):
        # non-integer in numeric list
        TargetToken('A one two')
    with pytest.raises(TargetingError):
        # malformed layer spec
        TargetToken('A L()')


def test_generate_label_match_regex_wildcards():
    # Directly test helper
    r = _generate_label_match_regex('X*Y')
    # Should match XANY, XfooYbar, etc.
    assert r.match('XANY')
    assert r.match('XfooYbar')
    assert r.match('XYX')
    assert not r.match('Y')
