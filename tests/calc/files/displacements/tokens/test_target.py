# test_target_token.py
import re

import pytest

from viperleed_jax.files.displacements.tokens.target import (
    TargetingError,
    TargetToken,
    _generate_label_match_regex,
)


@pytest.mark.parametrize(
    'raw, expected_regex, expected_nums, expected_layers',
    [
        # only label, no nums or layers
        ('A', re.compile(r'^A\w*'), None, None),
        # list of numbers
        ('B 1 3 5', re.compile(r'^B\w*'), [1, 3, 5], None),
        # numeric range
        ('C 2-4', re.compile(r'^C\w*'), [2, 3, 4], None),
        # single layer
        ('D L(1)', re.compile(r'^D\w*'), None, [1]),
        # layer range
        ('E L(2-3)', re.compile(r'^E\w*'), None, [2, 3]),
        # wildcard label
        ('* 1-2', re.compile(r'^\w*\w*'), [1, 2], None),
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
        # numberic without label
        '1 2 3',
        '2',
        '[1 0 0'
        # layer without label
        'L(1)',
        # wrongly formatted layer
    ])
def test_invalid_target_raises(raw):
    with pytest.raises(TargetingError) as exc:
        TargetToken(raw)


def test_empty_string_raises():
    with pytest.raises(TargetingError) as exc:
        TargetToken('')
    assert 'Target string is empty' in str(exc.value)


def test_repr():
    tok = TargetToken('A')
    repr_str = repr(tok)
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
    assert not (a == d)
    # comparing to different type
    assert not (a == 'A 1-2')


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
    with pytest.raises(ValueError):
        # non-integer in numeric list
        TargetToken('A one two')
    with pytest.raises(ValueError):
        # malformed layer spec
        TargetToken('A L()')


def test_generate_label_match_regex_wildcards():
    # Directly test helper
    r = _generate_label_match_regex('X*Y')
    # Should match XANY, XfooYbar, etc.
    assert r.match('XANY')
    assert r.match('XfooYbar')
    assert not r.match('XYX') == None  # ensure it matches "XYX"
    assert not r.match('Y')
