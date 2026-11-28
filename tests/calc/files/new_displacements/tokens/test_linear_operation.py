"""Tests for module viperleed.calc.files.new_displacements.tokens.linear_operation."""

__authors__ = ('Alexandra Mia Imre (@alexmiame)',)
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2025-04-11'
__license__ = 'GPLv3+'


import numpy as np
import pytest
from pytest_cases import parametrize

from viperleed.calc.files.new_displacements.tokens.linear_operation import (
    LinearOperationToken,
    LinearOperationTokenParserError,
)


@pytest.mark.parametrize(
    'input_str, expected_array',
    [
        # Scalars
        ('1', 1.0),
        ('0.5', 0.5),
        # (1x1)
        ('[0.7]', [[0.7]]),
        ('[[ 0.5] ]', [[0.5]]),
        # (2x2)
        ('[[1 0] [0 1]]', [[1.0, 0.0], [0.0, 1.0]]),
        ('[[1, 0], [0, 1]]', [[1.0, 0.0], [0.0, 1.0]]),
        ('[[1e-3, 2e-3], [3e-3, 4e-3]]', [[1e-3, 2e-3], [3e-3, 4e-3]]),
        ('[[1 0] [0 -1]]', [[1, 0], [0, -1]]),
        ('[[1 0] [0 .5]]', [[1, 0], [0, 0.5]]),
        # (3x3)
        ('[[1 2 3] [3 4 5] [6 7 8]]', [[1, 2, 3], [3, 4, 5], [6, 7, 8]]),
    ],
)
def test_linear_operation_token_parsing(input_str, expected_array):
    token = LinearOperationToken(input_str)
    np.testing.assert_allclose(token.arr, expected_array)


@pytest.mark.parametrize(
    'invalid_input',
    [
        '',  # empty string
        '[[1 2], [3]]',  # ragged array
        '[1, 2 3]',  # mixed comma + space
        '[[1, 2], [3, four]]',  # non-numeric
        '[[1, 2], 3]',  # mixed nesting levels
        '[1 2',  # unclosed bracket
        '[[1 2] [3 4]]]',  # extra closing bracket
    ],
)
def test_linear_operation_token_invalid_input(invalid_input):
    with pytest.raises(LinearOperationTokenParserError):
        LinearOperationToken(invalid_input)


@parametrize(bad=([1, 2, 3], np.eye(4), [[1, 2], [3]]))
def test_from_array_invalid(bad):
    with pytest.raises(LinearOperationTokenParserError):
        LinearOperationToken.from_array(bad)


def test_from_array_reshapes_scalar():
    token = LinearOperationToken.from_array(0.5)
    assert token.arr.shape == (1, 1)
