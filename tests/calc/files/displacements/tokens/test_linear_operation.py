__authors__ = ('Alexander M. Imre (@amimre)',)
__created__ = '2025-04-11'


import numpy as np
import pytest

from viperleed_jax.files.displacements.tokens.linear_operation import (
    LinearOperationToken,
    LinearOperationTokenParserError,
)


@pytest.mark.parametrize(
    'input_str, expected_array',
    [
        # Scalars
        ('1', 1.0),
        ('0.5', 0.5),
        # 1D Vectors
        ('[0.7]', [0.7]),
        ('[0.5, 0.75]', [0.5, 0.75]),
        ('[1 2 3]', [1.0, 2.0, 3.0]),
        ('1 2 3', [1.0, 2.0, 3.0]),
        ('[1, 2, 3]', [1.0, 2.0, 3.0]),
        ('[1.0 2.5 3.2]', [1.0, 2.5, 3.2]),
        # 2D Matrices
        ('[[0.5]]', [[0.5]]),
        ('[[1 0] [0 1]]', [[1.0, 0.0], [0.0, 1.0]]),
        ('[[1, 0], [0, 1]]', [[1.0, 0.0], [0.0, 1.0]]),
        ('[[1 2 3] [3 4 5] [6 7 8]]', [[1, 2, 3], [3, 4, 5], [6, 7, 8]]),
        ('[[1.0, 2.0, 2.0], [3.0, 4.0, 4.0]]', [[1.0, 2.0, 2], [3.0, 4.0, 4]]),
        ('[[1e-3 2e-3] [3e-3 4e-3]]', [[1e-3, 2e-3], [3e-3, 4e-3]]),
        ('[[1e-3, 2e-3], [3e-3, 4e-3]]', [[1e-3, 2e-3], [3e-3, 4e-3]]),
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
