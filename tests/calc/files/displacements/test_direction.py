import numpy as np
import pytest

from viperleed_jax.files.displacements.direction import Direction


# Test the Direction class
@pytest.mark.parametrize(
    'direction_str, expected_vectors, ucell, fractional, num_free_directions',
    [
        (
            'xyz',
            np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]]),
            np.eye(3),
            False,
            3,
        ),
        ('xy[1 0]', np.array([[1, 0, 0]]), np.eye(3), False, 1),
        ('x', np.array([[1, 0, 0]]), np.eye(3), False, 1),
        ('xy[0 1]', np.array([[0, 1, 0]]), np.eye(3), False, 1),
        (
            'abc',
            np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
            np.array([[3.5, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 5.0]]),
            True,
            3,
        ),
        (
            'ab[1 0]',
            np.array([[1.0, 0.0, 0.0]]),
            np.array([[3.5, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 5.0]]),
            True,
            1,
        ),
    ],
)
def test_direction(
    direction_str, expected_vectors, ucell, fractional, num_free_directions
):
    """Test different cases of Cartesian and fractional directions."""
    d = Direction(direction_str)
    assert d._fractional == fractional
    assert d.num_free_directions == num_free_directions
    np.testing.assert_allclose(d.get_cart_vectors(ucell), expected_vectors)
    if fractional:
        np.testing.assert_allclose(d.get_frac_vectors(ucell), expected_vectors)


@pytest.mark.parametrize(
    'invalid_direction',
    [
        'xab',  # Mixing fractional and Cartesian coordinates
        'xy[0 0]',  # Zero-length vector
        'invalid[1 0]',  # Invalid format
    ],
)
def test_invalid_directions(invalid_direction):
    """Test invalid direction formats."""
    with pytest.raises(ValueError):
        Direction(invalid_direction)
