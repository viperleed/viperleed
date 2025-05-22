"""Module utils of guilib.leedsim.classes.woods.

======================================
  ViPErLEED Graphical User Interface
======================================
 *** package guilib.leedsim.classes.woods ***

Author: Michele Riva (@michele-riva)
Created: 2024-03-01

Defines utility functions used by the Woods class. Refactored from the
former woods.py module.
"""

from collections import Counter

import numpy as np

from viperleed.guilib.helpers import is_integer_matrix
from viperleed.guilib.helpers import prime_numbers


def is_commensurate(matrix):
    """Return whether a matrix represent a commensurate structure.

    Parameters
    ----------
    matrix : Sequence
        Matrix to be tested

    Returns
    -------
    commensurate : bool
        True if matrix is commensurate
    """
    return (is_integer_matrix(matrix, 5e-3)
            and bool(round(np.linalg.det(matrix))))


def prime_factors(number):
    """Yield the prime factors of a number, with repetition."""
    for prime_factor in prime_numbers():
        if prime_factor*prime_factor > number:
            break
        while not number % prime_factor:
            yield prime_factor
            number //= prime_factor
    if number > 1:
        yield number


def square_to_prod_of_squares(number):
    """Decompose integer into two factors.

    Useful to decompose, e.g., sqrt(12) into 2sqrt(3).

    Parameters
    ----------
    number : int
        Integer to be decomposed

    Returns
    -------
    square : int
        The largest perfect square divisor of number
    remainder : int
        The remainder, i.e., number / squares
    """
    # Take number, find all prime factors, return a tuple:
    # first element is the product of all primes showing
    # up an even number of times, the second one the rest.
    factors = Counter(prime_factors(number))
    square, remainder = 1, 1
    for prime_factor, count in factors.items():
        pow2, rest_pow = (count // 2) * 2, count % 2
        square *= prime_factor ** pow2
        remainder *= prime_factor ** rest_pow
    return square, remainder
