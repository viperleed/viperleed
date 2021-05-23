"""Module woods of  guilib.leedsim.classes.

======================================
  ViPErLEED Graphical User Interface
======================================
 *** module guilib.leedsim.classes.woods ***

Defines the Woods class, used for conversion of a Woods' notation string into
a matrix and vice versa, as well as for formatting

Author: Michele Riva
Created: 2021-03-13
"""

import re

import numpy as np

from viperleed import guilib as gl

# TODO: allow also "rect" special Woods notation for hex lattices (ONLY??)

DEGREES = '\u00b0'


class Woods:
    """Class to represent a reconstruction in Wood's notation.

    Also allows conversion to/from a matrix representation.
    """

    common = {'p(1\u00d71)', 'p(1\u00d72)', 'p(2\u00d71)', 'p(2\u00d72)',
              'p(3\u00d71)', 'p(3\u00d72)', 'p(3\u00d73)', 'c(2\u00d72)',
              'c(4\u00d74)', 'c(6\u00d72)', 'c(8\u00d72)'}

    examples = {
        'Oblique': common,
        'Square': common | {'c(4\u00d72)',
                            'p(\u221a2\u00d7\u221a2)R45' + DEGREES,
                            'p(\u221a5\u00d7\u221a5)R26.6' + DEGREES,
                            'p(2\u221a2\u00d7\u221a2)R45' + DEGREES,
                            'c(3\u221a2\u00d7\u221a2)R45' + DEGREES,
                            'c(5\u221a2\u00d7\u221a2)R45' + DEGREES},
        'Rectangular': common,
        'Hexagonal': common | {'c(4\u00d72)',
                               'p(\u221a3\u00d7\u221a3)R30' + DEGREES,
                               'p(\u221a7\u00d7\u221a7)R19.1' + DEGREES,
                               'p(2\u221a3\u00d72\u221a3)R30' + DEGREES},
        'Rhombic': common,
        }

    def to_matrix(self, woods, bulk_basis):
        """Convert woods notation into matrix representation.

        Parameters
        ----------
        woods : str
            A string-representation of a Wood's notation.
        bulk_basis : iterable
            Basis vectors of the bulk. Shape (2, 2).

        Returns
        -------
        matrix : numpy.ndarray or None
            Matrix representation if commensurate, None otherwise
        """
        parsed = self.parse(woods)
        if parsed is None:
            return None

        prefix, gamma1, gamma2, alpha = parsed

        alpha = np.radians(alpha)
        norm1, norm2 = np.sqrt(bulk_basis.T[0]**2 + bulk_basis.T[1]**2)
        basis_ratio = norm2/norm1
        omega = np.arccos(np.dot(bulk_basis[0], bulk_basis[1])/(norm1*norm2))

        matrix = ((gamma1 * np.sin(omega - alpha),
                   gamma1 * np.sin(alpha)/basis_ratio),
                  (-gamma2 * basis_ratio * np.sin(alpha),
                   gamma2 * np.sin(omega + alpha)))
        matrix = np.array(matrix)/np.sin(omega)

        if prefix == 'c':
            matrix = np.dot([[1, 1], [-1, 1]], matrix)/2

        if not self.is_commensurate(matrix):
            return None
        # matrix = np.array([int(np.round(mij))
        #                    for mij in m.ravel()]).reshape(m.shape)
        return matrix.round()

    @staticmethod
    def parse(woods):
        """Parse a Wood's string into conventional blocks.

        Parameters
        ----------
        woods : str
            String to be parsed

        Returns
        -------
        tuple or None
            Returns None when input could not be parsed. Otherwise
            a tuple (prefix, gamma1, gamma2, alpha) where:

        prefix : {'p, 'c'}
            Primitive or centred
        gamma1, gamma2 : float
            Scaling factors along the two directions
        alpha : float
            Rotation angle in degrees
        """
        if not isinstance(woods, str):
            return None

        re_woods = re.compile(
            r'''^(?P<prefix>[pc])                # * primitive or centered
            \(                                   # * open parenthesis
            (?P<gamma1_int>\d+)?                 # * direction1, integer part
            (\u221a(?P<gamma1_sqrt>\d+))?        # * direction1, radical part
            \u00d7                               # * times
            (?P<gamma2_int>\d+)?                 # * direction2, integer part
            (\u221a(?P<gamma2_sqrt>\d+))?        # * direction2, radical part
            \)                                   # * close parenthesis
            ((R(?P<alpha>\d+(\.\d+)?))\u00b0)?$  # * rotation angle
            ''', re.VERBOSE)
        # notice that the ^ and $ anchors make sure that the string is
        # matched as a whole

        # check if it matches the full woods
        match = re_woods.match(woods)
        if not match:
            return None

        groups = match.groupdict()
        gamma1 = 1
        gamma2 = 1
        if groups['gamma1_int'] is not None:
            gamma1 *= float(groups['gamma1_int'])
        if groups['gamma1_sqrt'] is not None:
            gamma1 *= np.sqrt(float(groups['gamma1_sqrt']))
        if groups['gamma2_int'] is not None:
            gamma2 *= float(groups['gamma2_int'])
        if groups['gamma2_sqrt'] is not None:
            gamma2 *= np.sqrt(float(groups['gamma2_sqrt']))
        if groups['alpha']:
            alpha = float(groups['alpha'])
        else:
            alpha = 0
        return groups['prefix'], gamma1, gamma2, alpha

    def from_matrix(self, matrix, bulk_basis):
        """Construct woods string from matrix.

        Parameters
        ----------
        matrix : iterable
            Matrix to be converted to Wood's notation. Shape (2, 2).
        bulk_basis : iterable
            Basis vectors of the bulk lattice.

        Returns
        -------
        woods : str or None
            Woods representiation of the matrix, if
            representable, otherwise None.
        """
        matrix = np.asarray(matrix)

        representable = self.primitive_or_centered(matrix, bulk_basis)
        if not representable:
            return None

        prefix, *gammas, cos_alpha = representable

        woods = f"{prefix}({self.format_scaling_factors(gammas)})"

        if abs(cos_alpha) > 1e-3 and 1 - abs(cos_alpha) > 1e-3:
            # Angle not 0, 90 nor 180
            alpha = round(np.degrees(np.arccos(cos_alpha)),
                          digits=1)
            woods += f"R{alpha}{DEGREES}"
        return woods

    @staticmethod
    def format_scaling_factors(gammas):
        """Format scaling factors in Woods notation.

        Parameters
        ----------
        gammas : tuple
            Two element, both integers. Woods notation
            is ...(gammas[0] <times> gammas[1])...

        Returns
        -------
        str
            String formatted with the gammas passed. Each gamma
            is in the form <gamma_int><sqrt><gamma_radical>, where
            <sqrt> is the unicode symbol for square root '\u221a'.
            The two blocks are joined with a <times> ('\u00d7')
            character.
        """
        to_format = []

        for gamma in gammas:
            gamma_square = round(gamma**2)
            (gamma_int_squared,
             gamma_sqrt) = Woods.square_to_prod_of_squares(gamma_square)
            gamma_int = round(np.sqrt(gamma_int_squared))

            format_direction = ''  # Format the direction in here.
            if gamma_sqrt > 1:     # Root part
                format_direction = f"\u221a{gamma_sqrt}"

            if not format_direction:
                # If there is no root part, always
                # insert the integer part
                format_direction = str(gamma_int)
            elif gamma_int > 1:
                # Otherwise add it if it's not 1
                format_direction = str(gamma_int) + format_direction
            to_format.append(format_direction)

        return '\u00d7'.join(to_format)

    @staticmethod
    def is_commensurate(matrix, eps=1e-3):
        """Return whether a matrix represent a commensurate structure.

        Parameters
        ----------
        matrix : iterable
            Matrix to be tested
        eps : float, default=1e-3
            Relative tolerance for assuming thins equal

        Returns
        -------
        bool
            True if matrix is commensurate
        """
        if matrix is None:
            return False

        return np.all(np.mod(np.abs(matrix), 1.0) < eps)

        # OLD CODE FOLLOWS, to check if new ne returns correct results
        #
        # det = np.linalg.det(matrix)
        # if round(det) == 0:  # matrix is singular
        #     return False

        # if det % 1 > eps*abs(round(det)):  # determinant is not int
        #     return False

        # # now check whether any element is non-integer
        # for mij in matrix.ravel():
        #     if np.round(mij) == 0:
        #         if abs(mij)/np.sqrt(np.abs(mu)) > eps:
        #             return False
        #     elif np.abs(mij/np.round(mij) - 1) > eps:
        #         return False
        # return True

    @staticmethod
    def is_representable(matrix, basis):
        """Return whether a matrix is Woods-representable."""
        transformed_basis = np.dot(matrix, basis)

        basis_norm = np.linalg.norm(basis, axis=1)
        transf_norm = np.linalg.norm(transformed_basis, axis=1)
        gamma1, gamma2 = transf_norm/basis_norm

        det = np.abs(np.linalg.det(matrix))

        return abs(det/(gamma1*gamma2) - 1) < 1e-8

    def primitive_or_centered(self, matrix, basis):
        """Check if matrix is representable as primitive or centered.

        Parameters
        ----------
        matrix : iterable
            Matrix to be checked
        basis : iterable
            Basis vectors of the bulk lattice.

        Returns
        -------
        None or tuple
            None if not representable, otherwise

        prefix : {'p, 'c'}
            Primitive or centred
        gamma1, gamma2 : float
            Scaling factors in the two directions of basis
        cos_alpha : float
            Cosine of the rotation angle
        """
        centering_inv_transpose = (1, -1), (1, 1)

        primitive = self.is_representable(matrix, basis)
        centered = self.is_representable(
            np.dot(centering_inv_transpose, matrix), basis
            )
        if primitive:
            prefix = 'p'
        elif centered:
            prefix = 'c'
            matrix = np.dot(centering_inv_transpose, matrix)
        else:
            return None

        transformed_basis = np.dot(matrix, basis)

        basis_norm = np.linalg.norm(basis, axis=1)
        transformed_norm = np.linalg.norm(transformed_basis, axis=1)
        gamma1, gamma2 = transformed_norm/basis_norm

        cos_alpha = np.dot(transformed_basis[0], basis[0])
        cos_alpha *= 1/(transformed_norm[0]*basis_norm[0])

        return (prefix, gamma1, gamma2, cos_alpha)

    @staticmethod
    def square_to_prod_of_squares(number):
        """Decompose integer into two factors.

        Useful to decompose, e.g., sqrt(12) into 2sqrt(3).

        Parameters
        ----------
        number : int
            Integer to be decomposed

        Returns
        -------
        squares : int
        remainder : int

        Such that number = squares * remainder, where square
        is the largest perfect square factor of number
        """
        # takes number, finds all prime factors, returns a tuple, first
        # element is a product of all primes showing up an even number of
        # times, the second one the rest. Useful to turn, e.g., sqrt(12)
        # into 2sqrt(3)

        factors = list(Woods.prime_factors(number))
        if not factors:
            factors = [1]
        unique_factors = sorted(list(set(factors)))
        count_factors = [((factors.count(fac) // 2)*2, factors.count(fac) % 2)
                         for fac in unique_factors]
        pow2, rest_pow = zip(*count_factors)
        squares, remainders = zip(*[(fact**power, fact**rem)
                                    for (fact, power, rem)
                                    in zip(unique_factors, pow2, rest_pow)])

        return round(np.prod(squares)), round(np.prod(remainders))

    @staticmethod
    def prime_factors(number):
        """Yield the prime factors of a number, with repetition."""
        for prime_factor in gl.prime_numbers():
            if prime_factor*prime_factor > number:
                break
            while number % prime_factor == 0:
                yield prime_factor
                number //= prime_factor
        if number > 1:
            yield number
