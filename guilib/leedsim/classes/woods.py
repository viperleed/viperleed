"""Module woods of guilib.leedsim.classes.

======================================
  ViPErLEED Graphical User Interface
======================================
 *** module guilib.leedsim.classes.woods ***

Defines the Woods class, used for conversion of a Wood's notation
string into a matrix and vice versa, as well as for formatting

Author: Michele Riva
Created: 2021-03-13 (was part of leedsim.classes.py before)
"""

from collections import Counter
from dataclasses import dataclass, field
import re

import numpy as np

from viperleed.guilib.helpers import array_to_string
from viperleed.guilib.helpers import is_integer_matrix
from viperleed.guilib.helpers import prime_numbers
from viperleed.guilib.mathparse import MathParser, UnsupportedMathError

# TODO: allow also "rect" special Woods notation for hex lattices (ONLY??)

# Unicode symbols
DEGREES = '\u00b0'
TIMES = '\u00d7'
SQRT = '\u221a'

# In the following regular expression (spaces ignored everywhere):
# (1) <prefix>
#     Any letter sequence, provided it does not contain 's', 'q',
#     'r', 't' (also capital). Will be turned into 'p' (most of
#     the times) or 'c' (if there is a 'c'/'C')
# (2) optional open parenthesis
# (3) <gamma1>
#     Any character. Will be parsed with MathParser.
# (4) A single 'times' character ('x', 'X', '\u00d7')
# (5) <gamma2>
#     Any character, except for an 'R', used as delimiter for
#     the rotations block. Will be parsed with MathParser.
#     Also, it may 'eat up' the closing parenthesis. This is
#     checked for before parsing.
# (6) optional closed parenthesis
# (7) rotation block
#     structure is 'R'<alpha> (possibly with spaces in
#     between), optionally followed by a degrees designator
#     ('\u00b0', or a contraction of 'deg'), with <alpha>
#     a (possibly signed) integer or floating-point number.
MATCH_WOODS = re.compile(
    fr'''
    ^(?P<prefix>((?![sqrtSQRT])[\sa-zA-Z])*)
    \s*\(?
    (?P<gamma1>.*)
    [xX{TIMES}]
    (?P<gamma2>[^R]*)
    \)?\s*
    ([R]\s*(?P<alpha>-?\d+(\.\d+)?)\s*({DEGREES}|d[eg]+)?)?
    \s*$''',
    re.VERBOSE
    )


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


class MatrixIncommensurateError(Exception):
    """Matrix is incommensurate."""

    def __init__(self, matrix, message=''):
        if not message:
            message = (f'Matrix {array_to_string(matrix)} '
                       'gives an incommensurate lattice, i.e., '
                       'it is singular or has non-integer elements')
        super().__init__(message)


class WoodsNotRepresentableError(Exception):
    """Matrix is not Wood's-representable."""

    def __init__(self, matrix, message=''):
        if not message:
            message = (f'Matrix {array_to_string(matrix)} '
                       'is not Woods-representable')
        super().__init__(message)


class WoodsSyntaxError(Exception):
    """Exception raised in case a Wood's string has invalid syntax."""


@dataclass(frozen=True)
class _WoodsStyle:
    """Simple collection of style-dependent symbols for Wood's notations."""
    style: str
    degrees: str = field(init=False)
    times: str = field(init=False)
    sqrt: str = field(init=False)

    def __eq__(self, other):
        """Return whether self == other."""
        if isinstance(other, str):
            return self.style == other
        return super().__eq__(other)

    def __repr__(self):
        """Return a string version of this style."""
        return repr(self.style)

    def __post_init__(self):
        """Fill in the symbols given a style."""
        style = self.style
        if not isinstance(style, (str, _WoodsStyle)):
            raise TypeError(f'Woods: invalid style {type(style).__name__!r}. '
                            "Expected 'str'")
        if isinstance(style, _WoodsStyle):
            style = style.style

        if style[0].lower() not in 'ua':
            raise ValueError(f'Woods: invalid style {style!r}. '
                             "Expected 'unicode' or 'ascii'")
        style = 'unicode' if style[0].lower() == 'u' else 'ascii'
        setattr_ = object.__setattr__
        setattr_(self, 'style', style)
        setattr_(self, 'degrees', '' if style == 'ascii' else DEGREES)
        setattr_(self, 'times', 'x' if style == 'ascii' else TIMES)
        setattr_(self, 'sqrt', 'sqrt' if style == 'ascii' else SQRT)

    def startswith(self, string):
        """Return whether self.style starts with string."""
        return self.style.startswith(string)


class Woods:
    """Class to represent a reconstruction in Wood's notation.

    Also allows conversion to/from a matrix representation.
    """

    __common = {f'p(1{TIMES}1)', f'p(1{TIMES}2)', f'p(2{TIMES}1)',
                f'p(2{TIMES}2)', f'p(3{TIMES}1)', f'p(3{TIMES}2)',
                f'p(3{TIMES}3)', f'c(2{TIMES}2)', f'c(4{TIMES}4)',
                f'c(6{TIMES}2)', f'c(8{TIMES}2)'}

    __examples = {
        'Oblique': __common.copy(),
        'Square': __common | {f'c(4{TIMES}2)',
                              f'p({SQRT}2{TIMES}{SQRT}2)R45{DEGREES}',
                              f'p({SQRT}5{TIMES}{SQRT}5)R26.6{DEGREES}',
                              f'p(2{SQRT}2{TIMES}{SQRT}2)R45{DEGREES}',
                              f'c(3{SQRT}2{TIMES}{SQRT}2)R45{DEGREES}',
                              f'c(5{SQRT}2{TIMES}{SQRT}2)R45{DEGREES}'},
        'Rectangular': __common.copy(),
        'Hexagonal': __common | {f'c(4{SQRT}2)',
                                 f'p({SQRT}3{TIMES}{SQRT}3)R30{DEGREES}',
                                 f'p({SQRT}7{TIMES}{SQRT}7)R19.1{DEGREES}',
                                 f'p(2{SQRT}3{TIMES}2{SQRT}3)R30{DEGREES}'},
        'Rhombic': __common.copy(),
        }

    # Disable for 'style' due to pylint bug: optional
    # not allowed in 'list of allowed values'
    # pylint: disable-next=missing-type-doc
    def __init__(self, string='', bulk_basis=None,
                 matrix=None, style='unicode'):
        """Initialize instances.

        Parameters
        ----------
        string : str, optional
            String representation of Woods notation. Will
            be automatically formatted to the 'style' given.
        bulk_basis : sequence, optional
            Shape (2, 2). Real-space basis of the bulk lattice,
            a == basis[0], b == basis[1]. 'bulk_basis' is
            mandatory when 'matrix' is given. Default is None.
        matrix : sequence, optional
            Shape (2, 2). Matrix representation of Woods
            notation. If given, 'bulk_basis' becomes a
            mandatory argument. Default is None.
        style : {'unicode', 'ascii'}, optional
            The character set that will be used to represent
            the Wood's notation. Only the first letter passed
            is used (case insensitive). Default is 'unicode'.

        Raises
        ------
        MatrixIncommensurateError
            If matrix represents a non-commensurate lattice,
            i.e., it is singular, or its elements are not
            integer-representable.
        TypeError
            If matrix is given, but bulk_basis isn't
        TypeError
            If woods is not a string
        ValueError
            If 'style' is not acceptable
        ValueError
            If shape of matrix or bulk_basis is not (2, 2)
        ValueError:
            If both string and matrix are given, but they
            are inconsistent (i.e., they give different
            Wood's notations)
        ValueError
            If the parsed string contains scaling factors
            for the two directions whose square is not an
            integer
        WoodsNotRepresentableError
            If matrix is valid but not Wood's-representable.
        WoodsSyntaxError
            If the string does not match the structure of a
            Wood's notation, or could not be evaluated due
            to either unsupported, unmatched brackets or
            other any unexpected characters.
        """
        self.__style = _WoodsStyle(style)
        if matrix is not None and bulk_basis is None:
            raise TypeError('Woods: bulk_basis is mandatory when '
                            'a matrix is given')

        # Store only the string representation and the bulk
        # basis. If a matrix is given, it will be used to
        # construct the string. For the stored properties
        # defer checks and processing to the property setters
        self.__string = ''
        self.__bulk_basis = None
        self.bulk_basis = bulk_basis

        if matrix is None:
            # Only string given
            self.string = string
            return

        string = self.__fix_string(string)
        self.matrix = matrix   # Also sets self.string
        if string and string != self.string:
            # Both matrix and string given
            raise ValueError('Woods: matrix and string are inconsistent.')

    def __repr__(self):
        """Return a representation string of self."""
        txt = f'{self.string!r}'
        if self.bulk_basis is not None:
            txt += f', bulk_basis={array_to_string(self.bulk_basis)}'
        if self.style != 'unicode':  # Hide the default
            txt += f', style={self.style!r}'
        return f'Woods({txt})'

    def __str__(self):
        """Return a string representation of self.

        The style used is the one in self.style.
        """
        return self.string

    def __format__(self, format_spec=''):
        """Format self according to a given format_spec.

        Parameters
        ----------
        format_spec : {'u', 'a', 'su', 'sa', 's', ''}, default=''
            'u' for unicode formatting, 'a' for ascii, '' for native,
            i.e., using self.style. Prefixing with 's' (i.e., simple)
            removes the 'p' prefix for primitive cells.

        Returns
        -------
        formatted : str
            Formatted string representation of self

        Raises
        ------
        TypeError
            If format_spec is invalid
        """
        fmt_style = format_spec[-1:]  # 'u', 'a', or 's'
        if fmt_style == 's':
            # This means format_spec was 's'
            fmt_style = ''
        if (fmt_style and fmt_style not in 'ua') or (len(format_spec) > 2):
            raise TypeError('unsupported format string '
                            'passed to Woods.__format__')
        is_default_style = not fmt_style or self.style.startswith(fmt_style)
        fmt = (self.string if is_default_style
               else Woods(self.string, style=fmt_style).string)
        if fmt[0] == 'p' and 's' in format_spec:
            return fmt[1:]
        return fmt

    @staticmethod
    def parse(woods):
        """Parse a Wood's string into conventional blocks.

        A rather versatile regular expression is used, which requires
        essentially any two blocks separated by a <times> character
        (i.e., 'x', 'X', or '\u00d7'), with optional parentheses,
        optional rotation, and ignoring spaces. Each block can contain
        any arithmetic expression, also including square roots.

        Parameters
        ----------
        woods : str
            String to be parsed

        Returns
        -------
        prefix : {'p, 'c'}
            Primitive or centered
        gamma1, gamma2 : float
            Scaling factors along the two directions
        alpha : float
            Rotation angle in degrees

        Raises
        ------
        TypeError
            If woods is not a string
        ValueError
            If the parsed string contains scaling factors for the two
            directions that are not the square root of an integer
        WoodsSyntaxError
            If the string does not match the structure of a
            Wood's notation, or could not be evaluated due
            to either unsupported math or unmatched brackets
        """
        if not isinstance(woods, str):
            raise TypeError('Woods.parse: argument woods should '
                            f"be 'str', not {type(woods).__name__!r}")
        match = MATCH_WOODS.match(woods)
        if not match:
            raise WoodsSyntaxError(f'Woods: {woods} is not a '
                                   'valid Woods notation.')
        # Force prefix to be p or c
        prefix = 'c' if 'c' in match['prefix'].lower() else 'p'
        gammas = (Woods.__parse_and_check_gamma(match[key])
                  for key in ('gamma1', 'gamma2'))
        alpha = float(match['alpha']) if match['alpha'] else 0.
        return prefix, *gammas, alpha

    @staticmethod
    def __parse_and_check_gamma(gamma_str):
        """Return a float from the string of a direction-scaling factor.

        A scaling factor gamma is part of a Wood's notation, that
        takes the form <prefix>(<gamma1>x<gamma2>)R<angle>.

        Parameters
        ----------
        gamma_str : str
            The string representation of a scaling factor. Can
            contain only arithmetic expressions and square roots

        Returns
        -------
        gamma : float
            The parsed value of the string given

        Raises
        ------
        WoodsSyntaxError
            If gamma_str could not be evaluated due to either
            unsupported math or unmatched brackets
        ValueError
            If gamma_str can be parsed, but its value is not the
            square root of an integer
        """
        parser = MathParser()

        # Given the way the parser is, gamma2 may
        # eat up the closing parenthesis
        if gamma_str.count('(')  == gamma_str.count(')') - 1:
            # Get rid of the swallowed closing
            # parenthesis, i.e., the last ')'
            idx = gamma_str.rindex(')')
            gamma_str = gamma_str[:idx] + gamma_str[idx+1:]

        # parser needs at least 'rt' for square root.
        # Support also 'r' for Woods notations
        parser.expression = re.sub(r'r([^t])', r'rt\1', gamma_str)
        try:
            gamma = parser.evaluate()
        except SyntaxError as exc:
            raise WoodsSyntaxError(f'Woods: {gamma_str} could '
                                   'not be parsed, likely because '
                                   'of unmatched brackets') from exc
        except UnsupportedMathError as exc:
            raise WoodsSyntaxError(f'Woods: {gamma_str} could '
                                   'not be parsed as it contains '
                                   'unsupported math operations') from exc

        # Make sure that gamma**2 is close to int
        if abs(gamma**2 - round(gamma**2)) > 1e-3:
            raise ValueError(
                f'Woods: {gamma_str} (evaluated to {gamma}) cannot '
                'be represented as the square root of an integer'
                )
        return gamma

    @property
    def bulk_basis(self):
        """Return the bulk basis used for conversion to/from matrix.

        Returns
        -------
        numpy.ndarray
        """
        return self.__bulk_basis

    @bulk_basis.setter
    def bulk_basis(self, basis):
        """Set the bulk basis used to convert to/from matrix.

        Parameters
        ----------
        basis : sequence or None
            If a sequence: shape==(2, 2). a==basis[0], b==basis[1].

        Raises
        ------
        ValueError
            If shape of basis is not (2, 2)
        """
        if basis is None:
            self.__bulk_basis = basis
            return
        basis = np.asarray(basis)
        if basis.shape != (2, 2):
            raise ValueError('Woods: invalid bulk_basis. Should have '
                             f'shape==(2, 2). Found {basis.shape}')
        self.__bulk_basis = basis

    @property
    def matrix(self):
        """Return the matrix representation of self."""
        return self.to_matrix()

    @matrix.setter
    def matrix(self, matrix):
        """Convert matrix to string Wood's.

        Parameters
        ----------
        matrix : sequence
            shape == (2, 2), integer elements

        Raises
        ------
        ValueError
            If shape of matrix is not (2, 2)
        MatrixIncommensurateError
            If matrix represents a non-commensurate lattice, i.e., it
            is singular, or its elements are not integer-representable.
        WoodsNotRepresentableError
            If matrix is valid but not Wood's-representable.
        """
        self.from_matrix(matrix)

    @property
    def string(self):
        """Return a string representation of self.

        The style used is the one in self.style.
        """
        return self.__string

    @string.setter
    def string(self, woods):
        """Set the string representation of self.

        Parameters
        ----------
        woods : str
            String representation of Wood's notation. Will be
            preprocessed to conform to the style in self.style

        Raises
        ------
        TypeError
            If woods is not a string
        ValueError
            If the parsed string contains scaling factors for the
            two directions that are not the square root of an integer
        WoodsSyntaxError
            If the string does not match the structure of a
            Wood's notation, or could not be evaluated due
            to either unsupported math or unmatched brackets
        """
        self.__string = self.__fix_string(woods)

    @property
    def style(self):
        """Return the style of this instance."""
        return self.__style

    @classmethod
    def add_example(cls, woods, shape):
        """Add an example Wood's notation for a given cell shape.

        Parameters
        ----------
        woods : str or Woods
            The example to add
        shape : {'Oblique', 'Square', 'Rectangular',
                 'Hexagonal', 'Rhombic'}
            Shape of the surface unit cell

        Returns
        -------
        None.
        """
        if isinstance(woods, Woods):
            woods = woods.string
        woods = format(woods)
        cls.__examples[shape].add(woods)

    def from_matrix(self, matrix, bulk_basis=None):
        """Construct woods string from matrix.

        If the conversion is successful, the string
        representation of the Wood's notation is also
        stored into self.string.

        Parameters
        ----------
        matrix : Sequence
            Matrix to be converted to Wood's notation. Shape (2, 2).
        bulk_basis : Sequence or None, optional
            Basis vectors of the bulk lattice. If given and acceptable,
            it will be stored into self.bulk_basis. If not given or
            None, self.bulk_basis is used instead. Default is None.

        Returns
        -------
        woods : str
            Woods representation of the matrix.

        Raises
        ------
        MatrixIncommensurateError
            If matrix represents a non-commensurate lattice, i.e., it
            is singular, or its elements are not integer-representable.
        ValueError
            If shape of matrix or bulk_basis is not (2, 2)
        WoodsNotRepresentableError
            If matrix is valid but not Wood's-representable.
        """
        matrix = np.asarray(matrix)
        if matrix.shape != (2, 2):
            raise ValueError('Woods: invalid matrix. Should have '
                             f'shape==(2, 2). Found {matrix.shape}')
        if not is_commensurate(matrix):
            raise MatrixIncommensurateError(matrix)
        matrix = matrix.round().astype(int)

        if bulk_basis is None and self.bulk_basis is None:
            raise ValueError('Woods: cannot convert from matrix '
                             'without a bulk basis. Pass a bulk_basis '
                             'or set one via the bulk_basis property')
        if bulk_basis is not None:
            self.bulk_basis = bulk_basis

        # Next line may raise WoodsNotRepresentableError
        prefix, *gammas, alpha = self.__primitive_or_centered(matrix)

        # Don't use the setter of self.string to skip
        # reformatting, which is anyway already correct
        self.__string = self.__format_parsed(prefix, gammas, alpha)
        return self.string

    def get_examples(self, shape):
        """Return examples of Wood's for a given cell shape.

        The examples are formatted to the current style.

        Parameters
        ----------
        shape : {'Oblique', 'Square', 'Rectangular',
                 'Hexagonal', 'Rhombic'}
            Shape of the surface unit cell

        Returns
        -------
        examples : set of str
            Examples of Wood's notation for the given shape.
        """
        if self.style == 'unicode':
            return self.__examples[shape]
        cls = type(self)
        return {format(cls(ex)) for ex in self.__examples[shape]}

    def guess_correct_rotation(self, woods_txt=None):
        """Guess the rotation angle that gives an acceptable Woods.

        Typically, this would be used in case the user gives
        an angle that is off by (at most) one degree, and this
        results in a MatrixIncommensurateError. If the guess is
        successful, the corrected text is stored into self.string.

        Parameters
        ----------
        woods_txt : str
            The Woods string from which to guess the angle.

        Returns
        -------
        guessed_angle : float
            The new angle guessed, if successful

        Raises
        ------
        MatrixIncommensurateError
            If the angle guessed and the one passed are far off.
        """
        if woods_txt is None:
            woods_txt = self.string

        # Prepare a copy that will be used for messing around
        tmp_woods = Woods(string=woods_txt, bulk_basis=self.bulk_basis,
                          style=self.style)
        orig_prefix, *orig_gamma, orig_alpha = self.parse(tmp_woods.string)
        orig_matrix = tmp_woods.to_matrix(check_commensurate=False)

        if is_commensurate(orig_matrix):
            # Skip using the setter of self.string to avoid
            # reformatting, as the format is already correct
            self.__string = tmp_woods.string
            return self.string

        # Matrix is incommensurate. Try to round
        # it to int, and redo the parsing
        rounded_matrix = orig_matrix.round()
        rounded_txt = tmp_woods.from_matrix(rounded_matrix)
        round_prefix, *round_gamma, round_alpha = self.parse(rounded_txt)

        delta_gamma = np.abs(np.subtract(orig_gamma, round_gamma)/orig_gamma)
        # Prefix should stay the same; same for gammas
        # (the user rarely inputs the wrong scaling);
        # New angle is acceptable if it is within +-1
        # degrees from the original input
        if (round_prefix != orig_prefix
            or np.any(delta_gamma > 1e-3)
                or abs(round_alpha - orig_alpha) > 1):
            raise MatrixIncommensurateError(orig_matrix)

        # Correction was successful.  Skip using self.string
        # to avoid reformatting, as the format is already correct
        self.__string = rounded_txt
        return self.string

    def to_matrix(self, woods='', bulk_basis=None, check_commensurate=True):
        """Convert woods notation into matrix representation.

        Parameters
        ----------
        woods : str, optional
            A string-representation of a Wood's notation. Will be set
            to self.string if parsing is successful. If not given (or
            empty) self.string is used instead for the conversions.
            Default is an empty string.
        bulk_basis : Sequence or None, optional
            Basis vectors of the bulk. Shape (2, 2). If not given or
            None, self.bulk_basis is used, otherwise the value passed
            is stored into self.bulk_basis. Default is None.
        check_commensurate : bool, optional
            If True, raise MatrixIncommensurateError in case the
            resulting matrix is incommensurate. Default is True.

        Returns
        -------
        matrix : numpy.ndarray
            Matrix representation. If check_commensurate is True and
            matrix is commensurate, it returns a rounded matrix of
            integers, otherwise the matrix is returned as-is.

        Raises
        ------
        MatrixIncommensurateError
            If check_commensurate is True, and the matrix resulting
            from conversion gives an incommensurate lattice
        TypeError
            If woods is not a string
        ValueError
            It the parsed string contains scaling factors for the two
            directions that are not the square root of an integer
        WoodsSyntaxError
            If the string does not match the structure of a
            Wood's notation, or could not be evaluated due
            to either unsupported math or unmatched brackets
        """
        if not woods:
            if not self.string:
                raise WoodsSyntaxError('Woods: empty string is not a '
                                       'valid Woods notation. Pass a '
                                       'valid, non-empty string, or '
                                       'set self.string beforehand')
            woods = self.string
        if bulk_basis is not None:
            self.bulk_basis = bulk_basis
        if self.bulk_basis is None:
            raise ValueError(f'{type(self).__name__}: Cannot convert {self} '
                             'to_matrix() without a bulk_basis. Pass a valid '
                             'bulk_basis, or set self.bulk_basis beforehand')

        # The next line may raise
        prefix, gamma1, gamma2, alpha = self.parse(woods)
        if woods != self.string:
            self.string = woods

        alpha = np.radians(alpha)
        norm1, norm2 = np.linalg.norm(self.bulk_basis, axis=1)
        basis_ratio = norm2/norm1
        omega = np.arccos(np.dot(*self.bulk_basis)/(norm1*norm2))

        matrix = np.array(((gamma1 * np.sin(omega - alpha),
                            gamma1 * np.sin(alpha)/basis_ratio),
                           (-gamma2 * basis_ratio * np.sin(alpha),
                            gamma2 * np.sin(omega + alpha))))
        matrix /= np.sin(omega)

        if prefix == 'c':
            matrix = np.dot([[1, 1], [-1, 1]], matrix)/2

        if not check_commensurate:
            return matrix

        if not is_commensurate(matrix):
            raise MatrixIncommensurateError(matrix)
        return matrix.round().astype(int)

    def __fix_string(self, woods):
        """Reformat a woods string to standard style.

        A rather versatile regular expression is used, which requires
        essentially any two blocks separated by a <times> character
        (i.e., 'x', 'X', or '\u00d7'), with optional parentheses,
        optional rotation, and ignoring spaces. Each block can contain
        any arithmetic expression, also including square roots

        Parameters
        ----------
        woods : str
            The string to be fixed.

        Returns
        -------
        fixed_woods : str
            String reformatted according to self.style.

        Raises
        ------
        TypeError
            If woods is not a string
        ValueError
            If the parsed string contains scaling factors for the two
            directions that are not the square root of an integer
        WoodsSyntaxError
            If the string does not match the structure of a
            Wood's notation, or could not be evaluated due
            to either unsupported math or unmatched brackets
        """
        if not isinstance(woods, str):
            raise TypeError("Woods: argument should be 'str', "
                            f'not {type(woods).__name__!r}')
        if not woods:
            return woods
        prefix, *gammas, alpha = self.parse(woods)
        return self.__format_parsed(prefix, gammas, alpha)

    def __format_parsed(self, prefix, gammas, alpha):
        """Return a standardized Woods notation from parsed values."""
        woods = f'{prefix}({self.__format_scaling_factors(gammas)})'
        fmt_alpha = self.__format_rotation(alpha)
        if fmt_alpha:
            woods += f'R{fmt_alpha}{self.style.degrees}'
        return woods

    def __format_rotation(self, alpha):
        """Return a formatted version of the rotation angle in degrees."""
        if abs(alpha) < 0.25:  # degrees
            return ''
        abs_cos_alpha = abs(np.cos(np.radians(alpha)))
        # Limits below are different because cos(90 + x) ~ x, but
        # cos(x) ~ 1 - x**2/2. Both limits tolerate x ~ 0.25deg
        if abs_cos_alpha > 4e-3 and 1 - abs_cos_alpha > 8e-6:
            # Angle is neither 90, nor 0 or 180
            return f'{alpha:.1f}'
        return ''

    def __format_scaling_factors(self, gammas):
        """Format scaling factors in Woods notation.

        Parameters
        ----------
        gammas : sequence
            Two elements, both integers. Woods notation
            is ...(gammas[0] <times> gammas[1])...

        Returns
        -------
        str
            <gamma_1><times><gamma_2>. Each <gamma_i> is in the
            form <gamma_int><sqrt><gamma_radical>. <gamma_radical>
            is omitted if it is 1. <gamma_int> is omitted if it is
            1, unless also <gamma_radical> is omitted. <times> and
            <sqrt> have the style given by self.style.
        """
        to_format = []

        for gamma in gammas:
            gamma_square = round(gamma**2)
            (gamma_int_squared,
             gamma_sqrt) = square_to_prod_of_squares(gamma_square)
            gamma_int = round(np.sqrt(gamma_int_squared))

            format_direction = ''  # Format the direction in here.
            if gamma_sqrt > 1:     # Root part
                format_direction = f'{self.style.sqrt}{gamma_sqrt}'

            if not format_direction:
                # If there is no root part, always
                # insert the integer part
                format_direction = str(gamma_int)
            elif gamma_int > 1:
                # Otherwise add it if it's not 1
                format_direction = str(gamma_int) + format_direction
            to_format.append(format_direction)

        times = self.style.times
        if self.style == 'ascii' and any('sqrt' in g for g in to_format):
            # Add some spaces around the 'x' if there
            # are 'sqrt's to make it look less messy
            times = f' {times} '
        return times.join(to_format)

    def __is_representable(self, matrix):
        """Return whether a matrix is Woods-representable.

        Make sure self.bulk_basis is up to date before calling.

        Parameters
        ----------
        matrix : sequence
            The matrix to be checked. Shape (2, 2).

        Returns
        -------
        bool
            True if matrix is Wood's-representable
        """
        basis = self.bulk_basis
        transformed_basis = np.dot(matrix, basis)

        basis_norm = np.linalg.norm(basis, axis=1)
        transf_norm = np.linalg.norm(transformed_basis, axis=1)
        gamma1, gamma2 = transf_norm/basis_norm

        det = np.abs(np.linalg.det(matrix))

        return abs(det/(gamma1*gamma2) - 1) < 1e-8

    def __primitive_or_centered(self, matrix):
        """Check if matrix is representable as primitive or centered.

        Make sure self.bulk_basis is up to date before calling.

        Parameters
        ----------
        matrix : Sequence
            Matrix to be checked

        Returns
        -------
        prefix : {'p, 'c'}
            Primitive or centered
        gamma1, gamma2 : float
            Scaling factors in the two directions of basis
        alpha : float
            Rotation angle, in degrees.

        Raises
        ------
        WoodsNotRepresentableError
            If matrix not Wood's-representable, with neither
            a 'p' or a 'c' notation.
        """
        centering_inv_transpose = (1, -1), (1, 1)

        primitive = self.__is_representable(matrix)
        centered = self.__is_representable(
            np.dot(centering_inv_transpose, matrix)
            )
        if primitive:
            prefix = 'p'
        elif centered:
            prefix = 'c'
            matrix = np.dot(centering_inv_transpose, matrix)
        else:
            raise WoodsNotRepresentableError(matrix)

        basis = self.bulk_basis
        transformed_basis = np.dot(matrix, basis)

        basis_norm = np.linalg.norm(basis, axis=1)
        transformed_norm = np.linalg.norm(transformed_basis, axis=1)
        gamma1, gamma2 = transformed_norm/basis_norm

        # Here is the calculation of the rotation angle  from the
        # first lattice vectors. In theory, the same calculation
        # from the second lattice vectors should yield the same
        alpha = np.arctan2(np.cross(basis[0], transformed_basis[0]),
                           np.dot(basis[0], transformed_basis[0]))
        return prefix, gamma1, gamma2, np.degrees(alpha)
