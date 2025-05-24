"""Module _woods of guilib.leedsim.classes.woods.

======================================
  ViPErLEED Graphical User Interface
======================================
 *** package guilib.leedsim.classes.woods ***

Author: Michele Riva (@michele-riva)
Created: 2021-03-13 (was part of leedsim.classes.py before)
Refactored: 2024-03-01

Defines the Woods class, used for conversion of a Wood's notation
string into a matrix and vice versa, as well as for formatting.
"""

from dataclasses import dataclass, field
import re

import numpy as np

from viperleed.guilib.helpers import array_to_string
from viperleed.guilib.mathparse import MathParser
from viperleed.guilib.mathparse import UnsupportedMathError
from viperleed.guilib.mathparse import TooComplexMathError

from .errors import MatrixIncommensurateError
from .errors import WoodsInvalidForBasisError
from .errors import WoodsNotRepresentableError
from .errors import WoodsSyntaxError
from .utils import is_commensurate
from .utils import square_to_prod_of_squares

# TODO: allow also 'rect' special Woods notation for hex lattices (ONLY??)

# Unicode symbols
DEGREES = '\u00b0'
TIMES = '\u00d7'
SQRT = '\u221a'

# In the following regular expression (spaces ignored everywhere):
# (1) <prefix>, optional
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
# (7) rotation block, optional
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

# Some examples of common Wood's notations. In Unicode style.
_COMMON = {f'p(1{TIMES}1)', f'p(1{TIMES}2)', f'p(2{TIMES}1)', f'p(2{TIMES}2)',
           f'p(3{TIMES}1)', f'p(3{TIMES}2)', f'p(3{TIMES}3)', f'c(2{TIMES}2)',
           f'c(4{TIMES}4)', f'c(6{TIMES}2)', f'c(8{TIMES}2)'}
_EXAMPLES = {
        'Oblique': _COMMON,
        'Square': _COMMON | {f'c(4{TIMES}2)',
                             f'p({SQRT}2{TIMES}{SQRT}2)R45{DEGREES}',
                             f'p({SQRT}5{TIMES}{SQRT}5)R26.6{DEGREES}',
                             f'p(2{SQRT}2{TIMES}{SQRT}2)R45{DEGREES}',
                             f'c(3{SQRT}2{TIMES}{SQRT}2)R45{DEGREES}',
                             f'c(5{SQRT}2{TIMES}{SQRT}2)R45{DEGREES}'},
        'Rectangular': _COMMON,
        'Hexagonal': _COMMON | {f'c(4{TIMES}2)',
                                f'p({SQRT}3{TIMES}{SQRT}3)R30{DEGREES}',
                                f'p({SQRT}7{TIMES}{SQRT}7)R19.1{DEGREES}',
                                f'p(2{SQRT}3{TIMES}2{SQRT}3)R30{DEGREES}'},
        'Rhombic': _COMMON,
        }


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
        if isinstance(style, _WoodsStyle):
            style = style.style
        if not isinstance(style, str):
            raise TypeError(f'Woods: invalid style {type(style).__name__!r}. '
                            "Expected 'str'")
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
    __examples = {}  # Filled once via make_examples

    # Disable for 'style' due to pylint bug: optional
    # not allowed in 'list of allowed values'
    # pylint: disable-next=missing-type-doc
    def __init__(self, string='', bulk_basis=None,
                 matrix=None, style='unicode'):
        """Initialize instance.

        Parameters
        ----------
        string : str, optional
            String representation of Woods notation. Will
            be automatically formatted to the `style` given.
        bulk_basis : sequence, optional
            Shape (2, 2). Real-space basis of the bulk lattice,
            a == basis[0], b == basis[1]. `bulk_basis` is
            mandatory when `matrix` is given. Default is None.
        matrix : sequence, optional
            Shape (2, 2). Matrix representation of Woods
            notation. If given, `bulk_basis` becomes a
            mandatory argument. Default is None.
        style : {'unicode', 'ascii'}, optional
            The character set that will be used to represent
            the Wood's notation. Only the first letter passed
            is used (case insensitive). Default is 'unicode'.

        Raises
        ------
        MatrixIncommensurateError
            If `matrix` is singular, or its elements are not
            integer-representable.
        TypeError
            If matrix is given, but `bulk_basis` isn't.
        TypeError
            If `string` or `style` is not a string.
        ValueError
            If 'style' is not acceptable
        ValueError
            If shape of `matrix` or `bulk_basis` is not (2, 2)
        ValueError:
            If both `string` and `matrix` are given, but they are
            inconsistent (i.e., they give different Wood's notations).
        ValueError
            If the parsed string contains scaling factors for
            the two directions whose square is not an integer.
        WoodsNotRepresentableError
            If matrix is valid but not Wood's-representable.
        WoodsSyntaxError
            If `string` does not match the structure of a
            Wood's notation, or could not be evaluated due
            to either unsupported, unmatched brackets or
            other any unexpected characters.
        """
        self.__style = _WoodsStyle(style)
        if matrix is not None and bulk_basis is None:
            raise TypeError('Woods: bulk_basis is mandatory when '
                            'a matrix is given')
        # Here the private attributes. Checking and processing
        # is deferred to the property setters called below.
        self.__string = ''
        self.__bulk_basis = None
        self.__matrix = None
        self.bulk_basis = bulk_basis
        if matrix is None:  # Only string given
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
        # Use the private attribute not to trigger a conversion
        if self.__matrix is not None:
            txt += f', matrix={array_to_string(self.__matrix)}'
        if self.style != 'unicode':  # Hide the default
            txt += f', style={self.style!r}'
        return f'Woods({txt})'

    def __str__(self):
        """Return a self.style-d string version of self."""
        return self.string

    def __format__(self, format_spec=''):
        """Format self according to a given `format_spec`.

        Parameters
        ----------
        format_spec : {'u', 'a', 'su', 'sa', 's', ''}, default=''
            'u' for Unicode formatting, 'a' for ASCII, '' for native,
            i.e., using self.style. Prefixing with 's' (i.e., simple)
            removes the 'p' prefix for primitive cells.

        Returns
        -------
        str

        Raises
        ------
        TypeError
            If `format_spec` is invalid
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

    @classmethod
    def add_example(cls, woods, shape):
        """Add an example Wood's notation for a given cell shape.

        Parameters
        ----------
        woods : str or Woods
            The example to add.
        shape : {'Oblique', 'Square', 'Rectangular', 'Hexagonal', 'Rhombic'}
            Shape of the surface unit cell.

        Returns
        -------
        None.
        """
        # Store examples always in Unicode style
        if isinstance(woods, str):
            woods = cls(woods)
        if not woods.style.startswith('u'):
            woods = cls(woods.string)
        cls.__examples[shape].add(woods)

    @classmethod
    def make_examples(cls):
        """Populate the __examples class attribute."""
        if cls.__examples:
            return
        for shape, examples in _EXAMPLES.items():
            cls.__examples[shape] = {cls(example) for example in examples}

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
        return (prefix, *gammas, alpha)

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
            If `gamma_str` could not be evaluated due to either
            unsupported math or unmatched brackets
        ValueError
            If `gamma_str` can be parsed, but its value is not the
            square root of an integer
        """
        parser = MathParser()

        # Given the regular expression we use, gamma2 may eat up the
        # closing parenthesis potentially separating the rotation
        if gamma_str.count('(')  == gamma_str.count(')') - 1:
            # Get rid of the swallowed closing
            # parenthesis, i.e., the last ')'
            idx = gamma_str.rindex(')')
            gamma_str = gamma_str[:idx] + gamma_str[idx+1:]

        # parser needs at least 'rt' for square root.
        # Support also 'r' for Woods notations
        try:
            parser.expression = re.sub(r'r([^t])', r'rt\1', gamma_str)
        except TooComplexMathError as exc:
            raise WoodsSyntaxError(f'Woods: {gamma_str} could '
                                   'not be parsed as it is too '
                                   'complex.') from exc
        try:
            gamma = parser.evaluate()
        except SyntaxError as exc:
            raise WoodsSyntaxError(f'Woods: {gamma_str} could '
                                   'not be parsed, likely because '
                                   'of unmatched brackets') from exc
        except (UnsupportedMathError, TooComplexMathError) as exc:
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
        self.__set_bulk_basis(basis)

    @property
    def matrix(self):
        """Return the matrix representation of self."""
        if self.__matrix is not None:
            return self.__matrix
        # Both to_matrix() and guess_correct_rotation() below will
        # store a valid, commensurate matrix into self.__matrix
        try:
            return self.to_matrix()
        except MatrixIncommensurateError as exc:
            # See if we can fix the rotation angle
            try:
                self.guess_correct_rotation()
            except MatrixIncommensurateError:
                raise exc from None
        return self.__matrix

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
        if matrix is None:
            self.__matrix = matrix
            return
        self.from_matrix(matrix)  # Also stores matrix if representable

    @property
    def string(self):
        """Return a self.style-d string representation of self."""
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
            If `woods` is not a string
        ValueError
            If `woods` is inappropriate for the current basis or
            the parsed string contains scaling factors for the
            two directions whose square is not an integer
        WoodsInvalidForBasisError
            If `woods` is incompatible with the current bulk basis.
        WoodsSyntaxError
            If `woods` does not match the structure of a
            Wood's notation, or could not be evaluated due
            to either unsupported math or unmatched brackets
        """
        old_string = self.__string
        new_string = self.__fix_string(woods)
        if new_string == old_string:
            return
        if not new_string or self.bulk_basis is None:
            self.__matrix = None  # No matrix with no string or basis
            self.__string = new_string
            return
        self.__string = new_string
        self.__matrix, old_matrix = None, self.__matrix
        try:  # Make sure new_string is OK for the current basis
            _ = self.matrix
        except MatrixIncommensurateError as exc:
            self.__string = old_string
            self.__matrix = old_matrix
            raise WoodsInvalidForBasisError(
                f'{woods!r} is incompatible with bulk basis '
                f'{array_to_string(self.bulk_basis)}. Gives '
                f'matrix={exc.matrix}.'
                ) from None

    @property
    def style(self):
        """Return the style of this instance."""
        return self.__style

    def from_matrix(self, matrix, bulk_basis=None):
        """Construct a Wood's string from matrix.

        If the conversion is successful, `matrix` and the string
        representation of the Wood's notation are also stored.

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
            If `matrix` represents a non-commensurate lattice, i.e., it
            is singular, or its elements are not integer-representable.
        ValueError
            If shape of `matrix` or `bulk_basis` is not (2, 2)
        WoodsNotRepresentableError
            If `matrix` is valid but not Wood's-representable.
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
            self.__set_bulk_basis(bulk_basis, check_matrix=False)

        # Next line may raise WoodsNotRepresentableError
        prefix, *gammas, alpha = self.__primitive_or_centered(matrix)

        # Don't use the setter of self.string to skip
        # reformatting, which is anyway already correct
        self.__string = self.__format_parsed(prefix, gammas, alpha)
        self.__matrix = matrix
        return self.string

    def get_examples(self, shape):
        """Return examples of Wood's for a given cell shape.

        The examples are formatted to the current style.

        Parameters
        ----------
        shape : {'Oblique', 'Square', 'Rectangular', 'Hexagonal', 'Rhombic'}
            Shape of the surface unit cell

        Returns
        -------
        examples : set of Woods
            Examples of Wood's notation for the given shape.
        """
        if self.style == 'unicode':
            return self.__examples[shape].copy()
        cls = type(self)
        return {cls(woods.string, style=self.style)
                for woods in self.__examples[shape]}

    def guess_correct_rotation(self, woods_txt=None):
        """Guess the rotation angle that gives an acceptable Woods.

        Typically, this would be used in case the user gives
        an angle that is off by (at most) one degree, and this
        results in a MatrixIncommensurateError. If the guess is
        successful, the corrected text and matrix are stored.

        Parameters
        ----------
        woods_txt : str
            The Woods string from which to guess the angle.

        Returns
        -------
        fixed_string : float
            Wood's notation with the new angle guessed, if successful

        Raises
        ------
        MatrixIncommensurateError
            If the angle guessed and the one passed are far off.
        """
        if woods_txt is None:
            woods_txt = self.string

        # Prepare a copy that will be used for messing around.
        # Notice that we DO NOT pass the bulk_basis at this stage
        # not to trigger a consistency check for woods_txt. The
        # check happens when we see if the matrix is commensurate
        tmp_woods = Woods(string=woods_txt, style=self.style)
        orig_prefix, *orig_gamma, orig_alpha = self.parse(tmp_woods.string)
        orig_matrix = tmp_woods.to_matrix(bulk_basis=self.bulk_basis,
                                          check_commensurate=False)
        if is_commensurate(orig_matrix):
            # Don't use the setter of self.string to avoid reformatting
            self.__string = tmp_woods.string
            self.__matrix = orig_matrix.round().astype(int)
            return self.string

        # Matrix is incommensurate. Try to round
        # it to int, and redo the parsing
        rounded_matrix = orig_matrix.round().astype(int)
        rounded_txt = tmp_woods.from_matrix(rounded_matrix)
        round_prefix, *round_gamma, round_alpha = self.parse(rounded_txt)
        delta_gamma = abs(np.subtract(orig_gamma, round_gamma)/orig_gamma)
        # Prefix should stay the same; same for gammas (the user rarely
        # inputs the wrong scaling); new angle is acceptable if it is
        # within +-1 degrees of the original input
        if (round_prefix != orig_prefix
            or any(delta_gamma > 1e-3)
                or abs(round_alpha - orig_alpha) > 1):
            raise MatrixIncommensurateError(orig_matrix)

        # Correction was successful
        self.__string = rounded_txt  # Private attribute skip reformat
        self.__matrix = rounded_matrix
        return self.string

    def to_matrix(self, woods='', bulk_basis=None, check_commensurate=True):
        """Convert woods notation into matrix representation.

        Upon successful conversion, `woods` and the new matrix are stored.

        Parameters
        ----------
        woods : str, optional
            A string-representation of a Wood's notation. If not given
            (or empty) self.string is used instead for the conversions.
            Default is an empty string.
        bulk_basis : Sequence or None, optional
            Basis vectors of the bulk. Shape (2, 2). If not given or
            None, self.bulk_basis is used, otherwise the value passed
            is stored into self.bulk_basis. Default is None.
        check_commensurate : bool, optional
            If True, raise MatrixIncommensurateError in case the
            resulting matrix is incommensurate. Notice that only
            a commensurate matrix is stored internally, whatever
            the value of `check_commensurate` is. Default is True.

        Returns
        -------
        matrix : numpy.ndarray
            Matrix representation. If `matrix` is commensurate,
            a rounded version is returned. The matrix is instead
            returned as-is if `check_commensurate` is False.

        Raises
        ------
        MatrixIncommensurateError
            If `check_commensurate` is True, and the matrix resulting
            from conversion gives an incommensurate lattice
        TypeError
            If `woods` is not a string
        ValueError
            It the parsed string contains scaling factors for
            the two directions whose square is not an integer
        WoodsSyntaxError
            If `woods` does not match the structure of a
            Wood's notation, or could not be evaluated due
            to either unsupported math or unmatched brackets
        """
        woods = woods or self.string
        if not woods:
            raise WoodsSyntaxError('Woods: empty string is not a valid Woods '
                                   'notation. Pass a valid, non-empty string, '
                                   'or set self.string beforehand')
        if bulk_basis is not None:
            self.__set_bulk_basis(bulk_basis, check_matrix=False)
        if self.bulk_basis is None:
            raise ValueError(f'{type(self).__name__}: Cannot convert {self} '
                             'to_matrix() without a bulk_basis. Pass a valid '
                             'bulk_basis, or set self.bulk_basis beforehand')
        prefix, gamma1, gamma2, alpha = self.parse(woods)  # May raise!
        woods = self.__format_parsed(prefix, (gamma1, gamma2), alpha)
        if woods != self.string:
            self.__string = woods
            self.__matrix = None

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

        if is_commensurate(matrix):
            self.__matrix = matrix.round().astype(int)
            return self.__matrix
        if not check_commensurate:
            return matrix
        raise MatrixIncommensurateError(matrix)

    def __fix_string(self, woods):
        """Return a woods string formatted to standard style.

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
            If `woods` is not a string
        ValueError
            If the parsed string contains scaling factors for
            the two directions whose square is not an integer
        WoodsSyntaxError
            If `woods` does not match the structure of a
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
        if abs(alpha - round(alpha)) < 0.1:
            alpha = round(alpha)
            fmt = 'd'
        else:
            fmt = '.1f'
        return f'{alpha:{fmt}}'

    def __format_scaling_factors(self, gammas):
        """Format scaling factors in Woods notation.

        Parameters
        ----------
        gammas : sequence
            Two elements, both integers. Woods notation
            is ...(gammas[0] <times> gammas[1])...

        Returns
        -------
        formatted_factors : str
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

        Make sure self.bulk_basis is up to date before calling this.

        Parameters
        ----------
        matrix : sequence
            The matrix to be checked. Shape (2, 2).

        Returns
        -------
        bool
        """
        basis = self.bulk_basis
        transformed_basis = np.dot(matrix, basis)

        basis_norm = np.linalg.norm(basis, axis=1)
        transf_norm = np.linalg.norm(transformed_basis, axis=1)
        gammas = transf_norm/basis_norm
        for gamma in gammas:
            if abs(gamma**2 - round(gamma**2)) > 1e-3:
                return False
        det = abs(np.linalg.det(matrix))
        return abs(det/np.prod(gammas) - 1) < 1e-8

    def __primitive_or_centered(self, matrix):
        """Check if matrix is representable as primitive or centered.

        Make sure self.bulk_basis is up to date before calling this.

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
            If `matrix` not Wood's-representable with neither
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
        alpha = np.arctan2(np.cross(basis[0], transformed_basis[0]),
                           np.dot(basis[0], transformed_basis[0]))
        return prefix, gamma1, gamma2, np.degrees(alpha)

    def __set_bulk_basis(self, basis, check_matrix=True):
        """Assign a new bulk_basis attribute."""
        if basis is None:
            self.__bulk_basis = basis
            self.__matrix = None  # Cannot have a matrix without basis
            return
        basis = np.asarray(basis)
        if basis.shape != (2, 2):
            raise ValueError('Woods: invalid bulk_basis. Should have '
                             f'shape==(2, 2). Found {basis.shape}')
        old_basis, old_matrix = self.bulk_basis, self.__matrix
        if old_basis is not None and np.allclose(basis, old_basis):
            return
        self.__bulk_basis = basis
        self.__matrix = None  # Outdated
        if not self.string or not check_matrix:
            return
        try:
            _ = self.matrix
        except MatrixIncommensurateError as exc:
            # Revert changes
            self.__bulk_basis, self.__matrix = old_basis, old_matrix
            raise ValueError(f'Invalid basis={array_to_string(basis)} for '
                             f'{self.string}. Would give an incommensurate '
                             f'matrix={exc.matrix}.') from None

# A bit of a trick to assign examples as class attributes. This
# class method will only do something the first time it is called
Woods.make_examples()
