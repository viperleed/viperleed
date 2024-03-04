"""Module mathparse of viperleed.guilib.

======================================
  ViPErLEED Graphical User Interface
======================================

Defines a MathParser class that allows parsing simple mathematical
expressions. The parser is safe with respect to code injection, as
it only accepts numeric expressions and square root calculations.

Author: Michele Riva
Created: 2021-06-13
"""

# The code is based on https://stackoverflow.com/questions/20748202/
# See https://greentreesnakes.readthedocs.io/en/latest/index.html
# for a nice description of ast, as the official python docs are
# a bit tight...

__all__ = ('MathParser', 'TooComplexMathError', 'UnsupportedMathError')

import ast
import math
import operator
import re

BINARY_OPERATIONS = {
    ast.Add: operator.add,
    ast.Sub: operator.sub,
    ast.Mult: operator.mul,
    ast.Div: operator.truediv,
    ast.Mod: operator.mod,
    ast.Pow: operator.pow,
    }

UNARY_OPERATIONS = {
    ast.USub: operator.neg,
    ast.UAdd: operator.pos,
    }

MATH_OPERATIONS = {
    'sqrt': math.sqrt,
    }

MAX_LEN = 10000  # characters
MAX_POW = 1000   # Limit range of operands in pow


def _fix_expression(txt):
    """Fix a math expression such that it can be evaluated.

    The following fixes are considered:
    (1) remove all white-space characters
    (2) replace square and curly braces with round ones
    (3) add '*' whenever a digit and a letter are next to each other
        and around parentheses if there are no operators yet
    (4) replace 'rt' and '\u221a' with 'sqrt(...)', with brackets
        surrounding any group of digits that follows

    This means that, e.g., '2[3+2*rt 95] rt3 -9sqrt(7)/\u221a10'
    becomes '2*(3+2*sqrt(95))*sqrt(3)-9*sqrt(7)/sqrt(10)'

    Parameters
    ----------
    txt : str
        The text to be fixed

    Returns
    -------
    fixed_txt : str
        Expression suitable to be evaluated by MathParser
    """
    remove_spaces = r'\s*'
    txt = re.sub(remove_spaces, '', txt)

    txt = _brackets_to_parentheses(txt)
    txt = _fix_sqrt(txt)
    txt = _fix_multiplication(txt)

    return txt


def _brackets_to_parentheses(txt):
    """Replace '['/'{' [and closed counterparts] with '(' [and ')']."""
    replace_left_brackets = r'[\{\[]'
    txt = re.sub(replace_left_brackets, '(', txt)
    replace_right_brackets = r'[\}\]]'
    txt = re.sub(replace_right_brackets, ')', txt)

    return txt


def _fix_sqrt(txt):
    """Replace 'rt' and '\u221a' with 'sqrt' and fix parentheses."""
    # Replace 'rt' (not preceded by 'sq')
    # or '\u221a' with 'sqrt'
    replace_rt = r'((?<!sq)rt)|\u221a'
    txt = re.sub(replace_rt, r'sqrt', txt)

    # Add parentheses surrounding numbers
    # that follow a 'sqrt'
    add_parentheses = r'(sqrt)([\d]+)'
    txt = re.sub(add_parentheses, r'\1(\2)', txt)

    return txt


def _fix_multiplication(txt):
    """Add a '*' wherever it's missing.

    This means:
        (1) between a number and anything that is not
            an operator or another number
        (2) around parentheses, unless other operators
            or parentheses are already present

    Parameters
    ----------
    txt : str
        The string to be fixed

    Returns
    -------
    fixed_txt : str
        Text with multiplication signs added where needed.
    """
    # Add a star whenever a digit is followed by any
    # character except a digit, an operator, or a parenthesis
    insert_star = r'((\d)([^+\-()*/.\d]))'
    txt = re.sub(insert_star, r'\2*\3', txt)

    # And again stars, this time any time there is
    # an '(' not preceded by an operator or a letter,
    # and any time there is an ')' followed by
    # a letter or a number
    # add_star_open_parentheses = r'(.+)((?<!t|[+\-*/])[(])'
    add_star_open_parentheses = r'(.+)((?<![a-zA-Z+\-*/])[(])'
    txt = re.sub(add_star_open_parentheses, r'\1*\2', txt)
    add_star_close_parentheses = r'(\))([a-zA-Z\d])'
    txt = re.sub(add_star_close_parentheses, r'\1*\2', txt)

    return txt


class TooComplexMathError(ValueError):
    """Expression is too long or contains too large numbers."""


class UnsupportedMathError(Exception):
    """Unsupported expression in math."""

    def __init__(self, node):
        """Initialize exception.

        Parameters
        ----------
        node : object

        Returns
        -------
        None.
        """
        err = f'Unsupported ast node/operation {node}'
        super().__init__(err)


class MathParser:
    """A simple math-expression parser.

    Especially, any arithmetic operation is allowed,
    as well as sqrt. The parsed expressions are
    preprocessed, so that the only requirement for
    a successful parsing is having matching brackets
    to start with.
    """

    def __init__(self, expression=''):
        """Initialize instance.

        Parameters
        ----------
        expression : str, optional
            The mathematical expression to be parsed. Can be
            read/written via the .expression property. The
            expression is manipulated to bring it to a standard
            format that can then be evaluated.

        Raises
        ------
        TooComplexMathError
            If `expression` contains too many characters
            after preprocessing.
        TypeError
            If `expression` is not a string.
        """
        self._expression = ''
        self.expression = expression  # Use setter for processing

    @property
    def expression(self):
        """Return the expression that will be evaluated as a string."""
        return self._expression

    @expression.setter
    def expression(self, new_expression):
        """Set a new expression to be evaluated.

        Parameters
        ----------
        new_expression : str
            Expression to be evaluated. Call .evaluate() to
            get the numeric value of the expression.

        Raises
        ------
        TooComplexMathError
            If `new_expression` contains too many characters
            after preprocessing.
        TypeError
            If `new_expression` is not a string.
        """
        if not isinstance(new_expression, str):
            raise TypeError('MathParser: expression should be a string, '
                            f'not {type(new_expression).__name__!r}')
        fixed_expression = _fix_expression(new_expression)
        if len(fixed_expression) > MAX_LEN:
            raise TooComplexMathError('Too many characters '
                                      f'({len(fixed_expression)}>{MAX_LEN}) '
                                      'in expression after preprocessing')
        self._expression = fixed_expression

    def evaluate(self):
        """Evaluate self.expression.

        Returns
        -------
        number : float
            Numeric value of the expression.

        Raises
        ------
        SyntaxError
            If the expression passed contains errors,
            typically unmatched brackets.
        TooComplexMathError
            If the expression would evaluate to a too-large number.
        UnsupportedMathError
            If the expression contains calls to mathematical
            functions that are not supported.
        """
        node = ast.parse(self.expression, mode='eval')
        return self._eval(node)

    def _eval(self, node):
        """Recursively evaluate a node."""
        if isinstance(node, ast.Expression):
            ret = self._eval(node.body)
        elif isinstance(node, ast.Constant):  # Available since 3.6
            ret = node.value
        elif isinstance(node, ast.Str):       # Deprecated since 3.8
            ret = node.s
        elif isinstance(node, ast.Num):       # Deprecated since 3.8
            ret = node.n
        elif isinstance(node, ast.BinOp):
            ret = self._eval_binary_operation(node)
        elif isinstance(node, ast.UnaryOp):
            ret = self._eval_unary_operation(node)
        elif isinstance(node, ast.Call):
            ret = self._eval_math_function(node)
        else:
            raise UnsupportedMathError(node)
        return ret

    def _eval_math_function(self, node):
        """Evaluate a node containing the call to a function.

        Supported mathematical functions are:
         * sqrt(x)

        Parameters
        ----------
        node : ast.node
            The node to be evaluated. It is NOT checked to
            actually be a node describing a math function.

        Return
        ------
        number : float
            Result of evaluating the mathematical function

        Raises
        ------
        UnsupportedMathError
            If the mathematical function is not supported.
        """
        try:
            operation = MATH_OPERATIONS[node.func.id]
        except KeyError as exc:
            raise UnsupportedMathError(node.func.id) from exc
        operands = (self._eval(op) for op in node.args)
        return operation(*operands)

    def _eval_binary_operation(self, node):
        """Evaluate a node containing a binary operation.

        Supported binary operations are addition, subtraction,
        multiplication, division, modulo division, and raising
        to power.

        Parameters
        ----------
        node : ast.node
            The node to be evaluated. It is NOT checked to
            actually be a node describing a binary operation.

        Return
        ------
        number : float
            Numeric value of the binary operation

        Raises
        ------
        UnsupportedMathError
            If the binary operation is not one of the supported ones
        """
        try:
            operation = BINARY_OPERATIONS[type(node.op)]
        except KeyError as exc:
            raise UnsupportedMathError(node.op) from exc
        operands = self._eval(node.left), self._eval(node.right)
        if (operation is operator.pow  # Limit ranges for a**b
                and any(op > MAX_POW for op in operands)):
            raise TooComplexMathError('Result too large. Cannot compute '
                                      f'{operands[0]}**{operands[1]}.')
        return operation(*operands)

    def _eval_unary_operation(self, node):
        """Evaluate a node containing a unary operation.

        Supported unary operations are sign change (-x) and
        no-sign-change (+x).

        Parameters
        ----------
        node : ast.node
            The node to be evaluated. It is NOT checked to
            actually be a node describing a unary operation.

        Return
        ------
        number : float
            Numeric value of the unary operation

        Raises
        ------
        UnsupportedMathError
            If the unary operation is not one of the supported ones
        """
        try:
            operation = UNARY_OPERATIONS[type(node.op)]
        except KeyError as exc:
            raise UnsupportedMathError(node.op) from exc
        return operation(self._eval(node.operand))


if __name__ == '__main__':
    # Run some tests
    parser = MathParser('2*sqrt[2]')
    print(parser.evaluate())

    parser.expression = '2*2**0.5'
    print(parser.evaluate())

    parser.expression = '2[3+2*rt 95] rt3 -9sqrt(7) / \u221a3'
    print(parser.expression)
    print(parser.evaluate())

    parser.expression = '3.254*rt2'
    print(parser.expression)
    print(parser.evaluate())
