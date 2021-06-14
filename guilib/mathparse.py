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

__all__ = ['UnsupportedMathError', 'MathParser']

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
    """
    remove_spaces = r'\s*'
    txt = re.sub(remove_spaces, '', txt)

    txt = __brackets_to_parentheses(txt)
    txt = __fix_sqrt(txt)
    txt = __fix_multiplication(txt)

    return txt


def __brackets_to_parentheses(txt):
    """Replace '['/'{' [and closed counterparts] with '(' [and ')']."""
    replace_left_brackets = r'[\{\[]'
    txt = re.sub(replace_left_brackets, '(', txt)
    replace_right_brackets = r'[\}\]]'
    txt = re.sub(replace_right_brackets, ')', txt)

    return txt


def __fix_sqrt(txt):
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


def __fix_multiplication(txt):
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

        Returns
        -------
        None.
        """
        self.__expression = _fix_expression(expression)

    @property
    def expression(self):
        """Return the expression that will be evaluated as a string."""
        return self.__expression

    @expression.setter
    def expression(self, new_expression):
        """Set a new expression to be evaluated.

        Parameters
        ----------
        new_expression : str
        """
        if not isinstance(new_expression, str):
            raise TypeError("MathParser: expression should be a string, "
                            f"not {type(new_expression).__name__!r}")
        self.__expression = _fix_expression(new_expression)

    def evaluate(self):
        """Evaluate the expression given at instantiation.

        Returns
        -------
        number : float

        Raises
        ------
        SyntaxError
            If the expression passed contains errors,
            typically unmatched brackets
        UnsupportedMathError
            If the expression contains calls to mathematical
            functions that are not supported
        """
        node = ast.parse(self.expression, mode='eval')
        return self.__eval(node)

    def __eval(self, node):
        """Recursively evaluate a node."""
        if isinstance(node, ast.Expression):
            ret = self.__eval(node.body)

        elif isinstance(node, ast.Constant):  # Available since 3.6
            ret = node.value

        elif isinstance(node, ast.Str):       # Deprecated since 3.8
            ret = node.s

        elif isinstance(node, ast.Num):       # Deprecated since 3.8
            ret = node.n

        elif isinstance(node, ast.BinOp):
            try:
                operation = BINARY_OPERATIONS[type(node.op)]
            except KeyError as err:
                raise UnsupportedMathError(node.op) from err
            ret = operation(self.__eval(node.left), self.__eval(node.right))

        elif isinstance(node, ast.UnaryOp):
            try:
                operation = UNARY_OPERATIONS[type(node.op)]
            except KeyError as err:
                raise UnsupportedMathError(node.op) from err
            ret = operation(self.__eval(node.operand))

        elif isinstance(node, ast.Call):
            try:
                operation = MATH_OPERATIONS[node.func.id]
            except KeyError as err:
                raise UnsupportedMathError(node.func.id) from err

            # Allow only single-argument math operations
            ret = operation(self.__eval(node.args[0]))

        else:
            raise UnsupportedMathError(node)

        return ret


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
