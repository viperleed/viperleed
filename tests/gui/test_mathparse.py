"""Tests for module mathparse of viperleed.gui."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2024-03-04'
__license__ = 'GPLv3+'

import pytest
from pytest_cases import parametrize

from viperleed.gui import mathparse
from viperleed.gui.mathparse import MathParser
from viperleed.gui.mathparse import TooComplexMathError
from viperleed.gui.mathparse import UnsupportedMathError


class TestMathParser:
    """Collection of tests for the MathParser class."""

    _expr = {
        '2*sqrt(2)': 2 * (2 ** 0.5),  # Simple
        '-5': -5,                     # Unary
        'sqrt(9)': 3,                 # sqrt
        '2*2**0.5': 2 * (2 ** 0.5),   # pow
        '(2+3)*4': 20,                # round brackets
        '1.5+2.5': 4,                 # float
        '3.254*sqrt(2)': 3.254 * (2 ** 0.5),             # float & sqrt
        '2*(3+2*sqrt(95))*sqrt(3)-9*sqrt(7)/sqrt(3)': (  # complex
            2*(3+2*(95**0.5))*(3**0.5) - 9*(7**0.5)/(3**0.5)
            ),
        }

    @parametrize('expr,result', _expr.items(), ids=_expr)
    def test_evaluate(self, expr, result):
        """Check correct evaluation of a mathematical expression."""
        parser = MathParser(expr)
        assert parser.evaluate() == result

    _processed_expr = {
        '2[3+2*rt 95] rt3 -9sqrt(7) / \u221a3': (
            2*(3+2*(95**0.5))*(3**0.5) - 9*(7**0.5)/(3**0.5)
            ),
        '-1(1+2)(3-4)': 3,
        '-1((1+2))(3-4)': 3,
        }

    @parametrize('expr,result', _processed_expr.items(), ids=_processed_expr)
    def test_fix_and_evaluate(self, expr, result):
        """Check correct evaluation of an expression that needs processing."""
        self.test_evaluate(expr, result)


class TestMathParserRaises:
    """Collection of tests for complaints from a MathParser."""

    def test_no_expression(self):
        """Check complaints when no expression is defined."""
        parser = MathParser()
        with pytest.raises(SyntaxError):
            parser.evaluate()

    def test_not_a_string(self):
        """Check complaints for a non-string expression."""
        with pytest.raises(TypeError):
            MathParser(1)

    _complex = (
        '9**9**9',  # Number too large
        )

    @parametrize(expr=_complex)
    def test_too_complex(self, expr):
        """Check complaints for too-complex mathematical expressions."""
        parser = MathParser(expr)
        with pytest.raises(TooComplexMathError):
            parser.evaluate()

    def test_too_long(self):
        """Check complaints for a too long expression."""
        backup, mathparse.MAX_LEN = mathparse.MAX_LEN, 50
        with pytest.raises(TooComplexMathError):
            MathParser('2*sqrt(2)' * 10)  # Exceeds MAX_LEN
        mathparse.MAX_LEN = backup
    
    def test_max_recursion(self):
        """Check complaints for a too deep expression."""
        backup, mathparse.MAX_DEPTH = mathparse.MAX_DEPTH, 3
        expr = '_'
        for _ in range(mathparse.MAX_DEPTH-1):
            expr = expr.replace('_', 'sqrt(_)')
        expr = expr.replace('_', '1')
        parser = MathParser(expr)
        with pytest.raises(TooComplexMathError):
            parser.evaluate()
        mathparse.MAX_DEPTH = backup

    _unsupported = (
        'sin(0)',      # No support for trigonometric
        'abc',         # No support for arbitrary strings
        '1<2',         # No support for ordering
        'ge(1,2)',     # No support for ordering
        'truth(1)',    # No support for unary truth
        'abs(1)',      # No support for absolute value
        'not 1',       # No support for unary not
        '1 << 8',      # No support for bit shift
        '~1',          # No support for bit inversion
        '__module__',  # No support for underscore
        )

    @parametrize(expr=_unsupported)
    def test_unsupported_operation(self, expr):
        """Check complaints for unsupported mathematical functions."""
        parser = MathParser()
        with pytest.raises(UnsupportedMathError):
            parser.expression = expr
            parser.evaluate()
