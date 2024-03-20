"""Package parameters of viperleed.calc.files.

Contains definitions useful for reading, writing, editing, and
interpreting a PARAMETERS file. This package comes from a refactor
of the former parameters.py module.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__created__ = '2023-10-13'

from . import errors
from ._interpret import interpret, ParameterInterpreter
from ._read import read, update
from ._write import comment_out, modify
