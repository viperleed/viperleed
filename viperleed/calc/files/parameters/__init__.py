"""Package parameters of viperleed.calc.files.

Contains definitions useful for reading, writing, editing, and
interpreting a PARAMETERS file. This package comes from a refactor
of the former parameters.py module.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-10-13'
__license__ = 'GPLv3+'

from . import errors
from .interpret import interpret, ParameterInterpreter
from .read import read, update
from .write import comment_out, modify
