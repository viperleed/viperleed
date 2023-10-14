"""Package parameters of viperleed.tleedmlib.files.

Created on 2023-10-13

@author: Michele Riva (@michele-riva)

Contains definitions useful for reading, writing, editing, and
interpreting a PARAMETERS file. This package comes from a refactor
of the former parameters.py module.
"""

from . import errors
from ._interpret import interpretPARAMETERS, ParameterInterpreter
from ._read import read, update
from ._write import modify
