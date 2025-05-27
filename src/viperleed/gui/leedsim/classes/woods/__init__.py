"""Package woods of viperleed.gui.leedsim.classes.

Defines the Woods class, used for conversion of a Wood's notation
string into a matrix and vice versa, as well as for formatting. It
also defines exceptions related to Wood's notation interpretation.

This package was refactored from the former woods.py module, and
exposes the same API.

Modules
-------
_woods
    The Woods class for interpreting a Wood's notation.
errors
    Exceptions that may occur when interpreting a Wood's notation.
utils
    Functions used for processing and interpreting a Wood's notation.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2024-03-01'
__license__ = 'GPLv3+'

from ._woods import Woods
from .errors import MatrixIncommensurateError
from .errors import WoodsError
from .errors import WoodsInvalidForBasisError
from .errors import WoodsNotRepresentableError
from .errors import WoodsSyntaxError
